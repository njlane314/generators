#!/usr/bin/env python3
import argparse
import csv
import html
import math
import os
import sys


AXIS_COLORS = {
    "primary_mechanism": "#0072B2",
    "transport": "#009E73",
    "generator_configuration": "#D55E00",
    "generator_version": "#CC79A7",
    "generator_model": "#555555",
    "other": "#777777",
}


def parse_float(text, fallback=0.0):
    try:
        value = float(text)
    except (TypeError, ValueError):
        return fallback
    return value if math.isfinite(value) else fallback


def shorten(text, limit):
    text = " ".join(str(text).split())
    if len(text) <= limit:
        return text
    return text[: max(0, limit - 3)].rstrip() + "..."


def default_output(input_path):
    if input_path.endswith("_members.tsv"):
        return input_path[: -len("_members.tsv")] + "_ranked_impacts.svg"
    base, _ = os.path.splitext(input_path)
    return base + "_ranked_impacts.svg"


def compact_group(row):
    category = row.get("category", "")
    pieces = [piece for piece in row.get("group", "").split("|") if piece and piece != category]
    if pieces:
        return " ".join(pieces)
    pieces = [row.get("beam_mode", ""), row.get("beam_species", ""), row.get("interaction", "")]
    return " ".join(piece for piece in pieces if piece)


def read_records(input_path, metric):
    with open(input_path, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    records = []
    for row in rows:
        expected = parse_float(row.get("expected_events"))
        delta = parse_float(row.get("delta_from_central"))
        central = expected - delta
        stat_pull = parse_float(row.get("stat_pull"))
        frac = delta / central if abs(central) > 1.0e-12 else math.nan
        signed_pull = math.copysign(stat_pull, delta) if stat_pull != 0.0 else 0.0

        if metric == "frac":
            value = 100.0 * frac
        elif metric == "pull":
            value = signed_pull
        else:
            value = delta

        if not math.isfinite(value) or abs(value) <= 0.0:
            continue

        axis = row.get("model_axis", "") or "other"
        group = compact_group(row)
        label = row.get("label", "")
        if group:
            display = f"{group} | {label}"
        else:
            display = label

        records.append(
            {
                "axis": axis,
                "color": AXIS_COLORS.get(axis, AXIS_COLORS["other"]),
                "display": display,
                "full_label": display,
                "value": value,
                "frac_percent": 100.0 * frac if math.isfinite(frac) else math.nan,
                "delta": delta,
                "pull": stat_pull,
                "role": row.get("envelope_role", ""),
            }
        )

    records.sort(key=lambda item: abs(item["value"]), reverse=True)
    return records


def nice_limit(value):
    value = abs(value)
    if value <= 0.0 or not math.isfinite(value):
        return 1.0
    exponent = math.floor(math.log10(value))
    base = 10.0 ** exponent
    scaled = value / base
    for step in (1.0, 2.0, 5.0, 10.0):
        if scaled <= step:
            return step * base
    return 10.0 * base


def tick_values(limit):
    return [-limit, -0.5 * limit, 0.0, 0.5 * limit, limit]


def fmt_value(value, metric):
    if metric == "frac":
        return f"{value:+.1f}%"
    if metric == "pull":
        return f"{value:+.2f}"
    return f"{value:+.3g}"


def metric_axis_title(metric):
    if metric == "frac":
        return "Yield shift from central [%]"
    if metric == "pull":
        return "Signed statistical pull"
    return "Yield shift from central [expected events]"


def render_svg(records, args):
    width = args.width
    left = args.left_margin
    right = 190
    top = 118
    row_height = 30
    bottom = 74
    min_rows = max(1, len(records))
    height = top + bottom + row_height * min_rows
    plot_width = max(100, width - left - right)
    max_abs = max((abs(record["value"]) for record in records), default=1.0)
    limit = nice_limit(1.08 * max_abs)

    def x_pos(value):
        return left + (value + limit) / (2.0 * limit) * plot_width

    axis_title = metric_axis_title(args.metric)
    title = "Phenomenology sensitivity: largest yield shifts"
    subtitle = f"{os.path.basename(args.input)}; top {len(records)} rows by |{axis_title}|"

    lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        "<style>",
        "text{font-family:Arial,Helvetica,sans-serif;fill:#202020}",
        ".small{font-size:12px}.label{font-size:12px}.tick{font-size:11px;fill:#555}.title{font-size:22px;font-weight:bold}.subtitle{font-size:13px;fill:#555}",
        ".grid{stroke:#d8d8d8;stroke-width:1}.zero{stroke:#202020;stroke-width:1.4}.row{fill:#f7f7f7}.bar{shape-rendering:geometricPrecision}.marker{stroke:#202020;stroke-width:0.8}",
        "</style>",
        '<rect x="0" y="0" width="100%" height="100%" fill="#ffffff"/>',
        f'<text class="title" x="{left}" y="34">{html.escape(title)}</text>',
        f'<text class="subtitle" x="{left}" y="56">{html.escape(subtitle)}</text>',
        f'<text class="small" x="{left}" y="{height - 22}">{html.escape(axis_title)}</text>',
    ]

    legend_x = left
    legend_y = 80
    used_axes = []
    for record in records:
        if record["axis"] not in used_axes:
            used_axes.append(record["axis"])
    for axis in used_axes:
        color = AXIS_COLORS.get(axis, AXIS_COLORS["other"])
        lines.append(f'<rect x="{legend_x}" y="{legend_y - 10}" width="12" height="12" fill="{color}"/>')
        lines.append(f'<text class="small" x="{legend_x + 18}" y="{legend_y}">{html.escape(axis)}</text>')
        legend_x += 18 + 7 * len(axis) + 24

    plot_top = top - 12
    plot_bottom = top + row_height * min_rows - 5
    if records:
        for index in range(len(records)):
            if index % 2 == 1:
                y = top + index * row_height
                lines.append(f'<rect class="row" x="0" y="{y - 7}" width="{width}" height="{row_height}"/>')

    for tick in tick_values(limit):
        x = x_pos(tick)
        css = "zero" if abs(tick) < 1.0e-12 else "grid"
        lines.append(f'<line class="{css}" x1="{x:.1f}" y1="{plot_top}" x2="{x:.1f}" y2="{plot_bottom}"/>')
        lines.append(
            f'<text class="tick" text-anchor="middle" x="{x:.1f}" y="{plot_bottom + 18}">{html.escape(fmt_value(tick, args.metric))}</text>'
        )

    if not records:
        lines.append(f'<text class="label" x="{left}" y="{top + 20}">No non-zero yield shifts found.</text>')
    else:
        zero = x_pos(0.0)
        for index, record in enumerate(records):
            y = top + index * row_height
            center = y + 12

            label = shorten(record["display"], args.label_chars)
            lines.append(
                f'<text class="label" text-anchor="end" x="{left - 14}" y="{center + 4}">'
                f'<title>{html.escape(record["full_label"])}</title>{html.escape(label)}</text>'
            )

            x = x_pos(record["value"])
            bar_x = min(zero, x)
            bar_w = max(1.0, abs(x - zero))
            lines.append(
                f'<rect class="bar" x="{bar_x:.1f}" y="{center - 6}" width="{bar_w:.1f}" height="12" fill="{record["color"]}" opacity="0.86"/>'
            )
            marker_stroke = "#202020" if record["role"] in ("min", "max") else "#ffffff"
            marker_width = 2.0 if record["role"] in ("min", "max") else 0.8
            lines.append(
                f'<circle class="marker" cx="{x:.1f}" cy="{center}" r="4.2" fill="{record["color"]}" stroke="{marker_stroke}" stroke-width="{marker_width}"/>'
            )

            value_text = fmt_value(record["value"], args.metric)
            pull_text = f"pull {record['pull']:.2f}"
            if args.metric != "frac" and math.isfinite(record["frac_percent"]):
                pull_text += f", {record['frac_percent']:+.1f}%"
            if record["role"] in ("min", "max"):
                pull_text += f", {record['role']}"
            text_anchor = "start" if record["value"] >= 0.0 else "end"
            text_x = x + 8 if record["value"] >= 0.0 else x - 8
            lines.append(
                f'<text class="small" text-anchor="{text_anchor}" x="{text_x:.1f}" y="{center + 4}">{html.escape(value_text + " (" + pull_text + ")")}</text>'
            )

    lines.append("</svg>")
    return "\n".join(lines) + "\n"


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description="Plot ranked phenomenology sensitivity impacts from a *_members.tsv table."
    )
    parser.add_argument("--input", "-i", required=True, help="Sensitivity members TSV.")
    parser.add_argument("--output", "-o", help="Output SVG path. Defaults next to the input TSV.")
    parser.add_argument("--max-rows", type=int, default=24, help="Maximum rows to draw.")
    parser.add_argument("--metric", choices=("frac", "pull", "events"), default="frac")
    parser.add_argument("--width", type=int, default=1200)
    parser.add_argument("--left-margin", type=int, default=480)
    parser.add_argument("--label-chars", type=int, default=72)
    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv)
    if args.max_rows <= 0:
        raise SystemExit("--max-rows must be positive")
    if not os.path.exists(args.input):
        raise SystemExit(f"missing input: {args.input}")

    records = read_records(args.input, args.metric)[: args.max_rows]
    output = args.output or default_output(args.input)
    output_dir = os.path.dirname(output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    with open(output, "w", encoding="utf-8") as handle:
        handle.write(render_svg(records, args))
    print(output)


if __name__ == "__main__":
    main(sys.argv[1:])
