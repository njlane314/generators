BEGIN {
  FS = "\n"
  OFS = "\t"
  if (category == "") category = "detector_visible_Lambda_to_p_piminus"
  if (grouping == "") grouping = "beam"
  if (min_pull == "") min_pull = 2.0
  if (min_frac == "") min_frac = 0.20

  print "group", "beam_mode", "beam_species", "interaction", "category", \
        "n_members", "central_label", "central_expected_events", \
        "central_stat_uncertainty", "min_expected_events", "min_label", \
        "max_expected_events", "max_label", "envelope_range", \
        "envelope_frac_of_central", "max_abs_delta", "max_abs_delta_frac", \
        "max_stat_pull", "max_stat_pull_label", "sensitive_by_stat", \
        "sensitive_by_fraction" > summary

  print "group", "beam_mode", "beam_species", "interaction", "category", \
        "model_axis", "label", "generator", "version", "variation", "knob", \
        "fsi_state", "expected_events", "stat_uncertainty", \
        "delta_from_central", "ratio_to_central", "stat_pull", \
        "envelope_role", "leave_one_out_min", "leave_one_out_max", \
        "leave_one_out_range", "exclusion_contracts_range_by" > members

  print "group", "beam_mode", "beam_species", "interaction", "category", \
        "model_axis", "n_members", "min_expected_events", "min_label", \
        "max_expected_events", "max_label", "axis_range", \
        "max_abs_delta_vs_central", "max_stat_pull_vs_central" > axes
}

NR == 1 {
  nf = csvsplit($0, h)
  for (i = 1; i <= nf; ++i) col[h[i]] = i
  next
}

{
  nf = csvsplit($0, f)
  row_category = get("category")
  if (category != "all" && row_category != category) next

  id = ++n
  sample[id] = get("sample")
  gen[id] = get("generator")
  version[id] = get("version")
  variation[id] = get("variation")
  knob[id] = get("knob")
  bm[id] = get("beam_mode")
  bs[id] = get("beam_species")
  interaction[id] = get("interaction")
  fsi[id] = get("fsi_state")
  cat[id] = row_category
  y[id] = get("expected_events") + 0
  stat[id] = get("expected_stat_uncertainty") + 0
  label[id] = mklabel(id)
  axis[id] = model_axis(id)

  g = group_key(id)
  groups[g] = 1
  g_bm[g] = group_value("beam_mode", id)
  g_bs[g] = group_value("beam_species", id)
  g_int[g] = group_value("interaction", id)
  g_cat[g] = row_category
  group_n[g]++
  group_id[g, group_n[g]] = id

  a = g SUBSEP axis[id]
  axis_seen[a] = 1
  axis_group[a] = g
  axis_name[a] = axis[id]
  axis_n[a]++
  axis_id[a, axis_n[a]] = id
}

END {
  for (g in groups) write_group(g)
}

function csvsplit(s, a,    i, c, n, q, nextc, field) {
  for (i in a) delete a[i]
  n = 1
  q = 0
  field = ""
  for (i = 1; i <= length(s); ++i) {
    c = substr(s, i, 1)
    if (q) {
      if (c == "\"") {
        nextc = substr(s, i + 1, 1)
        if (nextc == "\"") {
          field = field "\""
          ++i
        } else {
          q = 0
        }
      } else {
        field = field c
      }
    } else if (c == "\"") {
      q = 1
    } else if (c == ",") {
      a[n++] = field
      field = ""
    } else {
      field = field c
    }
  }
  a[n] = field
  return n
}

function get(name) { return f[col[name]] }

function abs(x) { return x < 0 ? -x : x }

function clean(x) {
  gsub(/[\t\r\n]+/, " ", x)
  return x
}

function lower(x) {
  return tolower(x)
}

function mklabel(id,    out) {
  out = gen[id] " " version[id] " " variation[id]
  if (knob[id] != "" && index(variation[id], knob[id]) == 0) out = out " " knob[id]
  if (fsi[id] != "") out = out " " fsi[id]
  return clean(out)
}

function model_axis(id,    text, k) {
  text = lower(gen[id] " " version[id] " " variation[id] " " knob[id] " " fsi[id])
  k = lower(knob[id] " " variation[id])
  if (fsi[id] != "" && lower(fsi[id]) != "fsi_on") return "transport"
  if (text ~ /fsi_off|transport|cascade|formation|gibuu/) return "transport"
  if (k ~ /dis|hyp|hyperon|lambda|sigma|xi|omega|strange|kaon|charm/) return "primary_mechanism"
  if (text ~ /g18|ar23|tune|config/) return "generator_configuration"
  if (version[id] != "") return "generator_version"
  return "generator_model"
}

function group_value(field, id) {
  if (grouping == "all_beam" && (field == "beam_mode" || field == "beam_species")) return "all"
  if (grouping == "beam_envelope" && field == "beam_mode") return "all"
  if (grouping == "generator" && field == "generator") return gen[id]
  if (field == "beam_mode") return bm[id]
  if (field == "beam_species") return bs[id]
  if (field == "interaction") return interaction[id]
  return ""
}

function group_key(id,    key) {
  if (grouping == "all_beam") return "all|all|" interaction[id] "|" cat[id]
  if (grouping == "beam_envelope") return "all|" bs[id] "|" interaction[id] "|" cat[id]
  if (grouping == "generator") return gen[id] "|" bm[id] "|" bs[id] "|" interaction[id] "|" cat[id]
  return bm[id] "|" bs[id] "|" interaction[id] "|" cat[id]
}

function central_score(id,    text, score) {
  text = lower(label[id])
  score = 100000 + id
  if (text ~ /nominal|central|ar23_20i_00_000/) score = 0
  else if (text ~ /all_strange/) score = 100
  else if (text ~ /dis_only/) score = 200
  else if (text ~ /hyp_all/) score = 300
  if (text ~ /fsi_off/) score += 25
  return score + id / 1000000.0
}

function write_group(g,    i, id, c, best, minid, maxid, min2, max2, full, maxd, maxp, maxpid, d, p, den, frac, role, contraction) {
  if (group_n[g] < 1) return
  c = group_id[g, 1]
  best = central_score(c)
  minid = maxid = c

  for (i = 1; i <= group_n[g]; ++i) {
    id = group_id[g, i]
    if (central_score(id) < best) { c = id; best = central_score(id) }
    if (y[id] < y[minid]) minid = id
    if (y[id] > y[maxid]) maxid = id
  }

  full = y[maxid] - y[minid]
  maxd = maxp = 0
  maxpid = c
  for (i = 1; i <= group_n[g]; ++i) {
    id = group_id[g, i]
    d = y[id] - y[c]
    den = sqrt(stat[id] * stat[id] + stat[c] * stat[c])
    p = den > 0 ? abs(d) / den : 0
    if (abs(d) > abs(maxd)) maxd = d
    if (p > maxp) { maxp = p; maxpid = id }
  }

  frac = y[c] != 0 ? full / y[c] : 0
  print g, g_bm[g], g_bs[g], g_int[g], g_cat[g], group_n[g], label[c], y[c], stat[c], \
        y[minid], label[minid], y[maxid], label[maxid], full, frac, abs(maxd), \
        (y[c] != 0 ? abs(maxd) / y[c] : 0), maxp, label[maxpid], \
        (maxp >= min_pull ? "yes" : "no"), (frac >= min_frac ? "yes" : "no") >> summary

  for (i = 1; i <= group_n[g]; ++i) {
    id = group_id[g, i]
    d = y[id] - y[c]
    den = sqrt(stat[id] * stat[id] + stat[c] * stat[c])
    p = den > 0 ? abs(d) / den : 0
    role = id == minid ? "min" : (id == maxid ? "max" : "interior")
    leave_one_out(g, id)
    contraction = full - loo_range
    print g, g_bm[g], g_bs[g], g_int[g], g_cat[g], axis[id], label[id], \
          gen[id], version[id], variation[id], knob[id], fsi[id], y[id], stat[id], \
          d, (y[c] != 0 ? y[id] / y[c] : 0), p, role, loo_min, loo_max, \
          loo_range, contraction >> members
  }

  for (a in axis_seen) if (axis_group[a] == g) write_axis(g, a, c)
}

function leave_one_out(g, drop,    i, id, seen) {
  seen = 0
  loo_min = loo_max = loo_range = 0
  for (i = 1; i <= group_n[g]; ++i) {
    id = group_id[g, i]
    if (id == drop) continue
    if (!seen) {
      loo_min = loo_max = y[id]
      seen = 1
    } else {
      if (y[id] < loo_min) loo_min = y[id]
      if (y[id] > loo_max) loo_max = y[id]
    }
  }
  loo_range = seen ? loo_max - loo_min : 0
}

function write_axis(g, a, c,    i, id, minid, maxid, d, p, den, maxd, maxp) {
  if (axis_n[a] < 1) return
  minid = maxid = axis_id[a, 1]
  maxd = maxp = 0
  for (i = 1; i <= axis_n[a]; ++i) {
    id = axis_id[a, i]
    if (y[id] < y[minid]) minid = id
    if (y[id] > y[maxid]) maxid = id
    d = y[id] - y[c]
    den = sqrt(stat[id] * stat[id] + stat[c] * stat[c])
    p = den > 0 ? abs(d) / den : 0
    if (abs(d) > abs(maxd)) maxd = d
    if (p > maxp) maxp = p
  }
  print g, g_bm[g], g_bs[g], g_int[g], g_cat[g], axis_name[a], axis_n[a], \
        y[minid], label[minid], y[maxid], label[maxid], y[maxid] - y[minid], \
        abs(maxd), maxp >> axes
}
