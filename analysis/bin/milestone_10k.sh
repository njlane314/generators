#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
. "${script_dir}/../lib/common.sh"
exec "${script_dir}/pipeline.sh" "${ANA_SHARE_DIR}/config/milestone_10k.env"
