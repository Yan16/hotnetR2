#!/usr/bin/env bash
set -euo pipefail

analysis_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
config_file="$analysis_root/analysis.yaml"
dispatcher="$analysis_root/run_analysis.R"
stage="${1:-check}"
expected_version="${HOTNETR2_VERSION:-0.2.1}"
r_script="${R_SCRIPT:-Rscript}"

case "$stage" in
  install|check|all|annotations|classification|ldak|scores|networks|hotnet|exports) ;;
  -h|--help|help)
    echo "Usage: bash run_analysis.sh {install|check|all|annotations|classification|ldak|scores|networks|hotnet|exports} [tarball]"
    exit 0
    ;;
  *)
    echo "Unknown stage: $stage" >&2
    exit 2
    ;;
esac

[[ -f "$config_file" ]] || { echo "Missing $config_file" >&2; exit 2; }
[[ -f "$dispatcher" ]] || { echo "Missing $dispatcher" >&2; exit 2; }
command -v "$r_script" >/dev/null 2>&1 || {
  echo "Rscript executable not found: $r_script" >&2
  exit 127
}

if [[ "$stage" == "install" ]]; then
  tarball="${2:-${HOTNETR2_TARBALL:-}}"
  [[ -n "$tarball" && -f "$tarball" ]] || {
    echo "Supply a hotnetR2 source tarball." >&2
    exit 2
  }
  exec R CMD INSTALL "$tarball"
fi

if (( $# > 1 )); then
  echo "Only install accepts a tarball argument." >&2
  exit 2
fi

if [[ "$stage" == "check" ]]; then
  exec "$r_script" --vanilla "$dispatcher" "$config_file" "$stage" "$expected_version"
fi

mkdir -p "$analysis_root/logs"
log_file="$(mktemp "$analysis_root/logs/${stage}.XXXXXX.log")"
echo "Log: $log_file"
"$r_script" --vanilla "$dispatcher" "$config_file" "$stage" "$expected_version" 2>&1 | tee "$log_file"
