#!/usr/bin/env bash
# Run kneed's own KneeLocator, one release per behaviour in its own pixi environment, on the rank
# curves of public data.
#
#   bash oracles/kneedle/run.sh OUT_DIR
#
# Needs pixi and uv. The environments are built in $KNEEDLE_ORACLE_WORK (default
# $RUNNER_TEMP/kn or /tmp/kn). Writes OUT_DIR/<dataset>/<release>/{knees.tsv,versions.txt}
# (see oracle.py). Runs on linux-64 and osx-64; the committed outputs are linux-64's, as CI's.
set -euo pipefail

mkdir -p "$1"
out=$(cd "$1" && pwd)
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)
work="${KNEEDLE_ORACLE_WORK:-${RUNNER_TEMP:-/tmp}/kn}"
mkdir -p "$work"
cp "$here/pixi.toml" "$here/pixi.lock" "$work/"

datasets=(pbmc_1k_v3_3_0 pbmc4k_2_1)
# environment:release
runs=(kn060:0.6.0 kn085:0.8.5 kn086:0.8.6 kn090:0.9.0)

(cd "$work" && pixi install --all --locked)

for dataset in "${datasets[@]}"; do
  counts="$work/counts/$dataset.tsv"
  (cd "$repo" && uv run python oracles/make_counts.py "$dataset" "$counts")
  for run in "${runs[@]}"; do
    IFS=: read -r env release <<< "$run"
    (cd "$work" && pixi run -e "$env" python "$here/oracle.py" "$counts" "$out/$dataset/$release")
  done
done
