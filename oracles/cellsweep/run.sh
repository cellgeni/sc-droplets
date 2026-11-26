#!/usr/bin/env bash
# Run cellsweep's own infer_empty_droplets, one release (and, for 1.0.2, one kneed) per pixi
# environment, on public data.
#
#   bash oracles/cellsweep/run.sh OUT_DIR
#
# Needs pixi and uv. The environments are built in $CELLSWEEP_ORACLE_WORK (default
# $RUNNER_TEMP/cs or /tmp/cs). Writes OUT_DIR/<dataset>/<release>/{cutoffs.tsv,versions.txt}
# (see oracle.py). The committed outputs are linux-64's, as CI's.
set -euo pipefail

mkdir -p "$1"
out=$(cd "$1" && pwd)
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)
work="${CELLSWEEP_ORACLE_WORK:-${RUNNER_TEMP:-/tmp}/cs}"
mkdir -p "$work"
cp "$here/pixi.toml" "$here/pixi.lock" "$work/"

# dataset:expected cells
datasets=(pbmc_1k_v3_3_0:1000 pbmc4k_2_1:5000)
# environment:release (1.0.2 with each kneed it may run with)
runs=(cs101:1.0.1 cs102-kneed085:1.0.2-kneed0.8.5 cs102-kneed086:1.0.2-kneed0.8.6 cs102-kneed090:1.0.2-kneed0.9.0)

(cd "$work" && pixi install --all --locked)

for entry in "${datasets[@]}"; do
  IFS=: read -r dataset expected <<< "$entry"
  counts="$work/counts/$dataset.tsv"
  (cd "$repo" && uv run python oracles/make_counts.py "$dataset" "$counts" --all)
  for run in "${runs[@]}"; do
    IFS=: read -r env release <<< "$run"
    (cd "$work" && pixi run -e "$env" python "$here/oracle.py" "$counts" "$out/$dataset/$release" "$expected")
  done
done
