#!/usr/bin/env bash
# Run UMI-tools' own whitelist functions, one release per pixi environment, on public data.
#
#   bash oracles/umitools/run.sh OUT_DIR
#
# Needs pixi and uv. The environments are built in $UMITOOLS_ORACLE_WORK (default
# $RUNNER_TEMP/ut or /tmp/ut): old conda packages cannot be installed under a long path. Writes
# OUT_DIR/<dataset>/<release>/{<case>.txt.gz,versions.txt} (see oracle.py). Runs on linux-64 and
# osx-64 (Apple silicon through Rosetta); the committed outputs are linux-64's, as CI's.
set -euo pipefail

mkdir -p "$1"
out=$(cd "$1" && pwd)
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)
work="${UMITOOLS_ORACLE_WORK:-${RUNNER_TEMP:-/tmp}/ut}"
mkdir -p "$work"
cp "$here/pixi.toml" "$here/pixi.lock" "$work/"

# dataset:--expect-cells values for density (one that picks the default's valley, one that picks
# another valley or none)
datasets=(pbmc_1k_v3_3_0:1500,20 pbmc4k_2_1:5000,1500)
# environment:release
runs=(ut055:0.5.5 ut100:1.0.0 ut116:1.1.6)

(cd "$work" && pixi install --all --locked)

for entry in "${datasets[@]}"; do
  IFS=: read -r dataset expect <<< "$entry"
  counts="$work/counts/$dataset.tsv"
  (cd "$repo" && uv run python oracles/umitools/make_counts.py "$dataset" "$counts")
  for run in "${runs[@]}"; do
    IFS=: read -r env release <<< "$run"
    (cd "$work" && pixi run -e "$env" python "$here/oracle.py" "$counts" "$out/$dataset/$release" "$expect")
  done
  find "$out/$dataset" -name '*.txt' ! -name versions.txt -exec gzip -f {} +
done
