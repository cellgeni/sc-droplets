#!/usr/bin/env bash
# Run bustools' own whitelist (allowlist), one release per pixi environment, on BUS files written
# from public raw matrices (make_bus.py: one record per UMI, in bustools' barcode order).
#
#   bash oracles/bustools/run.sh OUT_DIR
#
# Needs pixi and uv. The environments and BUS files go to $BUSTOOLS_ORACLE_WORK (default
# $RUNNER_TEMP/bt or /tmp/bt). Writes OUT_DIR/<dataset>/<release>/{<case>.txt.gz,cases.tsv,
# versions.txt}: the barcodes written (sorted, once each), and per case the threshold bustools
# printed, the lines it wrote and the distinct barcodes among them. The committed outputs are
# linux-64's, as CI's.
set -euo pipefail

mkdir -p "$1"
out=$(cd "$1" && pwd)
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)
work="${BUSTOOLS_ORACLE_WORK:-${RUNNER_TEMP:-/tmp}/bt}"
mkdir -p "$work/bus"
cp "$here/pixi.toml" "$here/pixi.lock" "$work/"

datasets=(pbmc_1k_v3_3_0 pbmc4k_2_1)
# environment:release
releases=(bt0393:0.39.3 bt0410:0.41.0 bt0432:0.43.2 bt0451:0.45.1)
# case:extra arguments
cases=("default:" "threshold500:-f 500")

(cd "$work" && pixi install --all --locked)

for dataset in "${datasets[@]}"; do
  bus="$work/bus/$dataset.bus"
  (cd "$repo" && uv run python oracles/bustools/make_bus.py "$dataset" "$bus")
  for run in "${releases[@]}"; do
    IFS=: read -r env release <<< "$run"
    dir="$out/$dataset/$release"
    mkdir -p "$dir"
    : > "$dir/cases.tsv"
    for entry in "${cases[@]}"; do
      name=${entry%%:*}
      extra=${entry#*:}
      # shellcheck disable=SC2086  # extra holds separate arguments
      (cd "$work" && pixi run -e "$env" bustools whitelist $extra -o "$work/list.txt" "$bus") 2> "$work/log.txt"
      threshold=$(sed -n 's/.*with threshold \([0-9-]*\).*/\1/p' "$work/log.txt")
      lines=$(wc -l < "$work/list.txt" | tr -d ' ')
      sort -u "$work/list.txt" | gzip -n > "$dir/$name.txt.gz"
      unique=$(gzip -dc "$dir/$name.txt.gz" | wc -l | tr -d ' ')
      printf '%s\t%s\t%s\t%s\n' "$name" "$threshold" "$lines" "$unique" >> "$dir/cases.tsv"
      echo "$dataset $release $name: threshold $threshold, $lines lines, $unique barcodes"
    done
    (cd "$work" && pixi run -e "$env" bustools version) > "$dir/versions.txt" 2>&1 || true
    echo "platform $(uname -s) $(uname -m)" >> "$dir/versions.txt"
  done
done
