#!/usr/bin/env bash
# Build CellBender's own SingleCellRNACountsDataset, one release per run, from its tagged source,
# on public raw matrices.
#
#   bash oracles/cellbender/run.sh OUT_DIR
#
# Needs pixi, uv and git. The environments and the checkouts go to $CELLBENDER_ORACLE_WORK (default
# $RUNNER_TEMP/cb or /tmp/cb). Writes OUT_DIR/<dataset>/<release>/{priors.tsv,versions.txt}
# (see oracle.py). The committed outputs are linux-64's, as CI's.
set -euo pipefail

mkdir -p "$1"
out=$(cd "$1" && pwd)
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)
work="${CELLBENDER_ORACLE_WORK:-${RUNNER_TEMP:-/tmp}/cb}"
mkdir -p "$work/src"
cp "$here/pixi.toml" "$here/pixi.lock" "$work/"

# dataset:expected cells
datasets=(pbmc_1k_v3_3_0:1000 pbmc4k_2_1:5000)
# release:environment (0.2.2 stands for 0.2.0 – 0.2.2, whose prior code is the same)
releases=(0.1.0:py37 0.2.2:py37 0.3.0:py311 0.3.2:py311 0.4.0:py311)

(cd "$work" && pixi install --all --locked)

for entry in "${releases[@]}"; do
  release=${entry%%:*}
  src="$work/src/v$release"
  if [ ! -d "$src" ]; then
    git -c advice.detachedHead=false clone --quiet --depth 1 --branch "v$release" https://github.com/broadinstitute/CellBender.git "$src"
  fi
done

for entry in "${datasets[@]}"; do
  IFS=: read -r dataset expected <<< "$entry"
  raw=$(cd "$repo" && uv run python -c 'import sys; sys.path.insert(0, "tests"); from datasets import raw_path; print(raw_path(sys.argv[1]))' "$dataset")
  for run in "${releases[@]}"; do
    IFS=: read -r release env <<< "$run"
    src="$work/src/v$release"
    (cd "$work" && PYTHONPATH="$src" pixi run -e "$env" python "$here/oracle.py" "$release" "$raw" "$out/$dataset/$release" "$expected")
    echo "source $(git -C "$src" rev-parse HEAD)" >> "$out/$dataset/$release/versions.txt"
  done
done
