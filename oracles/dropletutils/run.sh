#!/usr/bin/env bash
# Run DropletUtils' cell calling, release by release, on public raw matrices and save what
# tests/test_validation.py compares sc-droplets with.
#
#     bash oracles/dropletutils/run.sh OUT_DIR
#
# Needs pixi. Creates one R environment per release from pixi.toml (bioconda builds), in a
# short folder: old conda packages cannot be installed under a long path. Writes
# OUT_DIR/<dataset>/<release>/{barcoderanks.txt,defaultdrops.txt,emptydrops.tsv,
# emptydropscellranger.tsv}. Runs on linux-64 and osx-64 (Apple silicon through Rosetta).
set -euo pipefail

mkdir -p "$1"
out=$(cd "$1" && pwd)
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)
work="${DROPLETUTILS_ORACLE_WORK:-${RUNNER_TEMP:-/tmp}/du}"
mkdir -p "$work"
cp "$here/pixi.toml" "$here/pixi.lock" "$work/"

datasets=(pbmc_1k_v3_3_0)
# environment:release[:extra] (extra arguments of oracle.R)
runs=(du100:1.0.3 du102:1.2.1:alphainf du104:1.4.2 du106:1.6.1 du108:1.8.0 du110:1.10.3 du112:1.12.1
      du114:1.14.2 du118:1.18.0 du120:1.20.0 du126:1.26.0 du130:1.30.0 du130:1.28:br128)

(cd "$work" && pixi install --all --locked)
for run in "${runs[@]}"; do
  env=${run%%:*}
  (cd "$work" && pixi run -e "$env" R CMD INSTALL "$here/stub/GenomeInfoDbData" > /dev/null)
done

for dataset in "${datasets[@]}"; do
  mtx="$work/mtx/$dataset"
  (cd "$repo" && uv run python oracles/starsolo/make_mtx.py "$dataset" "$mtx")
  for run in "${runs[@]}"; do
    IFS=: read -r env release extra <<< "$run"
    # shellcheck disable=SC2086  # extra is one optional word
    (cd "$work" && pixi run -e "$env" Rscript "$here/oracle.R" "$mtx" "$out/$dataset/$release" $extra)
  done
  gzip -f "$out/$dataset"/*/*.tsv "$out/$dataset"/*/defaultdrops.txt
done
