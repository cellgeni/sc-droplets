#!/usr/bin/env bash
# Run dropestr's own EstimateCellsNumber (dropestr/R/cell_number.R at dropEst's tag v0.8.6,
# dropestr 0.7.9) in R from conda-forge on the per-barcode UMI totals of public raw matrices.
#
#   bash oracles/dropest/run.sh OUT_DIR
#
# Needs pixi, uv and git. Work goes to $DROPEST_ORACLE_WORK (default $RUNNER_TEMP/de or /tmp/de).
# Writes OUT_DIR/<dataset>/0.7.9/{bounds.tsv,versions.txt} (see oracle.R). The committed outputs are
# linux-64's, as CI's.
set -euo pipefail

mkdir -p "$1"
out=$(cd "$1" && pwd)
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)
work="${DROPEST_ORACLE_WORK:-${RUNNER_TEMP:-/tmp}/de}"
mkdir -p "$work/counts"
cp "$here/pixi.toml" "$here/pixi.lock" "$work/"

# Check out PATHs of REPO at REF (a tag or commit) into DIR, once.
fetch() {
  local repo=$1 ref=$2 dir=$3
  shift 3
  if [ ! -d "$dir/.git" ]; then
    git init --quiet "$dir"
    git -C "$dir" remote add origin "$repo"
    git -C "$dir" fetch --quiet --depth 1 origin "$ref"
    git -C "$dir" -c advice.detachedHead=false checkout --quiet FETCH_HEAD -- "$@"
  fi
}

(cd "$work" && pixi install --locked)
fetch https://github.com/kharchenkolab/dropEst.git refs/tags/v0.8.6 "$work/dropEst" dropestr/R/cell_number.R

for dataset in pbmc_1k_v3_3_0 pbmc4k_2_1; do
  counts="$work/counts/$dataset.tsv"
  (cd "$repo" && uv run python oracles/make_counts.py "$dataset" "$counts")
  (cd "$work" && pixi run Rscript "$here/oracle.R" "$work/dropEst/dropestr/R/cell_number.R" "$counts" "$out/$dataset/0.7.9")
  echo "dropEst v0.8.6 $(git -C "$work/dropEst" rev-parse FETCH_HEAD)" >> "$out/$dataset/0.7.9/versions.txt"
done
