#!/usr/bin/env bash
# Run dropbead's own estimateCellNumber (R/Tools.R at commit 73fd121, version 0.3.1) in R from
# conda-forge on the per-barcode UMI totals of public raw matrices.
#
#   bash oracles/dropseq/run.sh OUT_DIR
#
# Needs pixi, uv and git. Work goes to $DROPSEQ_ORACLE_WORK (default $RUNNER_TEMP/ds or /tmp/ds).
# Writes OUT_DIR/<dataset>/0.3.1/{cells.tsv,versions.txt} (see oracle.R). The committed outputs are
# linux-64's, as CI's.
set -euo pipefail

mkdir -p "$1"
out=$(cd "$1" && pwd)
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)
work="${DROPSEQ_ORACLE_WORK:-${RUNNER_TEMP:-/tmp}/ds}"
mkdir -p "$work/counts"
cp "$here/pixi.toml" "$here/pixi.lock" "$work/"
commit=73fd12129c9705eee21aad109b52ac1e9c8d21f4

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
fetch https://github.com/rajewsky-lab/dropbead.git "$commit" "$work/dropbead" R/Tools.R

for dataset in pbmc_1k_v3_3_0 pbmc4k_2_1; do
  counts="$work/counts/$dataset.tsv"
  (cd "$repo" && uv run python oracles/make_counts.py "$dataset" "$counts")
  (cd "$work" && pixi run Rscript "$here/oracle.R" "$work/dropbead/R/Tools.R" "$counts" "$out/$dataset/0.3.1")
  echo "dropbead $commit" >> "$out/$dataset/0.3.1/versions.txt"
done
