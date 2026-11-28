#!/usr/bin/env bash
# Run zUMIs' own cellBC (barcodeIDFUN.R at tags zUMIs2.5.6b and 2.9.7) with inflection 1.3.5 and
# 1.3.7 (CRAN's sources, from the cran/inflection mirror) in R from conda-forge, on the per-barcode
# UMI totals of public raw matrices.
#
#   bash oracles/zumis/run.sh OUT_DIR
#
# Needs pixi, uv and git. Work goes to $ZUMIS_ORACLE_WORK (default $RUNNER_TEMP/zu or /tmp/zu).
# Writes OUT_DIR/<dataset>/<zUMIs>-inflection<version>/{cells.tsv,cells_*.txt.gz,versions.txt} (see
# oracle.R). The committed outputs are linux-64's, as CI's.
set -euo pipefail

mkdir -p "$1"
out=$(cd "$1" && pwd)
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)
work="${ZUMIS_ORACLE_WORK:-${RUNNER_TEMP:-/tmp}/zu}"
mkdir -p "$work/counts"
cp "$here/pixi.toml" "$here/pixi.lock" "$work/"

# zUMIs release:tag
zumis=(2.5.6b:zUMIs2.5.6b 2.9.7:2.9.7)
inflection=(1.3.5 1.3.7)

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
for entry in "${zumis[@]}"; do
  IFS=: read -r release tag <<< "$entry"
  fetch https://github.com/sdparekh/zUMIs.git "refs/tags/$tag" "$work/zumis-$release" barcodeIDFUN.R
done
for version in "${inflection[@]}"; do
  fetch https://github.com/cran/inflection.git "refs/tags/$version" "$work/inflection-$version" R
done

for dataset in pbmc_1k_v3_3_0 pbmc4k_2_1; do
  counts="$work/counts/$dataset.tsv"
  (cd "$repo" && uv run python oracles/make_counts.py "$dataset" "$counts")
  for entry in "${zumis[@]}"; do
    IFS=: read -r release tag <<< "$entry"
    for version in "${inflection[@]}"; do
      dir="$out/$dataset/$release-inflection$version"
      (cd "$work" && pixi run Rscript "$here/oracle.R" "$work/zumis-$release/barcodeIDFUN.R" "$work/inflection-$version" "$counts" "$dir")
      echo "zUMIs $tag $(git -C "$work/zumis-$release" rev-parse FETCH_HEAD)" >> "$dir/versions.txt"
      echo "inflection $version $(git -C "$work/inflection-$version" rev-parse FETCH_HEAD)" >> "$dir/versions.txt"
    done
  done
done
