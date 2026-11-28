#!/usr/bin/env bash
# Run DropletQC's own identify_empty_drops (R/identify_empty_drops.R at tag v1.0.0) with ks 1.14.2
# and 1.15.3 from conda-forge, on the (nuclear fraction, UMIs) tables of make_inputs.py.
#
#   bash oracles/dropletqc/run.sh OUT_DIR
#
# Needs pixi, uv and git. Work goes to $DROPLETQC_ORACLE_WORK (default $RUNNER_TEMP/dqc or
# /tmp/dqc). Writes OUT_DIR/inputs/*.tsv.gz and OUT_DIR/ks-<version>/{summary.tsv,*.status.gz,
# versions.txt}. The committed outputs are linux-64's, as CI's.
set -euo pipefail

mkdir -p "$1"
out=$(cd "$1" && pwd)
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)
work="${DROPLETQC_ORACLE_WORK:-${RUNNER_TEMP:-/tmp}/dqc}"
mkdir -p "$work"
cp "$here/pixi.toml" "$here/pixi.lock" "$work/"

tag=v1.0.0
src="$work/DropletQC"
if [ ! -d "$src/.git" ]; then
  git init --quiet "$src"
  git -C "$src" remote add origin https://github.com/powellgenomicslab/DropletQC.git
  git -C "$src" fetch --quiet --depth 1 origin "refs/tags/$tag"
  git -C "$src" -c advice.detachedHead=false checkout --quiet FETCH_HEAD -- R/identify_empty_drops.R
fi

(cd "$repo" && uv run python oracles/dropletqc/make_inputs.py "$out/inputs")
for env in ks1142 ks1153; do
  (cd "$work" && pixi install --locked -e "$env")
  version=$(cd "$work" && pixi run -e "$env" Rscript -e 'cat(as.character(packageVersion("ks")))')
  dir="$out/ks-$version"
  (cd "$work" && pixi run -e "$env" Rscript "$here/oracle.R" "$src/R/identify_empty_drops.R" "$out/inputs" "$dir")
  echo "DropletQC $tag $(git -C "$src" rev-parse FETCH_HEAD)" >> "$dir/versions.txt"
  echo "platform $(uname -s)-$(uname -m)" >> "$dir/versions.txt"
done
