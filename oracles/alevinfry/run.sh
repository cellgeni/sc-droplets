#!/usr/bin/env bash
# Build a program per alevin-fry release from its own source (harness.py) and run its permit-list
# rules on the per-barcode UMI totals of public raw matrices.
#
#   bash oracles/alevinfry/run.sh OUT_DIR
#
# Needs pixi (for Rust), uv and git. Checkouts and builds go to $ALEVINFRY_ORACLE_WORK (default
# $RUNNER_TEMP/af or /tmp/af). Writes OUT_DIR/<dataset>/<release>/{cases.tsv,versions.txt} (see
# oracle.py). The committed outputs are linux-64's, as CI's.
set -euo pipefail

mkdir -p "$1"
out=$(cd "$1" && pwd)
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)
work="${ALEVINFRY_ORACLE_WORK:-${RUNNER_TEMP:-/tmp}/af}"
mkdir -p "$work/src" "$work/build" "$work/counts"
cp "$here/pixi.toml" "$here/pixi.lock" "$work/"

# dataset:expected cells
datasets=(pbmc_1k_v3_3_0:1000 pbmc4k_2_1:5000)
releases=(0.1.0 0.11.2 0.17.1 0.19.0)

(cd "$work" && pixi install --locked)

for release in "${releases[@]}"; do
  src="$work/src/v$release"
  if [ ! -d "$src" ]; then
    git -c advice.detachedHead=false clone --quiet --depth 1 --branch "v$release" https://github.com/COMBINE-lab/alevin-fry.git "$src"
  fi
  (cd "$repo" && uv run python "$here/harness.py" "$src" "$work/build/$release")
  (cd "$work" && pixi run cargo build --release --quiet --manifest-path "$work/build/$release/Cargo.toml")
done

for entry in "${datasets[@]}"; do
  IFS=: read -r dataset expected <<< "$entry"
  counts="$work/counts/$dataset.tsv"
  (cd "$repo" && uv run python oracles/make_counts.py "$dataset" "$counts" --all)
  for release in "${releases[@]}"; do
    dir="$out/$dataset/$release"
    (cd "$repo" && uv run python "$here/oracle.py" "$work/build/$release/target/release/harness" "$release" "$counts" "$dir" "$expected")
    {
      echo "alevin-fry $release"
      echo "source $(git -C "$work/src/v$release" rev-parse HEAD)"
      (cd "$work" && pixi run rustc --version)
      echo "platform $(uname -s) $(uname -m)"
    } > "$dir/versions.txt"
  done
done
