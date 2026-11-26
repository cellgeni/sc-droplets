#!/usr/bin/env bash
# Build a program per salmon release from its own alevin source (harness.py) and run its initial
# whitelist on the per-barcode UMI totals of public raw matrices.
#
#   bash oracles/alevin/run.sh OUT_DIR
#
# Needs pixi (for the C++ compiler), uv and git. Checkouts and builds go to $ALEVIN_ORACLE_WORK
# (default $RUNNER_TEMP/al or /tmp/al). Writes OUT_DIR/<dataset>/<release>/{cases.tsv,versions.txt}
# (see oracle.py). The committed outputs are linux-64's, as CI's.
set -euo pipefail

mkdir -p "$1"
out=$(cd "$1" && pwd)
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)
work="${ALEVIN_ORACLE_WORK:-${RUNNER_TEMP:-/tmp}/al}"
mkdir -p "$work/src" "$work/build" "$work/counts"
cp "$here/pixi.toml" "$here/pixi.lock" "$work/"

# dataset:expected cells
datasets=(pbmc_1k_v3_3_0:1000 pbmc4k_2_1:5000)
releases=(0.11.3 0.13.1 0.14.0 1.1.0 1.10.3)

(cd "$work" && pixi install --locked)
# conda-forge's compiler: gcc with libstdc++ on Linux (as salmon's Linux builds), clang on macOS
bin="$work/.pixi/envs/default/bin"
cxx=$(find "$bin" -maxdepth 1 \( -name '*-conda-linux-gnu-c++' -o -name '*-apple-darwin*-clang++' \) | head -1)

for release in "${releases[@]}"; do
  src="$work/src/v$release"
  if [ ! -d "$src" ]; then
    git init --quiet "$src"
    git -C "$src" remote add origin https://github.com/COMBINE-lab/salmon.git
    git -C "$src" -c advice.detachedHead=false fetch --quiet --depth 1 origin "refs/tags/v$release"
    git -C "$src" -c advice.detachedHead=false checkout --quiet FETCH_HEAD -- src/Alevin.cpp include/AlevinOpts.hpp
  fi
  (cd "$repo" && uv run python "$here/harness.py" "$src" "$work/build/alevin-$release.cpp")
  (cd "$work" && pixi run "$cxx" -std=c++17 -O2 -o "$work/build/alevin-$release" "$work/build/alevin-$release.cpp")
done

for entry in "${datasets[@]}"; do
  IFS=: read -r dataset expected <<< "$entry"
  counts="$work/counts/$dataset.tsv"
  (cd "$repo" && uv run python oracles/make_counts.py "$dataset" "$counts")
  for release in "${releases[@]}"; do
    dir="$out/$dataset/$release"
    (cd "$repo" && uv run python "$here/oracle.py" "$work/build/alevin-$release" "$release" "$counts" "$dir" "$expected")
    {
      echo "salmon $release"
      echo "source $(git -C "$work/src/v$release" rev-parse FETCH_HEAD)"
      (cd "$work" && pixi run "$cxx" --version | head -1)
      echo "platform $(uname -s) $(uname -m)"
    } > "$dir/versions.txt"
  done
done
