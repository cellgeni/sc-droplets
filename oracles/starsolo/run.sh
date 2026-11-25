#!/usr/bin/env bash
# Run STAR's own cell filtering on public raw matrices and save the cell barcodes
# it keeps, for tests/test_validation_starsolo.py.
#
#     bash oracles/starsolo/run.sh OUT_DIR
#
# Uses STAR's `--runMode soloCellFiltering` (2.7.8a+; 2.7.9a crashes in this mode),
# with STAR's released static Linux binaries, so it runs on Linux x86-64 only
# (GitHub Actions: .github/workflows/oracles.yml). Writes
# OUT_DIR/<dataset>/<method>@<version>.txt.gz, one barcode per line.
set -euo pipefail

out=$(realpath -m "$1")
work=$(mktemp -d "${RUNNER_TEMP:-/tmp}/starsolo.XXXXXX")
here=$(cd "$(dirname "$0")" && pwd)

datasets=(pbmc_1k_v2_3_0 pbmc_1k_v3_3_0)
versions=(2.7.10a 2.7.11b)

star_binary() {
  local version=$1 dir="$work/star-$1"
  mkdir -p "$dir"
  case "$version" in
    2.7.8a|2.7.9a|2.7.10a)
      curl -fsSL -o "$dir/STAR" "https://github.com/alexdobin/STAR/raw/$version/bin/Linux_x86_64_static/STAR" ;;
    *)
      curl -fsSL -o "$dir/star.zip" "https://github.com/alexdobin/STAR/releases/download/$version/STAR_$version.zip"
      unzip -q -o "$dir/star.zip" -d "$dir"
      cp "$dir/STAR_$version/Linux_x86_64_static/STAR" "$dir/STAR" ;;
  esac
  chmod +x "$dir/STAR"
  echo "$dir/STAR"
}

# method id -> --soloCellFilter arguments (STAR's defaults)
declare -A filters=(
  [starsolo-cellranger2.2]="CellRanger2.2 3000 0.99 10"
  [starsolo-topcells]="TopCells 3000"
  [starsolo-emptydrops-cr]="EmptyDrops_CR 3000 0.99 10 45000 90000 500 0.01 20000 0.01 10000"
)

for dataset in "${datasets[@]}"; do
  uv run python "$here/make_mtx.py" "$dataset" "$work/$dataset/raw"
done

for version in "${versions[@]}"; do
  star=$(star_binary "$version")
  "$star" --version
  for dataset in "${datasets[@]}"; do
    mkdir -p "$out/$dataset"
    for method in "${!filters[@]}"; do
      prefix="$work/$dataset/$method-$version/"
      mkdir -p "$prefix"
      # shellcheck disable=SC2086  # the filter is several words on purpose
      "$star" --runMode soloCellFiltering "$work/$dataset/raw/" "$prefix" --soloCellFilter ${filters[$method]} \
        --outFileNamePrefix "$prefix" > "$prefix/stdout.txt"
      found=$(find "$prefix" -name barcodes.tsv | head -n 1)
      if [ -z "$found" ]; then
        echo "::error::STAR $version wrote no barcodes.tsv for $dataset $method"; ls -R "$prefix"; cat "$prefix/stdout.txt"; exit 1
      fi
      gzip -n -c "$found" > "$out/$dataset/$method@$version.txt.gz"
      echo "$dataset $method@$version: $(wc -l < "$found") cells"
    done
  done
done
