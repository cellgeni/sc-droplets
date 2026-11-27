#!/usr/bin/env bash
# Run the reference implementations of the generic rules on the per-barcode UMI totals of public
# raw matrices: scuttle's isOutlier (1.0.4, 1.22.0) for `mad`, scikit-image's threshold_otsu and
# threshold_multiotsu (0.20.0, 0.26.0) for `otsu` and `multiotsu`.
#
#   bash oracles/generic/run.sh OUT_DIR
#
# Needs pixi and uv. The environments are built in $GENERIC_ORACLE_WORK (default
# $RUNNER_TEMP/gen or /tmp/gen). Writes OUT_DIR/<dataset>/scuttle-<release>/{mad.tsv,versions.txt}
# (see mad.R) and OUT_DIR/<dataset>/skimage-<release>/{otsu.tsv,multiotsu.tsv,versions.txt} (see
# otsu.py). The committed outputs are linux-64's, as CI's: scikit-image's arm64 builds fuse a
# multiply and an add in multi-Otsu's moments, which can break a tie between two splits the
# other way.
set -euo pipefail

mkdir -p "$1"
out=$(cd "$1" && pwd)
here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)
work="${GENERIC_ORACLE_WORK:-${RUNNER_TEMP:-/tmp}/gen}"
mkdir -p "$work/counts"
cp "$here/pixi.toml" "$here/pixi.lock" "$work/"

(cd "$work" && pixi install --all --locked)
# scuttle imports GenomeInfoDb, which imports a data package conda fills by a download script; the
# DropletUtils oracle's empty stand-in satisfies the import without running it.
for env in scuttle104 scuttle122; do
  (cd "$work" && pixi run -e "$env" R CMD INSTALL "$here/../dropletutils/stub/GenomeInfoDbData" > /dev/null)
done

for dataset in pbmc_1k_v3_3_0 pbmc4k_2_1; do
  counts="$work/counts/$dataset.tsv"
  (cd "$repo" && uv run python oracles/make_counts.py "$dataset" "$counts")
  for run in scuttle104:1.0.4 scuttle122:1.22.0; do
    IFS=: read -r env release <<< "$run"
    (cd "$work" && pixi run -e "$env" Rscript "$here/mad.R" "$counts" "$out/$dataset/scuttle-$release")
  done
  for run in sk020:0.20.0 sk026:0.26.0; do
    IFS=: read -r env release <<< "$run"
    (cd "$work" && pixi run -e "$env" python "$here/otsu.py" "$counts" "$out/$dataset/skimage-$release")
  done
done
