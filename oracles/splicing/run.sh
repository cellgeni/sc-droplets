#!/usr/bin/env bash
# Spliced and unspliced counts for 10x's pbmc_1k_v3: kb count --workflow nac (kb-python 0.28.2 with
# kallisto 0.50.1 and bustools 0.43.2, from bioconda) with kallisto's prebuilt human nac index (made
# in November 2023; kallisto 0.52 cannot read it), on the sample's FASTQs, summed per barcode. The
# splicing tests and the comparison page read them.
#
#   bash oracles/splicing/run.sh OUT_DIR [READS]
#
# Needs pixi and curl, about 15 GB of disk and 16 GB of memory. READS, when given, keeps only the
# first READS reads of lane 1 (a quick check of the commands). Work goes to $SPLICING_WORK
# (default $RUNNER_TEMP/splicing or /tmp/splicing). Writes
# OUT_DIR/pbmc_1k_v3_3_0/kb-0.28.2/{counts.tsv.gz,versions.txt}; the committed output is linux-64's.
set -euo pipefail

mkdir -p "$1"
out=$(cd "$1" && pwd)
reads=${2:-}
here=$(cd "$(dirname "$0")" && pwd)
work="${SPLICING_WORK:-${RUNNER_TEMP:-/tmp}/splicing}"
mkdir -p "$work"
cp "$here/pixi.toml" "$here/pixi.lock" "$work/"
fastqs=https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar
index=https://github.com/pachterlab/kallisto-transcriptome-indices/releases/download/v1/human_index_nac.tar.xz
target="$out/pbmc_1k_v3_3_0/kb-0.28.2"

(cd "$work" && pixi install --locked)

# The FASTQs, extracted as they download (5.5 GB): no copy of the tarball on disk.
mkdir -p "$work/fastqs"
if [ -z "$reads" ]; then
  [ -f "$work/fastqs/done" ] || { curl -sSfL "$fastqs" | tar -x -C "$work/fastqs" --strip-components 1; touch "$work/fastqs/done"; }
  r1=("$work"/fastqs/*_R1_001.fastq.gz)
  r2=("$work"/fastqs/*_R2_001.fastq.gz)
  pairs=()
  for i in "${!r1[@]}"; do pairs+=("${r1[$i]}" "${r2[$i]}"); done
else
  lines=$((reads * 4))
  for read in R1 R2; do
    member="pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001_${read}_001.fastq.gz"
    [ -s "$work/fastqs/head_$read.fastq.gz" ] || { curl -sSfL "$fastqs" | tar -xO -f - "$member" | gzip -dc | head -n "$lines" | gzip > "$work/fastqs/head_$read.fastq.gz" || true; }
  done
  pairs=("$work/fastqs/head_R1.fastq.gz" "$work/fastqs/head_R2.fastq.gz")
fi

# kallisto's prebuilt human nac index (2 GB).
if [ ! -f "$work/index/index.idx" ]; then
  mkdir -p "$work/index"
  curl -sSfL "$index" | tar -xJ -C "$work/index"
fi
idx=$(find "$work/index" -name '*.idx' | head -n 1)
t2g=$(find "$work/index" -name 't2g.txt' | head -n 1)
c1=$(find "$work/index" -name 'cdna.txt' | head -n 1)
c2=$(find "$work/index" -name 'nascent.txt' | head -n 1)
threads=$(getconf _NPROCESSORS_ONLN)
kallisto=$(cd "$work" && pixi run which kallisto)
bustools=$(cd "$work" && pixi run which bustools)

rm -rf "$work/count"
(cd "$work" && pixi run kb count --kallisto "$kallisto" --bustools "$bustools" -i "$idx" -g "$t2g" -c1 "$c1" -c2 "$c2" \
  -x 10xv3 --workflow nac -t "$threads" -o "$work/count" "${pairs[@]}")
(cd "$work" && pixi run python "$here/summarise.py" "$work/count/counts_unfiltered" "$target/counts.tsv.gz")
{
  (cd "$work" && pixi list --explicit 2>/dev/null | grep -E "^(kb-python|kallisto|bustools) " | awk '{print $1, $2}')
  echo "index $index"
  echo "reads ${reads:-all}"
  echo "platform $(uname -s)-$(uname -m)"
} > "$target/versions.txt"
cat "$target/versions.txt"
