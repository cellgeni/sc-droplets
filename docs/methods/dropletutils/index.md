---
title: DropletUtils
parent: Methods
nav_order: 4
has_children: true
---

# DropletUtils
{: .no_toc }

Aaron Lun's Bioconductor package for droplet data, where EmptyDrops was first published. It
has four ways to call cells: a knee point on the barcode-rank curve, the EmptyDrops test,
STARsolo's version of Cell Ranger's test, and Cell Ranger 2's rule. The knee changed four times
and `emptyDrops` uses it to choose which barcodes it always keeps, so the same call with the same
arguments gives different cells in different releases.

1. TOC
{:toc}

## Which release does what

| DropletUtils | Bioconductor | R | What changed |
|---|---|---|---|
| 1.0 | 3.7 (May 2018) | 3.5 | `barcodeRanks`, `emptyDrops` (multinomial), `defaultDrops` |
| 1.2 | 3.8 | 3.5 | `emptyDrops` fits a Dirichlet-multinomial by default (`alpha` estimated) |
| 1.4 | 3.9 | 3.6 | Monte Carlo moves to one PCG32 stream per iteration, seeded from R |
| 1.6 | 3.10 | 3.6 | Boost's `discrete_distribution` (alias table) for the draws; non-integer matrices rounded |
| 1.8 | 3.11 | 4.0 | `round=`; a pseudo-count when Good-Turing leaves a feature at 0 |
| 1.10 | 3.12 | 4.0 | knee fit skips the first 50 barcodes (`exclude.from`), indexing bug fixed; `by.rank=` |
| 1.12 | 3.13 | 4.1 | totals rounded to integers; `test.ambient` excludes the ambient droplets from the FDR |
| 1.14 | 3.14 | 4.1 | `emptyDropsCellRanger` (with a bug in the candidates' minimum until 1.18.0) |
| 1.18.1 | 3.16 | 4.2 | `emptyDropsCellRanger` bug fixed |
| 1.20 | 3.17 | 4.3 | `known.empty=`; `by.rank` counts exactly that many barcodes, ties by order |
| 1.28 | 3.21 | 4.5 | knee by distance to a chord; `emptyDrops` back to the multinomial (`alpha = Inf`) |
| 1.30 – 1.34 | 3.22 – 3.24 | 4.5 – 4.6 | knee by tracing the curve with a sliding window |

Releases 1.16, 1.22 – 1.26 and 1.32 – 1.34 change nothing these methods compute.

## Methods

| Page | Rule |
|---|---|
| [barcodeRanks](barcoderanks.md) | the knee (or inflection) of the log-log barcode-rank curve |
| [emptyDrops](emptydrops.md) | every droplet above `lower` tested against the ambient profile; the knee and above always kept |
| [emptyDropsCellRanger](emptydropscellranger.md) | STARsolo's EmptyDrops_CR, rewritten in R |
| [defaultDrops](defaultdrops.md) | Cell Ranger 2: 10 % of the 99th percentile of the top 3,000 totals |

## How DropletUtils sees a matrix

- Every column of the matrix is a barcode, including barcodes with no counts; ranks, `by.rank`
  and `emptyDropsCellRanger`'s ambient ranks count them.
- Totals are column sums; from 1.12 on, rounded to the nearest integer (half to even), as
  `as.integer(round(colSums(m)))`.
- Features with no count in any barcode are dropped before the ambient profile is estimated.

sc-droplets reads a raw matrix with every barcode it lists, as `DropletUtils::read10xCounts` does.

## Random numbers

DropletUtils draws its Monte Carlo from R's random number generator, so its p-values depend on
the seed given to `set.seed` before the call. sc-droplets reproduces R's Mersenne-Twister and
every stream built on it (pages: [emptyDrops](emptydrops.md#random-numbers)), so with the same
seed the p-values are R's own. The default seed is 100, the one in DropletUtils' vignette; pass
`--seed` to match another run.
