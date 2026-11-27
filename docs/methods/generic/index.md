---
title: Generic rules
parent: Methods
nav_order: 1
has_children: true
---

# Generic rules

Rules that belong to no single tool. Many tools offer one of them as an option (Cell Ranger's
`--force-cells` is a top N; Drop-seq's `MIN_TRANSCRIPTS_PER_CELL` is a fixed threshold), each
with its own small differences, which their own pages describe. The generic versions here are
the plain rules.

| Page | Rule |
|---|---|
| [Top N barcodes](top-n.md) | the N barcodes with the most counts |
| [Fixed threshold](threshold.md) | every barcode with at least a given number of UMIs, reads or genes |
| [MAD outliers](mad.md) | above a floor, every barcode not more than 3 MADs below the median log total (scuttle's `isOutlier`) |
| [Otsu threshold](otsu.md) | above a floor, the upper class of Otsu's split of the log totals (scikit-image) |
| [Multi-Otsu thresholds](multiotsu.md) | the top of three classes of the log totals (scikit-image) |
| [Kernel-density valley](kde-valley.md) | above a floor, the barcodes above the lowest point between the two highest peaks of the log totals' density |

The rules that work on the distribution of totals need the floor (`--min-counts`): a raw matrix
holds hundreds of thousands of barcodes with one or two counts, which would otherwise make the
distribution. The default floor, 100, is DropletUtils' `lower`.

| On 10x's raw matrices | pbmc_1k_v3 | pbmc4k |
|---|---|---|
| Cell Ranger (3.0, 2.1) | 1,222 | 4,340 |
| `mad` | 1,160 | 4,541 |
| `otsu` | 1,143 | 4,489 |
| `multiotsu` | 1,317 | 4,845 |
| `kde-valley` | 1,164 | 4,768 |

Coming in 0.1.0: Gaussian mixtures on log counts (one and two dimensions).
