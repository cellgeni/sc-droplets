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

Coming in 0.1.0: MAD outliers, Otsu and multi-Otsu thresholds, a kernel-density valley, Gaussian
mixtures on log counts (one and two dimensions) and a negative-binomial mixture.
