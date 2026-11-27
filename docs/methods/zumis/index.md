---
title: zUMIs
parent: Methods
nav_order: 14
has_children: true
---

# zUMIs
{: .no_toc }

zUMIs (Parekh et al. 2018) processes reads with cell barcodes and UMIs from many protocols. When
it is not given the cells, it detects them from the reads of each barcode: first with a Gaussian
mixture, from late 2019 with the unit invariant knee (UIK) of the cumulative curve.

1. TOC
{:toc}

## Which release does what

| zUMIs | What changed |
|---|---|
| 0.0.1 (Sep 2017) | no automatic detection: the barcodes, or their number, had to be given |
| 0.0.2 – 0.0.6 (2017 – 2018) | a Gaussian mixture of $$\log_{10}$$ reads (`mclust::Mclust`) |
| 2.0 – 2.5.4 (2018 – 2019) | the same mixture, in `barcodeIDFUN.R` (`.FindBCcut_mclust`), after a floor of `nReadsperCell` (100) reads |
| 2.5.6b (Nov 2019) – 2.9.7 | `inflection::uik` on the cumulative reads (`.FindBCcut`) |

The mixture is not reproducible: from mclust 5.4 (Nov 2017), `Mclust` initialises on a random
subset of 2,000 values whenever there are more, and zUMIs sets no seed. sc-droplets refuses those
releases (a documented gap).

## Algorithms

| Page | Id | Rule |
|---|---|---|
| [cell barcodes](cell-barcodes.md) | `zumis-cell-barcodes` | the UIK of the cumulative reads; the cells are the barcodes with at least its count |
