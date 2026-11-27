---
title: Drop-seq
parent: Methods
nav_order: 12
has_children: true
---

# Drop-seq
{: .no_toc }

Drop-seq's computational cookbook (Macosko et al. 2015) picks the number of cells by eye, at the
inflection of the cumulative fraction of reads against barcodes ranked by reads. dropbead, the R
package of Alles et al. (2017) for exploring Drop-seq data, turned that into a rule:
`estimateCellNumber` takes the barcode farthest from the diagonal of that curve.

1. TOC
{:toc}

## Which release does what

| dropbead | What changed |
|---|---|
| 0.3.1 (Jan 2018) | `estimateCellNumber` added (its first commit read an undefined `df`; fixed a week later, still 0.3.1) |

dropbead has no tags; its `DESCRIPTION` says 0.3.1 from then on.

## Algorithms

| Page | Id | Rule |
|---|---|---|
| [dropbead](dropbead.md) | `dropseq-dropbead` | the barcode farthest from the diagonal of the cumulative fraction of reads, among the top 15,000 |

## Not implemented

The cookbook's own choice is made by eye; Drop-seq's `SelectCellsByNumTranscripts` keeps the
barcodes above a fixed count, which the generic [threshold](../generic/threshold.md) does.
