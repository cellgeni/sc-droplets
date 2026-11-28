---
title: UMI-tools
parent: Methods
nav_order: 5
has_children: true
---

# UMI-tools
{: .no_toc }

`umi_tools whitelist` picks the cell barcodes of a droplet run from the reads, before any
alignment: it counts reads (or distinct UMIs) per barcode and keeps the barcodes above a
threshold. Drop-seq and inDrop workflows, and pipelines built on UMI-tools' `extract` and
`count`, take their cells from it. It has two rules for the threshold.

1. TOC
{:toc}

## Which release does what

| UMI-tools | Default rule | What changed |
|---|---|---|
| 0.5.0 (Aug 2017) | density | `whitelist` added: the first valley of a kernel density of log counts |
| 1.0.0 (Feb 2019) | distance | the farthest point from a chord on the cumulative curve becomes the default; `--knee-method=density` keeps the old rule |
| 1.1.0 – 1.1.6 | distance | error correction of barcodes rewritten; the thresholds are unchanged |

Neither rule changed after it was added: `whitelist_methods.py` (1.0+) and `umi_methods.py`
(0.5) compute the same threshold in every release.

## Algorithms

| Page | Id | Rule |
|---|---|---|
| [distance](distance.md) | `umitools-distance` | farthest point from the chord of the cumulative curve, repeated on the first three times as many barcodes |
| [density](density.md) | `umitools-density` | a local minimum of a Gaussian kernel density of $$\log_{10}$$ counts |

## Counting

UMI-tools counts the reads of each barcode in the first 100 million reads
(`--subset-reads`), or the distinct UMIs with `--method=umis`. sc-droplets takes the counts
from the matrix: UMIs by default, or reads given with `--reads`. A barcode with no counts is
not on UMI-tools' curve (it never appears in the reads), so both rules leave it out.

With `--method=umis`, UMI-tools counts the distinct UMI sequences of a barcode over all its
reads; a matrix's total counts a UMI once per gene, so it can be a little larger.

`whitelist --error-correct-threshold` maps barcodes one edit away from a cell barcode onto it;
that changes which reads count towards a cell, not which barcodes are cells. The opt-in
`--ed-above-threshold` drops cell barcodes within an edit of a more abundant cell barcode; it
needs the barcode sequences and is not offered here.
