---
title: dropEst
parent: Methods
nav_order: 13
has_children: true
---

# dropEst
{: .no_toc }

dropEst (Petukhov et al. 2018) processes droplet reads into corrected counts. Its R package,
dropestr, estimates how many cells a library holds from the log-log rank curve of UMI counts. The
estimate gives three bands: barcodes surely cells (High), surely not (Low), and between
(Unknown). dropestr's quality classifier is then trained on the first two.

1. TOC
{:toc}

## Which release does what

| dropestr | What changed |
|---|---|
| 0.6.0 (Jul 2017) | `EstimateCellsNumber` added; unchanged since (0.7.9 in dropEst 0.8.5 and 0.8.6) |

## Algorithms

| Page | Id | Rule |
|---|---|---|
| [cell-number bounds](bounds.md) | `dropest-bounds` | lagged derivatives of log count against log rank: the expected number of cells and its bounds |

## Not implemented

dropestr's classifier of cell quality, trained on the High and Low bands, is planned for 0.2.
