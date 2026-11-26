---
title: cellsweep
parent: Methods
nav_order: 7
has_children: true
---

# cellsweep
{: .no_toc }

cellsweep (Caskey et al. 2026) removes ambient and "bulk" contamination from droplet counts with a
multinomial mixture fitted by EM. It is not a cell caller, but it needs to know which barcodes
are empty droplets, since they anchor its ambient profile, and when it is not told it decides
from the barcodes' totals. That step is a cell-vs-empty call like the others, and `kb sweep`
(kallisto | bustools) runs it.

1. TOC
{:toc}

## Which release does what

| cellsweep | Default rule | What changed |
|---|---|---|
| 0.1.0 (Mar 2026) – 1.0.1 | `threshold` | the cut is the most negative second derivative of the smoothed $$\log_{10}$$ rank curve |
| 1.0.2 (Nov 2026) | `kneed` | the cut is the knee of the rank curve by kneed's Kneedle (`online=True`); `threshold` becomes `manual` |

`mx_filter` (a two-component Gaussian mixture of $$\log(1 + \text{count})$$) is offered by every
release. A given UMI cutoff or number of expected cells always wins (`manual`).

The function was named `denoise_count_matrix` before 1.0.0 and `denoise` since; the empty-droplet
step (`infer_empty_droplets` in `cellsweep/utils/data_utils.py`) is the same in both.

## Algorithms

| Page | Id | Rule |
|---|---|---|
| [empty droplets](empty.md) | `cellsweep-empty` | `infer_empty_droplets`: a cutoff on the totals by the release's rule; below it, empty |
