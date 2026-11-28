---
title: CellBender
parent: Methods
nav_order: 8
has_children: true
---

# CellBender
{: .no_toc }

CellBender's `remove-background` (Fleming et al. 2023) fits a deep generative model of droplet
counts, ambient RNA and barcode swapping. It decides which droplets hold cells as part of that
fit, but it starts from rules that are not part of the model. These rules give prior guesses of
the counts in cells and in empty droplets. They also set how many droplets are probably cells,
how many droplets the model analyses, and which empty droplets it learns the ambient profile
from. Those starting rules are a cell-vs-empty call on the rank curve like the others, and they
changed twice.

1. TOC
{:toc}

## Which release does what

| CellBender | Priors from | What changed |
|---|---|---|
| 0.1.0 | mode and median of $$\log(1 + \text{count})$$ | first release: the empty prior is the mode of rounded log counts, the cell prior a median above five times it; 25,000 droplets analysed |
| 0.2.0 – 0.2.2 | the same | the low-count cut drops from 0.8 to 0.5 of the empty prior; with `--expected-cells`, the cell prior is the median of the top barcodes |
| 0.3.0 (Aug 2023), 0.3.2, 0.4.0 (Aug 2026) | Otsu's threshold and a kernel-density peak | new heuristics (`priors.py`) set the expected cells and the number of droplets analysed, and drop features with almost no ambient counts |

0.3.0_rc, the code of the Nature Methods paper, fitted a Gaussian mixture with pyro instead; it
was never released on PyPI and is not implemented (it is stochastic and needs torch). 0.3.1 was
withdrawn before it was tagged. 0.3.2 and 0.4.0 change no line of the prior code of 0.3.0.

## Algorithms

| Page | Id | Rule |
|---|---|---|
| [priors](prior.md) | `cellbender-prior` | the release's prior heuristics; probable cells, the other analysed droplets, the surely-empty droplets |

## Not implemented

The model itself, and so CellBender's own cells (posterior cell probability above 0.5), is a
variational autoencoder trained with pyro, documented on [CellBender model](model.md). Its result
depends on the training, the hardware and the random seed, not only on the matrix.
