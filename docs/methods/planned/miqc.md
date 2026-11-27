---
title: miQC
parent: Planned
grand_parent: Methods
nav_order: 7
method_ids: [miqc]
status: planned
---

# miQC: a mixture of regressions for compromised cells
{: .no_toc }

miQC's `mixtureModel` and `filterCells`. miQC does not call cells from empty droplets: it removes
compromised cells, those with a high mitochondrial fraction for the number of genes they detect,
from barcodes already called. A mixture of two linear regressions of the mitochondrial
percentage on genes detected separates intact cells from compromised ones. In sc-droplets it is
planned as a quality step after a caller. Planned for 0.2.

1. TOC
{:toc}

## At a glance

| Id | `miqc` (planned) |
|---|---|
| Tool | miQC 1.0.0 – 1.22.0 (Bioconductor 3.13 – 3.24) |
| Category | mixture |
| Counts | genes detected and the mitochondrial percentage of each cell (scater's `addPerCellQC`) |
| Needs expected cells | no; it needs the cells already called |
| Deterministic | no: flexmix starts from a random assignment and miQC sets no seed; under the caller's `set.seed` and a fixed flexmix a run repeats |
| Original licence | BSD-3-Clause |
| Extra inputs | the called cells; the mitochondrial genes |

## How it works

Releases 1.2.0 – 1.22.0. Each cell has `detected` (genes with a count) and
`subsets_mito_percent` (mitochondrial UMIs over all UMIs, times 100), read from `colData` under
exactly those names: `mixtureModel`'s arguments of the same names are not used in the formula.

1. **Fit.** `flexmix(subsets_mito_percent ~ detected, k = 2)`: a mixture of two Gaussian linear
   regressions, fitted by EM at flexmix's defaults. One start, from a random hard assignment of
   the cells (R's `sample`); at most 200 iterations; tolerance $$10^{-6}$$; `minprior = 0.05`, so
   a component whose prior falls below 0.05 is dropped, with a warning if only one is left.
   Other `model_type`s change the formula: `"spline"` uses `bs(detected)` (cubic B-spline, 3
   degrees of freedom), `"polynomial"` uses `poly(detected, 2)`, `"one_dimensional"` uses `~ 1`,
   a mixture of two Gaussians on the percentage alone.
2. **Compromised component.** The one with the larger intercept (the first parameter).
3. **Posterior.** `filterCells` keeps a cell when its posterior probability of the compromised
   component is at most `posterior_cutoff` (0.75). From 1.6.0, it stops with an error when that
   keeps every cell.
4. **Boundary** (`keep_all_below_boundary`). Every cell whose percentage lies below the intact
   component's fitted line at its `detected` is kept as well.
5. **Left cut-off** (`enforce_left_cutoff`). With $$m$$ the lowest percentage among the cells
   discarded so far and $$d$$ the `detected` of the first cell, in column order, whose
   percentage equals $$m$$: every cell with `detected` $$\le d$$ and percentage $$\ge m$$ is
   discarded.
6. **Result.** The kept cells. `get1DCutoff` (from 1.8.0) returns $$m$$ for a one-dimensional
   model, or stops when the kept and discarded percentages overlap.

The random start is the only source of randomness; miQC's own error message suggests "setting a
different random seed". Under `set.seed` the fit repeats for a given flexmix. Ties in step 5 go
to the first cell in column order.

On a whole raw matrix the two components would separate empty droplets from cells rather than
compromised cells from intact ones; nobody has tested that, so sc-droplets plans miQC after a
caller, on its cells.

## Versions

| miQC | Bioconductor | Change to the call |
|---|---|---|
| 1.0.0 | 3.13 (May 2021) | linear, spline and polynomial models; a cell is kept when its posterior of the compromised component is at most 0.75 |
| 1.2.0 | 3.14 (Oct 2021) | the one-dimensional model; `keep_all_below_boundary` and `enforce_left_cutoff`, both on by default, so the default result changes |
| 1.4.0 | 3.15 | none |
| 1.6.0 | 3.16 (Nov 2022) | an error when no cell fails the posterior cut (1.2 and 1.4 returned every cell) |
| 1.8.0 | 3.17 (Apr 2023) | `get1DCutoff` added; filtering unchanged |
| 1.10.0 – 1.22.0 | 3.18 – 3.24 | none |

Two variants for the filter: 1.0.0, the posterior alone; 1.2.0 – 1.22.0, with the boundary and
the left cut-off. A pipeline today runs 1.18 – 1.22.

## Parameters

| Original | Default | Meaning |
|---|---|---|
| `model_type` (`mixtureModel`) | `"linear"` | the regression: linear, spline, polynomial or one-dimensional |
| `model` (`filterCells`) | `NULL` | a fitted model; with none, the linear model is fitted, with a warning |
| `posterior_cutoff` | 0.75 | the largest posterior of the compromised component a kept cell may have |
| `keep_all_below_boundary` | `TRUE` (from 1.2.0) | keep every cell below the intact line |
| `enforce_left_cutoff` | `TRUE` (from 1.2.0) | discard cells with few genes and a high percentage |
| flexmix's `iter.max`, `minprior`, `tolerance` | 200, 0.05, $$10^{-6}$$ | not set by miQC: flexmix's own defaults |

## Status

Planned for 0.2. Implementing it needs:

- the per-cell measures: genes detected, and the mitochondrial share, which the package already
  computes for each barcode (`mito_fraction`, from genes named `MT-` or `mt-`);
- flexmix's EM for Gaussian regressions: a weighted least-squares M-step, the `minprior`
  pruning, and the stopping rule;
- R's `sample` for the starting labels: the package has R's Mersenne-Twister, but `sample`'s
  rejection sampling (R 3.6 and later) for any range is still to port;
- a place in the package for a step that filters a caller's cells rather than calling them. The
  generic [Gaussian mixture](../generic/gmm.md) and [MAD outliers](../generic/mad.md) are the
  nearest rules it has.

The oracle is miQC from Bioconductor (bioconda's `bioconductor-miqc` 1.18.0, or 1.22.0 from
Bioconductor 3.24) with flexmix, scater for the QC columns, and `set.seed` before
`mixtureModel`. The model is small, so an exact match is feasible. Failing that, the original's
fitted parameters can be a fixture, and steps 2 to 5 checked exactly.

## References

- Hippen, Falco, Weber et al. 2021, miQC: an adaptive probabilistic framework for quality
  control of single-cell RNA-sequencing data, *PLoS Computational Biology* 17(8):e1009290,
  [doi:10.1371/journal.pcbi.1009290](https://doi.org/10.1371/journal.pcbi.1009290).
- Source: [miQC on Bioconductor](https://bioconductor.org/packages/miQC/)
  ([git](https://git.bioconductor.org/packages/miQC), `R/mixtureModel.R`, `R/filterCells.R`,
  `R/get1DCutoff.R`; mirror [greenelab/miQC](https://github.com/greenelab/miQC)).
- Research notes: `docs/_research/planned-rna.md` §7.
