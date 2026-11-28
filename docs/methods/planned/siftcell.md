---
title: SiftCell
parent: Planned
grand_parent: Methods
nav_order: 6
method_ids: [siftcell]
status: planned
---

# SiftCell: a shuffled null, then boosted trees
{: .no_toc }

SiftCell's Shuffle and Boost steps. Shuffling the gene of every UMI across the whole matrix keeps
each barcode's total but turns its profile into a draw from the pooled profile: a set of
droplets that hold no cell. An XGBoost classifier, trained to tell those shuffled barcodes and
mitochondria-flagged barcodes from the barcodes least like the pool, labels the rest. Planned
for 0.2.

1. TOC
{:toc}

## At a glance

| Id | `siftcell` (planned) |
|---|---|
| Tool | SiftCell 1.0 (GitHub `jyxi7676/SiftCell` only; tag `v1.0.0`, `master` adaaec1), `SiftCellShuffle` then `SiftCellBoost` |
| Category | classifier |
| Counts | UMIs per gene per barcode |
| Needs expected cells | a guess: `expectedN` (1,000) sets how many barcodes train as cells |
| Deterministic | no: the default `seed = 0` sets no seed, so the shuffle and the random tie-breaks differ from run to run. Under the caller's `set.seed` and a fixed xgboost a run repeats |
| Original licence | GPL (≥ 2) |
| Extra inputs | none; the raw matrix as an uncompressed Matrix Market folder with 10x v2 file names |

## How it works

The input is `workingdir/DGE/` with `matrix.mtx`, `barcodes.tsv` and `genes.tsv` (gene names
from column 2, made unique).

**Shuffle** (`SiftCellShuffle`): the matrix is expanded into one (gene, barcode) pair per UMI,
$$T$$ pairs. A Fisher–Yates-style pass swaps the gene label at position $$i$$ with the one at
$$i + \operatorname{round}(U(0, T - i - 1))$$, $$U$$ from R's `runif` through Rcpp. The rounding
gives the two end positions half the weight of the others, so the permutation is not exactly
uniform. Barcode and gene totals are kept; the shuffled matrix goes to `workingdir/shuffleDGE/`.

**Boost** (`SiftCellBoost`):

1. **Barcodes.** Those with more than `threshold` (100) UMIs, in both matrices: the same barcodes,
   since the totals are kept.
2. **Flagged.** For each gene matching `^MT-|^mt-`, its fraction in each barcode; its cut is the
   median plus 3 standard deviations over barcodes. A barcode above the cut of any of them is
   flagged.
3. **Positives.** $$p_0$$, the pooled gene proportions over every barcode with a count, and
   weights $$w_g = \sqrt{p_{0g}/(1 - p_{0g})}$$. Each barcode's squat score: over its detected
   genes, a $$z_g$$ for how extreme $$x_g$$ is under Binomial(total, $$p_{0g}$$) conditioned on
   $$x_g > 0$$ (expectation-based, variance-adjusted, capped at $$\pm 10$$), combined as
   $$\sum_g w_g z_g / \sqrt{\sum_g w_g^2}$$. The positives are the barcodes within the top
   `expectedN` (1,000) both by UMIs and by score, ranked with `ties.method = "random"`.
4. **Variable genes.** Genes with more than `geneUMIthreshold` (50) UMIs; a squat test of each
   across barcodes (sizes the barcode totals, weights $$\sqrt{\text{size}}$$); the 1,000 with the
   largest $$-\log_{10} p$$. With `quantileReg`, those above a logistic quantile-regression curve
   at $$\tau = 0.8$$ instead.
5. **PCs.** $$\log(10^4 x / \text{total} + 1)$$ on those genes, the original and the shuffled
   barcodes stacked ($$2N$$ rows); `prcomp`, centred and scaled; 100 PCs (`N.PCs`).
6. **Training.** The positives are `Nonsoup`; the flagged barcodes and every shuffled barcode are
   `Soup`. Every other original barcode is predicted. (`labels`, a file, adds positives and
   negatives by hand.)
7. **Classifier.** In Python, `xgboost.XGBClassifier(objective = "binary:logistic")`, every other
   setting at the installed version's default, fitted on the PCs; `predict` gives the class, a
   probability above 0.5.
8. **Call.** The cells are the positives and the barcodes predicted `Nonsoup`, written to
   `CellContaining.txt`. Flagged barcodes and those with at most 100 UMIs are never cells.

Defects that change a run: `seed = 0` sets no seed (only `seed = NULL` seeds, from the clock); the
Python script runs twice, its failure silenced by `try(..., silent = TRUE)`; the string labels
need xgboost before 1.6 (1.5.2 encodes them, 1.6.0 and later stop with "Invalid classes
inferred"), so with a current xgboost no file is written; `s.strip(dataName + "Org_")` strips
characters rather than a prefix, harmless for the default "Sample" but corrupting barcodes when
`dataName` holds A, C, G or T; `workingdir` must end in `/`.

SiftCell-Mix, the package's third step, is not a caller: it needs a cell type for every barcode
above the threshold and fits each barcode's weights over the cell types' profiles (each blended
with 1 % of the ambient profile).

## Versions

| SiftCell | Change to the call |
|---|---|
| 1.0.0 (`v1.0.0`, Jan 2023) | Shuffle and Boost as above |
| `master` (Mar – Jul 2023, last adaaec1) | Shuffle reads the matrix with `readMM` instead of parsing lines; Boost only frees memory; the XGBoost script unchanged in substance; SiftCell-Mix rewritten and SiftCell-Bayes removed. The same call |

One variant, but the installed xgboost decides whether there is a call at all:

| xgboost | Effect |
|---|---|
| up to 1.5.2 | string labels encoded; cells written |
| 1.6.0 and later | the script fails, silently: no output |

## Parameters

| Original | Default | Meaning |
|---|---|---|
| `threshold` | 100 | barcodes with more UMIs are classified; the rest are never cells |
| `expectedN` | 1,000 | positives come from the top this many by UMIs and by score |
| `mitoGenes` | `TRUE` | flag barcodes by their mitochondrial genes |
| `geneUMIthreshold` | 50 | genes with more UMIs are tested for variability |
| `quantileReg` | `FALSE` | variable genes by quantile regression rather than the top 1,000 |
| `N.PCs` | 100 | PCs the classifier sees |
| `labels` | `FALSE` | a file of manual positives and negatives |
| `seed` | 0 | 0 sets no seed; `NULL` seeds from the clock |

## Status

Planned for 0.2. Implementing it needs:

- R's `runif` stream for the shuffle, with its rounding (the package has R's Mersenne-Twister and
  `unif_rand`);
- squat's binomial overdispersion test (squat v1.11, Apache-2.0), ported exactly, for the scores
  and the variable genes;
- R's `prcomp`, centred and scaled;
- xgboost itself: the port can call the same xgboost on the same PCs, from Python.

The oracle is SiftCell built from GitHub (Rcpp, RcppArmadillo; squat from GitHub; quantreg,
Rsolnp, reticulate), Python with pandas, scikit-learn and **xgboost 1.5.x** (conda-forge's
`py-xgboost` has 1.5.0 and 1.5.1, PyPI 1.5.2), and `set.seed` called by the harness before each
step. Runs then repeat for fixed versions, xgboost's default tree method and threads aside. The
blocker is xgboost: 1.6 and later write nothing. The 1,000 positives and the shuffled null are
the parts most worth validating exactly.

## References

- Xi, Park, Lee & Kang 2023, SiftCell: a robust framework to detect and isolate cell-containing
  droplets from single-cell RNA sequence reads, *Cell Systems* 14(7):620–628.e3,
  [doi:10.1016/j.cels.2023.06.002](https://doi.org/10.1016/j.cels.2023.06.002).
- Source: [jyxi7676/SiftCell at adaaec1](https://github.com/jyxi7676/SiftCell/tree/adaaec1)
  (`R/SiftCellShuffle.R`, `R/SiftCellBoost.R`, `src/dgeShuffle.cpp`, `inst/runXGboost.py`), and
  [at v1.0.0](https://github.com/jyxi7676/SiftCell/tree/v1.0.0); squat:
  [hyunminkang/squat](https://github.com/hyunminkang/squat) at 3da0ca4 (`R/squat.R`,
  `src/squat.cpp`).
- Research notes: `docs/_research/planned-rna.md` §6.
