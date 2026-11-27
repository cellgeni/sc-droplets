---
title: EmptyDropsMultiome
parent: Planned
grand_parent: Methods
nav_order: 13
method_ids: [emptydropsmultiome]
status: planned
---

# EmptyDropsMultiome: emptyDrops on both modalities, then lines in the plane
{: .no_toc }

Megas, Lorenzi & Marioni's caller for single-cell multiome data. Each modality gets its own
bounds from a three-component normal mixture and its own [emptyDrops](../dropletutils/emptydrops.md)
test; the two p-values are averaged and adjusted together. Barcodes below a line parallel to the
authors' approximation of [Cell Ranger ARC](cellranger-arc.md)'s k-means boundary are set aside,
and, in one of its two call modes, barcodes above a second line are accepted outright. Planned
for 0.3.

1. TOC
{:toc}

## At a glance

| Id | `emptydropsmultiome` (planned) |
|---|---|
| Tool | EmptyDropsMultiome, tag `first_release` (2024-02-21, the only release), `emptydrops_multiome` |
| Category | statistical-test |
| Counts | RNA UMIs per barcode; ATAC counts per barcode summed over peaks (cut sites, from Cell Ranger ARC's raw matrix) |
| Needs expected cells | no |
| Deterministic | seeded: each emptyDrops runs after `set.seed(42)`; the RNA mixture's start weights are random and follow the caller's R seed |
| Original licence | GPL (≥ 3); MIT before May 2023 |
| Extra inputs | both libraries: a gene × barcode RNA matrix and a peak × barcode ATAC matrix of the same barcodes |

## How it works

With $$R$$ and $$A$$ a barcode's RNA and ATAC totals, step 1 fits, per modality, a normal mixture
from these starts:

| | RNA | ATAC |
|---|---|---|
| Start means; SDs | 5, 500, 1,200; 5, 50, 800 | 5, 100, 10,000; 5, 50, 800 |
| Start weights | random: three uniform draws, normalised, redrawn until each is at least 0.05 | 6/7, 0.8/7, 0.1/7, renormalised to sum 1 |
| Fit rejected when | $$\mu_1 > 100$$ or $$\mu_2 > 700$$ or $$\mu_1 > \mu_2 - 15$$ | $$\mu_1 > 60$$ or $$\mu_2 > 600$$ or $$\mu_1 > \mu_2 - 15$$ |
| Defaults (rejected fit, or error) | $$\mu$$ 25, 200, 1,000; $$\sigma$$ 1, 60, 300; $$\lambda$$ 1, 0.10, 0.05, divided by 1.15 | $$\mu$$ 1, 40, 800; $$\sigma$$ 1, 60, 600; $$\lambda$$ 1, 0.10, 0.05 |
| `lower` $$l$$ | $$\mu_2 + 1.5\,\sigma_2$$ | $$\mu_2 + 2\,\sigma_2$$ |

1. **Bounds per modality**, unless both `lower` and `barhop` are given for it. mixtools'
   `normalmixEM` (tolerance $$10^{-9}$$, at most 1,000 iterations) fits three normals to the raw,
   unlogged totals of every barcode. Ordered by mean, the components are barcode swapping
   ("barhop"), ambient and cells. `barhop` $$b$$ is where the weighted densities of the first two
   cross: $$x = (-\kappa_1 + \sqrt{\kappa_1^2 - 4\kappa_2\kappa_0}) / (2\kappa_2)$$ with
   $$\kappa_2 = 1/\sigma_1^2 - 1/\sigma_2^2$$, $$\kappa_1 = -2(\mu_1/\sigma_1^2 - \mu_2/\sigma_2^2)$$,
   $$\kappa_0 = \mu_1^2/\sigma_1^2 - \mu_2^2/\sigma_2^2 + 2\log(\lambda_2\sigma_1 / (\lambda_1\sigma_2))$$,
   or $$-\kappa_0/\kappa_1$$ if $$\kappa_2 = 0$$. With a discriminant $$\le 0$$ no bound is returned
   and the later steps fail. A variance below $$10^{-8}$$ makes mixtools restart from random
   values; after 20 restarts it fails and the defaults apply.
2. **emptyDrops per modality.** Barcodes with total $$< b$$ are dropped (NA). After
   `set.seed(seed)` (42), DropletUtils' [emptyDrops](../dropletutils/emptydrops.md) runs with
   `lower` $$= l$$, `test.ambient = TRUE` and 10,000 iterations for RNA, 25,000 for ATAC, every
   other argument at its default. The ambient profile comes from the barcodes with
   $$b \le$$ total $$\le l$$; every remaining barcode with a positive total gets a p-value, ambient
   ones included. Only the p-values are used, so emptyDrops' knee and FDR do not matter.
3. **Combine.** The two tables are joined by barcode (a barcode missing from one gets NA).
   $$p = (p_{RNA} + p_{ATAC}) / 2$$, the arithmetic mean (the paper cites Vovk & Wang 2020), NA if
   either is; Benjamini-Hochberg over every non-NA $$p$$.
4. **k-means line.** Missing totals count as 0. One point per distinct $$(A, R)$$ with both above
   0 (ARC needs 2 and drops its ATAC exclusions; this does neither). Per axis, on the points' raw
   totals: $$q_1$$ is the 0.99 quantile (R's type 7); $$q_2$$ the 0.99 quantile of the values
   above $$q_1/10$$ (if there are none, the threshold is 0); the threshold is the 0.99 quantile
   of all values below $$q_2$$, divided by 10 ($$q_1/10$$ if none is below). Points above both thresholds start as
   cells. If that is every point, or none (a bug), every point is a cell and the line has slope and
   intercept 0. Otherwise, in $$(x, y) = (\log_{10}(A + 0.1), \log_{10}(R + 0.1))$$, R's `kmeans`
   (Hartigan–Wong, at most 10 iterations, one start at the two groups' means). The cell cluster is
   the one with the larger RNA centre; duplicates share their point's label; a barcode with a zero
   is not a k-means cell. The line $$y = \beta x + \alpha$$ is the perpendicular bisector of the
   two final centres. Hartigan–Wong's partition need not be the nearest-centre split, so the
   `k_means` column and the line can disagree near the line.
5. **Two parallel lines**, slope $$\beta$$. *Tip*: the barcodes whose $$y - \beta x$$ exceeds its
   0.99 quantile; its centre is their mean $$R$$ and mean $$A$$, then $$\log_{10}(\cdot + 0.1)$$.
   *Higher line*: $$\alpha_{hi} = \alpha + \tfrac{2}{3} d \sqrt{\beta^2 + 1}$$, with $$d$$ meant
   as the tip centre's distance to the k-means line. The code passes the point with its
   coordinates swapped, so
   $$d = \lvert \beta \log_{10}(\bar R + 0.1) - \log_{10}(\bar A + 0.1) + \alpha \rvert / \sqrt{\beta^2 + 1}$$.
   *Lower line*: through $$(\log_{10}(l_{ATAC} + 0.1), \log_{10}(l_{RNA} + 0.1))$$, its intercept
   $$\alpha_{lo}$$ capped at $$\alpha$$.
6. **FDR.** `FDR` is step 3's adjusted value, set to 1 for barhop barcodes of either modality,
   for barcodes below the lower line ($$y < \beta x + \alpha_{lo}$$) and where NA. `FDR_multi` is
   `FDR`, but 0 above the higher line ($$y - \beta x > \alpha_{hi}$$, column `above_ambiguous`).

The package returns both columns and has no switch, so there are two call modes. The README calls
`FDR_multi` below 0.001, accepting very large droplets. The paper calls `FDR` below 0.001: it
calls the acceptance optional and says its results did not use it, as its script confirms.

Step 2's Monte Carlo is seeded per modality. The RNA mixture's start weights are drawn from R's
`runif` before the package seeds anything, so they follow the caller's state (the paper's script
calls `set.seed(0)` first); mixtools' restarts are random too. The `set.seed(123)` before `kmeans`
has no effect, as its start is given.

## Versions

| EmptyDropsMultiome | Change to the call |
|---|---|
| `first_release` (`a1b1bbbd`, 2024-02-21) | the only release; `main` (2024-06-03) changes only the README and a figure |

Earlier states (Fisher's method before the mean, the higher line at 1/2 rather than 2/3, no
`seed`; May to August 2023) are not worth reproducing. The call also depends on DropletUtils:
with no `alpha` passed, emptyDrops tests a Dirichlet-multinomial in 1.2 – 1.26 and a multinomial
from 1.28. The README installs R 4.3.0 with DropletUtils 1.20 or 1.22 (which one is not known), so
the paper's results are Dirichlet-multinomial. A call is a pair (EmptyDropsMultiome, DropletUtils
before or from 1.28). The paper's analyses vendor mixtools 2.0.0.

## Parameters

| Original | Default | Meaning |
|---|---|---|
| `lower_rna`, `barhop_rna` | `NULL`: estimated | the RNA bounds; both must be given to skip the fit |
| `lower_atac`, `barhop_atac` | `NULL`: estimated | the ATAC bounds, likewise |
| `niter_rna` | 10,000 | emptyDrops iterations for RNA |
| `niter_atac` | 25,000 | emptyDrops iterations for ATAC |
| `seed` | 42 | set before each emptyDrops |
| `verbose` | `TRUE` | prints and draws plots (needs a graphics device) |
| FDR cut (not exposed) | 0.001 | applied by the user to `FDR` or `FDR_multi` |
| not exposed | — | the mixture starts, checks and defaults; 1.5 and 2 SDs; pseudocount 0.1; tip quantile 0.99; factor 2/3; 10 k-means iterations |

## Status

Planned for 0.3. Implementing it needs:

- two modalities at once: the RNA counts and the ATAC peak counts of the same barcodes, both in
  Cell Ranger ARC's `raw_feature_bc_matrix.h5` (features "Gene Expression" and "Peaks"). No
  method in the package takes two feature types together yet;
- ports of mixtools' `normalmixEM` with its random start and restarts (R's `runif`, from the
  package's port of R's random numbers) and of R's Hartigan–Wong `kmeans` (Fortran, AS 136);
- [emptyDrops](../dropletutils/emptydrops.md), already implemented, at the DropletUtils release
  the run used.

The oracle is the R original in a pixi environment (conda-forge R and r-mixtools, bioconda
DropletUtils), installed from commit `a1b1bbbd`, with a `pdf()` device for the plots and
`set.seed(0)` before the call, as the paper. Two environments: DropletUtils 1.22 (the paper's)
and 1.28 or later. Compare `lower_*`, `barhop_*`, every p-value, `FDR`, `FDR_multi`, `k_means`,
the three intercepts and both call sets, on ARC 2.0.0's pbmc_granulocyte_sorted_10k (the paper's
PBMC sample) and pbmc_granulocyte_sorted_3k. With the ports an exact match is plausible;
otherwise the Monte Carlo rule applies. Not a blocker, but the package is on neither CRAN,
Bioconductor nor conda: it installs from GitHub.

## References

- Megas, Lorenzi & Marioni 2024, EmptyDropsMultiome discriminates real cells from background in
  single-cell multiomics assays, *Genome Biology* 25:121,
  [doi:10.1186/s13059-024-03259-x](https://doi.org/10.1186/s13059-024-03259-x).
- Source (GPL-3): [MarioniLab/EmptyDropsMultiome](https://github.com/MarioniLab/EmptyDropsMultiome/tree/first_release)
  at tag `first_release` (`R/emptydrops_multiome.R` and the files it calls), archived as
  [doi:10.5281/zenodo.11051741](https://doi.org/10.5281/zenodo.11051741).
- The paper's analyses: [MarioniLab/EmptyDropsMultiome2023_final](https://github.com/MarioniLab/EmptyDropsMultiome2023_final)
  at tag `first_release_reprod_repo`.
- Research notes: docs/_research/planned-atac.md §3.
