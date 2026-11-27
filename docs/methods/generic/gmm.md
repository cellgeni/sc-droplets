---
title: Gaussian mixture
parent: Generic rules
grand_parent: Methods
nav_order: 7
method_ids: [gmm]
status: implemented
---

# Gaussian mixture
{: .no_toc }

A mixture of Gaussians fitted to the log totals of the barcodes above a floor; the component with
the highest mean is the cells. zUMIs did this before 2019 (with mclust), cellsweep's `mx_filter`
does it on $$\log(1 + x)$$ of every barcode, and quality-control tools fit one per metric. The fit
is scikit-learn's `GaussianMixture`, with its defaults unless told otherwise.

1. TOC
{:toc}

## At a glance

| Id | `gmm` |
|---|---|
| Tool | none (generic); scikit-learn's `GaussianMixture` |
| Category | mixture |
| Counts | UMIs, reads or genes detected |
| Needs expected cells | no |
| Deterministic | yes, for a given seed (k-means initialisation, seed 0) |
| Original licence | BSD-3-Clause (scikit-learn) |
| Validated against | scikit-learn 1.5.2 and 1.9.1, `GaussianMixture` (same parameters and calls) |

## How it works

1. Keep the barcodes whose value $$v_b$$ is at least `--min-counts` (100) and above zero. The
   others are empty.
2. $$x_b = \log_{10} v_b$$ (the C library's `log10`).
3. Fit a mixture of `--components` (2) Gaussians to $$x$$ by expectation-maximisation:
   scikit-learn's `GaussianMixture(n_components, tol, max_iter, random_state)`, initialised by
   k-means with the seed (`--seed`, default 0). It stops when the mean log-likelihood per barcode changes by less
   than `--tol` (scikit-learn's default, $$10^{-3}$$) or after `--max-iter` (100) iterations.
4. The cell component is the one with the highest mean, $$\mu_c$$ with standard deviation
   $$\sigma_c$$. With `--rule posterior` (default), cells are the candidates whose posterior
   probability of that component is above 0.5. With `--rule quantile`, as zUMIs did, cells have
   $$x_b \ge \mu_c + \sigma_c\,\Phi^{-1}(q)$$, the `--quantile` $$q$$ (0.01) of the cell component.

The threshold reported is the smallest total called, and for `quantile` the cut itself.

### The fit depends on when it stops

The log totals of a raw matrix are not two Gaussians: above the floor, the empty droplets make a
long tail rather than a peak. The likelihood has more than one maximum, and where EM stops
matters:

| On 10x's raw matrices | pbmc_1k_v3 | pbmc4k |
|---|---|---|
| Cell Ranger (3.0, 2.1) | 1,222 | 4,340 |
| `gmm` (scikit-learn's defaults: stops after 4 iterations) | 1,121 | 4,413 |
| `gmm --rule quantile` | 1,108 | 4,351 |
| `gmm --tol 1e-8 --max-iter 10000` (converged) | 1,117 | 3,969 |

On pbmc4k the converged fit moves its second component into the cells (means $$10^{2.97}$$ and
$$10^{3.60}$$) and calls fewer barcodes; the early stop is closer to Cell Ranger by luck of where it
started. R's mclust (`Mclust(x, G = 2, modelNames = "V")`, 6.1) reaches the same optimum on
pbmc4k but stops at a lower one on pbmc_1k_v3 (log-likelihood −545.17, against −544.46 for the
converged fit here). zUMIs chose the number of components by BIC; on these matrices BIC picks 3
and 8 components, and the 1 % quantile of the top one calls 383 and 3,275 barcodes.

## Versions

None: the rule does not follow the releases of a tool.

## Parameters

| Option | Default | Original | Meaning |
|---|---|---|---|
| `--components` | 2 | `n_components` | Gaussians in the mixture |
| `--rule` | `posterior` | | `posterior` (above 0.5 for the cell component) or `quantile` |
| `--quantile` | 0.01 | | for `--rule quantile`: the cell component's quantile that is the cut |
| `--min-counts` | 100 | | candidates have at least this value; the others are empty |
| `--tol` | 0.001 | `tol` | stop when the mean log-likelihood changes by less |
| `--max-iter` | 100 | `max_iter` | stop after this many EM iterations |
| `--unit` | `umis` | | `umis`, `reads` (needs `--reads`) or `genes` |

## Run it

```bash
sc-droplets call gmm raw_feature_bc_matrix.h5 --barcodes cells.tsv.gz
sc-droplets call gmm raw_feature_bc_matrix.h5 --rule quantile --tol 1e-8 --max-iter 10000 -o calls.h5ad
```

```python
import sc_droplets as scd

adata = scd.read("raw_feature_bc_matrix.h5")
result = scd.call(adata, "gmm", components=3)
result.thresholds["umis"], result.diagnostics["means"]
adata.obs["gmm_probability"]
```

Each barcode's posterior probability of the cell component is stored (`<key>_probability`; 0
below the floor), with the fitted means, standard deviations and weights in the diagnostics.
`sc-droplets plot diagnostics` draws the candidates' histogram with each component's density.

## Differences from the original

Not applicable: scikit-learn's estimator is applied as it is. zUMIs' own mixture (mclust, before
2.5.6b) is a documented gap on the [zUMIs page](../zumis/cell-barcodes.md).

## Validation

`oracles/generic/run.sh` fits `GaussianMixture` in scikit-learn 1.5.2 and 1.9.1 (conda-forge) to
the same log totals of 10x's pbmc_1k_v3 and pbmc4k raw matrices, with 2 and 3 components, the
default and the converged tolerance: the same means, standard deviations and weights (to
$$10^{-9}$$ relative) and the same calls. It also fits mclust 6.1 (`modelNames = "V"`): the
converged fit here is at least as likely as mclust's.

## References

- Pedregosa et al. 2011, Scikit-learn: machine learning in Python, *Journal of Machine Learning
  Research* 12:2825.
- Scrucca, Fop, Murphy & Raftery 2016, mclust 5: clustering, classification and density
  estimation using Gaussian finite mixture models, *The R Journal* 8:289,
  [doi:10.32614/RJ-2016-021](https://doi.org/10.32614/RJ-2016-021).
- Parekh, Ziegenhain, Vieth, Enard & Hellmann 2018, zUMIs - A fast and flexible pipeline to
  process RNA sequencing data with UMIs, *GigaScience* 7(6),
  [doi:10.1093/gigascience/giy059](https://doi.org/10.1093/gigascience/giy059).
