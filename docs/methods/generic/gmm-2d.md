---
title: Gaussian mixture in two dimensions
parent: Generic rules
grand_parent: Methods
nav_order: 8
method_ids: [gmm-2d]
status: implemented
---

# Gaussian mixture in two dimensions
{: .no_toc }

A mixture of Gaussians with full covariances fitted to two measures of each barcode above a
floor: $$\log_{10}$$ UMIs and $$\log_{10}$$ genes detected. The component with the highest mean
UMIs is the cells. A second measure separates droplets that one cannot: an empty droplet with
many UMIs from a few highly expressed ambient genes detects fewer genes than a cell with as many
UMIs. DropletQC fits such a mixture to (nuclear fraction, log UMIs) to find damaged cells; the
nuclear fraction comes with the splicing readers.

1. TOC
{:toc}

## At a glance

| Id | `gmm-2d` |
|---|---|
| Tool | none (generic); scikit-learn's `GaussianMixture` |
| Category | mixture |
| Counts | UMIs and genes detected |
| Needs expected cells | no |
| Deterministic | yes, for a given seed (k-means initialisation, seed 0) |
| Original licence | BSD-3-Clause (scikit-learn) |
| Validated against | scikit-learn 1.5.2 and 1.9.1, `GaussianMixture` (same parameters and calls) |

## How it works

1. Keep the barcodes with at least `--min-counts` (100) UMIs. The others are empty.
2. $$x_b = (\log_{10} u_b, \log_{10} g_b)$$, with $$u_b$$ the UMIs and $$g_b$$ the genes with a
   count (the C library's `log10`).
3. Fit a mixture of `--components` (2) Gaussians with full covariance matrices:
   scikit-learn's `GaussianMixture(n_components, covariance_type="full", tol, max_iter,
   random_state)`, initialised by k-means with the seed (`--seed`, default 0).
4. The cell component is the one with the highest mean $$\log_{10}$$ UMIs. Cells are the
   candidates whose posterior probability of that component is above 0.5.

The threshold reported is the smallest UMI total called.

| On 10x's raw matrices | pbmc_1k_v3 | pbmc4k |
|---|---|---|
| Cell Ranger (3.0, 2.1) | 1,222 | 4,340 |
| `gmm-2d` (scikit-learn's defaults) | 1,151 | 4,475 |
| `gmm-2d --tol 1e-8 --max-iter 10000` (converged) | 1,151 | 4,451 |
| [`gmm`](gmm.md), log UMIs alone | 1,121 | 4,413 |

The second dimension makes the fit steadier: converged or not, it calls nearly the same barcodes,
where the one-dimensional mixture moves by 444 on pbmc4k.

## Versions

None: the rule does not follow the releases of a tool.

## Parameters

| Option | Default | Original | Meaning |
|---|---|---|---|
| `--components` | 2 | `n_components` | Gaussians in the mixture |
| `--min-counts` | 100 | | candidates have at least this many UMIs; the others are empty |
| `--tol` | 0.001 | `tol` | stop when the mean log-likelihood changes by less |
| `--max-iter` | 100 | `max_iter` | stop after this many EM iterations |

## Run it

```bash
sc-droplets call gmm-2d raw_feature_bc_matrix.h5 --barcodes cells.tsv.gz
sc-droplets call gmm-2d raw_feature_bc_matrix.h5 --components 3 -o calls.h5ad
```

```python
import sc_droplets as scd

adata = scd.read("raw_feature_bc_matrix.h5")
result = scd.call(adata, "gmm-2d")
adata.obs["gmm_2d_probability"]
```

Each barcode's posterior probability of the cell component is stored (`<key>_probability`; 0
below the floor), with the fitted means, covariances and weights in the diagnostics.
`sc-droplets plot diagnostics` draws the candidates in the plane, coloured by call, with each
component's two-standard-deviation ellipse.

## Differences from the original

Not applicable: scikit-learn's estimator is applied as it is.

## Validation

`oracles/generic/run.sh` fits `GaussianMixture` in scikit-learn 1.5.2 and 1.9.1 (conda-forge) to
the same pairs from 10x's pbmc_1k_v3 and pbmc4k raw matrices, with 2 and 3 components and the
default and converged tolerances: the same means, covariances and weights (to $$10^{-9}$$
relative) and the same calls.

## References

- Pedregosa et al. 2011, Scikit-learn: machine learning in Python, *Journal of Machine Learning
  Research* 12:2825.
- Muskovic & Powell 2021, DropletQC: improved identification of empty droplets and damaged cells
  in single-cell RNA-seq data, *Genome Biology* 22:329,
  [doi:10.1186/s13059-021-02547-0](https://doi.org/10.1186/s13059-021-02547-0).
