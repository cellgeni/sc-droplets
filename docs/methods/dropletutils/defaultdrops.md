---
title: defaultDrops
parent: DropletUtils
grand_parent: Methods
nav_order: 4
method_ids: [dropletutils-defaultdrops]
status: implemented
---

# DropletUtils defaultDrops
{: .no_toc }

DropletUtils' version of Cell Ranger 2's cell calling: take the 99th percentile of the totals of
the top `expected` barcodes, divide by 10, and call every barcode above it. Unlike Cell Ranger it
uses R's interpolated quantile and a strict `>`, and has no bootstrap.

1. TOC
{:toc}

## At a glance

| Id | `dropletutils-defaultdrops` |
|---|---|
| Tool | DropletUtils 1.0 and later (`defaultDrops()`) |
| Category | rank-threshold |
| Counts | column sums of the matrix (unrounded) |
| Needs expected cells | yes: 3,000 unless told otherwise |
| Deterministic | yes |
| Original licence | GPL-3 |
| Validated against | `DropletUtils::defaultDrops` in R (exact) |

## How it works

1. Library sizes $$\ell_b$$ are the column sums, unrounded.
2. Take the `expected` largest (all of them, if there are fewer).
3. $$q$$ is their `upper.quant` quantile, R's type 7: with $$h = 1 + (n - 1) \cdot 0.99$$,
   $$q = (1 - f) x_{(\lfloor h \rfloor)} + f x_{(\lceil h \rceil)}$$, $$f = h - \lfloor h \rfloor$$,
   over the ascending order statistics; when the two order statistics are equal the quantile is
   that value, without the interpolation's rounding (R 3.5 to 4.5 all do this).
4. The threshold is $$q \cdot$$ `lower.prop` (0.1); cells are the barcodes with $$\ell_b >$$ it.

## Versions

| DropletUtils | Variant | `--tool-version` | Behaviour |
|---|---|---|---|
| 1.0 – 1.34 | `v1.0` | `1.0` | as above; only the function's packaging changed (an S4 generic from 1.12) |

## Parameters

| Option | Default | Original | Meaning |
|---|---|---|---|
| `--expected-cells` | 3000 | `expected` | how many top barcodes the quantile is taken over |
| `--upper-quant` | 0.99 | `upper.quant` | the quantile |
| `--lower-prop` | 0.1 | `lower.prop` | fraction of the quantile that is the threshold |

## Run it

```bash
sc-droplets call dropletutils-defaultdrops raw_feature_bc_matrix.h5 --expected-cells 5000 --barcodes cells.tsv
```

```python
import sc_droplets as scd

adata = scd.read("raw_feature_bc_matrix.h5")
result = scd.call(adata, "dropletutils-defaultdrops", expected_cells=5000)
result.thresholds["umis"]
```

## Differences from the original

None known.

## Validation

Compared with `DropletUtils::defaultDrops` from bioconda (twelve releases, 1.0.3 under R 3.5.1 to
1.30.0 under R 4.5.3) on 10x's pbmc_1k_v3 raw matrix and a simulated matrix: identical calls.

## References

- Zheng et al. 2017, Massively parallel digital transcriptional profiling of single cells,
  *Nat Commun* 8:14049, [doi:10.1038/ncomms14049](https://doi.org/10.1038/ncomms14049).
- Source: [`R/defaultDrops.R`](https://code.bioconductor.org/browse/DropletUtils/blob/RELEASE_3_22/R/defaultDrops.R).
