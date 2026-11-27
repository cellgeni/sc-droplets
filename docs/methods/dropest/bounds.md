---
title: cell-number bounds
parent: dropEst
grand_parent: Methods
nav_order: 1
method_ids: [dropest-bounds]
status: implemented
---

# dropEst cell-number bounds
{: .no_toc }

dropestr's `EstimateCellsNumber`: lagged derivatives of $$\log$$ UMI count against $$\log$$ rank
give the expected number of cells and two bounds. Each barcode gets one of three classes, as
dropestr labels them for its classifier: `high` (within the lower bound), `unknown`, and `low`
(from the upper bound on). The cells are the expected number from the top.

1. TOC
{:toc}

## At a glance

| Id | `dropest-bounds` |
|---|---|
| Tool | dropestr 0.6.0 and later (`EstimateCellsNumber`; dropEst's reports and its quality classifier) |
| Category | knee |
| Counts | UMIs per barcode |
| Needs expected cells | no |
| Deterministic | yes |
| Original licence | GPL-3.0 |
| Validated against | dropestr 0.7.9 (dropEst 0.8.6), its own function in R (exact) |

## How it works

The barcodes with a count, sorted in decreasing order: $$u_1 \ge \dots \ge u_n$$. With
$$\lambda = \operatorname{round}(0.05\,n)$$ (R's `round`, halves to even):

1. $$x_i = \log i$$, $$y_i = (\log u_{i+\lambda} - \log u_i) / (\log(i + \lambda) - \log i)$$:
   the slope over $$\lambda$$ ranks, placed at rank $$i + \lambda$$.
2. The slope of that slope over $$\lambda$$ points, $$y'$$, placed at rank $$i + 2\lambda$$.
3. The upper bound: the rank where the longest run of $$y' > 0$$ starts (the first such run if
   several are longest), $$\text{max} = \operatorname{round}(e^{x'})$$ at that point.
4. The expected number of cells: the steepest slope before the upper bound,
   $$\text{expected} = \operatorname{round}(e^{x_j} - \lambda / 2)$$ for $$j$$ = the first
   minimum of $$y$$ over its first $$\text{max} - \lambda$$ values.
5. The lower bound: $$\text{min} = \operatorname{round}(0.75 \cdot \text{expected})$$.

The classes are `high` for ranks $$\le$$ min, `low` for ranks $$\ge$$ max, `unknown` between; the
cells are the top `expected` barcodes.

## Versions

| dropestr | Variant | `--tool-version` | Behaviour |
|---|---|---|---|
| 0.6.0 – 0.7.9 | `default` | | the rule above |

## Parameters

| Option | Default | Original | Meaning |
|---|---|---|---|
| `--lag` | 0.05 | `lag` | the derivatives' lag, as a fraction of the barcodes |

## Run it

```bash
sc-droplets call dropest-bounds raw_feature_bc_matrix.h5 --json
sc-droplets call dropest-bounds raw_feature_bc_matrix.h5 --lag 0.1 -o calls.h5ad
```

```python
import sc_droplets as scd

adata = scd.read("raw_feature_bc_matrix.h5")
result = scd.call(adata, "dropest-bounds")
result.thresholds["expected"], result.thresholds["min"], result.thresholds["max"]
adata.obs["dropest_bounds_class"].value_counts()
```

The thresholds hold `expected`, `min`, `max` and `lag`. `sc-droplets plot diagnostics` draws the
rank curve with the cut.

## Differences from the original

- dropestr takes the UMI counts it is given; dropEst's reports pass the counts of the barcodes it
  kept. Here every barcode with a count enters, since `log 0` has no place on the curve.
- Among equal counts at the expected number, the cells are the first in the matrix's order;
  dropestr gives only the numbers.

## Validation

`oracles/dropest/run.sh` sources dropestr's `R/cell_number.R` from dropEst's tag v0.8.6 in R from
conda-forge and calls `EstimateCellsNumber` on the UMI totals of 10x's pbmc_1k_v3 and pbmc4k raw
matrices, with `lag` 0.05 and 0.1. Identical bounds and expected numbers.

## References

- Petukhov, Guo, Baryawno, Severe, Scadden, Samsonova & Kharchenko 2018, dropEst: pipeline for
  accurate estimation of molecular counts in droplet-based single-cell RNA-seq experiments,
  *Genome Biology* 19:78, [doi:10.1186/s13059-018-1449-6](https://doi.org/10.1186/s13059-018-1449-6).
- Source: [`dropestr/R/cell_number.R`](https://github.com/kharchenkolab/dropEst/blob/v0.8.6/dropestr/R/cell_number.R)
  (`EstimateCellsNumber`, `ArrayDerivative`, `GetLongestTrue`).
