---
title: MAD outliers
parent: Generic rules
grand_parent: Methods
nav_order: 3
method_ids: [mad]
status: implemented
---

# MAD outliers
{: .no_toc }

Among the barcodes above a floor, the ones whose log total is more than three median absolute
deviations (MADs) below the median are empty; the rest are cells. This is the low-count rule of
quality control, `scuttle::isOutlier(log = TRUE, type = "lower")`, applied to cell calling.

1. TOC
{:toc}

## At a glance

| Id | `mad` |
|---|---|
| Tool | none (generic); the rule of scuttle's `isOutlier` (scater before Bioconductor 3.12) |
| Category | density |
| Counts | UMIs, reads or genes detected |
| Needs expected cells | no, but needs a floor that leaves out most empty droplets |
| Deterministic | yes |
| Original licence | GPL-3.0 (scuttle) |
| Validated against | scuttle 1.0.4 and 1.22.0, `isOutlier` in R (exact) |

## How it works

1. Keep the barcodes whose value $$v_b$$ (total UMIs by default) is at least `--min-counts`
   (100: DropletUtils' `lower`, below which emptyDrops takes a barcode to be empty). They are the
   candidates; every other barcode is empty. With `--log`, a barcode with no count is never a
   candidate, whatever the floor: its log is $$-\infty$$.
2. $$m_b = \log_2 v_b$$ (the C library's `log2`, as R calls it).
3. The median $$\tilde m$$ and the MAD, scaled to estimate a standard deviation under normality:
   $$\operatorname{MAD} = 1.4826 \cdot \operatorname{median}_b |m_b - \tilde m|$$ (R's `median` and
   `mad`; the median of an even number of values is the mean of the middle two).
4. The lower limit $$L = \tilde m - \max(d, k \cdot \operatorname{MAD})$$, with $$k$$ = `--nmads`
   (3) and $$d$$ = `--min-diff` (none by default, as scuttle's `min.diff = NA`).
5. A candidate is a cell unless $$m_b < L$$. High outliers are kept: they are large cells or
   doublets, not empty droplets.

The threshold reported is $$2^{L}$$, as `isOutlier` reports it.

The rule assumes the candidates are mostly cells, with empty droplets in the lower tail. On a raw
matrix that depends on the floor: with every barcode that has a count, the median lands among
the empty droplets and nothing is an outlier. With the default floor of 100 UMIs it calls 1,160
cells on 10x's pbmc_1k_v3 (Cell Ranger 3.0: 1,222) and 4,541 on pbmc4k (Cell Ranger 2.1: 4,340).

## Versions

None: the rule does not follow the releases of a tool. `isOutlier` moved from scater to scuttle
in Bioconductor 3.12 (2020) with the same arithmetic.

## Parameters

| Option | Default | Original | Meaning |
|---|---|---|---|
| `--nmads` | 3 | `nmads` | how many MADs below the median a barcode may be |
| `--min-counts` | 100 | | candidates have at least this value; the others are empty |
| `--log` / `--no-log` | `--log` | `log` | work on $$\log_2$$ values |
| `--min-diff` | none | `min.diff` | the least distance from the median to the limit |
| `--unit` | `umis` | | `umis`, `reads` (needs `--reads`) or `genes` |

## Run it

```bash
sc-droplets call mad raw_feature_bc_matrix.h5 --barcodes cells.tsv.gz
sc-droplets call mad raw_feature_bc_matrix.h5 --nmads 2 --min-counts 200 -o calls.h5ad
```

```python
import sc_droplets as scd

adata = scd.read("raw_feature_bc_matrix.h5")
result = scd.call(adata, "mad", nmads=3)
result.thresholds["umis"], result.thresholds["median"], result.thresholds["mad"]
```

The thresholds hold the limit in counts (`umis`, `reads` or `genes`), the candidates' `median`
and their `mad` (in $$\log_2$$ units when `--log`). `sc-droplets plot diagnostics` draws the
candidates' distribution with the median and the limit.

## Differences from the original

`isOutlier` is a quality-control rule for barcodes already called cells; here it is the whole
call, on the candidates above the floor. `batch` and `subset` are not offered.

## Validation

`oracles/generic/run.sh` calls `scuttle::isOutlier(x, nmads, type = "lower", log = TRUE)` (scuttle
1.0.4 in R 4.0 and 1.22.0 in R 4.5, from bioconda) on the candidates of 10x's pbmc_1k_v3 and pbmc4k
raw matrices (UMI totals from the matrix), with floors 50, 100 and 200 and 2 or 3 MADs: identical
limits and calls in all 24 cases. With a floor of 50 on pbmc4k, the candidates' median falls among
the empty droplets and every candidate is kept: the floor matters.

## References

- McCarthy, Campbell, Lun & Wills 2017, Scater: pre-processing, quality control, normalization and
  visualization of single-cell RNA-seq data in R, *Bioinformatics* 33:1179,
  [doi:10.1093/bioinformatics/btw777](https://doi.org/10.1093/bioinformatics/btw777).
- Amezquita et al. 2020, Orchestrating single-cell analysis with Bioconductor, *Nature Methods*
  17:137, [doi:10.1038/s41592-019-0654-x](https://doi.org/10.1038/s41592-019-0654-x).
- Source: [`scuttle/R/isOutlier.R`](https://github.com/LTLA/scuttle/blob/master/R/isOutlier.R).
