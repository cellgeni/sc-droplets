---
title: Kernel-density valley
parent: Generic rules
grand_parent: Methods
nav_order: 6
method_ids: [kde-valley]
status: implemented
---

# Kernel-density valley
{: .no_toc }

The distribution of log totals above a floor has two peaks: the ambient droplets and the cells.
The lowest point of a kernel density estimate between the two highest peaks is the cut. The same
idea, with rules of their own for which valley to take, is in UMI-tools' `density`
([page](../umitools/density.md)), alevin's Gauss correction ([page](../alevin/initial-whitelist.md)),
BAP, DropSift and DropletQC.

1. TOC
{:toc}

## At a glance

| Id | `kde-valley` |
|---|---|
| Tool | none (generic) |
| Category | density |
| Counts | UMIs, reads or genes detected |
| Needs expected cells | no |
| Deterministic | yes |
| Validated against | simulated two-population matrices with a known valley; SciPy's `gaussian_kde` |

## How it works

1. Keep the barcodes whose value $$v_b$$ is at least `--min-counts` (100). The others are empty.
2. $$x_b = \log_{10} v_b$$ (the C library's `log10`). A barcode with no count is never a
   candidate, whatever the floor.
3. A Gaussian kernel density of $$x$$ (SciPy's `gaussian_kde`) with bandwidth factor
   `--bandwidth`: `scott` (default, $$n^{-1/5}$$), `silverman`, or a number, which multiplies the
   standard deviation of $$x$$ (UMI-tools uses 0.1).
4. The density on 10,000 equally spaced points from the smallest $$x$$ to the largest, as
   UMI-tools evaluates it.
5. The peaks are the grid points higher than both neighbours. The two highest (the lower grid
   point first when two are equally high) bound the search; the valley is the first grid point
   between them with the lowest density.
6. Cells have $$v_b > 10^{x^*}$$, where $$x^*$$ is the valley's grid point.

With fewer than two peaks there is no valley, and the call fails.

On 10x's pbmc_1k_v3 and pbmc4k raw matrices the default calls 1,164 and 4,768 cells (Cell Ranger:
1,222 and 4,340). With every barcode that has a count, the highest peaks are at one or two UMIs
and the valley between them is meaningless: the rule needs the floor.

## Versions

None: the rule does not follow the releases of a tool.

## Parameters

| Option | Default | Original | Meaning |
|---|---|---|---|
| `--min-counts` | 100 | | candidates have at least this value; the others are empty |
| `--bandwidth` | `scott` | `bw_method` | `scott`, `silverman`, or a factor such as 0.1 |
| `--unit` | `umis` | | `umis`, `reads` (needs `--reads`) or `genes` |

## Run it

```bash
sc-droplets call kde-valley raw_feature_bc_matrix.h5 --barcodes cells.tsv.gz
sc-droplets call kde-valley raw_feature_bc_matrix.h5 --bandwidth 0.1 --min-counts 10 -o calls.h5ad
```

```python
import sc_droplets as scd

adata = scd.read("raw_feature_bc_matrix.h5")
result = scd.call(adata, "kde-valley", bandwidth="silverman")
result.thresholds["umis"], result.thresholds["log10"]
```

The thresholds hold the cut in counts and its $$\log_{10}$$ (`log10`). `sc-droplets plot
diagnostics` draws the density with its two peaks and the valley.

## Differences from the original

Not applicable: the rule is this package's. The tools that use a valley choose it differently:
UMI-tools takes the first valley (from the top) that passes its tests, alevin the first with a
barcode count in a range, DropletQC the first minimum of the density of the nuclear fraction.
Their own pages reproduce them.

## Validation

Unit tests: on matrices simulated from two log-normal populations, the valley falls within a
grid step of the minimum of the true mixture density; the density is SciPy's on the same grid;
with one population the call fails.

## References

- Smith, Heger & Sudbery 2017, UMI-tools: modeling sequencing errors in Unique Molecular
  Identifiers to improve quantification accuracy, *Genome Research* 27:491,
  [doi:10.1101/gr.209601.116](https://doi.org/10.1101/gr.209601.116).
- Muskovic & Powell 2021, DropletQC: improved identification of empty droplets and damaged cells
  in single-cell RNA-seq data, *Genome Biology* 22:329,
  [doi:10.1186/s13059-021-02547-0](https://doi.org/10.1186/s13059-021-02547-0).
- Scott 1992, *Multivariate Density Estimation*, Wiley; Silverman 1986, *Density Estimation for
  Statistics and Data Analysis*, Chapman and Hall.
