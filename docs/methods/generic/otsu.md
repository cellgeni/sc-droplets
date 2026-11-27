---
title: Otsu threshold
parent: Generic rules
grand_parent: Methods
nav_order: 4
method_ids: [otsu]
status: implemented
---

# Otsu threshold
{: .no_toc }

Otsu's method splits a histogram in two where the variance between the two classes is largest.
On the log totals of the barcodes above a floor, the upper class is the cells. CellBender 0.3
uses the same cut to separate cells from empty droplets before it trains (see its
[priors](../cellbender/prior.md)); dropkick uses its multi-class form
([multi-Otsu](multiotsu.md)).

1. TOC
{:toc}

## At a glance

| Id | `otsu` |
|---|---|
| Tool | none (generic); scikit-image's `threshold_otsu` |
| Category | density |
| Counts | UMIs, reads or genes detected |
| Needs expected cells | no |
| Deterministic | yes |
| Original licence | BSD-3-Clause (scikit-image) |
| Validated against | scikit-image 0.20.0 and 0.26.0, `threshold_otsu` (exact) |

## How it works

1. Keep the barcodes whose value $$v_b$$ is at least `--min-counts` (100). The others are empty.
2. $$x_b = \log_{10} v_b$$ (the C library's `log10`); `--transform arcsinh` takes
   $$\operatorname{arcsinh} v_b$$ instead, `none` the values themselves.
3. A histogram of $$x$$ in `--bins` (256) equal bins from its minimum to its maximum
   (`numpy.histogram`), with bin centres $$c_i$$ and counts $$h_i$$ held as 32-bit floats.
4. For each split after bin $$i$$: the weights $$w_1(i) = \sum_{j \le i} h_j$$ and
   $$w_2(i) = \sum_{j > i} h_j$$, the means $$\mu_1(i)$$ and $$\mu_2(i)$$ of the bin centres they
   hold, and the variance between the classes, up to a constant,
   $$\sigma^2_B(i) = w_1(i)\, w_2(i)\, (\mu_1(i) - \mu_2(i))^2$$.
5. The threshold is the bin centre $$t = c_{i^*}$$ at the first $$i^*$$ that maximises
   $$\sigma^2_B$$. Cells are the candidates with $$x_b > t$$, as scikit-image's documentation
   applies it.

When every candidate has the same value, that value is the threshold and nothing is a cell.

On 10x's pbmc_1k_v3 and pbmc4k raw matrices the default calls 1,143 and 4,489 cells (Cell Ranger:
1,222 and 4,340). With every barcode that has a count, the split falls between barcodes with one
or two UMIs and the rest (105,520 "cells" on pbmc_1k_v3): Otsu needs the floor.

## Versions

None: the rule does not follow the releases of a tool. scikit-image 0.20 (2022) started holding
the histogram counts as 32-bit floats, which this follows; before, they were 64-bit, and a
maximum shared by two splits could fall on the other one.

## Parameters

| Option | Default | Original | Meaning |
|---|---|---|---|
| `--min-counts` | 100 | | candidates have at least this value; the others are empty |
| `--transform` | `log10` | | `log10`, `arcsinh` or `none` |
| `--bins` | 256 | `nbins` | histogram bins |
| `--unit` | `umis` | | `umis`, `reads` (needs `--reads`) or `genes` |

## Run it

```bash
sc-droplets call otsu raw_feature_bc_matrix.h5 --barcodes cells.tsv.gz
sc-droplets call otsu raw_feature_bc_matrix.h5 --unit genes --transform arcsinh --min-counts 10 -o calls.h5ad
```

```python
import sc_droplets as scd

adata = scd.read("raw_feature_bc_matrix.h5")
result = scd.call(adata, "otsu")
result.thresholds["umis"], result.thresholds["transformed"]
```

The thresholds hold the cut in counts and on the transformed scale (`transformed`).
`sc-droplets plot diagnostics` draws the histogram with the cut.

## Differences from the original

Not applicable: scikit-image's function is applied as it is, to the candidates' transformed
values.

## Validation

`oracles/generic/run.sh` calls `skimage.filters.threshold_otsu` (scikit-image 0.20.0 and 0.26.0
from conda-forge) on the same transformed values of 10x's pbmc_1k_v3 and pbmc4k raw matrices,
with floors 1, 10, 100 and 200, on the $$\log_{10}$$ and $$\operatorname{arcsinh}$$ scales:
identical thresholds in all 32 cases.

## References

- Otsu 1979, A threshold selection method from gray-level histograms, *IEEE Transactions on
  Systems, Man, and Cybernetics* 9:62,
  [doi:10.1109/TSMC.1979.4310076](https://doi.org/10.1109/TSMC.1979.4310076).
- van der Walt et al. 2014, scikit-image: image processing in Python, *PeerJ* 2:e453,
  [doi:10.7717/peerj.453](https://doi.org/10.7717/peerj.453).
- Source: [`skimage/filters/thresholding.py`](https://github.com/scikit-image/scikit-image/blob/v0.26.0/src/skimage/filters/thresholding.py)
  (`threshold_otsu`).
