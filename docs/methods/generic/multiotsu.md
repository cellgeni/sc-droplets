---
title: Multi-Otsu thresholds
parent: Generic rules
grand_parent: Methods
nav_order: 5
method_ids: [multiotsu]
status: implemented
---

# Multi-Otsu thresholds
{: .no_toc }

Otsu's method with more than two classes: the thresholds that maximise the variance between the
classes. On the log totals of every barcode with a count, three classes separate barcodes with
almost nothing, the ambient plateau, and the cells; the top class is the cells. dropkick labels
its training droplets this way, on $$\operatorname{arcsinh}$$ of the genes detected.

1. TOC
{:toc}

## At a glance

| Id | `multiotsu` |
|---|---|
| Tool | none (generic); scikit-image's `threshold_multiotsu` |
| Category | density |
| Counts | UMIs, reads or genes detected |
| Needs expected cells | no |
| Deterministic | yes |
| Original licence | BSD-3-Clause (scikit-image) |
| Validated against | scikit-image 0.20.0 and 0.26.0, `threshold_multiotsu` (exact) |

## How it works

1. Keep the barcodes whose value $$v_b$$ is at least `--min-counts` (1: every barcode with a
   count). The others are empty.
2. $$x_b = \log_{10} v_b$$ (the C library's `log10`), or $$\operatorname{arcsinh} v_b$$, or
   $$v_b$$ (`--transform`).
3. A histogram of $$x$$ in `--bins` (256) equal bins, normalised to probabilities $$p_i$$ and held
   as 32-bit floats, as everything below is.
4. Cumulative moments $$P_i = \sum_{j \le i} p_j$$ and $$S_i = p_0 + \sum_{1 \le j \le i} j\,p_j$$.
   The first term of $$S$$ is $$p_0$$, not $$0 \cdot p_0$$: scikit-image's code starts it so.
5. For each range of bins $$[u, v]$$, $$H(u, v) = (S_v - S_{u-1})^2 / (P_v - P_{u-1})$$, or 0 when the
   range is empty ($$H(0, 0)$$ is left at 0).
6. For $$C$$ classes, every increasing choice of $$C - 1$$ split bins $$t_1 < \dots < t_{C-1}$$
   scores $$\sigma = H(0, t_1) + H(t_{C-1} + 1, n - 1) + \sum_k H(t_k + 1, t_{k+1})$$ (added in that
   order). The first choice in lexicographic order with the largest score wins; the thresholds
   are the bin centres at the split bins.
7. Each candidate's class is the number of thresholds at or below $$x_b$$ (`numpy.digitize`). The
   top class is the cells.

When the histogram has exactly $$C$$ nonzero bins, the thresholds are the centres of the first
$$C - 1$$ of them; with fewer, the call fails.

On 10x's pbmc_1k_v3 and pbmc4k raw matrices the default calls 1,317 and 4,845 cells (Cell Ranger:
1,222 and 4,340).

## Versions

None: the rule does not follow the releases of a tool. scikit-image added `threshold_multiotsu`
in 0.16 (2019); 0.20 moved its arithmetic from 64-bit to 32-bit floats, which this follows.

## Parameters

| Option | Default | Original | Meaning |
|---|---|---|---|
| `--classes` | 3 | `classes` | how many classes, 2 to 4 |
| `--min-counts` | 1 | | candidates have at least this value; the others are empty |
| `--transform` | `log10` | | `log10`, `arcsinh` or `none` |
| `--bins` | 256 | `nbins` | histogram bins |
| `--unit` | `umis` | | `umis`, `reads` (needs `--reads`) or `genes` |

## Run it

```bash
sc-droplets call multiotsu raw_feature_bc_matrix.h5 --barcodes cells.tsv.gz
sc-droplets call multiotsu raw_feature_bc_matrix.h5 --unit genes --transform arcsinh -o calls.h5ad
```

```python
import sc_droplets as scd

adata = scd.read("raw_feature_bc_matrix.h5")
result = scd.call(adata, "multiotsu", classes=3)
result.thresholds["umis"], result.thresholds["threshold_1"]
adata.obs["multiotsu_class"].value_counts()
```

The thresholds hold the top class's cut in counts (`umis`, `reads` or `genes`) and every
threshold on the transformed scale (`threshold_1`, …). Each barcode's class is stored (`class 1`
is the lowest; barcodes below the floor are `below floor`). `sc-droplets plot diagnostics` draws
the histogram with the thresholds.

## Differences from the original

scikit-image's function is applied as it is, with two remarks:

- Its search is exhaustive, so more than 4 classes on 256 bins would take too long and are
  refused.
- Its builds for arm64 (Apple silicon, Linux on ARM) let the compiler fuse the multiply and the
  add of the first moments, $$S_{i-1} + i\,p_i$$, into one rounding. When two choices of splits
  score within a rounding of each other, those builds can pick the other one. This package
  rounds twice, as scikit-image's x86-64 builds do; one case in 300 random histograms differed.

## Validation

`oracles/generic/run.sh` calls `skimage.filters.threshold_multiotsu` (scikit-image 0.20.0 and
0.26.0 from conda-forge, x86-64) on the same transformed values of 10x's pbmc_1k_v3 and pbmc4k raw
matrices, with 3 and 4 classes, floors 1, 10, 100 and 200, on the $$\log_{10}$$ and
$$\operatorname{arcsinh}$$ scales: identical thresholds in all 64 cases.

## References

- Liao, Chen & Chung 2001, A fast algorithm for multilevel thresholding, *Journal of Information
  Science and Engineering* 17:713, [doi:10.6688/JISE.2001.17.5.1](https://doi.org/10.6688/JISE.2001.17.5.1).
- Heiser, Wang, Chen, Hughes & Lau 2021, Automated quality control and cell identification of
  droplet-based single-cell data using dropkick, *Genome Research* 31:1742,
  [doi:10.1101/gr.271908.120](https://doi.org/10.1101/gr.271908.120).
- Source: [`skimage/filters/_multiotsu.pyx`](https://github.com/scikit-image/scikit-image/blob/v0.26.0/src/skimage/filters/_multiotsu.pyx)
  (`_get_multiotsu_thresh_indices_lut`).
