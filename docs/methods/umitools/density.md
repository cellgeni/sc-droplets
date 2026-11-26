---
title: density
parent: UMI-tools
grand_parent: Methods
nav_order: 2
method_ids: [umitools-density]
status: implemented
---

# UMI-tools whitelist, density valley
{: .no_toc }

`umi_tools whitelist --knee-method=density`, the only rule in 0.5 and the default until 1.0.0.
It smooths the distribution of $$\log_{10}$$ counts with a Gaussian kernel and cuts at a valley
between the barcodes of empty droplets and those of cells. Which valley is chosen is decided by
fixed rules, or by an expected number of cells.

1. TOC
{:toc}

## At a glance

| Id | `umitools-density` |
|---|---|
| Tool | UMI-tools 0.5.0 and later (`whitelist --knee-method=density`) |
| Category | density |
| Counts | reads per barcode in UMI-tools; here UMIs, or reads with `--unit reads` |
| Needs expected cells | no (`--expect-cells` optional) |
| Deterministic | yes |
| Original licence | MIT |
| Validated against | `umi_tools` 0.5.5, 1.0.0 and 1.1.6 (exact) |

## How it works

Barcodes with a count of 0 are left out. With $$M$$ the largest count:

1. **Kept for the density:** the barcodes with count $$c > 0.001\,M$$ (strictly);
   $$x_j = \log_{10} c_j$$ for them.
2. **Density:** a Gaussian kernel density of the $$x_j$$ (scipy's `gaussian_kde` with
   `bw_method=0.1`: bandwidth $$h = 0.1\,\hat\sigma$$, with $$\hat\sigma$$ the standard deviation
   of the $$x_j$$ with denominator $$n - 1$$):
   $$f(t) = \frac{1}{n h} \sum_j \varphi\!\left(\frac{t - x_j}{h}\right)$$
3. **Grid:** $$t_0, \dots, t_{9999}$$, 10,000 evenly spaced points from $$\min x$$ to $$\max x$$.
4. **Valleys:** the grid points $$t_i$$, not at either end, with $$f(t_i) < f(t_{i-1})$$ and
   $$f(t_i) < f(t_{i+1})$$ (strictly: a valley with a flat bottom has none).
5. **Choice:** the valleys are tried from the highest $$t_i$$ down. For each, $$n_i$$ is the number
   of barcodes (all of them, not only those of step 1) with $$c > 10^{t_i}$$. The first valley that
   passes is chosen:
   - with `--expect-cells E`: $$0.1\,E < n_i \le E$$;
   - otherwise: $$i \ge 2000$$ (at least a fifth of the way along the grid), and
     $$\max x - t_i > 0.5$$ or $$t_i < \max x / 2$$.
6. **Cells:** the barcodes with $$c > 10^{t_i}$$ for the chosen valley. If no valley passes,
   UMI-tools stops with an error ("No local minima was accepted"); with
   `--allow-threshold-error` it reports no cells, and so does this.

With `--n-cells N` (`--set-cell-number`), the density is still computed, but the cells are the
barcodes with a count strictly above the $$(N+1)$$-th largest count, as in
[distance](distance.md).

### Arithmetic

The density, the logarithms and the powers of ten are computed with the same numpy and scipy
calls as UMI-tools (`np.log10`, `scipy.stats.gaussian_kde`, `np.linspace`, `np.power`), so given
the same numpy and scipy the grid and the valleys are the same to the last bit. Between scipy
releases the density's last bits can differ; a valley moves only if two neighbouring grid values
are within those bits of each other.

## Versions

| UMI-tools | Variant | `--tool-version` | Behaviour |
|---|---|---|---|
| 0.5.0 – 1.1.6 | `v0.5` | `0.5.0` | as above (the default rule until 1.0.0) |

## Parameters

| Option | Default | Original | Meaning |
|---|---|---|---|
| `--expect-cells` | none | `--expect-cells` | choose the first valley with between a tenth of this and this many barcodes above it |
| `--n-cells` | none | `--set-cell-number` | call the barcodes above the (N+1)-th count instead of a valley |
| `--allow-threshold-error` | off | `--allow-threshold-error` | call no cells, instead of failing, when no valley passes |
| `--unit` | `umis` | `--method` | count `umis` (the matrix's totals) or `reads` (needs `--reads`) |

## Run it

```bash
sc-droplets call umitools-density raw_feature_bc_matrix.h5 --barcodes whitelist.tsv
sc-droplets call umitools-density raw_feature_bc_matrix.h5 --expect-cells 1000 --json
```

```python
import sc_droplets as scd

adata = scd.read("raw_feature_bc_matrix.h5")
result = scd.call(adata, "umitools-density", expect_cells=1000)
result.thresholds["umis"]
```

The thresholds hold the cut ($$10^{t_i}$$, in the unit used; cells are strictly above it).
`sc-droplets plot diagnostics` draws the rank curve with the cut, and the density with every
valley, the chosen one marked.

## Differences from the original

- UMI-tools counts reads (or UMIs) itself from the first 100 million reads; this takes the
  counts from the matrix. Given the same counts, the calls are the same.
- UMI-tools 0.5 ran on Python 2 and 3; nothing in this rule depends on which.

## Validation

Compared with UMI-tools' own function (`umi_methods.getKneeEstimate` in 0.5.5,
`whitelist_methods.getCellWhitelist(knee_method="density")` in 1.0.0 and 1.1.6, each in its own
environment from bioconda, with the numpy and scipy of its time), given the per-barcode UMI
totals of 10x's pbmc_1k_v3 and pbmc4k raw matrices: without options, with `--expect-cells`, and
with `--set-cell-number`. Identical barcodes.

## References

- Smith, Heger & Sudbery 2017, UMI-tools, *Genome Res* 27:491–499,
  [doi:10.1101/gr.209601.116](https://doi.org/10.1101/gr.209601.116).
- Tom Smith, "Estimating the number of true cell barcodes in single cell RNA-Seq" (CGAT blog,
  2017), which introduced the rule.
- Source: [`umi_tools/whitelist_methods.py`](https://github.com/CGATOxford/UMI-tools/blob/1.1.6/umi_tools/whitelist_methods.py)
  (`getKneeEstimateDensity`, 1.1.6);
  [`umi_tools/umi_methods.py`](https://github.com/CGATOxford/UMI-tools/blob/0.5.0/umi_tools/umi_methods.py)
  (`getKneeEstimate`, 0.5.0).
