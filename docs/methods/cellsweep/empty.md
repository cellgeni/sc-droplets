---
title: empty droplets
parent: cellsweep
grand_parent: Methods
nav_order: 1
method_ids: [cellsweep-empty]
status: implemented
---

# cellsweep empty droplets
{: .no_toc }

`infer_empty_droplets`, the step that decides which barcodes cellsweep treats as empty droplets
when the AnnData does not say. Each rule finds a cutoff on the barcodes' totals; a barcode below
the cutoff is empty, every other barcode is a cell.

1. TOC
{:toc}

## At a glance

| Id | `cellsweep-empty` |
|---|---|
| Tool | cellsweep 0.1.0 and later (`denoise(..., empty_droplet_method=...)`; `kb sweep`) |
| Category | knee |
| Counts | the matrix's totals over every feature (here, the selected feature types: Gene Expression by default) |
| Needs expected cells | no (optional) |
| Deterministic | yes (the Gaussian mixture has a fixed seed) |
| Original licence | BSD-2-Clause |
| Validated against | cellsweep 1.0.1 and 1.0.2 (with kneed 0.8.5, 0.8.6 and 0.9.0) (exact) |

## How it works

The totals are cellsweep's `adata.X.sum(axis=1)`: in the matrix's dtype, float32 for the 10x and
kb files scanpy and anndata read. sc-droplets sums in float64 and converts (exact for totals below
$$2^{24}$$), then computes in that dtype as cellsweep does (`--matrix-dtype`).

**Manual** (with `--umi-cutoff` or `--expected-cells`, in every release): the cutoff is the given
one, or the `--expected-cells`-th largest total.

**threshold** (the automatic rule before 1.0.2):

1. The totals above 10, sorted in decreasing order: $$c_1 \ge c_2 \ge \dots$$, ranks $$1, 2, \dots$$.
2. $$y = \log_{10} c$$, smoothed with a Gaussian filter of $$\sigma = 2$$ points
   (`scipy.ndimage.gaussian_filter1d`, reflecting at the ends).
3. The second derivative by central differences, twice (`numpy.gradient` against the rank).
4. The cutoff is $$c$$ at the most negative second derivative (the first, if tied).

**kneed** (the automatic rule from 1.0.2):

1. The distinct totals above 10, in decreasing order; each is placed at the last rank holding it
   (the cumulative number of barcodes with at least that total).
2. kneed's `KneeLocator` on (rank, total), convex and decreasing, `online=True`: the
   [Kneedle](../kneedle/kneedle.md) search that keeps the last knee it finds, $$S = 1$$, in the
   behaviour of the installed kneed (`--kneed-version`: 0.8.5, 0.8.6 or 0.9.0 give different
   knees here, since the online search changed in 0.8.6 and 0.9.0).
3. The cutoff is the total at the knee (`knee_y`). No knee, or a single distinct total: cellsweep
   stops with an error, and so does this.

**mx_filter** (any release):

1. All totals, sorted in decreasing order (zeros included); $$v = \log(1 + c)$$.
2. scikit-learn's `GaussianMixture(n_components=2, random_state=42)` fitted to $$v$$.
3. The cutoff is the total of the barcode whose two posterior probabilities have the highest
   entropy (the first, if tied): the barcode the mixture is least sure about.

The cells are the barcodes with total $$\ge$$ the cutoff.

## Versions

| cellsweep | Variant | `--tool-version` | Behaviour |
|---|---|---|---|
| 0.1.0 – 1.0.1 | `v0.1` | `0.1.0` | default `threshold`; `mx_filter` on request |
| 1.0.2 | `v1.0.2` | `1.0.2` | default `kneed`; `mx_filter` on request |

## Parameters

| Option | Default | Original | Meaning |
|---|---|---|---|
| `--empty-method` | `threshold` (0.1.0 – 1.0.1), `kneed` (1.0.2) | `empty_droplet_method` | the automatic rule: `threshold`, `kneed` or `mx_filter` (each where its release has it) |
| `--umi-cutoff` | none | `umi_cutoff` | the cutoff itself (wins over everything) |
| `--expected-cells` | none | `expected_cells` | the cutoff is the total of the barcode ranked this |
| `--kneed-version` | `0.9.0` | the installed `kneed` | the kneed whose behaviour `kneed` follows: `0.8.5`, `0.8.6` or `0.9.0` |
| `--matrix-dtype` | `float32` | the dtype of `adata.X` | the dtype cellsweep sums and computes in |

## Run it

```bash
sc-droplets call cellsweep-empty raw_feature_bc_matrix.h5 --json
sc-droplets call cellsweep-empty raw_feature_bc_matrix.h5 --tool-version 1.0.1 --barcodes cells.tsv
sc-droplets call cellsweep-empty raw_feature_bc_matrix.h5 --empty-method mx_filter
```

```python
import sc_droplets as scd

adata = scd.read("raw_feature_bc_matrix.h5")
result = scd.call(adata, "cellsweep-empty", kneed_version="0.8.5")
result.thresholds["umis"]
```

The thresholds hold the cutoff (`umis`). `sc-droplets plot diagnostics` draws the rank curve with
the cut, and the rule's own view: the smoothed curve and its second derivative, the Kneedle
difference curve, or the mixture's components and entropy.

## Differences from the original

- cellsweep sums every feature in the AnnData it is given; sc-droplets counts the selected
  feature types (Gene Expression by default; `--feature-type` for others).
- The automatic rules use numpy's `log10` and `log1p` in float32 and scikit-learn's mixture, as
  cellsweep does: with the same numpy and scikit-learn the cut is the same; others can move it.

## Validation

Compared with cellsweep's own `infer_empty_droplets` (1.0.1, and 1.0.2 with kneed 0.8.5, 0.8.6
and 0.9.0, each in its own environment) on float32 AnnData of the per-barcode UMI totals of 10x's
pbmc_1k_v3 and pbmc4k raw matrices: the default rule, `mx_filter` and `expected_cells`. Identical
cutoffs.

## References

- Caskey, Rich, Weber, Mortazavi, Pachter & Hallgrímsdóttir 2026, Single-cell genomics
  decontamination with CellSweep, *bioRxiv*,
  [doi:10.64898/2026.03.04.709349](https://doi.org/10.64898/2026.03.04.709349).
- Source: [`cellsweep/utils/data_utils.py`](https://github.com/pachterlab/cellsweep/blob/c81ab360/cellsweep/utils/data_utils.py)
  (`infer_empty_droplets`, `get_umi_cutoff_from_kneed`, `gmm`; 1.0.2) and
  [1.0.1](https://github.com/pachterlab/cellsweep/blob/ed8083a8/cellsweep/utils/data_utils.py)
  (`automatic_umi_cutoff_detection`).
