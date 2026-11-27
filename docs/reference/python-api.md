---
title: Python API
parent: Reference
nav_order: 2
---

# Python API
{: .no_toc }

```python
import sc_droplets as scd
```

1. TOC
{:toc}

## Reading

### `scd.read(path, *, solo_feature="Gene", genome=None)`

Any input in [Inputs](inputs.md) as an AnnData, barcodes × features, raw counts in `X`. Spliced
and unspliced counts, where the input has them, go to `layers["spliced"]` and
`layers["unspliced"]`.

## Calling

### `scd.call(adata, method, *, tool_version=None, key_added=None, layer=None, feature_types=None, genome=None, reads=None, random_state=None, inplace=True, **params)`

Runs one method and returns a `CallResult`. `method` is a spec, `id[@tool_version][:name=value,...]`;
`**params` are the method's parameters by name (`n_cells=3000`), and take precedence over those in
the spec.

| Argument | Meaning |
|---|---|
| `tool_version` | Release of the original tool to reproduce; default the latest known |
| `key_added` | Key in `obs` and `uns["sc_droplets"]`; default the method id (plus version and spec parameters) |
| `layer` | Count this layer instead of `X` |
| `feature_types` | A type or list of types from `var["feature_types"]`, or `"all"`; default what the method counts |
| `genome` | Count only this genome's features |
| `reads` | An obs column of read counts per barcode, for methods that rank reads |
| `random_state` | Seed; default the seed the original tool uses |
| `inplace` | Store the call on `adata` (default `True`) |

```python
result = scd.call(adata, "top-n", n_cells=3000)
result = scd.call(adata, "threshold:min_counts=200,unit=genes", key_added="genes200")
```

### `scd.compare(adata, methods, **common)`

Runs every spec in `methods` with the same `common` options (`layer`, `feature_types`, ...) and
returns a table with one row per call: key, method, tool version, variant, number of cells and
smallest cell total. Every call is stored on `adata`. When `obs` holds `unspliced_fraction` or
`mito_fraction` (`add_metrics`), each row also has `separation`'s columns.

```python
scd.compare(adata, ["top-n:n_cells=3000", "threshold:min_counts=500"])
```

### `scd.add_metrics(adata, splicing=None)`

Stores each barcode's measures in `obs` and returns the columns written: `unspliced_fraction`,
U / (S + U) over gene-expression features (from `splicing`, matched by barcode, or from the
AnnData's own `spliced` and `unspliced` layers), and `mito_fraction`, the share of counts in genes
named `MT-` or `mt-` (when there are any). `splicing` is anything `unspliced_fraction` takes, or a
Series of fractions indexed by barcode. Barcodes match by name, or else by the sequence alone
(velocyto's `sample:` prefix and `x` suffix and a `-1` suffix dropped); unmatched barcodes get NaN.

```python
scd.add_metrics(adata, splicing="velocyto/sample.loom")
```

### `scd.unspliced_fraction(source)`

Each barcode's U / (S + U) as a Series indexed by barcode (NaN where it has neither). `source` is
an AnnData with `spliced` and `unspliced` layers, or a path: a velocyto `.loom`, STARsolo's
`Velocyto/raw/` or `Solo.out/` (with `Gene` and `GeneFull` but no Velocyto, the estimate
(GeneFull − Gene) / GeneFull), kb `--workflow nac`, alevin-fry USA mode, an `.h5ad`.

### `scd.separation(adata, keys=None, metric="unspliced_fraction", min_counts=10)`

For each stored call: the median of `obs[metric]` in its cells (`<metric>_cells`) and in the
barcodes it left out with at least `min_counts` (`<metric>_others`), and the AUROC between them
(`<metric>_auroc`: the chance that a cell's value exceeds an outsider's, ties counting half; 0.5
is no separation). Unsupervised indicators, not ground truth.

```python
scd.separation(adata, ["cellranger", "otsu"])
```

### `scd.agreement(adata, keys=None)`

The pairwise Jaccard index of the cell sets of stored calls (all of them by default).

### `scd.stored(adata)`

The keys of the calls stored on an AnnData.

### `scd.list_methods()`

A table of every method: tool, title, category, counts, whether deterministic, variants, status
and page.

## `CallResult`

| Attribute | Meaning |
|---|---|
| `method`, `title`, `variant`, `tool_version` | What ran |
| `params` | Every parameter, resolved |
| `barcodes`, `totals`, `is_cell` | Per barcode, in input order |
| `classes` | Per-barcode categories, for methods with more than two states (else `None`) |
| `columns` | Per-barcode values: p-values, FDR, scores |
| `thresholds` | Named cut-offs the method used |
| `diagnostics` | Small arrays the plots use |
| `notes` | Warnings about the run |
| `n_cells`, `min_cell_total` | Number of cells; smallest total among them |
| `table` | A DataFrame, one row per barcode |
| `summary()`, `to_json()` | The call without per-barcode data |
| `store(adata, key)` | Write the call into `adata` |
| `key` | Where it was stored (`None` if not) |

## Plots

### `scd.pl.knee(adata, calls=None, *, layer=None, feature_types=None, rotated=False, cumulative=False, ax=None, title=None, legend=True)`

The barcode-rank curve with the calls on it; returns the matplotlib Axes. `calls` takes stored
keys or method specs (a spec not yet stored is run first); `None` draws every stored call. With
one call its cells are marked on the curve; with several, each gets a dashed line at its smallest
cell total.

```python
scd.pl.knee(adata, ["top_n", "threshold"])
scd.pl.knee(adata, "top-n:n_cells=3000", rotated=True)
```

### `scd.pl.density(adata, calls=None, *, layer=None, feature_types=None, min_counts=1, bins=100, weighted=False, ax=None, title=None, legend=True)`

The histogram of $$\log_{10}$$ totals (barcodes with at least `min_counts`), on a logarithmic
y axis, with a dashed line at each call's smallest cell total; returns the Axes. `calls` as for
`knee`. `weighted` weights each barcode by its total, so a bar shows the share of all counts its
barcodes hold.

```python
scd.pl.density(adata, ["otsu", "kde-valley"])
scd.pl.density(adata, "mad", weighted=True)
```

### `scd.pl.diagnostics(adata, call, *, fig=None)`

One call's own view, as `plot diagnostics` draws it; returns the Figure. `call` is a stored key or
a method spec.

```python
scd.pl.diagnostics(adata, "starsolo-emptydrops-cr")
```

### `scd.pl.knee_text(totals, *, width=60, height=14, marks=None)`

The knee as lines of Unicode braille, as `inspect` prints it.

## Data for trying things

### `scd.simulate_droplets(n_cells=1500, n_empty=8000, n_low=30000, n_genes=1000, ..., seed=0, splicing=False)`

A raw matrix with known cells: `obs["truth"]` is `cell`, `empty` (ambient RNA only) or `low`
(a handful of counts). With `splicing=True`, `spliced` and `unspliced` layers.

```python
adata = scd.simulate_droplets(seed=1)
```

## Errors

| Exception | Raised when |
|---|---|
| `scd.UnknownMethod` | no method has that id |
| `scd.UnsupportedVersion` | the tool release is in a documented gap, predates the method, or the method follows no tool |
