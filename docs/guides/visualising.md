---
title: Visualising calls
parent: Guides
nav_order: 4
---

# Visualising calls
{: .no_toc }

Every plot takes an AnnData of raw counts and draws the calls stored on it, or runs the method
specs it is given first. The Python functions return the matplotlib Axes or Figure; each has a
command that writes a `.png`, `.pdf` or `.svg`. The figures below are 10x's pbmc_1k_v3 raw
matrix with five calls (Cell Ranger 3.0 as 10x ran it, DropletUtils' emptyDrops, Otsu, a
Gaussian mixture and zUMIs), drawn by `docs/_scripts/figures.py`.

1. TOC
{:toc}

## The calls to draw

```python
import sc_droplets as scd

adata = scd.read("raw_feature_bc_matrix.h5")
for spec in ["cellranger", "dropletutils-emptydrops", "otsu", "gmm", "zumis-cell-barcodes"]:
    scd.call(adata, spec)
adata.write_h5ad("calls.h5ad")
```

```bash
sc-droplets compare raw_feature_bc_matrix.h5 -m cellranger -m dropletutils-emptydrops -m otsu -m gmm -o calls.h5ad
```

A plot given `calls.h5ad` with neither `-m` nor `--key` draws every stored call; `--key` picks
some, `-m` runs more.

## The knee

The barcode-rank curve: each barcode's total against its rank, both on log scales. Cells make the
high plateau, empty droplets the long tail; each call is a dashed line at its smallest cell total
with a dot where its cells end.

```bash
sc-droplets plot knee calls.h5ad -o knee.png
sc-droplets plot knee calls.h5ad --key cellranger --rotated -o knee-rotated.png
```

![The knee with five calls]({{ '/assets/guide/knee.png' | relative_url }})

`--rotated` puts totals on x and ranks on y; `--cumulative` plots the share of all counts against
rank, the view Drop-seq and UMI-tools' distance method work on.

![The rotated knee]({{ '/assets/guide/knee-rotated.png' | relative_url }})

### Colour and tracks

`--color` colours each barcode by a measure, and `--track` adds a panel under the curve with the
measure's median in log-spaced bins of rank. With the unspliced fraction (see
[Inputs](../reference/inputs.md#splicing)) the change in composition at the knee shows: cells hold
more unspliced RNA than the cytoplasmic ambient RNA of empty droplets. On PBMC 1k (the figure),
just below the knee the mitochondrial fraction peaks and the unspliced fraction drops to near 0.
Those barcodes are mostly mitochondrial RNA, which has no introns.

```bash
sc-droplets plot knee raw_feature_bc_matrix.h5 -m cellranger --splicing counts_unfiltered --color unspliced_fraction --track unspliced_fraction --track mito_fraction -o knee-tracks.png
```

![The knee coloured by the unspliced fraction, with tracks]({{ '/assets/guide/knee-tracks.png' | relative_url }})

## The density of totals

The histogram of log10 totals on a log axis, with each call's cut. Empty droplets make the left
mass, cells the right peak; rules such as Otsu, the kernel-density valley and the mixtures cut
between them. `--weighted` weights each barcode by its total, so the cells' peak stands out.

```bash
sc-droplets plot density calls.h5ad -o density.png
```

![The density of log10 totals with five calls]({{ '/assets/guide/density.png' | relative_url }})

## Several calls side by side

The curve with a strip per call under it, shaded by the share of each rank bin the call took, so
where the calls part is plain; the cells per call; and the Jaccard index of their cell sets.

```bash
sc-droplets plot compare calls.h5ad -o compare.png
```

![Five calls side by side]({{ '/assets/guide/compare.png' | relative_url }})

## Quality measures

Genes detected and the mitochondrial fraction against count, each barcode coloured by how many
calls took it; then the barcodes the calls disagree on (some took them) beside those every call
took and those none took in the band of totals where the calls disagree.

```bash
sc-droplets plot qc calls.h5ad -o qc.png
```

![Genes, mitochondrial fraction and the disputed barcodes]({{ '/assets/guide/qc.png' | relative_url }})

## The unspliced fraction

The unspliced fraction against count by agreement; in each call's cells and the barcodes it left
out, with the AUROC between them; and where the calls disagree. Whole-cell samples: empty droplets
sit low, cells higher, damaged cells high. Nuclei samples: the nuclei sit high. On PBMC 1k the
cells sit near 0.31 and the barcodes no call took near 0.14. Below the cells, the barcodes only
some calls took spread out: some near 0, mostly mitochondrial RNA, and some above 0.4.

```bash
sc-droplets plot splicing calls.h5ad --splicing counts_unfiltered -o splicing.png
```

![The unspliced fraction against five calls]({{ '/assets/guide/splicing.png' | relative_url }})

`sc-droplets compare` with `--splicing` adds the same scores to its table: each call's median
fraction in its cells and in the barcodes it left out, and the AUROC.

## One call's own view

Each method draws what it computed: the test methods their classes, likelihoods and p-values;
the distribution rules their histogram or density with the cut; the mixtures their components.

```bash
sc-droplets plot diagnostics calls.h5ad --key emptydrops -o diagnostics.png
sc-droplets plot diagnostics raw_feature_bc_matrix.h5 -m gmm -o diagnostics.png
```

![emptyDrops: classes, likelihoods and p-values]({{ '/assets/guide/diagnostics-emptydrops.png' | relative_url }})

![The Gaussian mixture's components]({{ '/assets/guide/diagnostics-gmm.png' | relative_url }})

## The report

One self-contained HTML page with every plot above, the calls and their separation scores, and
each call's parameters, notes and diagnostic view. It loads nothing from elsewhere, and names the
input file but not its path, so it can be mailed or archived.

```bash
sc-droplets report calls.h5ad --splicing counts_unfiltered -o report.html --title "PBMC 1k"
```

```python
scd.add_metrics(adata, splicing="counts_unfiltered")
scd.report(adata, path="report.html", title="PBMC 1k")
```

## In Python

Every plot is a function of `scd.pl` taking the same calls; `knee` and `density` return Axes and
accept `ax=` to draw into a figure of your own; `compare`, `qc`, `splicing` and `diagnostics`
return a Figure. See the [Python API](../reference/python-api.md#plots).

```python
import matplotlib.pyplot as plt

fig, (left, right) = plt.subplots(1, 2, figsize=(12, 4.8))
scd.pl.knee(adata, ["cellranger", "otsu"], ax=left)
scd.pl.density(adata, ["cellranger", "otsu"], ax=right)
fig.savefig("knee-and-density.png")
```
