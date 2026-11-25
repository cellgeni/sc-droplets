---
title: Getting started
nav_order: 2
---

# Getting started
{: .no_toc }

1. TOC
{:toc}

## Install

```bash
pip install sc-droplets
sc-droplets --version
```

## What to give it

Cell calling needs the **raw** matrix: every barcode, including the empty droplets. A filtered
matrix holds only the cells some tool already called, and the knee methods have nothing to work
with; `sc-droplets` warns when an input looks filtered.

Any of these work as input (details in [Inputs](reference/inputs.md)):

| Input | Example |
|---|---|
| A Cell Ranger output folder | `outs/` (count or multi) |
| A 10x HDF5 matrix | `raw_feature_bc_matrix.h5`, `raw_gene_bc_matrices_h5.h5` |
| An MTX folder | `raw_feature_bc_matrix/`, STARsolo `Solo.out/Gene/raw/` |
| A STARsolo output | `Solo.out/`, with `--solo-feature` to pick Gene, GeneFull or Velocyto |
| kb or alevin-fry output | `counts_unfiltered/`, `af_quant/alevin/` |
| AnnData | `.h5ad` or `.zarr`, raw counts in `X` or a layer |

## Look at the sample first

```bash
sc-droplets inspect raw_feature_bc_matrix.h5
```

`inspect` prints the number of barcodes and features, the counts at a few ranks and whether
spliced/unspliced counts are present, then draws the knee in the terminal.

## Call cells

```bash
sc-droplets methods
sc-droplets methods top-n
sc-droplets call top-n raw_feature_bc_matrix.h5 --n-cells 3000 -o calls.h5ad --barcodes cells.tsv.gz
```

Each method is a subcommand of `call` with its own options; `sc-droplets call top-n --help` lists
them. A call can be kept in any of these forms:

| Option | Writes |
|---|---|
| `-o calls.h5ad` | the input with the call added (`--cells-only` keeps only the cells) |
| `--barcodes cells.tsv.gz` | the cell barcodes, one per line |
| `--table calls.csv` | one row per barcode: total, call, and the method's scores |
| `--json` | the call's summary on stdout: parameters, thresholds, number of cells |

## Compare methods

```bash
sc-droplets compare raw_feature_bc_matrix.h5 -m top-n:n_cells=3000 -m threshold:min_counts=500 --table calls.csv --plot knee.png
```

Each `-m` is a method spec, `id[@tool_version][:name=value,...]`. `compare` prints how many cells
each call found and writes one column per call.

## Plot

```bash
sc-droplets plot knee raw_feature_bc_matrix.h5 -m threshold:min_counts=500 -o knee.png
sc-droplets plot knee calls.h5ad -o knee.pdf --rotated
```

Given an `.h5ad` written by `call -o` or `compare -o`, `plot knee` draws the calls stored in it.

## From Python

```python
import sc_droplets as scd

adata = scd.read("raw_feature_bc_matrix.h5")
result = scd.call(adata, "top-n", n_cells=3000)
result.n_cells, result.thresholds
adata.obs["top_n"]

scd.compare(adata, ["top-n:n_cells=3000", "threshold:min_counts=500"])
scd.agreement(adata)
scd.pl.knee(adata)
```

`scd.call` stores each call in `adata.obs` and `adata.uns["sc_droplets"]`; [Outputs](reference/outputs.md)
describes the layout. To try things without data, `scd.simulate_droplets()` makes a raw matrix
with known cells.
