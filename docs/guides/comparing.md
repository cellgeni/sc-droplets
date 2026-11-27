---
title: Comparing methods on a sample
parent: Guides
nav_order: 3
---

# Comparing methods on a sample
{: .no_toc }

When no pipeline's call is the reference, run several methods and look at where they disagree.
The barcodes every method calls are cells. The barcodes no method calls are empty. What decides
the sample is the band in between.

1. TOC
{:toc}

## Look at the sample

```bash
sc-droplets inspect raw_feature_bc_matrix.h5
```

`inspect` draws the knee in the terminal, says whether the matrix looks filtered (cell calling
needs the raw one), and whether spliced and unspliced counts are present.

## Run one method from each family

Methods of one [family](../overview/families.md) tend to agree with each other, so a panel of one
per family shows more than many methods of one kind:

```bash
sc-droplets compare raw_feature_bc_matrix.h5 -m cellranger -m dropletutils-emptydrops -m umitools-distance -m otsu -m kde-valley -m gmm -o calls.h5ad --table calls.csv --json
```

| Family | A representative |
|---|---|
| test against the ambient profile | `cellranger`, `dropletutils-emptydrops` |
| farthest point from a chord | `umitools-distance` |
| distribution of log counts | `otsu`, `kde-valley` |
| mixture | `gmm` |

`compare` prints each call's number of cells and smallest cell total. `--json` adds the Jaccard
index of every pair of cell sets. Method specs take parameters and releases:
`-m cellranger@3.0.0:expect_cells=1000`, `-m dropletutils-emptydrops:niters=100000`.

## Add the unspliced and mitochondrial fractions

Two per-barcode measures help judge the disputed barcodes:

- the **unspliced fraction**, from spliced and unspliced counts (STARsolo Velocyto, velocyto,
  kb nac, alevin-fry USA; [Inputs](../reference/inputs.md#splicing));
- the **mitochondrial fraction**, when genes are named `MT-` or `mt-`.

```bash
sc-droplets compare raw_feature_bc_matrix.h5 -m cellranger -m otsu -m gmm --splicing Solo.out --table calls.csv --json
```

With them, `compare` scores each call. The score is the median of each measure in the call's
cells and in the barcodes it left out (those with at least 10 gene-expression counts), and the
AUROC between the two. These are indicators, not truth:

- In whole cells, empty droplets hold cytoplasmic, spliced RNA, so cells have the higher unspliced
  fraction.
- Barcodes made mostly of mitochondrial RNA have a low unspliced fraction, because mitochondrial
  transcripts have no introns. A call that takes them in raises its mitochondrial score and lowers
  its unspliced one. On PBMC 1k the barcodes the deeper calls add are mostly of this kind
  ([comparison](../overview/comparison.md#what-it-shows)).
- In nuclei the mitochondrial score flips: nuclei carry less mitochondrial RNA than the ambient
  pool, so a lower score is the cleaner call.
- A call that leaves out almost nothing above the floor has a score that means little.

## Look at the disputed barcodes

```bash
sc-droplets plot compare calls.h5ad -o compare.png
sc-droplets plot qc calls.h5ad -o qc.png
sc-droplets plot splicing calls.h5ad --splicing Solo.out -o splicing.png
sc-droplets report calls.h5ad -o report.html
```

- `plot compare` shows where along the rank curve the calls part, and their Jaccard indices.
- `plot qc` sets genes and the mitochondrial fraction against count, and the barcodes only some
  calls took beside those every call took.
- `plot splicing` does the same for the unspliced fraction.
- `report` puts everything on one page to share.

[Visualising calls](visualising.md) shows each plot on PBMC 1k.

## Decide

Read the disputed band as cells or not from what is in it, not from which method took it:

- **Genes per count like the agreed cells, unspliced fraction like theirs:** small cells. A call
  that keeps them is right for this sample.
- **Mostly mitochondrial RNA, few genes:** damaged cells or debris. A call that leaves them out
  saves a filter later; one that keeps them needs it.
- **Like the empty droplets:** a call that takes them in is too deep for this sample.

Then store the chosen call (`--key`) and keep its record: the method, the release and the
parameters are in `uns["sc_droplets"]`, so the choice can be repeated.
