---
title: Inputs
parent: Reference
nav_order: 3
---

# Inputs
{: .no_toc }

Cell calling needs the **raw** matrix, with every barcode the pipeline saw. A folder is searched
for the raw matrix inside it; a folder with only a filtered matrix is refused, with that reason.

1. TOC
{:toc}

## Files

| Input | Read as |
|---|---|
| `.h5ad` | AnnData, as written |
| `.zarr` | AnnData, as written |
| 10x HDF5, v3 (`raw_feature_bc_matrix.h5`, Cell Ranger 3+, also CellBender's output) | counts in `X`; `var` has `gene_ids`, `feature_types`, `genome`, indexed by gene symbol made unique |
| 10x HDF5, v2 (`raw_gene_bc_matrices_h5.h5`, Cell Ranger 2) | every genome group side by side, `var["genome"]` saying which; `--genome` reads one |

Counts stay integers, as the pipeline wrote them.

## Matrix Market folders

| Layout | Files | Notes |
|---|---|---|
| Cell Ranger 3+ | `matrix.mtx.gz`, `barcodes.tsv.gz`, `features.tsv.gz` | |
| Cell Ranger 2 | `matrix.mtx`, `barcodes.tsv`, `genes.tsv` | |
| STARsolo | `matrix.mtx`, `barcodes.tsv`, `features.tsv` | `Solo.out/<feature>/raw/` |
| STARsolo Velocyto | `spliced.mtx`, `unspliced.mtx`, `ambiguous.mtx` | counts are their sum; each is also a layer |
| kb | `cells_x_genes.mtx`, `.barcodes.txt`, `.genes.txt` (`.genes.names.txt` for symbols) | barcodes are rows |
| kb `--workflow nac` | `cells_x_genes.mature.mtx`, `.nascent.mtx`, `.ambiguous.mtx` | counts are their sum; mature and nascent become `spliced` and `unspliced` |
| alevin-fry | `quants_mat.mtx`, `quants_mat_rows.txt`, `quants_mat_cols.txt` | barcodes are rows; in USA mode the three column blocks (S, U, A) are summed and kept as layers |

## Pipeline output folders

Given a folder, the raw matrix is looked for in this order:

1. the folder itself, if it is an MTX folder;
2. STARsolo: `<folder>/<feature>/raw`, `<folder>/Solo.out/<feature>/raw`, or `<folder>/raw` when
   the folder is the feature folder (`--solo-feature`, default `Gene`);
3. Cell Ranger: `raw_feature_bc_matrix.h5`, `raw_gene_bc_matrices_h5.h5`,
   `raw_feature_bc_matrix/`, the same under `outs/`, and `outs/multi/count/raw_feature_bc_matrix.h5`
   for `cellranger multi`;
4. kb: `counts_unfiltered/`;
5. alevin-fry and simpleaf: `alevin/`, `af_quant/alevin/`.

## Read counts

Some methods rank barcodes by reads rather than UMIs (UMI-tools, alevin, alevin-fry, bustools).
Give them with `--reads`: the name of an obs column of an `.h5ad`, or a CSV/TSV with a barcode and
a count on each line (a header line is skipped). In Python, pass `reads=` an obs column.
