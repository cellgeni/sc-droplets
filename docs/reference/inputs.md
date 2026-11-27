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
| velocyto `.loom` | `spliced`, `unspliced` and `ambiguous` layers (velocyto's genes × cells, read a block of genes at a time); counts are their sum; barcodes as velocyto writes them (`sample:BARCODEx`) |

Counts stay integers, as the pipeline wrote them.

## Matrix Market folders

| Layout | Files | Notes |
|---|---|---|
| Cell Ranger 3+ | `matrix.mtx.gz`, `barcodes.tsv.gz`, `features.tsv.gz` | |
| Cell Ranger 2 | `matrix.mtx`, `barcodes.tsv`, `genes.tsv` | |
| STARsolo | `matrix.mtx`, `barcodes.tsv`, `features.tsv` | `Solo.out/<feature>/raw/` |
| STARsolo Velocyto | `spliced.mtx`, `unspliced.mtx`, `ambiguous.mtx` | counts are their sum; each is also a layer |
| kb | `cells_x_genes.mtx`, `.barcodes.txt`, `.genes.txt` (`.genes.names.txt` for symbols) | barcodes are rows |
| kb `--workflow nac` | `cells_x_genes.mature.mtx`, `.nascent.mtx`, `.ambiguous.mtx` | counts are their sum; mature + ambiguous become `spliced`, nascent `unspliced` |
| alevin-fry | `quants_mat.mtx`, `quants_mat_rows.txt`, `quants_mat_cols.txt` | barcodes are rows; in USA mode the three column blocks (S, U, A) are summed, S + A become `spliced`, U `unspliced` |

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

## Splicing

The unspliced fraction U / (S + U) of each barcode (DropletQC's nuclear fraction) comes from
spliced and unspliced counts, whichever of these holds them:

| Source | Spliced | Unspliced |
|---|---|---|
| the input's own `spliced` and `unspliced` layers (an `.h5ad`, or any layout below read as the input) | `spliced` | `unspliced` |
| velocyto `.loom` | `spliced` | `unspliced` |
| STARsolo `Solo.out/Velocyto/raw/` | `spliced.mtx` | `unspliced.mtx` |
| STARsolo `Solo.out/` with `Gene/raw` and `GeneFull/raw` but no Velocyto | Gene | GeneFull − Gene: the reads on introns, an estimate |
| kb `--workflow nac` | `mature` + `ambiguous` | `nascent` |
| alevin-fry USA mode | S + A blocks | U block |

Reads within exons count as spliced, whichever tool counted them. velocyto and STARsolo count them
as spliced themselves; their ambiguous molecules (exonic in one transcript model, intronic in
another) are left out of the fraction. kb's nac `ambiguous` counts and alevin-fry's A block are the
reads that fit both the spliced and the unspliced sequence, mostly reads within exons: they are
added to the spliced counts, as kb's `--sum cell` and pyroe's velocity output add them. Left out,
they would double the fraction: on 10x's PBMC 1k, kb's counts give the cells 0.63 and the
barcodes left out 0.50 without them, 0.32 and 0.14 with them. Only gene-expression features count.

`--splicing SRC` (and `splicing=` in `scd.add_metrics`) reads another output for the same sample
and matches it by barcode: by name, or else by the barcode sequence alone, so velocyto's
`sample:ACGT...x` matches Cell Ranger's `ACGT...-1`. Barcodes the source lacks get no value. When two barcodes share a sequence
(two samples in one loom, two GEM groups in one matrix), matching by sequence would be a guess,
and it is refused: give a source whose barcode names match the input's. Cell
Ranger's own outputs have no spliced/unspliced matrix (the intron flag lives in the BAM's `RE`
tag); run velocyto or STARsolo with `--soloFeatures Gene GeneFull Velocyto` on the same reads.

The mitochondrial fraction needs no extra input: the share of each barcode's counts in genes named
`MT-` (human) or `mt-` (mouse).

## Read counts

Some methods rank barcodes by reads rather than UMIs (UMI-tools, alevin, alevin-fry, bustools).
Give them with `--reads`: the name of an obs column of an `.h5ad`, or a CSV/TSV with a barcode and
a count on each line (a header line is skipped). In Python, pass `reads=` an obs column.
