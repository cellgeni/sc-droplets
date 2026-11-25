---
title: Outputs
parent: Reference
nav_order: 4
---

# Outputs
{: .no_toc }

1. TOC
{:toc}

## In AnnData

A call stored under the key `k` (by `scd.call`, or `call -o` / `compare -o`) writes:

| Where | What |
|---|---|
| `obs[k]` | bool: the barcode is a cell |
| `obs[k + "_class"]` | categorical, for methods with more than two states (for example cell / transition / empty) |
| `obs[k + "_" + name]` | per-barcode values the method computes (p-values, FDR, scores) |
| `uns["sc_droplets"][k]` | the record below |

The record holds:

| Key | Meaning |
|---|---|
| `method`, `title` | the method id and name |
| `tool_version` | the release asked for (absent when none was) |
| `variant` | the behaviour it resolved to |
| `params` | every parameter, resolved |
| `n_barcodes`, `n_cells`, `min_cell_total` | counts, and the smallest total among the cells |
| `thresholds` | the named cut-offs the method used |
| `counting` | what was counted: unit, layer, feature types, genome, number of features, read column, seed |
| `notes` | warnings from the run |
| `columns` | the obs columns the call wrote (beside `obs[k]`) |
| `diagnostics` | small arrays the plots use |
| `package_version` | the sc-droplets version |

Storing a call again under the same key replaces the old one, including its columns. Everything
in the record is plain (strings, numbers, arrays and dicts), so it survives `write_h5ad`, and a
plot can be redrawn from the file.

The default key is the method id with `-` and `.` turned into `_`, followed by the tool version
and any parameters given in the spec: `top-n:n_cells=3000` is stored as `top_n_n_cells_3000`,
`cellranger-emptydrops@9.0` as `cellranger_emptydrops_9_0`.

## Files

| Option | Format |
|---|---|
| `--barcodes cells.tsv.gz` | one cell barcode per line, in input order; gzipped when the name ends in `.gz` |
| `--table calls.csv` | `call`: barcode, total, `is_cell`, `class` and every per-barcode value; `compare`: barcode, total and one column per call. TSV when the name ends in `.tsv` or `.txt`; gzipped for `.gz` |
| `--json` | `call`: the record above without diagnostics; `compare`: `calls` (one entry per key) and `jaccard` |
| `-o calls.h5ad` | the input with the calls added (`--cells-only` keeps the cells) |

Every file is written beside its final name and renamed into place once complete, so an
interrupted run leaves no partial file.
