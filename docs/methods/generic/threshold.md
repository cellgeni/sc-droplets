---
title: Fixed threshold
parent: Generic rules
grand_parent: Methods
nav_order: 2
method_ids: [threshold]
status: implemented
---

# Fixed threshold
{: .no_toc }

Every barcode with at least a given number of UMIs, reads or genes detected is a cell. Simple
and transparent; the number has to come from somewhere else (a look at the knee plot, a
protocol's typical yield, or another method).

1. TOC
{:toc}

## At a glance

| | |
|---|---|
| Id | `threshold` |
| Tool | none (generic) |
| Category | rank-threshold |
| Counts | UMIs, reads or genes detected |
| Needs expected cells | no, but needs the threshold |
| Deterministic | yes |
| Validated against | hand-computed cases |

## How it works

A barcode is a cell when its value is at least the threshold, or, with `--strict`, more than it:

$$\text{cell}_b = \left[\, v_b \ge t \,\right] \quad\text{or}\quad \left[\, v_b > t \,\right]$$

where \(v_b\) is the barcode's total UMIs, reads, or number of genes with a nonzero count.

Tools that use a fixed threshold:

| Tool | Setting | Unit |
|---|---|---|
| Drop-seq tools | `SelectCellsByNumTranscripts MIN_TRANSCRIPTS_PER_CELL` | UMIs |
| alevin-fry | `--unfiltered-pl` with `--min-reads 10` | reads |
| zUMIs | `nReadsperCell` (100) before its knee | reads |
| scanpy and Seurat tutorials | `min_genes=200` | genes |

## Versions

None: the rule does not follow the releases of a tool.

## Parameters

| Option | Default | Original | Meaning |
|---|---|---|---|
| `--min-counts` | required | Drop-seq `MIN_TRANSCRIPTS_PER_CELL` | the threshold |
| `--unit` | `umis` | | `umis`, `reads` (needs `--reads`) or `genes` |
| `--strict` / `--no-strict` | `--no-strict` | | require more than the threshold |

## Run it

```bash
sc-droplets call threshold raw_feature_bc_matrix.h5 --min-counts 500 --barcodes cells.tsv.gz
sc-droplets call threshold raw_feature_bc_matrix.h5 --min-counts 200 --unit genes -o calls.h5ad
```

```python
import sc_droplets as scd

adata = scd.read("raw_feature_bc_matrix.h5")
scd.call(adata, "threshold", min_counts=200, unit="genes")
```

## Differences from the original

Not applicable.

## Validation

Unit tests check exact calls on small hand-made matrices, for each unit and with and without
`--strict`.

## References

- Macosko et al. 2015, Drop-seq, *Cell* 161:1202,
  [doi:10.1016/j.cell.2015.05.002](https://doi.org/10.1016/j.cell.2015.05.002).
