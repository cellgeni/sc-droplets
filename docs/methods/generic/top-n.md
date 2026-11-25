---
title: Top N barcodes
parent: Generic rules
grand_parent: Methods
nav_order: 1
method_ids: [top-n]
status: implemented
---

# Top N barcodes
{: .no_toc }

The N barcodes with the most counts are cells. It is the rule to use when the number of cells is
known, and the one most tools fall back on when told it (Cell Ranger `--force-cells`, STARsolo
`TopCells`, UMI-tools `--set-cell-number`, alevin `--forceCells`).

1. TOC
{:toc}

## At a glance

| Id | `top-n` |
|---|---|
| Tool | none (generic) |
| Category | rank-threshold |
| Counts | UMIs, reads or genes detected |
| Needs expected cells | yes: N is the input |
| Deterministic | yes |
| Validated against | hand-computed cases (ties, N = 0, N above the number of barcodes) |

## How it works

1. Rank the barcodes by their total, largest first.
2. Call the first N.

The only choice is what to do with barcodes that tie with the N-th:

- `--ties first` (default) takes exactly N, keeping the earlier barcode of a tie in input order;
- `--ties include` takes every barcode whose total is at least the N-th total, so it can call
  more than N.

The tools that offer this rule differ exactly here; their pages say which they do.

## Versions

None: the rule does not follow the releases of a tool.

## Parameters

| Option | Default | Original | Meaning |
|---|---|---|---|
| `--n-cells` | required | Cell Ranger `--force-cells` | how many barcodes to call |
| `--ties` | `first` | | `first` or `include` |
| `--unit` | `umis` | | rank by `umis`, `reads` (needs `--reads`) or `genes` |

## Run it

```bash
sc-droplets call top-n raw_feature_bc_matrix.h5 --n-cells 3000 --barcodes cells.tsv.gz
sc-droplets call top-n raw_feature_bc_matrix.h5 --n-cells 3000 --ties include --json
```

```python
import sc_droplets as scd

adata = scd.read("raw_feature_bc_matrix.h5")
scd.call(adata, "top-n", n_cells=3000)
```

The call's threshold is the N-th total, under the unit's name (`thresholds["umis"]`).

## Differences from the original

Not applicable. Each tool's own top-N rule, with its tie handling and off-by-one quirks, will
be implemented on that tool's page.

## Validation

Unit tests check exact calls on small hand-made matrices: ties broken by input order, ties kept
with `include`, N = 0, and N larger than the number of barcodes.

## References

- Cell Ranger `--force-cells`, STARsolo `--soloCellFilter TopCells`, UMI-tools
  `--set-cell-number`, alevin `--forceCells`, alevin-fry `--force-cells`.
