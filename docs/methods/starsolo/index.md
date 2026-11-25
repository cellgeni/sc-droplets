---
title: STARsolo
parent: Methods
nav_order: 3
has_children: true
---

# STARsolo
{: .no_toc }

STAR's single-cell mode filters cells with `--soloCellFilter`, after counting. It has three rules
of its own, plus a script shipped beside it. None is a copy of Cell Ranger's code: the one named
after Cell Ranger 2.2 drops its bootstrap, and EmptyDrops_CR has its own Good-Turing estimate,
random numbers and tie order.

1. TOC
{:toc}

## Which release does what

| STAR | Default filter | What changed |
|---|---|---|
| 2.7.0a – 2.7.2d | none (raw output only) | `extras/scripts/soloBasicCellFilter.awk` added in 2.7.2a |
| 2.7.3a | CellRanger2.2 3000 0.99 10 | `--soloCellFilter` added: CellRanger2.2, TopCells, None |
| 2.7.5a | CellRanger2.2 | the cut is floored at 1 (before, a cut of 0 kept every barcode) |
| 2.7.6a | CellRanger2.2 | TopCells fixed (it used to keep every barcode) |
| 2.7.8a | CellRanger2.2 | EmptyDrops_CR added; CellRanger2.2's parameters become optional; rounding with C `round` |
| 2.7.10a – 2.7.11b | CellRanger2.2 | `--runMode soloCellFiltering` usable (it crashed in 2.7.9a); EmptyDrops_CR unchanged |

STAR reports "Estimated Number of Cells" in `Summary.csv` as the size of the filtered set; there
is no separate knee anywhere in STAR. Velocyto counts reuse the Gene feature's cells.

## Methods

| Page | Rule |
|---|---|
| [CellRanger2.2](cellranger2-2.md) | the 31st-largest total ÷ 10 |
| [TopCells](topcells.md) | the top N + 1 barcodes, plus ties |
| [EmptyDrops_CR](emptydrops-cr.md) | CellRanger2.2, then STAR's test of each droplet against the ambient profile |
| [soloBasicCellFilter.awk](basic-awk.md) | the awk script: the same cut, unrounded |

## How STAR sees a matrix

- Only **detected** barcodes count: those with at least one count, in the order of the barcode
  list. STAR calls their number nCB.
- When STAR loads a matrix (`--runMode soloCellFiltering`), it rounds each entry to the nearest
  integer, halves away from zero.
- Totals are sorted in decreasing order; ranks are 0-based.

sc-droplets does the same, so a raw matrix gives the same calls whether STAR made it or not.
