---
title: Cell Ranger
parent: Methods
nav_order: 2
has_children: true
---

# Cell Ranger
{: .no_toc }

10x Genomics' pipeline has called cells in at least six different ways. The two parts have stayed
recognisable: **OrdMag**, a cut at a tenth of a robust maximum, and, from 3.0, an
**EmptyDrops-style test** of the barcodes below it against the ambient profile. Nearly everything
around them has changed: rounding, the random numbers, the ambient ranks, the candidates' minimum,
the FDR, the number of simulations, and how expected cells are found when not given.

1. TOC
{:toc}

## Which release does what

| Cell Ranger | Default calling | What changed | sc-droplets |
|---|---|---|---|
| 1.0 – 1.3 | OrdMag | the paper's rule; no source published | gap |
| 2.0 – 2.2 | OrdMag | 100 bootstraps; per genome and GEM group | `v2` |
| 3.0 – 3.1 | OrdMag + EmptyDrops | ambient ranks 45,000–90,000, top 20,000 candidates, 10,000 simulations, FDR 0.01 | `v3` |
| 4.0 – 4.x | OrdMag + EmptyDrops | the simulation seeds its own stream; targeted calling for hybrid-capture panels | `v4` (from 10x's outputs) |
| 5.0 – 5.x | | antibody aggregates removed by their counts first | `v5` (from 10x's outputs) |
| 6.0 – 7.1 | | LT ambient range (by 6.0), species fix (6.0), per-species EmptyDrops (6.1), expected cells estimated (7.0) | gap: no source, and no rule matches 10x's outputs exactly |
| 7.2 | OrdMag (cells estimated) + EmptyDrops | chemistry-dependent ambient ranges; antibody aggregates removed first; high-occupancy GEM filter for multiplexed Flex | `v7.2` |
| 8.0 | | candidates must lie above the ambient range; GEM-X ranges and FDR 0.001; global-minimum and mito filters | `v8.0`, `v8.0.1` |
| 9.0 | | 100,000 simulations, new simulator; FDR 0.001 for more chemistries | `v9.0` |
| 10.0 – 10.1 | | Flex v2 ranges; features sorted by id in the simulation; LT range removed (10.1) | `v10.0`, `v10.1` |

## Methods

| Page | Id | What it reproduces |
|---|---|---|
| [Cell Ranger cell calling](cellranger.md) | `cellranger` | `cellranger count` / `multi` as the release ran it: OrdMag, then (3.0+) the EmptyDrops step; `--force-cells` |
| [OrdMag](ordmag.md) | `cellranger-ordmag` | the OrdMag step alone |

## Reproducing a run

Cell Ranger's calls depend on the options of the run, which the output files do not record:

- **`--expect-cells`**: changes OrdMag's robust maximum. 10x's public datasets were run with it
  (1,000 for the "1k" datasets, 5,000 for pbmc4k; from 4.0 the dataset pages say which); without
  it the calls differ by a few dozen barcodes.
- **`--force-cells`**: replaces calling by the top N.
- **The chemistry** (7.2+): sets the ambient ranks and the FDR.

Cell Ranger calls cells on Gene Expression features only, per genome for OrdMag, and per GEM group
(the number after the barcode's `-`). From 7.2 the chemistry matters too; Cell Ranger writes it
into its HDF5 files, and `sc-droplets` reads it from there.

Targeted (gradient) calling, which 7.2–9.0 used for hybrid-capture panels, is not implemented:
no public dataset made with those releases runs it.

## Source and licence

Cell Ranger 3.0.2's source was released by 10x under the MIT licence; the 2.x and 3.x methods here
are ports of it (its OrdMag is the 2.x algorithm). Cell Ranger 2.0–2.1 was released under the
AGPL, and 7.2 onwards under a 10x licence that does not allow copying: those releases are
implemented from the specification on these pages, not from their code. 1.x and 4.0–7.1 published
no source: 4.0 and 5.x are the 3.x port with the changes 10x's own outputs of those releases show
(see [Cell Ranger cell calling](cellranger.md#validation)); for 6.0–7.1 no rule reproduces them
exactly, and the same page records how close each comes.
