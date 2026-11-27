---
title: CeleScope
parent: Documented only
grand_parent: Methods
nav_order: 5
method_ids: [celescope]
status: docs-only
---

# CeleScope: STARsolo's EmptyDrops_CR, and rules of its own before
{: .no_toc }

The cell call of Singleron's CeleScope `rna` pipeline. Since 2.0 it is STARsolo 2.7.11a's
[EmptyDrops_CR](../starsolo/emptydrops-cr.md) at an FDR of 0.001, which the package runs. Before,
CeleScope called cells itself: with a port of Cell Ranger 3.0.2's EmptyDrops (the default in
1.9 – 1.17) or with a count threshold of its own (the default in 1.0 – 1.8). Documented, not
implemented as a method of its own.

1. TOC
{:toc}

## At a glance

| Id | `celescope` (documented only) |
|---|---|
| Tool | CeleScope 1.0.1 – 2.15.6 (`celescope rna`, `multi_rna`; PyPI, conda, Docker) |
| Category | statistical-test |
| Counts | UMIs per barcode: STARsolo's `GeneFull_Ex50pAS` matrix from 2.0 (intronic reads counted); CeleScope's own featureCounts-based matrix in 1.x |
| Needs expected cells | uses 3,000 unless told otherwise |
| Deterministic | yes: STARsolo seeds each simulation (2.x); the 1.x port seeds Python's `random` and NumPy with 0 |
| Original licence | MIT |
| Extra inputs | none; the raw matrix CeleScope writes |

## How it works

### 2.0 and later (the default)

The `starsolo` step runs STAR 2.7.11a with
`--soloCellFilter EmptyDrops_CR 3000 0.99 10 45000 90000 500 0.01 20000 0.001 10000`: STARsolo's
EmptyDrops_CR with STAR's defaults except the FDR, 0.001 instead of 0.01. It counts the features
`GeneFull_Ex50pAS` and `Gene`, and reports `GeneFull_Ex50pAS`, whose filtered matrix is the call.
The steps are those of [STARsolo EmptyDrops_CR](../starsolo/emptydrops-cr.md): the CellRanger2.2
cut at 3,000 expected cells; an ambient profile from ranks 45,000 – 89,999; up to 20,000
candidates with at least 500 UMIs, or 1 % of the middle initial cell's; 10,000 simulations;
Benjamini–Hochberg at 0.001. The 2.0.0 changelog also adds a filter of UMIs supported by fewer
reads than `read_threshold`, which changes the counts before calling (not read further).

The optional `celescope rna cells` step, which `multi_rna` does not run, re-filters the raw matrix:
with a list of barcodes (`--barcode`), the top $$N$$ barcodes by UMI (`--force_cells`), or another
`--soloCellFilter` string (through `STAR --runMode soloCellFiltering`); then optionally with a
largest mitochondrial fraction (`--max_mito`) and a fewest genes (`--min_gene`). The original call
is kept as `default_filtered`.

### 1.9 – 1.17: `EmptyDrops_CR` (1.0 – 1.8: `cellranger3`)

A port of Cell Ranger 3.0.2's calling (the 3.x rule on the [Cell Ranger](../cellranger/cellranger.md)
page), on CeleScope's own UMI matrix:

1. The initial cells are Cell Ranger's bootstrapped OrdMag: `recovered_cells` =
   `--expected_cell_num` (3,000), baseline index $$\mathrm{round}(0.01 \times 3000)$$, 100
   bootstrap resamples.
2. The ambient profile is simple Good-Turing on the barcodes at decreasing ranks 45,000 to 89,999
   with nonzero counts; features with no ambient count share $$p_0$$ equally.
3. The candidates are the next barcodes by count, up to 20,000, with at least
   $$\max(500,\ \lceil 0.01 \times \text{median of the initial cells} \rceil)$$ counts.
4. Multinomial log-likelihoods, 10,000 simulations, Monte Carlo p-values and Benjamini–Hochberg;
   a cell if the adjusted $$p \le 0.01$$.

Python's `random` and NumPy's global generator are seeded with 0 at import. The module is marked
as a modified copy of Cell Ranger's `sgt.py` and `stats.py`; how it differs from 3.0.2's was not
checked.

### 1.0 – 1.8: `auto`

Sort the barcodes by UMI count, decreasing. With $$E$$ the expected cells (3,000),
$$i = \lfloor 0.01\,E \rfloor$$ and the threshold is $$\max(1, \lfloor 0.1\,u_{(i)} \rfloor)$$,
where $$u_{(i)}$$ is the $$i$$-th largest count (from 1). Cells are the barcodes with at least the
threshold: with the defaults, a tenth of the 30th-largest count, truncated. Ties in the sort do
not matter, as the threshold is a count. This is an order statistic, not Cell Ranger 2's
[OrdMag](../cellranger/ordmag.md), which takes the 99th percentile of the top $$E$$ counts.

`--force_cell_num` $$N$$ (1.x) picks a count threshold. From 1.7.2 it is the count whose number of
barcodes at or above it is closest to $$N$$, looking only at the two counts around the first that
reaches $$N$$. Before 1.7.2 it is, among the counts whose cumulative numbers span $$N \pm 10\,\%$$,
the one with the most barcodes. Ties make the number of cells differ from $$N$$.

## Versions

| CeleScope | Default `rna` cell call |
|---|---|
| 1.0.1 – 1.8.x (2020-08 – 2022-03) | `auto`, expected cells 3,000 (`--cells 3000` in 1.0.1); options `cellranger3`, the Cell Ranger 3.0.2 port, and, to 1.5.x, `inflection` (an R script, `rescue.R`, removed in 1.6.0) |
| 1.7.2 | `--force_cell_num`'s rule changed, as above |
| 1.9.0 – 1.17.x (2022-04 – 2023-08) | `EmptyDrops_CR`, the same port renamed; FDR 0.01, expected cells 3,000; `auto` still offered |
| 2.0.0 – 2.15.6 (2023-10 – 2026-11) | STARsolo 2.7.11a's EmptyDrops_CR at FDR 0.001 on `GeneFull_Ex50pAS` (STAR 2.6.1 before 2.0); CeleScope's own calling (`tools/count.py`, `--expected_cell_num` now 8,000) is in no assay's steps |

## Parameters

| Original | Default | Meaning |
|---|---|---|
| `--soloCellFilter` (2.x) | `EmptyDrops_CR 3000 0.99 10 45000 90000 500 0.01 20000 0.001 10000` | STARsolo's filter and its parameters, in STAR's order |
| `--soloFeatures` (2.x) | `GeneFull_Ex50pAS Gene` | the matrices STARsolo counts |
| `--report_soloFeature` (2.x) | `GeneFull_Ex50pAS` | the matrix whose filter is the call; `Gene` leaves out intronic reads |
| `--force_cells`, `--barcode`, `--max_mito`, `--min_gene` (`cells` step) | none | re-filter the raw matrix |
| `--cell_calling_method` (1.x) | `auto` (to 1.8), `EmptyDrops_CR` (from 1.9) | the rule |
| `--expected_cell_num` (1.x) | 3000 | the expected cells of `auto` and of the port |
| `--force_cell_num` (1.x) | none | the number of cells to force; also column 4 of the `multi` mapfile |

## Status

Documented only, because since 2.0 it is not a method of its own: CeleScope's default call is
STARsolo 2.7.11a's EmptyDrops_CR at FDR 0.001, which the package runs as
[`starsolo-emptydrops-cr`](../starsolo/emptydrops-cr.md):

```bash
sc-droplets call starsolo-emptydrops-cr raw_feature_bc_matrix.h5 --tool-version 2.7.11a --fdr 0.001
```

The input is the raw `GeneFull_Ex50pAS` matrix of the same sample: CeleScope's own `outs/raw`, an
MTX folder, can be given in place of the `.h5`. Its UMIs have already passed the `read_threshold`
filter. This route has not been checked end to end against a CeleScope run. The 1.x rules are historical; `auto` would be simple to
implement if wanted.

CeleScope's own call sets beside the package's directly. `outs/filtered` holds the called
barcodes, and the counts file the `cells` step writes marks each barcode `CB` (called) or `UB` in
its `mark` column; in 1.x, `*_counts.txt` has the same column, beside `*_matrix_10X` raw and
filtered. Mark the filtered barcodes in a boolean `obs` column of `outs/raw` read with `scd.read`,
and `scd.agreement`, given that column's name among the keys, sets it beside the stored calls.

## References

- Kaminow, Yunusov & Dobin 2021, STARsolo, bioRxiv,
  [doi:10.1101/2021.05.05.442755](https://doi.org/10.1101/2021.05.05.442755).
- Source: [singleron-RD/CeleScope](https://github.com/singleron-RD/CeleScope/tree/v2.15.6) at
  v2.15.6 (`celescope/tools/starsolo.py`, `celescope/tools/cells.py`, `celescope/rna/__init__.py`,
  `doc/CHANGELOG.md`, `LICENSE`); at
  [v1.17.0](https://github.com/singleron-RD/CeleScope/tree/v1.17.0) and
  [v1.9.0](https://github.com/singleron-RD/CeleScope/tree/v1.9.0) for `celescope/tools/count.py`
  and `celescope/tools/emptydrop_cr/cell_calling_3.py`; at
  [v1.0.1](https://github.com/singleron-RD/CeleScope/tree/v1.0.1) for `tools/count.py`. No method
  paper.
- Research notes: docs/_research/docs-only.md §6.
