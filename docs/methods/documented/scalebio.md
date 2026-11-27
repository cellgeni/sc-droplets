---
title: ScaleBio
parent: Documented only
grand_parent: Methods
nav_order: 4
method_ids: [scalebio]
status: docs-only
---

# ScaleBio: a count threshold, then an EmptyDrops-like test
{: .no_toc }

The cell calling of ScaleBio's ScaleRna workflow. A threshold of a tenth of the 99th percentile of
counts (TopCell) calls the clear cells. From 2.0, CellFinder then tests the barcodes between 100
counts and that threshold against an ambient profile with a Dirichlet-multinomial likelihood, as
[emptyDrops](../dropletutils/emptydrops.md) does, and rescues those that differ. ScaleBio's
licence lets anyone read the source but not port it. Documented, not implemented.

1. TOC
{:toc}

## At a glance

| Id | `scalebio` (documented only) |
|---|---|
| Tool | ScaleRna 1.4.1 and earlier – 2.1.2 (Nextflow, GitHub `ScaleBio/ScaleRna`) |
| Category | statistical-test |
| Counts | unique transcript counts per barcode from STARsolo (`GeneFull_Ex50pAS`, multimappers `PropUnique`, so fractional) |
| Needs expected cells | no (`expectedCells` optional, per sample) |
| Deterministic | yes: one NumPy generator seeded 42 makes every draw |
| Original licence | ScaleBio's end-user licence: source available, no modification or derivative works |
| Extra inputs | for the ambient-bead filter (QuantumScale), each barcode's RT barcode and bead |

## How it works

2.1.x, per sample and per sub-library, on `counts`, the unique transcript count per barcode:

1. **Ambient.** Barcodes with fewer than `minUTC` (100) counts.
2. **Threshold.** `UTC` if it is set above 0. Otherwise TopCell: the preliminary set is the top
   `expectedCells` barcodes when $$0 <$$ `expectedCells` $$<$$ the number of barcodes, else every
   non-ambient barcode; the threshold is $$\max\big(\mathrm{round}(P_{99} / 10),\ 100\big)$$, with
   $$P_{99}$$ NumPy's linear `percentile` of that set's counts (`topCellPercent` 99,
   `minCellRatio` 10, floor `minUTC`).
3. **Without CellFinder** (barnyard samples, or `cellFinder` off): a cell if its count is at least
   the threshold. With `fixedCells`, the top `expectedCells` by count instead (NumPy's `argsort`,
   so ties follow its sort).
4. **CellFinder** (the default). A count at or above the threshold is a cell. A barcode is
   ambiguous if $$100 \le$$ count $$<$$ threshold and its count is at least `medianFraction` (0.05)
   times the median count of the threshold cells. Each ambiguous barcode is tested:
   - Genes with no count in the sample are dropped. The ambient profile is simple Good-Turing
     (Gale–Sampson) on the per-gene sums over the ambient barcodes, with a pseudo-gene of count 1
     appended and removed; genes unseen among them share the unseen mass in proportion to their
     totals over all barcodes. The switch from Turing to smoothed estimates is at 0.95 standard
     deviations, where Gale and Sampson, and edgeR's `goodTuring` that DropletUtils uses, take
     1.96.
   - The overdispersion $$\alpha$$, unless `--alpha` is given: 1,000 ambient barcodes with at
     least 10 counts, drawn without replacement; $$\alpha$$ maximises the 0.1 %-trimmed mean of
     their Dirichlet-multinomial log-likelihoods under $$\alpha \times$$ the profile, by SciPy's
     bounded scalar search on $$[1, 10000]$$ (`xatol` 0.001); 3,000 if the search fails.
   - Counts are rounded to integers, and each ambiguous barcode's Dirichlet-multinomial
     log-likelihood under $$\alpha \times$$ the profile is computed.
   - Monte Carlo: 10,000 probability vectors from the Dirichlet, each drawn multinomially to the
     smallest ambiguous total, then extended one count at a time, the gene drawn from that vector,
     up to the largest. After each count the log-likelihood gains
     $$\log\frac{n}{n + \alpha - 1} + \log\frac{x_g + \alpha_g - 1}{x_g}$$, with $$n$$ and $$x_g$$
     the new totals. A barcode with total $$n$$ gets $$p = (R + 1)/10001$$, $$R$$ the number of
     simulations at total $$n$$ whose log-likelihood is at most its own.
   - Benjamini–Hochberg (`scipy.stats.false_discovery_control`) on the p-values with $$k$$ zeros
     prepended, one per cell already called. A barcode is rescued, flagged `cellFinder`, if its
     adjusted $$p \le$$ `cellFinderFdr` (0.001).
5. **Ambient beads** (QuantumScale, `filterAmbientBeads`). Per bead, the base-2 KL divergence
   between its counts over RT barcodes and the library's, divided by $$\log_2$$ of the number of
   RT barcodes, over beads with at least two RT barcodes of 100 counts or more and libraries with
   at least 8 RT barcodes. Passing barcodes on beads below `minBeadDivergence` (0.05) fail, flagged
   `ambient_bead`.
6. **Outliers.** Log reads beyond 5 MADs either way, a passing-read fraction 8 MADs below, or a
   mitochondrial fraction 5 MADs above are flagged, and removed only with `filterOutliers`.

**A reading of the source, not a confirmed bug.** `call_cell_barcodes` passes `len(cell_bcs)` as
the number of cells already called, the $$k$$ zeros above. `cell_bcs` is a Boolean array over
every barcode in the sample, so $$k$$ would be the number of barcodes, not the number above the
threshold. With $$k$$ much larger than the number of ambiguous barcodes, the adjustment almost
vanishes: an ambiguous barcode's adjusted p-value is close to its raw one. How much this changes
the calls was not measured. In 1.5 – 2.0, Benjamini–Hochberg ran on the ambiguous barcodes alone.

Compared with DropletUtils' emptyDrops, from the code and not tested: the upper bound is the TopCell
threshold, not the knee; ambient means fewer than 100 counts, not at most `lower`; $$\alpha$$ is
always estimated, on 1,000 barcodes with a trimmed likelihood (DropletUtils 1.28 and later use the
plain multinomial by default); there is the `medianFraction` floor and the Benjamini–Hochberg set
above; fractional counts are rounded; and the draws are NumPy's, not R's.

A module-level `numpy.random.default_rng(42)` (PCG64) makes the ambient sample, the Dirichlet and
multinomial draws and the increments, so the call is the same for the same input and NumPy
release. In barnyard runs the generator's state carries over from one species to the next.

## Versions

| ScaleRna | Change to the call |
|---|---|
| 1.4.1 (2024-02-09) and earlier | TopCell only (`getSampleMetrics.py`): the larger of $$P_{99}/10$$ and `minReads` (100), not rounded; or a fixed number of cells (`useSTARthreshold` takes STARsolo's estimate) |
| 1.5.0 (2024-03-19) | CellFinder added, off by default; `minReads` renamed `minUTC`; `UTC`, `fixedCells` and `expectedCells` reorganised; the threshold rounded |
| 2.0.0, 2.0.1 | CellFinder on by default (QuantumScale RNA becomes the default kit); barnyard samples called per species, without CellFinder |
| 2.1.0 (2025-07-01) | calling per sub-library; the `medianFraction` floor (0.05); Benjamini–Hochberg with $$k$$ zeros prepended; the ambient-bead filter on by default |
| 2.1.1 – 2.1.2 (2026-08-21) | no change to the call |

So TopCell is the default up to 1.6, and CellFinder from 2.0.

## Parameters

| Original | Default | Meaning |
|---|---|---|
| `minUTC` | 100 | below it, ambient; the floor of the threshold and of the ambiguous range |
| `UTC` | 0 (off) | a fixed threshold instead of TopCell |
| `topCellPercent` | 99 | the percentile of the preliminary set |
| `minCellRatio` | 10 | the percentile is divided by it |
| `expectedCells` | 0 (per sample, in `samples.csv`) | the size of the preliminary set; with `fixedCells`, the number of cells |
| `fixedCells` | false | call exactly the top `expectedCells` |
| `cellFinder` | true | test the ambiguous barcodes |
| `cellFinderFdr` | 0.001 | the largest adjusted p-value rescued |
| `medianFraction` | 0.05 | the share of the threshold cells' median count an ambiguous barcode needs (the script's own default is 0; the workflow passes 0.05) |
| `--alpha` (script only) | 0 | the overdispersion; 0 estimates it |
| `filterAmbientBeads` | true | fail barcodes on beads with an ambient-like RT distribution |
| `minBeadDivergence` | 0.05 | the normalised divergence below which a bead is ambient |
| `filterOutliers` | false | remove the MAD outliers rather than flag them |
| `starFeature` | `GeneFull_Ex50pAS` | the STARsolo counts |
| `starMulti` | `PropUnique` (`Unique` for barnyard) | how multimapping reads are counted |

## Status

Documented only. ScaleBio's licence (`LICENSE.md`, an end-user licence agreement) grants use
"solely in connection with the use of ScaleBio Product(s)" and forbids modifying the software,
creating derivative works or combining it with other software. So nothing of it may be ported, and
the documentation is too thin to write an independent version that could be checked.

The workflow's own outputs are easy to set beside the package's calls. `<sample>_allBarcodes.parquet`
holds every barcode's `counts`, `pass` and `flags` (`cellFinder`, `ambient_bead`, `low_reads`,
`high_mito`; `species` in barnyard runs); `<sample>_allCells.csv` the passing barcodes; and
`<sample>_cellcalling_stats.csv` the ambiguous range, the estimated $$\alpha$$ and the smallest
p-value. Run the package's methods on STARsolo's raw matrix of the same sample, put `pass` in a
boolean `obs` column, and `scd.agreement`, given that column's name among the keys, sets it beside the stored calls.

## References

- Lun et al. 2019, EmptyDrops: distinguishing cells from empty droplets in droplet-based
  single-cell RNA sequencing data, *Genome Biology* 20:63,
  [doi:10.1186/s13059-019-1662-y](https://doi.org/10.1186/s13059-019-1662-y) (cited by
  CellFinder's docstring).
- Gale & Sampson 1995, Good-Turing frequency estimation without tears, *J Quant Linguist* 2:217,
  [doi:10.1080/09296179508590051](https://doi.org/10.1080/09296179508590051).
- Source: [ScaleBio/ScaleRna](https://github.com/ScaleBio/ScaleRna/tree/v2.1.2) at v2.1.2
  (`bin/call_cells.py`, `bin/cell_finder.py`, `bin/filter_beads.py`, `docs/cellCalling.md`,
  `changelog.md`, `LICENSE.md`); v1.5.0 (`bin/callCells.py`, `bin/cellFinder.py`); v1.4.1
  (`bin/getSampleMetrics.py`).
- Research notes: docs/_research/docs-only.md §5.
