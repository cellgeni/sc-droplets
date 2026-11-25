---
title: Cell Ranger cell calling
parent: Cell Ranger
grand_parent: Methods
nav_order: 1
method_ids: [cellranger]
status: implemented
---

# Cell Ranger cell calling
{: .no_toc }

What `cellranger count` (and `multi`) put in the filtered matrix. Cell Ranger 2 ran
[OrdMag](ordmag.md) alone. From 3.0 it runs OrdMag, then tests the barcodes below its cut against
the ambient profile of a range of low-ranked droplets, in the manner of DropletUtils' emptyDrops,
and rescues those that differ. Each release since has changed the details.

1. TOC
{:toc}

## At a glance

| | |
|---|---|
| Id | `cellranger` |
| Tool | Cell Ranger 2.0 – 3.1 and 7.2 – 10.1 |
| Category | statistical-test |
| Counts | UMIs, Gene Expression features (Antibody Capture in an antibody-only run) |
| Needs expected cells | 2.x–3.x: 3,000 per GEM group unless given; 7.2+: estimated unless given |
| Deterministic | yes: Cell Ranger's random streams, seeded as Cell Ranger seeds them |
| Original licence | MIT (3.0.2 source); 10x licence (7.2+: implemented from this page) |
| Validated against | 10x's filtered matrices for 2.1, 3.0, 7.2 (Flex, 4 samples), 8.0, 9.0 (exact) and 10.0 (all but 15 of 15,257, at the FDR cut) |

## How it works

### 1. Before calling (7.2+)

- **Antibody aggregates.** With Antibody Capture features, barcodes that are bright in most
  antibodies at once are removed from the matrix before anything else. Of the antibodies with at
  least 1,000 UMIs in total (at least 5 are needed), each nominates its top 25 barcodes (× the
  number of probe barcodes); a barcode among the top 25 by total over these antibodies is an
  aggregate if at least \(\operatorname{round}(n \cdot f)\) antibodies nominate it, with
  \(f = 1.1 - 0.02\,n\) (0.6 above 26 antibodies). Cell Ranger also removes barcodes with a high
  share of corrected reads, which needs read counts; give its `aggregate_barcodes.csv` to remove
  exactly its list.

### 2. OrdMag, per genome and GEM group

See [OrdMag](ordmag.md). From 7.2 the expected number of cells, if not given, is estimated from
the data; the estimate's search is capped at the low end of the ambient range below.
With `--force-cells N`, the top N barcodes are called instead and nothing more is done.

### 3. The ambient range and the FDR (by release and chemistry)

The droplets ranked `low` to `high` (by totals, decreasing; ties in numpy's quicksort order) are
taken as empty. With N partitions, `(low, high) = (N // 2, N)`:

| Kit | 3.x | 7.2 | 8.0.0 | 8.0.1 | 9.0 | 10.0 | 10.1 |
|---|---|---|---|---|---|---|---|
| 3′ v3, 5′ v1/v2 | 90,000 / 0.01 | 90,000 / 0.01 | 90,000 / 0.01 | 90,000 / 0.01 | 90,000 / **0.001** | 90,000 / 0.001 | 90,000 / 0.001 |
| 3′ v3 LT | 90,000 / 0.01 | 9,000 / 0.01 | 9,000 | 9,000 | 9,000 | 9,000 | 90,000 |
| 3′ v3 HT, 5′ HT | 90,000 / 0.01 | 90,000 / 0.01 | 90,000 / 0.01 | **160,000** / 0.01 | 160,000 / 0.001 | 160,000 / 0.001 | 160,000 / 0.001 |
| 3′ v4 (and OH/OCM), 5′ v3 (and R2 OH/OCM) | — | — | 160,000 / 0.001 | 160,000 / 0.001 | 160,000 / 0.001 | 160,000 / 0.001 | 160,000 / 0.001 |
| Flex v1 | 90,000 / 0.01 | 90,000 / 0.01 | 90,000 / 0.01 | 90,000 / 0.01 | 90,000 / 0.01 | 90,000 / 0.01 | 90,000 / 0.01 |
| Flex v2 | — | — | — | — | — | 160,000 / 0.01 | 160,000 / 0.01 |
| anything else | 90,000 / 0.01 | 90,000 / 0.01 | 90,000 / 0.01 | 90,000 / 0.01 | 90,000 / 0.01 | 90,000 / 0.01 | 90,000 / 0.01 |

For a sample of a multiplexed run (Flex probe barcodes, OCM tags; `--n-probe-barcodes m`), N
changes: 7.2–8.x multiply 45,000 (80,000 for v4) by m when m > 1; 9.0+ use 90,000·m (Flex v1),
40,000·m (3′ v4, 5′ v3, HT) and 160,000·m (Flex v2). In 8.0.1 the HT range applies only when the
chemistry was set by hand, since Cell Ranger detects HT kits as their standard version.

The chemistry is read from the input's `chemistry_description` (Cell Ranger writes it into
its HDF5 files) unless `--chemistry` names one, by name (`SC3Pv4`) or description.

### 4. The EmptyDrops step, per GEM group

1. **Ambient droplets**: the barcodes at ranks low to high − 1 with at least one count.
2. **Ambient profile**: their summed counts over the features with any count, smoothed by Cell
   Ranger's Simple Good-Turing; the unseen probability \(P_0\) is shared by the features with no
   ambient count (if there are none: 3.x fails, 7.2+ renormalises). Up to 8.x the step is skipped
   when the log-log slope is above −1; 9.0+ sets such a slope to −1 and goes on. It is skipped
   with fewer than 10 distinct ambient counts.
3. **Candidates**, not called by OrdMag, with at least:
   - 3.x: \(\max(500, \lceil 0.01 \times \text{median OrdMag cell} \rceil)\) UMIs, the 20,000 with the most;
   - 7.2: \(\max(m, \lceil 0.01 \times \text{median} \rceil)\) with \(m\) = `--emptydrops-minimum-umis`
     (500), no limit on their number;
   - 8.0+: \(\max(m, 1 + \text{the largest ambient total})\): above the ambient range.
4. **Observed likelihood** \(\log n! - \sum_g \log x_g! + \sum_g x_g \log p_g\).
5. **Simulations**, a droplet grown from the profile through every candidate total:
   - 3.x – 8.x: 10,000, from numpy's legacy stream seeded 0 (3.x: once, before OrdMag; 7.2+: at
     each genome and GEM group); a block of 1,000,000 pre-drawn features; gaps of 1,000 UMIs or
     more drawn by one multinomial. From 7.2, gaps of 20–999 refill the block early (its unused
     tail plus fresh draws) when they would reach its end.
   - 9.0+: 100,000, from a PCG64 stream seeded 42 once per run; every feature drawn by inverting
     the profile's cumulative sum; 10.0+ orders the features by id for this.
6. **p-value** \((1 + \#\{\text{simulations below the observed likelihood}\})/(1 + \text{simulations})\),
   Benjamini-Hochberg, called if adjusted \(p \le\) the FDR.

With several genomes, 3.x tests all genomes' features together; 7.2+ tests each genome on its own
features (7.2–8.x one genome after another, each treating the earlier genomes' rescues as cells;
9.0+ GEM group by GEM group).

### 5. After calling

- **8.0+**: cells with fewer than `--global-minimum-umis` UMIs (0: none) and with a mitochondrial
  share above `--max-mito-percent` (100: none) are removed. Mitochondrial genes are the 13
  protein-coding human and mouse genes, by Ensembl id.
- **Multiplexed Flex**: after all samples are merged, barcodes in GEMs holding more cells than a
  simulation of random GEM occupancy allows are removed. Call each sample on its own raw matrix,
  then pass the union of their cells to `sc_droplets.methods.cellranger.high_occupancy_gems` (see
  [Run it](#run-it)).

Each barcode's class is `initial`, `rescued`, `candidate` (tested, kept empty), `removed`
(aggregate, or taken out by an 8.0+ filter) or `background`; candidates also get `pvalue`,
`fdr` (adjusted) and `loglik`.

## Versions

| Cell Ranger | Variant | `--tool-version` | Behaviour |
|---|---|---|---|
| 1.0 – 1.3 | — | refused | no source published |
| 2.0 – 2.2 | `v2` | `2.1.0` | OrdMag only |
| 3.0 – 3.1 | `v3` | `3.0.0`, `3.0.2`, `3.1.0` | OrdMag, then EmptyDrops (3.x column of the table) |
| 4.0 – 7.1 | — | refused | no source published (see [Cell Ranger](index.md)) |
| 7.2 | `v7.2` | `7.2.0` | aggregates; OrdMag with estimated cells; chemistry ranges; candidates above 1 % of the median, unlimited; batched simulation |
| 8.0.0 | `v8.0` | `8.0.0` | candidates above the ambient range; GEM-X ranges, FDR 0.001; global-minimum and mito filters |
| 8.0.1 | `v8.0.1` | `8.0.1` | HT range for HT chemistries set by hand |
| 9.0 | `v9.0` | `9.0.0`, `9.0.1` | 100,000 simulations, PCG64; FDR 0.001 for 3′ v3, 5′, HT; slope clamped |
| 10.0 | `v10.0` | `10.0.0` | Flex v2 range; feature-id order in the simulation; genomes sorted |
| 10.1 | `v10.1` | `10.1.0` | LT range removed |

## Parameters

| Option | Default | Original | Meaning |
|---|---|---|---|
| `--expect-cells` | 3,000 per GEM group (2.x/3.x); estimated (7.2+) | `--expect-cells` | expected number of cells |
| `--force-cells` | none | `--force-cells` | call exactly the top N per GEM group instead |
| `--chemistry` | the input file's | `--chemistry` | chemistry name or description (7.2+) |
| `--n-probe-barcodes` | none | `probe_barcode_ids` | probe barcodes or OCM tags of this sample of a multiplexed run (7.2+) |
| `--emptydrops-minimum-umis` | 500 | `emptydrops-minimum-umis` | fewest UMIs a candidate may have (7.2+) |
| `--global-minimum-umis` | 0 | `global-minimum-umis` | remove called cells with fewer UMIs (8.0+) |
| `--max-mito-percent` | 100 | `max-mito-percent` | remove called cells with a higher mitochondrial % (8.0+) |
| `--filter-aggregates` / `--no-filter-aggregates` | on | `filter-aggregates` | remove antibody aggregates before calling (7.2+) |
| `--aggregate-barcodes` | none | `aggregate_barcodes.csv` | Cell Ranger's list of aggregates, removed instead of the count rule (7.2+) |

## Run it

```bash
sc-droplets call cellranger raw_feature_bc_matrix.h5 --tool-version 9.0.0 -o calls.h5ad
sc-droplets call cellranger raw_feature_bc_matrix.h5 --tool-version 3.0.0 --expect-cells 1000 --barcodes cells.tsv.gz
sc-droplets plot diagnostics calls.h5ad --key cellranger_9_0_0 -o diagnostics.png
```

```python
import sc_droplets as scd
from sc_droplets.methods.cellranger import high_occupancy_gems

adata = scd.read("raw_feature_bc_matrix.h5")
scd.call(adata, "cellranger@9.0.0")
adata.obs["cellranger_9_0_0_class"].value_counts()

# A multiplexed Flex run: call each sample, then remove crowded GEMs from the union.
samples = {name: scd.read(f"{name}/count/sample_raw_feature_bc_matrix.h5") for name in ["A", "B", "C", "D"]}
calls = {name: scd.call(a, "cellranger@7.2.0", inplace=False) for name, a in samples.items()}
union = set().union(*(set(r.barcodes[r.is_cell]) for r in calls.values()))
removed, info = high_occupancy_gems(union, "7.2.0")
```

## Differences from the original

- **Tie order.** Cell Ranger picks its ambient droplets, and 2.x/3.x its last OrdMag cells, by
  position in numpy's default sort, and many barcodes tie there. This method reproduces numpy's
  classic quicksort. Recent numpy sorts with SIMD code on some processors (AVX-512) and leaves
  ties elsewhere, so a Cell Ranger run on such a machine can differ at the FDR boundary; 10x's
  10.0 dataset below shows this.
- **Read-based filters** (antibody aggregates by corrected reads, 3.x's aggregate filter) need
  read counts a matrix does not hold; give `--aggregate-barcodes` for 7.2+.
- **Not reproduced**: the antigen (BEAM) aggregate filter; the removal of cells with no UMI in
  the target genes of hybrid-capture and Flex probe-set runs (a cell with no target UMI is
  vanishingly rare); targeted (gradient) calling, which 7.2–9.0 selected for hybrid-capture
  panels.
- 3.0.2 crashes when every feature has an ambient count; this skips the EmptyDrops step instead.

## Validation

`tests/test_validation.py` (marker `integration`) runs the method on 10x's raw matrices and
compares with 10x's filtered barcodes:

| Dataset | Cell Ranger | Options | Cells | Agreement |
|---|---|---|---|---|
| pbmc4k | 2.1.0 | `--expect-cells 5000` | 4,340 | identical |
| pbmc_1k_v2 | 3.0.0 | `--expect-cells 1000` | 996 (34 rescued) | identical |
| pbmc_1k_v3 | 3.0.0 | `--expect-cells 1000` | 1,222 (114 rescued) | identical |
| 4-plex mouse lymph node / spleen Flex, each sample's raw matrix | 7.2.0 | high-occupancy filter on the union | 7,579 / 6,002 / 7,810 / 7,554 | identical, every sample |
| 5k PBMC TotalSeq-B, 3′ v3.1 | 8.0.0 | none (5 antibody aggregates removed) | 5,078 | identical |
| 5k PBMC, 3′ v4 GEM-X | 9.0.0 | none | 5,710 | identical |
| mouse spinal cord nuclei, 3′ v4 | 10.0.0 | none | 15,257 | 12 extra, 3 missing; OrdMag's 5,924 identical; all 15 have adjusted p within 0.0007–0.0011 of the 0.001 cut (tie order, above) |

10x does not publish the `--expect-cells` it used for 2.x/3.x; the values above are the ones that
reproduce its calls.

## References

- Zheng et al. 2017, *Nat Commun* 8:14049, [doi:10.1038/ncomms14049](https://doi.org/10.1038/ncomms14049).
- Lun et al. 2019, EmptyDrops, *Genome Biol* 20:63,
  [doi:10.1186/s13059-019-1662-y](https://doi.org/10.1186/s13059-019-1662-y).
- 10x Genomics, [Cell Ranger cell calling](https://www.10xgenomics.com/support/software/cell-ranger/latest/algorithms-overview/cr-gex-algorithm).
- Source (MIT, 3.0.2): [`cell_calling.py`](https://github.com/lfperales/cellranger/blob/3.0.2/lib/python/cellranger/cell_calling.py),
  [`sgt.py`](https://github.com/lfperales/cellranger/blob/3.0.2/lib/python/cellranger/sgt.py),
  [`stats.py`](https://github.com/lfperales/cellranger/blob/3.0.2/lib/python/cellranger/stats.py).
  7.2 onwards: implemented from this page and a behavioural description of
  [Cell Ranger's public source](https://github.com/10XGenomics/cellranger), not from its code.
