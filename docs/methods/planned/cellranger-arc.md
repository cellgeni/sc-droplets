---
title: Cell Ranger ARC
parent: Planned
grand_parent: Methods
nav_order: 12
method_ids: [cellranger-arc]
status: planned
---

# Cell Ranger ARC: k-means on ATAC and gene-expression counts
{: .no_toc }

`cellranger-arc count`'s joint caller for 10x Single Cell Multiome. Each barcode is a point: its
ATAC cut sites in peaks and its gene-expression UMIs. Identical points are merged, Cell Ranger's
[OrdMag](../cellranger/ordmag.md) on each axis gives a first split, and two-cluster k-means in log
space settles it. 10x publishes no source; the rule on this page, written from 10x's
documentation, reproduces 10x's calls exactly on three published datasets. Planned for 0.3.

1. TOC
{:toc}

## At a glance

| Id | `cellranger-arc` (planned) |
|---|---|
| Tool | Cell Ranger ARC 1.0.0 – 2.2.0 (`cellranger-arc count`) |
| Category | mixture |
| Counts | ATAC cut sites (transposition events) in peaks, two per fragment inside a peak; GEX UMIs, exonic and intronic |
| Needs expected cells | no |
| Deterministic | yes: k-means from one fixed start; no randomness is documented |
| Original licence | binaries under 10x's EULA; no source. The command-line wrapper in Cell Ranger 10.1 is under 10x's licence (internal use, no redistribution). Never copied |
| Extra inputs | both libraries, and the ATAC exclusions (from the fragments file and peaks, or 10x's `per_barcode_metrics.csv`) |

## How it works

Per barcode, $$A$$ is the number of ATAC transposition events (cut sites) inside peaks
(`atac_peak_region_cutsites`) and $$G$$ the gene-expression UMIs, exonic plus intronic (introns
are counted by default). Cell Ranger ATAC counts fragments in peaks; ARC counts cut sites.

1. **Exclusions.** A barcode is set aside if it is a low-targeting ATAC barcode (fraction of
   fragments in peaks below the fraction of the genome within 2,000 bp of a peak), the minor
   member of an ATAC gel-bead doublet (both as in [Cell Ranger ATAC](cellranger-atac.md), steps 2
   and 3), or has at most 1 count in either library: it needs $$A \ge 2$$ and $$G \ge 2$$. The
   first two are written to `excluded_reason` (1 doublet, 2 low-targeting); the third is not.
   ARC 1.0 and 2.0 outputs carry no reason 3: there is no barcode-multiplet filter.
2. **De-duplication.** The eligible barcodes become points $$(A, G)$$; identical points merge into
   one, so the many low-count background barcodes weigh one point each.
3. **Start.** On the points, a threshold per axis "10-fold less than the maximum value after
   removing outliers" ([OrdMag](../cellranger/ordmag.md)). Points above both thresholds start as
   cells, the others as background. 10x does not say which robust maximum, how many expected cells
   or whether it bootstraps.
4. **k-means, K = 2.** From the two groups' centroids, one start, Euclidean distance on
   $$(\log A, \log G)$$. The base of the logarithm scales both axes alike and does not change the
   assignment. The cells are the points nearer the cell centroid at convergence.
5. **Map back.** Each eligible barcode takes its point's label.
6. **Cap.** Fewer than 20,000 cells (the mechanism is not documented), up to 2.0.2.

With `--min-atac-count X --min-gex-count Y` the caller is skipped: eligible barcodes (step 1
still applies) with $$A \ge X$$ or $$G \ge Y$$ are cells in 1.0, with $$A \ge X$$ and
$$G \ge Y$$ from 2.0. The wrapper refuses one option without the other.

### How sharp the rule is

Taking $$A$$, $$G$$ and `excluded_reason` from 10x's `per_barcode_metrics.csv`, steps 1–5 with
Lloyd iterations converge in 3–4 iterations to exactly 10x's cells from any of three starts
(per-axis 99th percentile ÷ 10 over all points, over the top 3,000, or maximum ÷ 10):

| Dataset (release) | Eligible points | 10x cells | Mismatches |
|---|---|---|---|
| pbmc_granulocyte_sorted_10k (1.0.0) | 25,420 | 11,909 | 0 |
| pbmc_granulocyte_sorted_10k (2.0.0) | 29,162 | 11,898 (11,895 points) | 0 |
| pbmc_granulocyte_sorted_3k (2.0.0) | 6,520 | 2,711 | 0 |

Small changes break it: $$\log(1 + x)$$ gives 2–6 mismatches; requiring 1 count per library
instead of 2 gives 12, requiring 3 gives 6; keeping the excluded barcodes gives 107. With
fragments in peaks, or all fragments,
instead of cut sites, no line in log or linear space separates 10x's cells from the rest. The
OrdMag start cannot be pinned down from these data, since every start reaches the same fixed
point; it would matter only on data where k-means has more than one.

## Versions

| Cell Ranger ARC | Change to the call |
|---|---|
| 1.0.0 | first release; the override calls barcodes meeting **either** threshold |
| 1.0.1 | none |
| 2.0.0 | new peak caller and barcode-aware duplicate marking (change $$A$$); barcode correction at Hamming distance 1 (was 2); the override needs **both** thresholds; `gex_umis_count` is all zeros in `per_barcode_metrics.csv` (an output bug) |
| 2.0.1, 2.0.2 | `gex_umis_count` fixed; nothing else |
| 2.1.0 | "ATAC filtering thresholds have been adjusted to avoid removing real cells", not described further (perhaps Cell Ranger ATAC 2.1's 250 bp padding, inferred); 20,000-cell cap removed; fragment strand used by the doublet and contamination filters |
| 2.2.0 | GEM-X Epi Multiome chemistry; no cell-calling change listed |

No release has public source. The rule is pinned by 10x's documentation and by its 1.0.0 and
2.0.0 outputs; 2.1 and 2.2 are unchecked, as no dataset of theirs was found.

## Parameters

| Original | Default | Meaning |
|---|---|---|
| `--min-atac-count` | unset | with `--min-gex-count`, skip the caller: the ATAC count a cell needs |
| `--min-gex-count` | unset | with `--min-atac-count`: the GEX count a cell needs |
| `--peaks` | the pipeline's | custom peaks, which change $$A$$ and the low-targeting fraction |
| `--gex-exclude-introns` | off (2.0+) | presumably $$G$$ = exonic UMIs only (not checked) |
| `reanalyze --barcodes` | unset | take the cells from a given list |
| internal | 2 per library; K = 2; 20,000 cells (up to 2.0.2) | eligibility; clusters; cap |

## Status

Planned for 0.3. Implementing it needs:

- inputs the package does not read: per barcode $$A$$, $$G$$ and the ATAC exclusions. A reader for
  10x's `per_barcode_metrics.csv` gives all three (in 2.0.0, $$G$$ is `gex_exonic_umis` plus
  `gex_intronic_umis`). The raw matrix gives $$G$$; if its ATAC rows (cut sites per peak) sum to
  $$A$$, which is not checked, it gives $$A$$ too. The exclusions need the fragments file, as for
  [Cell Ranger ATAC](cellranger-atac.md), with ARC's ATAC barcode list (believed to be
  `737K-arc-v1`, not checked);
- a choice of OrdMag start, since 10x documents none exactly; [OrdMag](../cellranger/ordmag.md)
  itself is implemented;
- Variants for 1.0 and 2.0 (the override); 2.1 and 2.2 would be Gaps until a dataset shows what
  changed.

The oracle is 10x's published `per_barcode_metrics.csv` (at `cf.10xgenomics.com/samples/cell-arc/`):
pbmc_granulocyte_sorted_10k and pbmc_granulocyte_sorted_3k (1.0.0 and 2.0.0),
e18_mouse_brain_fresh_5k (1.0.0 and 2.0.0), lymph_node_lymphoma_14k and human_brain_3k (2.0.0).
The comparison is exact. Blockers: the undocumented OrdMag start, 2.1's filter change and the
cap's mechanism.

## References

- 10x Genomics, Cell Ranger ARC, *Calling cell barcodes*, one page per release, e.g.
  [2.0](https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/2.0/algorithms/cell-calling),
  and the [release notes](https://www.10xgenomics.com/support/software/cell-ranger-arc/latest/release-notes/release-notes).
- Zheng et al. 2017, Massively parallel digital transcriptional profiling of single cells,
  *Nat Commun* 8:14049, [doi:10.1038/ncomms14049](https://doi.org/10.1038/ncomms14049), for OrdMag.
- Wrapper only (10x licence; the options, not the caller): [10XGenomics/cellranger](https://github.com/10XGenomics/cellranger/tree/cellranger-10.1.0/lib/rust/cr_wrap/src/arc)
  at tag `cellranger-10.1.0`, `lib/rust/cr_wrap/src/arc/`.
- Research notes: docs/_research/planned-atac.md §2.
