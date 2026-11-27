---
title: dropEst classifier
parent: dropEst
grand_parent: Methods
nav_order: 2
method_ids: [dropest-classifier]
status: planned
---

# dropEst classifier: kernel densities on quality measures
{: .no_toc }

dropestr's `ScorePipelineCells`. The [cell-number bounds](bounds.md) label the barcodes surely
cells (High) and surely not (Low). Quality measures of each barcode (reads per UMI, UMIs per
gene, the share of genes seen once, the intergenic and mitochondrial fractions) are reduced to
robust sparse principal components, and a kernel density per class on them gives each barcode
its probability of being High. dropEst's report counts a score above 0.9 as a cell. Planned for
0.2.

1. TOC
{:toc}

## At a glance

| Id | `dropest-classifier` (planned) |
|---|---|
| Tool | dropestr 0.7.9 (in dropEst 0.8.5 and 0.8.6), `ScorePipelineCells` |
| Category | classifier |
| Counts | UMIs and genes per barcode, with aligned and intergenic reads per barcode |
| Needs expected cells | no: the bounds estimate it |
| Deterministic | yes: no random numbers |
| Original licence | GPL-3 |
| Extra inputs | aligned reads and intergenic reads per barcode; optionally all reads per barcode (droptag's) and the mitochondrial genes or chromosome |

## How it works

The input is dropEst's `cell.counts.rds`: `cm_raw`, the UMI counts of barcodes with at least
`min_genes_before_merge` genes (20 in the shipped configurations); `cm`, the counts after barcode
merging, with at least 100 genes; aligned reads and UMIs per barcode; reads per chromosome and
barcode, exonic and intergenic; the merge targets; and optionally droptag's `reads_per_cb`, every
read of each barcode, aligned or not.

1. **Labels.** The barcodes of `cm_raw` sorted by UMIs; the [bounds](bounds.md) make ranks up to
   `min` High, ranks from `max` on Low, and the rest Unknown.
2. **Features** (`PrepareLqCellsData`):
   - `ReadsPerUmi`: aligned reads over UMIs;
   - `UmiPerGene`: UMIs over genes;
   - `LowExpressedGenesFrac`: genes with exactly 1 UMI over genes;
   - `IntergenicFrac`: intergenic reads over intergenic plus aligned reads;
   - with `reads_per_cb`, `NotAlignedUmisFrac`: $$\max(\text{total} - \text{aligned}, 0)$$ over
     total reads, divided by `ReadsPerUmi` (totals summed over merged barcodes);
   - with mitochondrial genes or the chromosome's name, `MitochondrionFraction`: the UMI fraction
     in those genes, or the exonic-read fraction on that chromosome.

   Features that are about 0 everywhere are dropped. Each is scaled to [0, 1] as
   $$(x - \min x)/\max(x - \min x)$$.
3. **Robust sparse PCA** (`GetOptimalPcs`). `pcaPP::sPCAgrid` with as many components as features,
   at pcaPP's defaults: `method = "mad"`, `lambda = 1`, `maxiter = 10`, `splitcircle = 25`,
   L1-median centring, a deterministic grid search. The fewest PCs that explain 98 % of the
   variance, at most `max.pcs.number` (3). A feature is used when its normalised absolute loading
   on a kept PC exceeds 0.075. The scores are rescaled to [0, 1].
4. **Outliers.** For the mitochondrial fraction (when given) and the intergenic fraction (when
   present): barcodes above `mean(x, trim = 0.2) + 4 * mad(x)`, or a threshold given. If the
   feature is used, they are relabelled Low before training; if not, their final score is set to
   the smallest score.
5. **Training** (`TrainClassifier`). High is class 1, Low class 0. When there are more than
   1.5 times as many Low as High barcodes, only the round(1.5 × High) Low barcodes with the most
   UMIs are kept. Priors 0.5 and 0.5. Per class, `ks::Hns`, the normal-scale bandwidth matrix
   $$H = \big(4/(n(d + 2))\big)^{2/(d+4)}\,\Sigma$$, with $$n$$ the class's barcodes, $$d$$ the
   PCs and $$\Sigma$$ their covariance.
6. **Score** (`PredictKDE`). For the barcodes of `cm` (every barcode of `cm_raw` with
   `predict.all = TRUE`): each class's density `ks::kde(H = H * kde.bandwidth.mult)` at the
   barcode, binned on ks's default grid when the class has more than 500 points in 2-D or 1,000
   in 3-D and then interpolated, exact otherwise, and clamped at 0. The score is
   $$f_1/(f_1 + f_0)$$, or 0.5 where $$f_0 + f_1 < 10^{-10}$$.
7. **Call.** The function returns scores only. dropEst's report (`scripts/report.Rmd`) and the
   vignette count a score above 0.9 as a cell.

`ScoreQualityData` runs the same classifier on a feature table the caller gives: scaled, reduced
by PCA only when `max.pcs` is below the number of features, without the outlier filters.

Nothing draws random numbers: `sPCAgrid`'s grid search, `Hns` and the KDE are deterministic. The
scores depend on the ks version, through its binning grid and the interpolation of `predict`.
R's `round` takes halves to even in step 5.

## Versions

| dropestr | Change to the call |
|---|---|
| 72e6428 (Jan 2018, untagged) | the current scoring scheme |
| 0.7.7 (eaf0411, Mar 2018) | negative KDE values clamped to 0; the intergenic filter skipped for pseudoaligners |
| d994488 (Feb 2019) | `kde.bandwidth.mult` takes effect; its default, 1, changes nothing |
| 0.7.9 (dropEst 0.8.5, Apr 2019; 0.8.6, Aug 2019) | the tagged releases, the same classifier; `master` changes only documentation after them |

One variant, dropestr 0.7.9. The labels come from the bounds at the same tag.

## Parameters

| Original | Default | Meaning |
|---|---|---|
| `mitochondrion.genes` | `NULL` | the mitochondrial genes, for `MitochondrionFraction` |
| `mit.chromosome.name` | `NULL` | or the mitochondrial chromosome |
| `tags.data` | `NULL` | droptag's reads per barcode, for `NotAlignedUmisFrac` |
| `filter.mitochondrial` | `TRUE` when mitochondrial genes or chromosome are given | the mitochondrial outlier filter |
| `filter.intergenic` | `TRUE` | the intergenic outlier filter |
| `mit.fraction.threshold`, `intergenic.fraction.threshold` | `NULL` (trimmed mean + 4 MAD) | the outlier cuts |
| `max.pcs.number` | 3 | the most PCs kept |
| `predict.all` | `FALSE` | score every barcode of `cm_raw`, not only `cm`'s |
| `kde.bandwidth.mult` | 1 | multiplies both bandwidth matrices |
| the report's cut | 0.9 | a score above it is a cell |

## Status

Planned for 0.2. Implementing it needs:

- the [cell-number bounds](bounds.md), already implemented, for the labels;
- reads per barcode, as `--reads` already takes them ([Inputs](../../reference/inputs.md#read-counts)):
  aligned reads (from Cell Ranger's `molecule_info.h5`, its reads per molecule summed per barcode,
  or STARsolo's per-barcode read statistics), and intergenic reads (from a BAM, or STARsolo's
  `CellReads.stats`). Without reads only `UmiPerGene` and `LowExpressedGenesFrac` can be
  computed. Passing UMIs as reads would make `ReadsPerUmi` constant: it is not dropped (only
  features about 0 everywhere are), and its [0, 1] scaling divides by zero. Such input goes
  through `ScoreQualityData` with the two features;
- pcaPP's `sPCAgrid` (a C++ grid search), ks's `Hns`, and ks's binned KDE with its grid and
  interpolation, ported.

The oracle is dropestr from bioconda's `dropest` 0.8.6 (Linux), or its R files sourced from the
`v0.8.6` tag as the oracle of the [bounds](bounds.md) does, with conda-forge's `r-ks` and
`r-pcapp`. It is deterministic. The practical route is staged: the PCs given as input, then the
KDE scores checked exactly.

## References

- Petukhov, Guo, Baryawno, Severe, Scadden, Samsonova & Kharchenko 2018, dropEst: pipeline for
  accurate estimation of molecular counts in droplet-based single-cell RNA-seq experiments,
  *Genome Biology* 19:78, [doi:10.1186/s13059-018-1449-6](https://doi.org/10.1186/s13059-018-1449-6).
- Source: [`dropestr/R/low_quality_cells.R`](https://github.com/kharchenkolab/dropEst/blob/v0.8.6/dropestr/R/low_quality_cells.R)
  and [`dropestr/R/kde_classifier.R`](https://github.com/kharchenkolab/dropEst/blob/v0.8.6/dropestr/R/kde_classifier.R)
  at v0.8.6 (`ScorePipelineCells`, `PrepareLqCellsData`, `GetOptimalPcs`, `TrainClassifier`,
  `PredictKDE`), and the report's cut in
  [`scripts/report.Rmd`](https://github.com/kharchenkolab/dropEst/blob/v0.8.6/scripts/report.Rmd).
- Research notes: `docs/_research/planned-rna.md` §8.
