---
title: DropSift
parent: Planned
grand_parent: Methods
nav_order: 4
method_ids: [dropsift]
status: planned
---

# DropSift: SVMs trained on exemplar nuclei, empties and debris
{: .no_toc }

The Broad Institute's DropSift, for single-nucleus data from Drop-seq or 10x. UMIs and the share
of reads in introns separate nuclei from empty droplets: DropSift draws boxes around exemplar
empty droplets, nuclei and debris on those two measures, adds gene-module scores, and trains
radial-kernel SVMs on the exemplars to classify every barcode. There is no paper, and the
repository has no tags, so its versions are commits. Planned for 0.2.

1. TOC
{:toc}

## At a glance

| Id | `dropsift` (planned) |
|---|---|
| Tool | DropSift (GitHub `broadinstitute/DropSift` only; DESCRIPTION 0.99.2, no tags), `SvmNucleusCaller` and its wrapper `runIntronicSVM`; commits f502f47 (Dec 2024) to ccaefb7 (Oct 2026) |
| Category | classifier |
| Counts | UMIs per barcode; UMIs per gene for the module scores |
| Needs expected cells | no |
| Deterministic | no: R's random numbers pick the silhouette subsample, order the pseudobulk and drive libsvm's probability cross-validation, and `random.seed` is unset by default |
| Original licence | MIT |
| Extra inputs | intronic and mitochondrial read fractions per barcode (from Optimus or Drop-seq tools); optionally CellBender's `frac_contamination` |

## How it works

Commit ccaefb7. Each barcode has `num_transcripts` (UMIs), `pct_intronic` and `pct_mt`
(fractions, despite the names): from Optimus, (intronic + antisense intronic reads) / uniquely
mapped reads and mitochondrial reads / uniquely mapped reads, for barcodes with at least 20 UMIs.
Both are read fractions from the aligner, not UMI fractions. Barcodes with no count are dropped.

**Exemplars from the density** (without CellBender, or `useCBRBInitialization = FALSE`):

1. Two runs, on $$x = \log_{10}(\text{UMI} + 1)$$ of the barcodes at or above a floor of 50
   (`maxUmisEmpty`) or 20.
2. **Threshold.** For 12 adjustments $$a = b^m$$, $$b = 0.25, 0.5, \dots, 3$$: R's
   `density(x, adjust = a, n = 200)` (Gaussian, `bw.nrd0`), turning points by
   `pastecs::turnpoints`, and the first pit after the highest peak. For each of the 12, the median
   and MAD of the thresholds within 2 neighbours on either side; the threshold is the local median
   with the smallest MAD relative to it.
3. **Boxes.** For smoothing multipliers $$m = 1.1^0, \dots, 1.1^{14}$$: a threshold on all of
   $$x$$; provisional boxes; barcodes above the empty box's UMI bound with `pct_intronic` < 0.25
   (likely debris) dropped; the threshold again. The empty box: below the threshold, the 75 %
   highest-density region of $$\log_{10}$$ UMI (`hdrcde::hdr`, bandwidth `hdrbw` times 1, 2, 8,
   … until at most 2 intervals), the interval holding the density's maximum, then the 95 % region
   of `pct_intronic` inside it; UMI bound at most the threshold. The nucleus box: above the
   threshold with `pct_intronic` at or above a floor (two candidates: 0.25, or 0.9 times the
   empty box's lower intronic bound, at most 0.25), 75 % regions of `pct_intronic` and then of
   UMI; UMI bounds from the threshold to the 99.5th percentile.
4. **Choice.** Each candidate (multiplier × floor) gets the mean silhouette (`cluster::silhouette`,
   Euclidean, on min-max rescaled $$\log_{10}$$ UMI and `pct_intronic`) of a random 10 % of its
   exemplars. Its score: the median silhouette over 2 neighbours on each side, minus
   $$0.02(1 - \text{stability})$$, the stability being the smaller of the empty and nucleus sets'
   mean overlaps with their neighbours (the smaller of containment and Jaccard). Of the two floors,
   the higher score among those with more than 0.2 empty exemplars per nucleus, if any. If its
   silhouette is at most 0.5, the two-peak solution (the lowest pit between the two highest
   peaks) is computed too, and taken when its silhouette is at least twice as large, or positive
   where the first is not.
5. **Debris box.** Barcodes above the threshold with `pct_intronic` below the floor, those at or
   above their median UMI; 95 % regions of `pct_intronic` then UMI, clipped to UMI at or above the
   threshold and intronic at most the floor; at least 10 barcodes.
6. **Labels.** Box membership on $$\log_{10}$$ UMI (without the + 1 the boxes were built on);
   debris overrides empty, which overrides nucleus.

**Exemplars from CellBender** (the default when `frac_contamination` is given), no debris class.
The empty box: 75 % regions of $$\log_{10}$$ UMI and `pct_intronic` over barcodes with fraction 1
and more than 20 UMIs. The nucleus box: 95 % regions over barcodes with fraction below 1 and at
or below its 25th percentile, joined with the 95 % box of those more than 1 above the empty box's
lower $$\log_{10}$$ UMI bound, in the nucleus intronic range, with fraction $$\le$$ 0.1.
Exemplars: empty in the empty box with fraction 1, nuclei in the nucleus box with fraction
$$\le$$ 0.1.

**Gene modules** (against empty droplets, then debris; each needs 100 exemplars per class): the
negative exemplars' counts, pooled in random order into pseudo-barcodes sized like the nuclei's
shuffled UMI totals; Seurat's `LogNormalize` (10,000); `FindMarkers` (nuclei against pseudobulk,
Wilcoxon, `min.pct = 0.25`, `logfc.threshold = 1`), adjusted $$p \le 0.05$$; the 100 genes most
down in nuclei. Seurat's `AddModuleScore` on them scores every barcode, calling `set.seed(1)`
itself. A module that fails is dropped with a warning.

**SVMs.** `e1071::svm` (C-classification, radial kernel, `probability = TRUE`, cost 1,
$$\gamma = 1/\text{features}$$), trained on the exemplars after $$\log_{10}$$ UMI and z-scores of
every feature over all barcodes: nucleus against empty (UMI, intronic, mito, CellBender's
fraction if used, the empty module), nucleus against debris (the debris module instead), debris
against empty (both modules). $$P_{\text{cell}}$$ is the smaller of P(nucleus) against empty and
against debris; $$P_{\text{debris}}$$ the smaller of P(debris) against empty and 1 − P(nucleus)
against debris.

**Call.** `nucleus` when $$P_{\text{cell}} \ge 0.5$$ (or `cellProbabilityThreshold`), `debris`
when $$P_{\text{debris}} \ge 0.5$$, the larger when both (ties to nucleus), else
`empty_or_other`. Barcodes whose $$\log_{10}(\text{UMI} + 1)$$ is below the empty box's upper UMI
bound get no class. The versions of Seurat (its fold change and Wilcoxon changed in Seurat 5;
presto is used if installed), e1071, hdrcde and R itself all reach the call.

## Versions

| Commits | Change to the call |
|---|---|
| f502f47 (Dec 2024) – be3decb (Apr 2026) | one SVM, nucleus against empty; the call is its predicted class (or probability $$\ge$$ `cellProbabilityThreshold`). 4b1bef3 (Apr 2025) adds `useCBRBInitialization`, 28d30d0 (Jun 2025) `random.seed`, 66f4aaf (Nov 2025) a debris fix |
| c0794aa (#35, May 2026), ed93bff (#36, Jun 2026) | thresholds learnt after filtering debris; the stability score; an experimental 2-D refinement, off by default |
| 6cbee01 (#38, Jun 2026) – 38bfe17 | the debris class, its gene module and two more SVMs; classes from probabilities $$\ge$$ 0.5 |
| f059e3d (#43, Sep 2026) – ccaefb7 (Oct 2026) | the automatic fallback to the two-peak solution (`twoClusterFallbackRatio = 2`); current |

Meanwhile DESCRIPTION went 1.0, 0.99.0, 0.99.1, 0.99.2, out of step with these changes.

## Parameters

| Original | Default | Meaning |
|---|---|---|
| `cellProbabilityThreshold` | none (0.5) | a nucleus has at least this $$P_{\text{cell}}$$ |
| `maxUmisEmpty` | 50 | the higher of the two UMI floors |
| `useCBRBFeatures` | `TRUE` | use CellBender's `frac_contamination` when present |
| `useCBRBInitialization` | `useCBRBFeatures` | exemplars from CellBender rather than the density |
| `forceTwoClusterSolution` | `FALSE` | always the two-peak threshold |
| `twoClusterFallbackRatio` | 2 | how much better the two-peak silhouette must be |
| `use2DTrainingRefinement` | `FALSE` | the experimental refinement |
| `random.seed` (`runIntronicSVM`) | `NA` | no seed |

## Status

Planned for 0.2. Implementing it needs:

- intronic and mitochondrial **read** fractions per barcode, from an Optimus h5ad, Drop-seq tools'
  metrics, a BAM, or perhaps STARsolo's `CellReads.stats` (whether its columns suffice is not
  checked). The package's unspliced and mitochondrial fractions count UMIs;
- R's `density` (`bw.nrd0`) and `sample` (the package has R's Mersenne-Twister), pastecs'
  `turnpoints`, hdrcde's `hdr` and `cluster::silhouette`;
- Seurat's `LogNormalize`, `FindMarkers` and `AddModuleScore`;
- libsvm as e1071 builds it: training and Platt scaling, whose cross-validation draws from R's
  random numbers.

The oracle is R 4.3 or later with DropSift at a pinned commit, Seurat 5, e1071, hdrcde 3.5 and
pastecs, with `random.seed = 1`. An exact match from the matrix would need all of the above to the
last bit, so the comparison is of the calls within a tolerance, or stage by stage: the exemplar
boxes for a fixed seed, then the SVM probabilities given the original's exemplars and features.

## References

- No paper. The authors are James Nemesh and Alec Wysoker (McCarroll lab, Broad Institute).
- Source: [broadinstitute/DropSift at ccaefb7](https://github.com/broadinstitute/DropSift/tree/ccaefb7)
  (`R/SvmNucleusCaller.R`, `R/svmIntronicCellSelection.R`, `R/findTrainingDataBounds.R`,
  `R/CellCallingUtils.R`, `R/geneModuleScore.R`, `R/parseOptimus.R`).
- Research notes: `docs/_research/planned-rna.md` §4.
