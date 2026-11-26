---
title: Overview
nav_order: 3
has_children: true
---

# Overview: how pipelines tell cells from empty droplets
{: .no_toc }

1. TOC
{:toc}

## The problem

A droplet run sequences far more barcodes than there are cells. Most droplets hold no cell, only
ambient RNA released by cells that lysed in the suspension; many barcodes are sequencing errors
holding a count or two. Ranked by their counts, barcodes trace a characteristic curve: a plateau
of cells, a drop (the "knee"), a shoulder of empty droplets, and a long tail. Every method
decides where the cells end, and they differ in what they look at:

- the **shape of the curve** alone (a cut at some rank or count);
- the **distribution of counts** (a valley between two modes, a fitted mixture);
- the **content of each droplet**: does it differ from the ambient profile?
- **other features** of each droplet (genes detected, mitochondrial or unspliced fraction), fed
  to a classifier.

## Families

| Family | Idea | Examples |
|---|---|---|
| Rank and threshold | A fixed rank or count, or a count derived from the top of the curve | top N, fixed UMI threshold, Cell Ranger 2 OrdMag ("99th percentile of the top N, ÷10"), STARsolo CellRanger2.2, alevin-fry `--expect-cells`, bustools allowlist |
| Knee geometry | A geometric feature of the rank curve: maximum distance to a chord, maximum curvature, a tangent | DropletUtils barcodeRanks, UMI-tools distance, alevin-fry knee, zUMIs, dropbead, Kneedle, cellsweep, dropEst, Cell Ranger targeted gradient |
| Density | A valley or cut in the distribution of log counts | UMI-tools density, alevin's Gauss correction, Otsu, CellBender 0.3's priors, BAP |
| Mixture | A fitted mixture whose components are cells and background | Gaussian mixtures on log counts (zUMIs before 2019, cellsweep mx_filter), Cell Ranger ATAC's negative-binomial mixture |
| Test against ambient | Each droplet's counts tested against the ambient profile estimated from surely-empty droplets | emptyDrops, Cell Ranger 3+, STARsolo EmptyDrops_CR, emptyDropsCellRanger, CB2 |
| Classifier | A model trained on droplet features or counts | alevin's naive Bayes whitelist, dropEst, DIEM, EmptyNN, dropkick, valiDrops, DropSift, SiftCell, CellBender's VAE |

Several methods with different names compute the same thing:

- **Maximum distance to a chord** on the (cumulative) rank curve: UMI-tools `distance`,
  alevin-fry `--knee-distance`, zUMIs, dropbead; Kneedle refines it.
- **A valley in a kernel density** of log counts: UMI-tools `density`, alevin's Gauss
  correction, BAP, DropSift, DropletQC.
- **A robust maximum divided by ten**: Cell Ranger 2 OrdMag, STARsolo CellRanger2.2,
  alevin `--expectCells`, DropletUtils `defaultDrops`, alevin-fry `--expect-cells` (with a
  twist: it takes the 0.99·N-th barcode, not the 0.01·N-th).

## Every method, and where it stands

**Implemented** have a page and run in the package; **planned** are scheduled for a release;
**docs only** will be described but not implemented (they need a heavy model, or their code is
closed).

| Tool | Methods | Status |
|---|---|---|
| Generic | [top N](../methods/generic/top-n.md), [fixed threshold](../methods/generic/threshold.md) | implemented |
| Generic | MAD, Otsu, multi-Otsu, KDE valley, Gaussian mixtures (1-D, 2-D), negative-binomial mixture | planned, 0.1.0 |
| Cell Ranger | [cell calling](../methods/cellranger/cellranger.md) (2.x, 3.x, 7.2, 8.0, 8.0.1, 9.0, 10.0, 10.1, with the antibody-aggregate, 8.0+ and high-occupancy filters), [OrdMag](../methods/cellranger/ordmag.md) | implemented; 1.x and 4.0–7.1 are documented gaps (no public source); targeted (gradient) calling planned |
| STARsolo | [CellRanger2.2](../methods/starsolo/cellranger2-2.md), [TopCells](../methods/starsolo/topcells.md), [EmptyDrops_CR](../methods/starsolo/emptydrops-cr.md), [`soloBasicCellFilter.awk`](../methods/starsolo/basic-awk.md) | implemented |
| DropletUtils | [barcodeRanks](../methods/dropletutils/barcoderanks.md) (four knee definitions), [emptyDrops](../methods/dropletutils/emptydrops.md) (ten behaviours), [emptyDropsCellRanger](../methods/dropletutils/emptydropscellranger.md), [defaultDrops](../methods/dropletutils/defaultdrops.md) | implemented; 1.0–1.2's Dirichlet-multinomial (R's `rgamma`) refused |
| CellBender | [priors](../methods/cellbender/prior.md): probable cells, analysed and surely-empty droplets, from mode and median (0.1, 0.2) or Otsu and a kernel density (0.3–0.4) | implemented; the VAE is docs only |
| cellsweep | [empty-droplet step](../methods/cellsweep/empty.md): second derivative (to 1.0.1), kneed (1.0.2), mx_filter | implemented |
| UMI-tools | whitelist [distance](../methods/umitools/distance.md) and [density](../methods/umitools/density.md) | implemented |
| alevin (salmon 1.x) | knee with Gauss correction, expect cells, force cells | planned, 0.1.0 |
| alevin-fry | [permit list](../methods/alevinfry/permit-list.md): knee distance, expect cells, force cells, unfiltered permit list; 0.12–0.17's own rules for Flex | implemented |
| Drop-seq | dropbead knee | planned, 0.1.0 |
| bustools / kb | [allowlist](../methods/bustools/allowlist.md) (`kb --filter bustools`) | implemented |
| dropEst | cell-number bounds | planned, 0.1.0 |
| zUMIs | UIK knee; Gaussian mixture (before Nov 2019) | planned, 0.1.0 |
| Kneedle | [kneed's KneeLocator](../methods/kneedle/kneedle.md) on the rank curve (0.5.1 – 0.9.0), with yascp's sensitivity search | implemented |
| CB2, DIEM, DropletQC, DropSift, dropEst classifier, alevin naive Bayes, valiDrops, SiftCell | statistical tests and classifiers | planned, 0.2 |
| Cell Ranger ATAC and ARC, EmptyDropsMultiome, BAP, SADE | ATAC and multiome | planned, 0.3 |
| CellBender VAE, EmptyNN, dropkick, BD Rhapsody, ScaleBio, CeleScope, Ambimux, AmbientMapper | | docs only |

## Versions matter

The same matrix gives different cells under different releases of one tool. A few of the changes
the methods' pages describe:

| Tool | Release | Change |
|---|---|---|
| Cell Ranger | 3.0 | an EmptyDrops-style test added after OrdMag |
| Cell Ranger | 7.x | expected cells estimated when not given; ambient range depends on chemistry |
| Cell Ranger | 8.0 | candidates must lie above the ambient range; FDR 0.001 for some chemistries |
| Cell Ranger | 9.0 | 100,000 simulations and a new simulator |
| DropletUtils | 1.2 | emptyDrops fits a Dirichlet-multinomial by default |
| DropletUtils | 1.10, 1.28, 1.30 | three new definitions of the knee, which sets emptyDrops' `retain` |
| DropletUtils | 1.28 | emptyDrops back to a plain multinomial |
| STARsolo | 2.7.8a | EmptyDrops_CR added; CellRanger2.2 rounding changed |
| CellBender | 0.3.0 | priors from Otsu thresholding and a density peak instead of mode and median |
| salmon | 2.0 | alevin removed; its rules are historical |
| kneed | 0.8.6, 0.9.0 | knee detection changed, which moves cellsweep's call |

Previous comparisons of these methods are collected in [Reviews and benchmarks](reviews.md).
