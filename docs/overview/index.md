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
- **A test against the ambient profile**: DropletUtils emptyDrops, Cell Ranger 3+, STARsolo
  EmptyDrops_CR, DropletUtils `emptyDropsCellRanger`.

[Same idea, different tools](families.md) sets each family's members side by side on the same
matrices. [Reads, UMIs, ties and rounding](arithmetic.md) collects the details that move a call
within a family: what is counted, which barcodes enter, ties, rounding and random numbers.

## Every method, and where it stands

**Implemented** have a page and run in the package; **planned** are scheduled for a release;
**docs only** will be described but not implemented (they need a heavy model, or their code is
closed).

| Tool | Methods | Status |
|---|---|---|
| Generic | [top N](../methods/generic/top-n.md), [fixed threshold](../methods/generic/threshold.md), [MAD outliers](../methods/generic/mad.md), [Otsu](../methods/generic/otsu.md), [multi-Otsu](../methods/generic/multiotsu.md), [kernel-density valley](../methods/generic/kde-valley.md), [Gaussian mixtures](../methods/generic/gmm.md) (one and [two dimensions](../methods/generic/gmm-2d.md)) | implemented |
| Generic | negative-binomial mixture (Cell Ranger ATAC's, on RNA) | planned, 0.3, with Cell Ranger ATAC |
| Cell Ranger | [cell calling](../methods/cellranger/cellranger.md) (2.x, 3.x, 4.0, 5.x, 7.2, 8.0, 8.0.1, 9.0, 10.0, 10.1, with the antibody-aggregate, 8.0+ and high-occupancy filters), [OrdMag](../methods/cellranger/ordmag.md) | implemented; 4.0 and 5.x from 10x's outputs (no public source); 1.x and 6.0–7.1 are documented gaps; targeted (gradient) calling planned |
| STARsolo | [CellRanger2.2](../methods/starsolo/cellranger2-2.md), [TopCells](../methods/starsolo/topcells.md), [EmptyDrops_CR](../methods/starsolo/emptydrops-cr.md), [`soloBasicCellFilter.awk`](../methods/starsolo/basic-awk.md) | implemented |
| DropletUtils | [barcodeRanks](../methods/dropletutils/barcoderanks.md) (four knee definitions), [emptyDrops](../methods/dropletutils/emptydrops.md) (ten behaviours), [emptyDropsCellRanger](../methods/dropletutils/emptydropscellranger.md), [defaultDrops](../methods/dropletutils/defaultdrops.md) | implemented; 1.0–1.2's Dirichlet-multinomial (R's `rgamma`) refused |
| CellBender | [priors](../methods/cellbender/prior.md): probable cells, analysed and surely-empty droplets, from mode and median (0.1, 0.2) or Otsu and a kernel density (0.3–0.4) | implemented; the [model](../methods/cellbender/model.md) is docs only |
| cellsweep | [empty-droplet step](../methods/cellsweep/empty.md): second derivative (to 1.0.1), kneed (1.0.2), mx_filter | implemented |
| UMI-tools | whitelist [distance](../methods/umitools/distance.md) and [density](../methods/umitools/density.md) | implemented |
| alevin (salmon 0.11 – 1.10) | [initial whitelist](../methods/alevin/initial-whitelist.md): tangent knee with a kernel-density correction, expected cells, forced cells; the low-confidence band | implemented; the [naive Bayes classifier](../methods/alevin/naive-bayes.md) is planned |
| alevin-fry | [permit list](../methods/alevinfry/permit-list.md): knee distance, expect cells, force cells, unfiltered permit list; 0.12–0.17's own rules for Flex | implemented |
| Drop-seq | [dropbead](../methods/dropseq/dropbead.md): farthest from the diagonal of the cumulative fraction | implemented |
| bustools / kb | [allowlist](../methods/bustools/allowlist.md) (`kb --filter bustools`) | implemented |
| dropEst | [cell-number bounds](../methods/dropest/bounds.md): lagged log-log derivatives, High / Unknown / Low | implemented; the [classifier](../methods/dropest/classifier.md) is planned |
| zUMIs | [cell barcodes](../methods/zumis/cell-barcodes.md): the UIK of the cumulative reads (2.5.6b on) | implemented; the earlier Gaussian mixture is a gap (unseeded) |
| Kneedle | [kneed's KneeLocator](../methods/kneedle/kneedle.md) on the rank curve (0.5.1 – 0.9.0), with yascp's sensitivity search | implemented |
| [CB2](../methods/planned/cb2.md), [DIEM](../methods/planned/diem.md), [DropletQC](../methods/planned/dropletqc.md), [DropSift](../methods/planned/dropsift.md), [valiDrops](../methods/planned/validrops.md), [SiftCell](../methods/planned/siftcell.md), [miQC](../methods/planned/miqc.md) | statistical tests, classifiers and quality mixtures | [planned](../methods/planned/index.md), 0.2 |
| [Cell Ranger ATAC](../methods/planned/cellranger-atac.md) and [ARC](../methods/planned/cellranger-arc.md), [EmptyDropsMultiome](../methods/planned/emptydropsmultiome.md), [BAP](../methods/planned/bap.md), [SADE](../methods/planned/sade.md) | ATAC and multiome | [planned](../methods/planned/index.md), 0.3 |
| [EmptyNN](../methods/documented/emptynn.md), [dropkick](../methods/documented/dropkick.md), [BD Rhapsody](../methods/documented/bd-rhapsody.md), [ScaleBio](../methods/documented/scalebio.md), [CeleScope](../methods/documented/celescope.md), [Ambimux](../methods/documented/ambimux.md), [AmbientMapper](../methods/documented/ambientmapper.md) | classifiers, vendor pipelines, genotype-based calls | [docs only](../methods/documented/index.md) |

## Versions matter

The same matrix gives different cells under different releases of one tool. A few of the changes
the methods' pages describe:

| Tool | Release | Change |
|---|---|---|
| Cell Ranger | 3.0 | an EmptyDrops-style test added after OrdMag |
| Cell Ranger | 4.0 | the test's simulation draws from a random stream of its own |
| Cell Ranger | 5.0 | antibody aggregates removed before calling |
| Cell Ranger | 7.2 | expected cells estimated when not given; ambient range depends on chemistry |
| Cell Ranger | 8.0 | candidates must lie above the ambient range; FDR 0.001 for some chemistries |
| Cell Ranger | 9.0 | 100,000 simulations and a new simulator |
| DropletUtils | 1.2 | emptyDrops fits a Dirichlet-multinomial by default |
| DropletUtils | 1.10, 1.28, 1.30 | three new definitions of the knee, which sets emptyDrops' `retain` |
| DropletUtils | 1.28 | emptyDrops back to a plain multinomial |
| STARsolo | 2.7.8a | EmptyDrops_CR added; CellRanger2.2 rounding changed |
| CellBender | 0.3.0 | priors from Otsu thresholding and a density peak instead of mode and median |
| salmon | 1.11 | alevin removed; its rules are historical |
| kneed | 0.8.6, 0.9.0 | knee detection changed, which moves cellsweep's call |

Previous comparisons of these methods are collected in [Reviews and benchmarks](reviews.md); this
package's own, every method on five public samples, is the [Comparison on public data](comparison.md).
