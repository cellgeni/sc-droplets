---
title: Choosing a method
parent: Guides
nav_order: 1
---

# Choosing a method
{: .no_toc }

No method is best on every sample. Three questions narrow the choice: whether the cells must be
a pipeline's, what kind of sample it is, and whether the number of cells is known. The numbers
below come from the [comparison on public data](../overview/comparison.md) and
[Same idea, different tools](../overview/families.md).

1. TOC
{:toc}

## Must the cells be a pipeline's?

To match the cells of a published object, or to rerun a paper's analysis, run the pipeline's own
method at the release it used: see [Reproducing a pipeline's calls](reproducing.md). The rest of
this guide is for choosing afresh.

## Whole cells or nuclei?

**Whole cells** with a sharp knee leave little to choose. On the three PBMC samples of the
comparison, 17 of 26 methods call a cell set within Jaccard 0.9 of 10x's. They agree on the
large cells and part at the bottom of the curve:

- A **test against the ambient profile** keeps small barcodes whose counts differ from the ambient
  RNA: [Cell Ranger](../methods/cellranger/cellranger.md) 3.0 and later, STARsolo's
  [EmptyDrops_CR](../methods/starsolo/emptydrops-cr.md), DropletUtils'
  [emptyDropsCellRanger](../methods/dropletutils/emptydropscellranger.md) and
  [emptyDrops](../methods/dropletutils/emptydrops.md).
- A **knee or density rule** stops above them: [OrdMag](../methods/cellranger/ordmag.md),
  [UMI-tools](../methods/umitools/index.md), [Otsu](../methods/generic/otsu.md), the
  [kernel-density valley](../methods/generic/kde-valley.md).

The barcodes the tests add below the knee are often damaged cells or debris. On PBMC 1k they are
mostly mitochondrial RNA. Keep them and filter afterwards, or leave them out; either way, look at
them (see [Comparing methods on a sample](comparing.md)).

**Nuclei** have a gradual bend, and their empty droplets carry much ambient RNA. On the spinal cord
nuclei of the comparison the methods split three ways:

- **Too few.** The knee rules stop in the bend: 1,140 to 6,028 cells, where 10x called 15,257.
- **About right.** The tests against the ambient profile, cellsweep, bustools, zUMIs and the
  distribution rules with a floor (Otsu, the kernel-density valley, the mixtures) find 13,500
  to 20,000.
- **Empty droplets taken in.** alevin-fry's knee, [MAD](../methods/generic/mad.md) and
  [multi-Otsu](../methods/generic/multiotsu.md) call about 155,000.

`dropletutils-emptydrops` at its defaults cannot resolve so many candidates with 10,000
simulations, and keeps only the barcodes above its knee. Raise `--niters`, or check its `limited`
column.

## Is the number of cells known?

Several rules take an expected number of cells: OrdMag before 7.2, STARsolo CellRanger2.2,
`defaultDrops`, alevin's and alevin-fry's `--expect-cells`, `top-n`. Their call moves with it.
alevin-fry's rule reads the count at rank $$0.99 E$$, so with too large an $$E$$ it keeps the empty
droplets: 119,362 barcodes on PBMC 1k with $$E = 3{,}000$$. Cell Ranger 7.2 and later estimate the
number when it is not given. Without a number you trust, prefer a rule that needs none.

## A short list

| Sample | Start with | Check against |
|---|---|---|
| Whole cells, 10x 3′ or 5′ | `cellranger` (the latest release's rule) | `dropletutils-emptydrops`, `otsu` |
| Nuclei | `cellranger` | `kde-valley`, `gmm`, and the unspliced fraction |
| Flex, OCM, or another multiplexed run | `cellranger` on each sample, at the release that made it | [Multiplexed and Flex runs](multiplexed.md) |
| Reads only, no UMIs | `umitools-distance` or `dropseq-dropbead` with `--unit reads` | [Read-based methods](reads.md) |
| Anything unusual | several families at once | [Comparing methods on a sample](comparing.md) |

```bash
sc-droplets compare raw_feature_bc_matrix.h5 -m cellranger -m dropletutils-emptydrops -m otsu -m kde-valley --table calls.csv --plot knee.png
```

The calls are a starting point. Quality control after calling (genes, mitochondrial fraction,
doublets) removes what any of them let through.
