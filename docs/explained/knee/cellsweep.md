---
title: cellsweep's bend
parent: Knee geometry
grand_parent: Explained
nav_order: 6
explains: [cellsweep-empty]
---

# cellsweep's bend
{: .no_toc }

cellsweep removes ambient RNA from droplet data. To estimate the ambient profile it needs to know
which droplets are empty, and when it is not told, its `infer_empty_droplets` finds a cutoff on the
totals. It has had two automatic rules, and they make a useful pair for this family:

- **to 1.0.1**, the most negative second derivative of the smoothed log curve. That is a
  curvature rule, and on a raw matrix it calls two cells.
- **from 1.0.2**, Kneedle on one point per distinct total, keeping the last knee it finds. That
  is a chord rule, and it lands at the foot of the cliff.

(A third rule, `mx_filter`, fits a two-component Gaussian mixture and belongs to the mixture
family.)

1. TOC
{:toc}

## The rule

**To 1.0.1** (`threshold`):

1. Sort the totals above 10 in decreasing order and take $$y = \log_{10} c$$ against rank
   $$1, 2, \dots$$
2. Smooth $$y$$ with a Gaussian filter of $$\sigma = 2$$ ranks.
3. Differentiate twice by central differences against the rank.
4. The cutoff is the total at the most negative second derivative.

**From 1.0.2** (`kneed`):

1. Take each distinct total above 10, in decreasing order, and place it at the last rank that
   holds it.
2. Run kneed's Kneedle on that curve (convex, decreasing) with $$S = 1$$ and `online=True`: walk
   the whole curve, and keep the last knee found.
3. The cutoff is the total at the knee.

The cells are the barcodes at or above the cutoff.

## Why it lands where it does

![The second derivative of the smoothed log curve against rank on PBMC 1k and the nuclei, and 1.0.2's curve of distinct totals with its knees]({{ '/assets/explained/cellsweep-rules.png' | relative_url }})

**The second derivative, against linear rank.** The [family page](index.md#curvature-and-slopes)
shows that a second derivative needs smoothing, and that it depends on how the axes are drawn.
This rule meets both problems at once:

- **The axes.** On a linear rank axis, the top barcodes stand one rank apart but several percent
  apart in count. A change of slope between ranks 1 and 3 is a change per rank squared. The
  cliff's bend, on the other hand, spreads over hundreds of ranks.
- **The smoothing.** A Gaussian of two ranks smooths almost nothing.

On PBMC 1k the second derivative reaches −0.0116 at rank 2. The cliff's dip, at rank 1,169, right
where 10x's cells end, is −0.0026: four times shallower. The cutoff falls at rank 2, and cellsweep
calls two cells. On the nuclei, with no cliff, the minimum is again at rank 2. A log rank axis
would not rescue it: there the tail's long runs of tied counts make the sharpest steps, and the
minimum falls at rank 11,583, among barcodes with 12 UMIs. Without real smoothing, a second
derivative finds the noise that the drawing magnifies most.

**Kneedle, online, on distinct totals.** Here every point is a distinct total, so the flat runs of
tied barcodes become single points. With $$S = 1$$ almost any maximum of the difference curve is a
knee, and `online` keeps walking to the last one. The script behind these pages checks that the
knee is the last local maximum of the difference curve, on both samples.

As on [Kneedle's page](kneedle.md#why-it-lands-where-it-does), the difference curve rises while
the curve falls faster than its chord, and falls otherwise. The last maximum is therefore the last
point at which the curve is steeper than its chord. Below the cliff the shoulder of empty droplets
is flatter than the chord all the way down, so that point is the foot of the cliff:

| | PBMC 1k | spinal cord nuclei |
|---|---|---|
| distinct totals above 10 | 1,364 | 7,478 |
| cutoff | 122 UMIs | 744 UMIs |
| cells | 1,345 | 14,181 |
| 10x's cells | 1,222 | 15,257 |

## Where it goes wrong

- **To 1.0.1, on raw matrices.** The cut is at rank 2 on both samples here: two cells, with
  52,133 and 62,026 UMIs.
- **From 1.0.2, a little generous.** The foot of the cliff includes the barcodes on the cliff
  itself: 1,345 on PBMC 1k, 123 more than 10x. On the nuclei it stops at 744 UMIs, 1,076 short of
  10x, inside the bend.
- **The installed kneed decides.** kneed 0.8.6 and 0.9.0 changed the online walk, so 1.0.2's
  cutoff can depend on which kneed is installed (`--kneed-version`). On these two samples all
  three releases agree.

## Details that matter

Each is on the [method page](../../methods/cellsweep/empty.md).

- **The matrix's dtype.** cellsweep sums totals in the AnnData's dtype, float32 for files read by
  scanpy, and so does this package (`--matrix-dtype`).
- **Manual cutoffs.** `--umi-cutoff` and `--expected-cells` replace the automatic rules in every
  release.
- **`mx_filter`** fits a two-component Gaussian mixture to $$\log(1 + c)$$ of every total, zeros
  included, and cuts at the barcode the mixture is least sure about. It is a mixture rule, for
  that family's pages.

## Further reading

- Caskey, Rich, Weber, Mortazavi, Pachter & Hallgrímsdóttir 2026, Single-cell genomics
  decontamination with CellSweep, *bioRxiv*,
  [doi:10.64898/2026.03.04.709349](https://doi.org/10.64898/2026.03.04.709349).
- Source: [`cellsweep/utils/data_utils.py`](https://github.com/pachterlab/cellsweep/blob/c81ab360/cellsweep/utils/data_utils.py)
  (1.0.2).
- The reference page: [cellsweep empty droplets](../../methods/cellsweep/empty.md).
- [Kneedle's difference curve](kneedle.md): the same walk, run offline with a searched
  sensitivity.
