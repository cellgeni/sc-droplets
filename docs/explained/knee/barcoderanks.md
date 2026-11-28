---
title: barcodeRanks' four knees
parent: Knee geometry
grand_parent: Explained
nav_order: 7
explains: [dropletutils-barcoderanks]
---

# barcodeRanks' four knees
{: .no_toc }

DropletUtils' `barcodeRanks()` reports two points on the log-log rank curve: the **knee**, where
the plateau of cells turns into the cliff, and the **inflection**, the steepest point of the
cliff. `emptyDrops` keeps every barcode above the knee without testing it, and many pipelines cut
at one point or the other. Across twelve releases the curve has stayed the same, but the knee
has had four definitions. Each is a different idea from the
[family page](index.md), applied to one drawing:

| Releases | Knee |
|---|---|
| 1.0 – 1.8 | least curvature of a smoothing spline, read off one run too early |
| 1.10 – 1.26 | the same, correctly indexed |
| 1.28 | farthest point above a chord |
| 1.30 on | the tightest bend of a window traced along the curve |

1. TOC
{:toc}

## The rule

Every release draws the same curve, and up to 1.28 brackets it the same way:

1. Sort the totals and collapse each run of tied totals into one point at its mid-rank. Keep the
   runs above 100 UMIs (`lower`), with $$x = \log_{10}$$ rank and $$y = \log_{10}$$ total.
2. The **right edge** is the steepest point, the least slope between neighbouring runs; this is
   also the inflection. The **left edge** is the flattest point before it, the largest slope. From
   1.10, the first 50 barcodes are skipped (`exclude.from`).
3. The knee, between the edges:
   - **to 1.26**: fit `smooth.spline` with 20 degrees of freedom, and take the least curvature;
   - **1.28**: draw the chord between the edges, and take the run farthest above it;
   - **1.30**: see below, with no edges at all.

![PBMC 1k and the spinal cord nuclei with the knees of four releases and the inflections of two]({{ '/assets/explained/barcoderanks-knees.png' | relative_url }})

## Why it lands where it does

**Log-log axes, scaled alike.** Both axes are $$\log_{10}$$, so they share a unit, and curvature
is well defined. Changing the base of the logarithm scales both axes alike, which never moves a
curvature extremum ([family page](index.md#curvature-and-slopes)). The plateau of cells is a gentle
power law, and the cliff is nearly vertical, so every definition finds the corner between them.
On PBMC 1k they agree to within 115 cells:

| | PBMC 1k | spinal cord nuclei |
|---|---|---|
| 1.0 – 1.8: curvature, mis-indexed | 876 | 34,089 |
| 1.10 – 1.26: spline curvature | 960 | 96,857 |
| 1.28: chord | 845 | 2,461 |
| 1.30: traced window | 873 | 2,290 |
| inflection, to 1.28 | 1,171 | 157,825 |
| inflection, 1.30 | 1,264 | 11,091 |
| 10x | 1,222 | 15,257 |

- **The spline** smooths the runs enough for a second derivative to mean something: from the raw
  points, the curvature on PBMC 1k swings between −2,309 and 3,090
  ([family page](index.md#curvature-and-slopes)).
- **The index bug of 1.0 – 1.8.** The least curvature's position is counted among the fitted
  runs, then used among *all* the runs, so the knee lands as many runs too early as the left edge
  is from the start. Those runs are on the plateau, a higher total: 876 cells instead of 960.
- **The chord of 1.28** sits where the curve's local exponent equals the chord's average exponent
  ([family page](index.md#linear-and-log-log-drawings)), at the top of the bend: 845.
- **The inflection** marks the middle of the cliff, not its top, and lands nearer 10x's call:
  1,171 and 1,264 on PBMC 1k.

**1.30's traced window.** Instead of a global construction, 1.30 slides a window of fixed length
along the curve. Length here is arc length on the log-log plot: one unit is a decade. For each
window it draws the chord between the window's two ends and asks two things. Is the window's
midpoint above the chord, so that the window straddles a concave bend? And how far apart are its
ends? A window of arc length 1 whose ends are close together has bent sharply. The knee is the
midpoint of the tightest such window, before the first *elbow* (a window bending the other way
on a slope steeper than −1). The inflection is the midpoint of the steepest window.

![1.30's knee window and inflection window on PBMC 1k's log-log curve]({{ '/assets/explained/barcoderanks-trace.png' | relative_url }})

The window has a length of its own, so 1.30 does not depend on edges found from neighbouring runs.
On PBMC 1k the window that bends most straddles the corner between plateau and cliff: 873 cells,
from 1,063 windows.

## Where it goes wrong

The nuclei separate the releases.

- **To 1.26, the edges are wrong.** The steepest drop above 100 UMIs is the second cliff, where
  the ambient plateau ends: the inflection is at 104 UMIs, 157,825 barcodes down. So the edges
  bracket the ambient plateau and the second cliff, and the spline's least curvature is a bend on
  the plateau, at 212 UMIs: 96,857 barcodes (34,089 with 1.0's mis-indexing). That is a knee of
  the empty droplets.
- **1.28 and 1.30 find the nuclei's upper bend**, 2,461 and 2,290: a sixth of 10x's call.
  1.30's inflection stays among the nuclei (11,091), because its search stops at the first elbow.
- **emptyDrops inherits the knee.** Every barcode above the knee is retained without a test. A new
  knee definition therefore changes emptyDrops' cells by release. That belongs to the
  test-against-ambient family, to come.

## Details that matter

Each is on the [method page](../../methods/dropletutils/barcoderanks.md).

- **Totals are rounded** to integers from 1.12. Before that, unrounded column sums can split a
  run of ties.
- **`exclude.from`** (50) skips the first barcodes from 1.10, so a few huge barcodes cannot be
  an edge.
- **`smooth.spline(df = 20)`** is ported operation for operation from R, down to its golden-section
  search; R built with OpenBLAS can differ in the last bits.
- **The knee is $$10^y$$**, computed in floating point, so it is often not a total. A barcode with
  1,852 UMIs sits above a knee of 1851.9999999999993
  ([differences](../../methods/dropletutils/barcoderanks.md#differences-from-the-original)).
- **`fit.bounds`** replaces the edges with fixed totals, up to 1.28.

## Further reading

- Lun et al. 2019, EmptyDrops: distinguishing cells from empty droplets in droplet-based
  single-cell RNA sequencing data, *Genome Biol* 20:63,
  [doi:10.1186/s13059-019-1662-y](https://doi.org/10.1186/s13059-019-1662-y).
- Source: [`R/barcodeRanks.R`](https://code.bioconductor.org/browse/DropletUtils/blob/RELEASE_3_22/R/barcodeRanks.R)
  (1.30), and the earlier releases' branches.
- The reference page: [DropletUtils barcodeRanks](../../methods/dropletutils/barcoderanks.md),
  with every release's options.
