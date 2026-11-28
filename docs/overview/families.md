---
title: Same idea, different tools
parent: Overview
nav_order: 1
---

# Same idea, different tools
{: .no_toc }

Many tools that look different compute the same thing. They differ in the details: which barcodes
enter, which curve they look at, how they round, which of two equal maxima they take. This page
sets the members of four families side by side and runs them on the same two matrices, to show
how far those details move the call.

1. TOC
{:toc}

## How the tables were made

Every method runs on 10x's raw matrices of PBMC 1k (whole cells, v3 chemistry) and mouse spinal
cord nuclei. Where a method takes a number of expected cells, it gets 3,000, the default of
STARsolo, DropletUtils and Cell Ranger 2; the others run with their defaults. *Smallest* is the
smallest count among the cells. *Jaccard* is the overlap of the cell set with the family's first
row (1: the same barcodes). The tables are made by `docs/_scripts/families.py`; the
[comparison on public data](comparison.md) runs every method with its own defaults instead.

## A robust maximum divided by ten

Take a count near the top of the curve that a few outliers cannot move: the count at rank 1 % of
the expected cells (the 99th percentile of the top $$N$$). Call every barcode with more than a
tenth of it. This is the rule of the Chromium paper (Zheng et al. 2017) and Cell Ranger 1 and 2.

| Tool | The maximum | The cut | The cells |
|---|---|---|---|
| Cell Ranger 2.x – 5.x OrdMag | the $$b$$-th largest total, $$b = \operatorname{round}(0.01 N)$$ (halves up), in each of 100 bootstrap samples (seed 0) | $$\max(1, \operatorname{round}(m/10))$$ in each sample | the top $$\operatorname{round}(\overline{n})$$, $$n$$ the count above the cut per sample |
| Cell Ranger 7.2+ OrdMag | the $$(b+1)$$-th largest, $$b$$ rounded half to even; $$N$$ estimated when not given | as 2.x | as 2.x, then one barcode past the run tied with the last |
| STARsolo CellRanger2.2 | the total at 0-based rank $$\operatorname{round}(0.01 N)$$, no bootstrap | $$\max(1, \operatorname{round}(m/10))$$ | every barcode $$\ge$$ the cut |
| DropletUtils defaultDrops | R's type-7 quantile 0.99 of the top $$N$$ totals | $$0.1\,m$$, unrounded | every barcode $$>$$ the cut |
| alevin `--expectCells` | the count at 0-based rank $$\max(1, \lfloor 0.01 E \rfloor)$$ | $$\max(1, 0.1\,m)$$ | the barcodes before the first below the cut (one more from 1.2.0) |
| alevin-fry `--expect-cells` | the count at rank $$\operatorname{round}(0.99 E)$$, not $$0.01 E$$ | $$\max(1, \operatorname{round}(m/10))$$ | every barcode $$\ge$$ the cut |

<!-- families:robust-maximum:start -->

| Method | Spec | PBMC 1k: cells | smallest | Jaccard | spinal cord nuclei: cells | smallest | Jaccard |
|---|---|---|---|---|---|---|---|
| Cell Ranger 2.x – 5.x OrdMag | `cellranger-ordmag@2.1.0:expect_cells=3000` | 1,133 | 1893 | 1.00 | 5,197 | 3655 | 1.00 |
| Cell Ranger 7.2+ OrdMag | `cellranger-ordmag@10.1.0:expect_cells=3000` | 1,135 | 1868 | 1.00 | 5,217 | 3630 | 1.00 |
| STARsolo CellRanger2.2 | `starsolo-cellranger2.2:expected_cells=3000` | 1,135 | 1868 | 1.00 | 5,228 | 3623 | 0.99 |
| DropletUtils defaultDrops | `dropletutils-defaultdrops:expected_cells=3000` | 1,135 | 1868 | 1.00 | 5,228 | 3623 | 0.99 |
| alevin `--expectCells` | `alevin-initial-whitelist:expect_cells=3000` | 1,136 | 1797 | 1.00 | 5,229 | 3620 | 0.99 |
| alevin-fry `--expect-cells` | `alevinfry-permit-list:expect_cells=3000` | 119,362 | 2 | 0.01 | 14,593 | 680 | 0.36 |

<!-- families:robust-maximum:end -->

Five of the six agree to within a few barcodes on both samples. Bootstrapping (Cell Ranger) or
not (STARsolo), and rounding or not (`defaultDrops`), move only the barcodes at the cut.
alevin-fry is the exception. It takes the count at rank $$0.99 E$$, not $$0.01 E$$. With
$$E = 3{,}000$$ on PBMC 1k, which has about 1,200 cells, that rank is deep among the empty
droplets: its count is 19, the cut is 2, and 119,362 barcodes pass. Given about the right number,
it keeps the cells and some barcodes below them (1,256 barcodes with $$E = 1{,}000$$, 1,421 with
$$E = 1{,}222$$); given too many, it keeps the empty droplets.

## The farthest point from a chord

Draw the curve of the ranked barcodes, and the straight line between two of its points; the knee
is the point of the curve farthest from that line. The tools differ in which curve, which points,
which line, and which of two equally far points.

| Tool | Curve | Points | Line | Knee |
|---|---|---|---|---|
| UMI-tools distance | cumulative counts against rank | every barcode with a count | first to last point | the first farthest point; refit on the first $$3k$$ points until $$k$$ stays |
| alevin-fry `--knee-distance` | cumulative counts, both axes scaled to end at 1 | every barcode with a count | $$(0, y_0)$$ to $$(1, 1)$$ | the last farthest point; refit on the first $$5k$$ until it stays |
| zUMIs | cumulative reads (here UMIs) against rank | barcodes with at least 100 | first to last point | inflection's `uik`: the first point farthest above the line |
| dropbead | cumulative fraction against rank / 15,000 | the top 15,000 | a slope that mixes the two scales | the first farthest point |
| DropletUtils barcodeRanks 1.28 | $$\log_{10}$$ total against $$\log_{10}$$ rank, one point per run of tied totals | runs above 100, between the plateau and the cliff | first to last of those | the farthest among the runs on or above the line |
| Kneedle (yascp) | counts against rank, both scaled to $$[0, 1]$$ | barcodes with at least 1,000 | the diagonal | the first maximum of the difference that the walk accepts; the sensitivity is lowered from 5,000 until one is |

<!-- families:chord:start -->

| Method | Spec | PBMC 1k: cells | smallest | Jaccard | spinal cord nuclei: cells | smallest | Jaccard |
|---|---|---|---|---|---|---|---|
| UMI-tools distance | `umitools-distance` | 1,089 | 2868 | 1.00 | 3,554 | 5683 | 1.00 |
| alevin-fry `--knee-distance` | `alevinfry-permit-list` | 1,144 | 1659 | 0.95 | 154,250 | 128 | 0.02 |
| zUMIs (inflection's `uik`) | `zumis-cell-barcodes` | 594 | 6772 | 0.55 | 15,092 | 615 | 0.24 |
| dropbead | `dropseq-dropbead` | 1,217 | 645 | 0.89 | 4,506 | 4326 | 0.79 |
| DropletUtils barcodeRanks 1.28 | `dropletutils-barcoderanks@1.28` | 845 | 5386 | 0.78 | 2,461 | 7990 | 0.69 |
| Kneedle with yascp's search | `kneedle` | 321 | 8717 | 0.29 | 4,389 | 4442 | 0.81 |

<!-- families:chord:end -->

The same idea gives 321 to 1,217 cells on PBMC 1k and 2,461 to 154,250 on the nuclei. The curve
decides more than the formula does.

- On a cumulative curve the farthest point sits on the cliff, where the counts stop adding up:
  UMI-tools, alevin-fry and dropbead find 1,089 to 1,217 cells on PBMC 1k.
- On a log-log curve (barcodeRanks) or a linear one scaled to a square (Kneedle), it sits at the
  top of the bend: 845 and 321.
- zUMIs draws its line once, over the barcodes with at least 100: it finds 594 cells on PBMC 1k,
  and on the nuclei the cliff (15,092).
- alevin-fry's knee, on the nuclei, lands past the first cliff, at the second one where the empty
  droplets end: 154,250.

## A valley in the density of log counts

The log counts of cells and of empty droplets make two humps; cut at the valley between them.

| Tool | Values | Bandwidth | Which valley |
|---|---|---|---|
| UMI-tools density | $$\log_{10}$$ counts above $$0.001 \times$$ the largest | $$0.1\,\hat\sigma$$ | from the highest down, the first that passes UMI-tools' rules (a fifth of the way along the grid, and far enough from the top) |
| alevin's density correction | $$\log_{10}$$ counts above $$0.001 \times$$ the largest | $$0.1\,\hat\sigma$$ (the same, written as a variance) | after its tangent knee, from the highest down: the first valley whose count lies between a tenth of the knee's and the knee's |
| kernel-density valley | $$\log_{10}$$ counts of at least 100 | Scott's, $$n^{-1/5}\hat\sigma$$, or a factor | the lowest point between the two highest peaks |

Each evaluates the density on a grid of 10,000 points. alevin's correction is UMI-tools' rule
written again in C++, with its own way of choosing the valley.
[DropletQC](../methods/dropletqc/empty-drops.md) cuts at the valley of the nuclear fraction's
density, not of the counts'; DropSift and BAP use the idea too, and are
[planned](index.md#every-method-and-where-it-stands).

<!-- families:density:start -->

| Method | Spec | PBMC 1k: cells | smallest | Jaccard | spinal cord nuclei: cells | smallest | Jaccard |
|---|---|---|---|---|---|---|---|
| UMI-tools density | `umitools-density` | 1,169 | 1182 | 1.00 | 2,925 | 6910 | 1.00 |
| alevin's knee and density correction | `alevin-initial-whitelist` | 1,169 | 1182 | 1.00 | 1,140 | 12113 | 0.39 |
| kernel-density valley, Scott's bandwidth | `kde-valley` | 1,164 | 1255 | 1.00 | 14,277 | 728 | 0.20 |
| kernel-density valley, UMI-tools' bandwidth | `kde-valley:bandwidth=0.1` | 1,171 | 1154 | 1.00 | 14,251 | 732 | 0.21 |

<!-- families:density:end -->

On PBMC 1k the four agree (1,164 to 1,171): there is one clear valley. On the nuclei they part:

- **UMI-tools: 2,925.** The first valley from the top that passes its rules lies at about 6,900
  counts, inside the nuclei.
- **alevin: 1,140.** Its density finds the same valley (2,924 cells), but the correction can only
  lower the number its tangent knee gives, never raise it. The knee gives 1,140, so the knee
  stands. On PBMC 1k the knee gives 1,362, and the valley lowers it to 1,169, UMI-tools' number.
- **The kernel-density valley: about 14,250 to 14,280,** the valley between the nuclei and the
  empty droplets.

The bandwidth barely matters: Scott's and UMI-tools' width differ by 26 cells. Which valley a
tool accepts is what decides.

## A test against the ambient profile

Estimate the ambient RNA from barcodes that surely hold no cell; a barcode whose counts are
unlikely under that profile holds a cell. Lun et al. (2019) proposed it as emptyDrops; Cell
Ranger 3.0 took it up beside OrdMag, and STARsolo and DropletUtils copied Cell Ranger's version.

| Tool | Ambient droplets | Tested | Simulations | Cells |
|---|---|---|---|---|
| DropletUtils emptyDrops | totals up to 100 | every barcode above 100 | 10,000, PCG32 seeded from R | FDR $$\le 0.001$$, and every barcode above the knee |
| Cell Ranger 3.x | ranks 45,000 to 90,000 | up to 20,000 below the OrdMag cells, with at least $$\max(500, 1\%$$ of the median cell$$)$$ | 10,000 | the OrdMag cells, and FDR $$\le 0.01$$ |
| STARsolo EmptyDrops_CR | ranks 45,000 to 89,999 of the detected barcodes | as Cell Ranger 3, stopping at the first below the minimum | 10,000, `std::mt19937` | as Cell Ranger 3; Benjamini-Hochberg not capped at 1 |
| DropletUtils emptyDropsCellRanger | positions 45,000 to 90,000 of every barcode | every total above the minimum, up to 20,000 | 10,000, PCG32 | as Cell Ranger 3 |

Each also estimates the ambient profile its own way: Cell Ranger's Simple Good-Turing, STAR's,
or edgeR's with a pseudo-count. Cell Ranger 7.2 and later change the ambient range by chemistry,
the candidates and the number of simulations ([Cell Ranger cell calling](../methods/cellranger/cellranger.md)).

<!-- families:ambient-test:start -->

| Method | Spec | PBMC 1k: cells | smallest | Jaccard | spinal cord nuclei: cells | smallest | Jaccard |
|---|---|---|---|---|---|---|---|
| Cell Ranger 3.x | `cellranger@3.0.0:expect_cells=3000` | 1,227 | 503 | 1.00 | 15,380 | 500 | 1.00 |
| STARsolo EmptyDrops_CR | `starsolo-emptydrops-cr:expected_cells=3000` | 1,226 | 503 | 1.00 | 15,375 | 500 | 1.00 |
| DropletUtils emptyDropsCellRanger | `dropletutils-emptydropscellranger:expected_cells=3000` | 1,226 | 503 | 1.00 | 15,378 | 501 | 1.00 |
| DropletUtils emptyDrops | `dropletutils-emptydrops` | 1,166 | 101 | 0.91 | 2,290 | 8374 | 0.15 |

<!-- families:ambient-test:end -->

The three copies of Cell Ranger 3's test agree with it (Jaccard 1.00 on both samples), although
each draws other random numbers and estimates the profile with another Good-Turing. With
10,000 simulations, a few borderline barcodes move; the cells do not.

emptyDrops asks the same question of other barcodes:

- It takes every barcode up to 100 as ambient and tests every barcode above it.
- It cuts at FDR 0.001, against Cell Ranger's 0.01.
- It keeps every barcode above its knee without a test.

On PBMC 1k that gives 1,166 cells. On the nuclei it tests 158,080 barcodes, which its 10,000
simulations cannot resolve at FDR 0.001. Only the 2,290 barcodes above its knee remain
([comparison](comparison.md#what-it-shows)).

## What moves the call

Within a family:

- **Small things move a few barcodes at the cut:** rounding, the order of tied counts, the random
  numbers, the estimator of the ambient profile.
- **Large things move the call wholesale:**
  - which curve the rule looks at: cumulative, log-log or linear;
  - which barcodes enter: a floor, the top 15,000, every barcode;
  - which of several candidate points the rule accepts: the first or last maximum, which valley;
  - the number of cells it is told to expect.
- **Whole cells hide these differences.** With a sharp knee, most members of a family land within
  a few percent of each other. Nuclei, with a gradual bend and a heavy ambient plateau, expose
  them.

A tool's name does not say which of these choices it made; its release does. Each method's page
records them, release by release.

The arithmetic behind these differences (ties, rounding, which barcodes count) has a page of its
own: [Reads, UMIs, ties and rounding](arithmetic.md).
