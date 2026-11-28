---
title: The barcode-rank curve
parent: Explained
nav_order: 1
---

# The barcode-rank curve
{: .no_toc }

A droplet run reports hundreds of thousands of barcodes, and a few thousand of them are cells.
To find those few thousand, nearly every method starts from the same picture: each barcode's
total, sorted from largest to smallest. This page is about that picture. It covers why it has
the shape it has, what each way of drawing it shows and hides, and how a sample of nuclei bends
it out of shape. Every family of methods reads this curve or something made from it, so the
pages after this one assume it.

1. TOC
{:toc}

## Three kinds of barcode

![The rank curve of PBMC 1k on log-log axes, with the plateau, the cliff, the shoulder and the tail labelled]({{ '/assets/explained/rank-anatomy.png' | relative_url }})

This is PBMC 1k. The top 1,222 barcodes, blue, are the ones 10x called cells. The curve has four
parts, and each comes from a different kind of barcode.

- **Cells.** A droplet that caught a cell holds thousands of its molecules: here a median of
  6,628 UMIs. Cell size and capture efficiency vary by factors, not by increments, so the totals
  spread roughly log-normally. On log axes that spread is the gently falling **plateau** at the
  top.
- **Empty droplets.** Most droplets catch a bead but no cell. They still hold RNA: transcripts
  that cells released into the suspension before it was loaded, the *ambient* RNA. That is tens
  of UMIs each (here a median of 13, at ranks two to ten times the number of cells), from tens of thousands of droplets:
  the **shoulder**.
- **Error barcodes.** Sequencing and synthesis errors create barcodes that match no bead well
  enough to be corrected. Each holds one or two UMIs, and there are many of them: the flat
  **tail**.

Between cells and empty droplets the curve falls through two orders of magnitude within a few
hundred ranks: the **cliff**. Where the cliff starts, the curve bends: the **knee**.

A simulation makes the parts explicit. The package's simulator (`scd.simulate_droplets`) draws
log-normal totals for 1,500 cells, 8,000 empty droplets and 30,000 error barcodes:

![A simulated rank curve coloured by kind of barcode, and the histogram of the same totals]({{ '/assets/explained/rank-toy.png' | relative_url }})

Each kind of barcode owns a stretch of the curve, and each is a hump in the histogram on the
right. Move the sliders to change the populations. The cliff only exists when cells hold far more
than empty droplets: bring the two medians within a factor of ten of each other, and the cliff
becomes a bend.

<div class="scd-widget" data-widget="toy"></div>

## A survival function on its side

The curve has a simpler description when you turn it on its side. Read it as "how many barcodes
have at least this many UMIs":

$$r(c) = \#\{\, b : t_b \ge c \,\}$$

This is the empirical survival function of the totals, scaled by the number of barcodes. The
rank curve is the same function with its axes swapped. With three populations it is a sum of
three survival functions:

$$r(c) = N_{\text{cell}}\, S_{\text{cell}}(c) + N_{\text{empty}}\, S_{\text{empty}}(c) + N_{\text{err}}\, S_{\text{err}}(c)$$

Several features of the curve follow from this sum:

- **Above the ambient range** only cells contribute, so $$r(c) = N_{\text{cell}} S_{\text{cell}}(c)$$.
  This is the plateau, read sideways.
- **Between the two populations** few barcodes have counts. There, $$r(c)$$ barely changes: it
  holds near $$N_{\text{cell}}$$ while $$c$$ falls by a factor of a hundred. Drawn the usual way
  round, that flat stretch of $$r(c)$$ is the vertical cliff, and it stands at the rank equal to
  the number of cells. The cliff is vertical when no barcode lies between the populations. It
  leans when barcodes fill the gap: low-quality cells, cell fragments, or droplets rich in ambient
  RNA.
- **Within a population**, $$r(c)$$ falls fastest where that population's counts are densest.

![The PBMC 1k rank curve turned on its side above the histogram of its log totals]({{ '/assets/explained/rank-survival.png' | relative_url }})

The last point is the bridge to the histogram. Write $$u = \log_{10} c$$. The barcodes whose log
total lies between $$u$$ and $$u + du$$ number $$r(10^u) - r(10^{u + du})$$. So the histogram of
log totals is $$-\,dr/du$$, the slope of the sideways curve. Each steep fall in the top panel
is a hump in the bottom one.

The rank curve and the histogram hold the same information: one is the integral, the other the
derivative. The families differ in which one they read:

- **Knee geometry** reads the integral. It is smooth, so rules on it need little smoothing.
- **Density rules** read the derivative. A derivative of data is noisy, so they choose a
  bandwidth.

## Four ways to draw it

The same totals can be drawn in several ways, and each way shows something different.

![PBMC 1k drawn four ways: linear, log-log, cumulative share against log rank, and the histogram of log totals]({{ '/assets/explained/rank-views.png' | relative_url }})

- **Linear axes** show only the largest barcodes. The top 2,000 ranks fill the panel. Below a few
  thousand UMIs everything sits on the axis, and the 233,000 other barcodes are a line at zero.
- **Log-log axes** give each order of magnitude the same room, so cells, empty droplets and the
  tail all show. A stretch where counts fall as a power of rank, $$c \propto r^{-a}$$, is straight,
  with slope $$-a$$. The plateau has $$a$$ well below 1; the cliff has $$a$$ above 10.
- **The cumulative share** $$v_i / v_{n-1}$$ is the fraction of all UMIs held by the top
  $$i + 1$$ barcodes. It rises steeply through the cells (10x's 1,222 hold 90.7 % of PBMC 1k's
  UMIs) and then flattens. Its slope is the count itself: $$v_i - v_{i-1} = c_i$$. So the curve
  is concave, and where it stops rising is where the counts become small.
- **The histogram of log totals** shows the populations as humps, with a valley between cells and
  empty droplets.

Every method in the knee family picks one of these drawings. On the
[next page](knee/index.md), the same rule applied to different drawings lands on different
barcodes.

## Whole cells against nuclei

![The rank curves of PBMC 1k and of mouse spinal cord nuclei side by side]({{ '/assets/explained/rank-nuclei.png' | relative_url }})

The two samples differ in how much a cell holds relative to the ambient RNA:

| | PBMC 1k | spinal cord nuclei |
|---|---|---|
| UMIs per called cell (median) | 6,628 | 2,237 |
| UMIs per empty droplet (median at ranks 2–10 × the cells) | 13 | 215 |
| ratio | about 500 | about 10 |

A nucleus holds a fraction of a cell's RNA, and dissociating tissue to nuclei releases far more
ambient RNA. With a ratio of 10 the two humps of the histogram touch, and the cliff softens into
a gradual bend. 10x's 15,257 cells run well into it.

The nuclei's empty droplets hold so much ambient RNA that they form a plateau of their own, about
200 UMIs high, which ends in a second cliff near rank 155,000. Below it are the error barcodes. A
rule that looks only at the shape of the curve now sees two bends, a gradual one and a sharp one,
and nothing in the shape says which of them is the cells. Most of the differences between knee
rules on the next pages appear on this sample.

## Where each family looks

| Family | What it reads |
|---|---|
| Rank and threshold | one point of the curve: a rank, or a count scaled from the top |
| [Knee geometry](knee/index.md) | the shape of the curve: where it bends |
| Density | its derivative: a valley in the histogram of log totals |
| Mixture | the derivative again, modelled as a sum of components |
| Test against ambient | not the curve but each droplet's content; the curve only chooses which droplets are surely empty |
| Classifier | features beyond the total: genes detected, the mitochondrial or unspliced fraction |

## Further reading

- Macosko et al. 2015, Highly parallel genome-wide expression profiling of individual cells using
  nanoliter droplets, *Cell* 161:1202–1214,
  [doi:10.1016/j.cell.2015.05.002](https://doi.org/10.1016/j.cell.2015.05.002). Drop-seq, and
  the cumulative curve whose bend sets the number of cells.
- Lun et al. 2019, EmptyDrops: distinguishing cells from empty droplets in droplet-based
  single-cell RNA sequencing data, *Genome Biol* 20:63,
  [doi:10.1186/s13059-019-1662-y](https://doi.org/10.1186/s13059-019-1662-y). The barcode-rank
  plot, and why the ambient RNA is information.
- Young & Behjati 2020, SoupX removes ambient RNA contamination from droplet-based single-cell RNA
  sequencing data, *GigaScience* 9:giaa151,
  [doi:10.1093/gigascience/giaa151](https://doi.org/10.1093/gigascience/giaa151). Where the
  ambient RNA comes from, and how much of each cell it makes up.
- [Visualising calls](../guides/visualising.md): the package's own knee and density plots.

<script src="{{ '/assets/explained/widgets.js' | relative_url }}" defer></script>
