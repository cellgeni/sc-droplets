---
title: Knee geometry
parent: Explained
nav_order: 2
has_children: true
---

# Knee geometry
{: .no_toc }

Look at a [barcode-rank curve](../rank-curve.md) and you see at once where the cells end: the
curve bends and then falls off a cliff. To turn that glance into an algorithm, you have to say
what a bend *is*. Nine tools say it in nine ways. This page builds what they share. One
construction, the point farthest from a chord, explains most of them. One of its properties
explains why they disagree: the construction ignores the scale of the axes, but not the shape of
the curve.

1. TOC
{:toc}

## A knee is not a point of the data

Take PBMC 1k and draw its totals four ways. In each drawing, join the curve's two ends with a
straight line and mark the point of the curve farthest from it:

![The farthest point from the chord on four drawings of PBMC 1k, at ranks 1,450, 839, 1,565 and 1,212]({{ '/assets/explained/knee-views.png' | relative_url }})

The rule is the same each time, and the data are the same, but the knee moves: rank 1,450, 839,
1,565 or 1,212. 10x called 1,222. On the spinal cord nuclei the same four drawings give 15,842,
2,578, 158,878 and 404 against 10x's 15,257: a factor of 400 between the extremes.

A knee is a property of a drawing, not of the data. So a knee rule is two choices: how to draw
the curve, and how to find the bend in the drawing. Most of the disagreement in this family
comes from the first choice.

## The chord and the parallel tangent

Fix two points $$A$$ and $$B$$ on a curve and draw the chord between them. The distance of a
point $$P$$ from the chord is

$$d(P) = \frac{\lvert \det(B - A,\; P - A) \rvert}{\lVert B - A \rVert}$$

The numerator is twice the area of the triangle $$ABP$$, and the denominator is its base, which is
the same for every $$P$$. So **the farthest point is the one that makes the largest triangle.**

![A curve with its chord, the farthest point, and the tangent there parallel to the chord; the same after stretching the axes]({{ '/assets/explained/knee-parallel.png' | relative_url }})

Two facts follow.

**The farthest point is where the curve runs parallel to the chord.** For a smooth curve
$$y(x)$$ and a chord of slope $$g$$, the vertical gap to the chord is $$h(x) = y(x) - y(x_A) - g\,(x - x_A)$$.
The perpendicular distance is $$h$$ times a constant, so it peaks where $$h'(x) = 0$$, that is
where $$y'(x) = g$$. The mean value theorem says such a point exists. A knee by chord distance is
the point whose local slope equals the average slope between the ends.

**It does not care about the scale of the axes.** Apply any invertible linear map $$M$$ to the
drawing: stretch the $$x$$ axis, squash the $$y$$ axis, normalise both to $$[0, 1]$$. Every
triangle's area is multiplied by the same $$\lvert \det M \rvert$$:

$$\det\big(M(B - A),\; M(P - A)\big) = \det M \cdot \det(B - A,\; P - A)$$

So the largest triangle stays the largest. In the right panel the distances, the triangle and
the direction of the perpendicular all change, but $$P$$ does not move.

Several tools normalise their axes first: Kneedle and alevin-fry to the unit square, and dropbead
divides rank by 15,000. That normalisation cannot move their knees. What does move a knee is
anything that bends the curve itself: a logarithm, a cumulative sum, or a different pair of
endpoints.

## On the cumulative curve, a chord is a mean

The most common drawing in this family is the cumulative one, with points $$P_i = (i, v_i)$$ and
$$v_i = c_0 + \dots + c_i$$. Take a chord from $$P_a$$ to $$P_b$$. Its slope is

$$g = \frac{v_b - v_a}{b - a} = \frac{c_{a+1} + \dots + c_b}{b - a}$$

That is the **mean count of the barcodes in the window**. The gap above the chord is
$$h_i = (v_i - v_a) - g\,(i - a)$$, and from one barcode to the next it changes by

$$h_i - h_{i-1} = c_i - g$$

The counts are sorted in decreasing order. So $$h$$ rises while $$c_i > g$$, and falls from the
first barcode with $$c_i \le g$$. It peaks at the last barcode above $$g$$:

> On a cumulative curve, the point farthest from a chord is **the last barcode whose count
> exceeds the mean count of the window**. A chord rule on this curve is a threshold in disguise.

Here is that identity at work in dropbead's rule, on PBMC 1k:

![dropbead's cumulative curve of the top 15,000 barcodes; its height above the diagonal; and the counts against the window's mean]({{ '/assets/explained/knee-mean.png' | relative_url }})

Every tool that draws the cumulative curve therefore computes a mean, and its choice of window
decides which mean. The identity holds exactly, pass by pass, against each tool's own code (to
the last barcode, on both samples):

| Tool | Window | PBMC 1k: mean | cells | nuclei: mean | cells |
|---|---|---|---|---|---|
| [dropbead](dropbead.md) | the top 15,000 | 633 | 1,217 | 4,325 | 4,506 |
| [UMI-tools](umitools.md), [alevin-fry](alevinfry.md), first pass | every barcode with a count | 43 | 1,565 | 84 | 158,878 |
| [zUMIs](zumis.md) | barcodes with at least 100 | 6,761 | 594 | 614 | 15,092 |

dropbead's line starts at the origin rather than at the first barcode, so its mean includes the
first barcode; UMI-tools and zUMIs leave it out.

Each row is easy to read once you know it is a mean:

- **dropbead's 15,000** is mostly empty droplets on PBMC 1k, so the mean is low and the cut lands
  at the cliff. On the nuclei the 15,000 are nearly all nuclei. The mean is that of an average
  nucleus, and the rule keeps only the nuclei above average: 4,506 of 15,257.
- **zUMIs' floor of 100** removes nearly all of PBMC 1k's empty droplets. The window is almost
  only cells, and the rule keeps the cells above the average cell: 594, about half. On the
  nuclei the floor lets the whole ambient plateau in (158,149 barcodes), and the mean of 614
  lands near 10x's call.
- **Every barcode** puts the error tail in the window. On PBMC 1k the mean of 43 still cuts at
  the cliff. On the nuclei, the mean of 84 sits below the ambient plateau, so the cut lands at
  the second cliff.

Try it below. "cumulative" and "every barcode" on the nuclei put the knee at the second cliff;
"≥ 100 UMIs" moves it to the first. Stretching the $$y$$ axis never moves the point, while
switching to the log-log drawing does. In the cumulative drawing the readout gives the window's
mean and the last barcode above it. That rank always equals the farthest point. The widget's chord
joins the window's first and last barcodes, as UMI-tools' does. On the top 15,000 of the nuclei
it therefore finds 4,511, where dropbead, whose line starts at the origin, finds 4,506.

<div class="scd-widget" data-widget="chord" data-pbmc="{{ '/assets/explained/curves/pbmc.json' | relative_url }}" data-nuclei="{{ '/assets/explained/curves/nuclei.json' | relative_url }}"></div>

## Moving the window: fixed points

One pass over every barcode puts the error tail in the window. UMI-tools and alevin-fry therefore
refit: find the knee $$k$$, redraw the curve over the first $$3k$$ barcodes (alevin-fry: $$5k$$),
find the knee again, and repeat until it stops moving. Write $$f(m)$$ for the knee of the first
$$m$$ barcodes. The loop is the iteration

$$k_{t+1} = f(c\,k_t), \qquad c = 3 \text{ or } 5$$

and it stops at a fixed point $$k^\ast = f(c\,k^\ast)$$. By the identity of the last section,
that is a self-consistent cell set: **the cells are exactly the barcodes above the mean count of
the top $$c$$ times as many barcodes.**

To see where an iteration goes, plot the map $$k \mapsto f(ck)$$ against the diagonal and bounce
between them: a cobweb diagram. Where the map crosses the diagonal from above there is a stable
fixed point; where it crosses from below, an unstable one.

![Cobweb diagrams of the refit loop on the spinal cord nuclei, with multipliers 3 and 5]({{ '/assets/explained/knee-cobweb.png' | relative_url }})

On the nuclei both loops start at the same first knee, 158,878, at the second cliff.

- **×3 (UMI-tools).** The map has one fixed point, near 3,570. From 158,878 the knee falls pass
  by pass, through the first cliff and into the nuclei, and settles after 24 passes at 3,554
  cells.
- **×5 (alevin-fry).** The map crosses the diagonal three times. There are stable points near
  10,450 and 154,450, with an unstable one near 91,000 between them. 158,878 lies beyond the
  unstable point, so the loop stays at the second cliff and ends at 154,250. Had it started
  anywhere below about 91,000, it would have settled at 10,390. The dotted path starts at 50,000.

A model curve explains the difference. Suppose that over a stretch the counts fall as a power of
rank, $$c_i \propto i^{-a}$$ with $$0 < a < 1$$. The mean of the top $$m$$ is then about
$$c_m / (1 - a)$$, because $$\sum_{i \le m} i^{-a} \approx m^{1-a}/(1-a)$$. A barcode is above
that mean when $$i < m\,(1 - a)^{1/a}$$. So

$$f(m) \approx \rho(a)\, m, \qquad \rho(a) = (1 - a)^{1/a}$$

and the loop multiplies the knee by $$c\,\rho(a)$$ at every pass. It moves down the curve where
$$c\,\rho(a) > 1$$, and back up where $$c\,\rho(a) < 1$$.

![The share of a window above its mean, against the local slope of the curve, with the critical slopes of multipliers 3 and 5]({{ '/assets/explained/knee-powerlaw.png' | relative_url }})

Three consequences follow.

- **On a gentle bend the multiplier must exceed $$e$$.** $$\rho(a)$$ is below $$1/e$$ for every
  slope $$a < 1$$. So with $$c \le e \approx 2.72$$ the knee shrinks along any stretch flatter than
  $$a = 1$$, and only a cliff can hold it. On the nuclei the loop collapses to a handful of
  barcodes for every multiplier below 2.70. A step holds any multiplier above 1: with $$K$$ cells
  of $$C$$ UMIs and nothing below them, the mean of the top $$cK$$ is $$C/c$$. PBMC 1k's cliff is
  close to a step, and there the loop holds down to 2.05. UMI-tools' 3 is just above $$e$$.
- **A larger multiplier needs a steeper bend to stop it.** A stable fixed point sits where the
  curve steepens past the critical slope $$a^\ast$$, where $$c\,\rho(a^\ast) = 1$$:
  - with ×3, $$a^\ast = 0.17$$, so the nuclei's gentle bend is enough;
  - with ×5, $$a^\ast = 0.65$$, which the bend reaches only further down (near 10,450) and the
    second cliff easily.
- **On a sharp cliff the multiplier hardly matters.** PBMC 1k's cliff is so steep that the map is
  almost flat. Any window from 3,000 barcodes to all of them gives a knee between 1,081 and
  1,565. Every whole multiplier from 3 to 8 settles between 1,089 and 1,174 cells. A gradual bend
  is the opposite: the map runs along the diagonal, and the multiplier decides. On the nuclei the
  same multipliers end anywhere from 3,554 to 158,878.

The power law is a local model; the real window's mean includes every barcode above it. The
maps in the figures are computed on the real curves and show the same structure.

Change the multiplier below and watch fixed points appear, merge and vanish. Drag the first knee
to start the loop elsewhere. [alevin-fry's page](alevinfry.md#a-second-fixed-point) draws every
multiplier at once.

<div class="scd-widget" data-widget="cobweb" data-sample="nuclei" data-pbmc="{{ '/assets/explained/curves/pbmc.json' | relative_url }}" data-nuclei="{{ '/assets/explained/curves/nuclei.json' | relative_url }}"></div>

The widget uses UMI-tools' form of the loop for every multiplier. alevin-fry's differs in two
small ways: its mean divides by the window's length rather than one less, and ties at the mean
count as cells. That moves a knee by a barcode or two.

## Linear and log-log drawings

The chord's parallel tangent explains the other two drawings as well.

![Left: the linear curve of barcodes with at least 1,000 UMIs, its chord, and the farthest point at rank 69. Right: barcodeRanks' chord on the log-log curve, with its knee at 845 cells]({{ '/assets/explained/knee-linear-loglog.png' | relative_url }})

- **On linear axes** ([Kneedle](kneedle.md)), the chord over the barcodes with at least 1,000
  UMIs falls by about 48 UMIs per rank. The curve falls that fast only among its first few dozen
  barcodes, so the farthest point is rank 69, at the very top of the bend. Kneedle's own walk
  then picks rank 321, still far above the cliff. Linear axes give the largest barcodes the most
  room, and the knee lands among them.
- **On log-log axes** ([barcodeRanks](barcoderanks.md) 1.28), the tangency condition reads
  $$d \log c / d \log r = g$$. The knee is where the local exponent of the curve first matches
  the chord's average exponent between the plateau and the cliff. That is near the top of the
  cliff: 845 cells on PBMC 1k.

So the drawings sort themselves. Log-log and linear knees sit at the top of the bend. Cumulative
knees sit at the foot of the cliff, where the counts fall below a mean.

## Curvature and slopes

The textbook knee is the point of greatest curvature:

$$\kappa = \frac{y''}{(1 + y'^2)^{3/2}}$$

which is most negative where a falling curve bends hardest. It has two problems the chord does
not have.

**It depends on the aspect ratio.** Stretch $$y$$ by $$\lambda$$. Then $$y'$$ and $$y''$$ are
multiplied by $$\lambda$$, and

$$\kappa_\lambda = \frac{\lambda\, y''}{(1 + \lambda^2 y'^2)^{3/2}}$$

The factor $$\lambda^2$$ inside the denominator moves the extremum. Scaling both axes by the
same factor is harmless: it divides $$\kappa$$ by the factor everywhere, which is why the base of
the logarithm does not matter. Stretching one axis is not harmless. A curvature knee is a
property of the plot's shape, not only of the curve.

**It amplifies noise.** A second difference $$(y_{i+1} - 2y_i + y_{i-1}) / h^2$$ of points with
noise $$\sigma$$ has a standard deviation of $$\sqrt{6}\,\sigma / h^2$$. On PBMC 1k's log-log
curve, the curvature taken straight from the points ranges from −2,309 to 3,090. Taken from a
smoothing spline it ranges from −8.2 to 3.7, and only then does the bend stand out:

![barcodeRanks' smoothing spline on PBMC 1k, and the curvature from the raw points against the curvature from the spline]({{ '/assets/explained/knee-curvature.png' | relative_url }})

Every curvature rule is therefore also a smoothing rule, and the smoothing is part of the
answer:

- **[barcodeRanks](barcoderanks.md) 1.0–1.26** fits a smoothing spline with 20 degrees of
  freedom to the log-log curve between the plateau and the cliff, and takes its least curvature:
  960 cells on PBMC 1k. On the nuclei, the steepest point it brackets is the second cliff, so
  the bend it finds is that one's: 96,857.
- **[cellsweep](cellsweep.md) to 1.0.1** smooths log counts with a Gaussian of two ranks, against
  *linear* rank, and takes the most negative second derivative. On a linear rank axis the
  sharpest change of slope lies among the first few barcodes, so the cut falls at rank 2 on both
  samples.
- **[dropEst](dropest.md)** takes slopes over a lag of 5 % of all barcodes. On a raw matrix that
  is 11,730 ranks for PBMC 1k, nine times the number of cells, and its estimate is 125,227.

A related point is the **inflection**: the steepest point of the cliff, where the log-log slope
is least. It marks the middle of the cliff, not its top. barcodeRanks reports it beside its knee
(1,264 cells on PBMC 1k with 1.30's definition).

## A tangent from the origin

[alevin](alevin.md) draws the log of the cumulative counts of its top barcodes, summed from the
smallest up. It takes the line from the origin that touches this curve, so every later point
lies on or above the line. At the touching point the curve's slope equals its average slope from
the origin, $$y'(j) = y(j)/j$$: the same parallel-tangent idea, with the origin as one end of the
chord. It finds 1,362 cells on PBMC 1k and 1,140 on the nuclei. alevin then corrects the count
with a valley in the density of log counts, which belongs to the density family.

## The family in one table

| Tool | Drawing | Barcodes on it | Quantity | Which point | Refit |
|---|---|---|---|---|---|
| [dropbead](dropbead.md) | cumulative share against rank / 15,000 | the top 15,000 | height above the diagonal | the first maximum | no |
| [UMI-tools](umitools.md) | cumulative counts against rank | every barcode with a count | distance from the chord | the first maximum | $$\times 3$$ until it settles |
| [alevin-fry](alevinfry.md) | cumulative counts, both axes scaled to end at 1 | every barcode with a count | distance from the chord | the last maximum | $$\times 5$$ until it settles |
| [zUMIs](zumis.md) | cumulative reads against rank | at least 100 reads | height above the chord (inflection's UIK) | the first extremum | no |
| [Kneedle](kneedle.md) | counts against rank, scaled to the unit square | at least 1,000 | difference from the diagonal | a maximum its walk accepts | the sensitivity is searched instead |
| [cellsweep](cellsweep.md) | log counts against rank (to 1.0.1); counts against rank (1.0.2) | more than 10 | second derivative; Kneedle | the most negative; the walk's last | no |
| [barcodeRanks](barcoderanks.md) | log counts against log rank, one point per run of ties | more than 100 | curvature (to 1.26), chord (1.28), a window traced along the curve (1.30) | least curvature; farthest above the chord; smallest gap | no |
| [dropEst](dropest.md) | log counts against log rank | every barcode | slopes over a lag | the steepest before the longest convex run | no |
| [alevin](alevin.md) | log cumulative counts, from the smallest | the top barcodes | the tangent from the origin | the first tangent | a density correction follows |

![Where each rule stops on PBMC 1k and on the spinal cord nuclei, against 10x's call]({{ '/assets/explained/knee-landing.png' | relative_url }})

On PBMC 1k most rules stop within a factor of two of 10x's call. The exceptions are rules whose
drawing gives the top barcodes the most room: Kneedle (321), zUMIs (594) and cellsweep's second
derivative (2). dropEst's lag is the other exception, and overshoots. On the nuclei the rules
spread over five orders of magnitude, from 2 to 181,638 cells. The
[comparison of the families](../../overview/families.md#the-farthest-point-from-a-chord) gives
the same numbers as tables.

## What to take away

- **A knee belongs to a drawing.** Choose the drawing (linear, log-log, cumulative) and you have
  mostly chosen the knee.
- **The farthest point from a chord is where the curve runs parallel to it.** Rescaling the axes
  never moves it; logarithms, cumulative sums and other endpoints do.
- **On the cumulative curve, the rule is a threshold:** the last barcode above the window's mean.
  The window (the top 15,000, at least 100, or every barcode) decides where it lands.
- **Refitting turns the rule into a map, and its fixed points are the answer.** On a sharp cliff
  the map is flat and every refit agrees. On a gradual bend the map runs along the diagonal, and
  the multiplier decides.
- **Curvature depends on the aspect ratio and amplifies noise.** Every curvature rule is also a
  smoothing rule.

## Further reading

- Satopää, Albrecht, Irwin & Raghavan 2011, Finding a "Kneedle" in a haystack: detecting knee
  points in system behavior, *ICDCS Workshops* 166–171,
  [doi:10.1109/ICDCSW.2011.20](https://doi.org/10.1109/ICDCSW.2011.20).
- Christopoulos 2016, Introducing Unit Invariant Knee (UIK) as an objective choice for elbow
  point in multivariate data analysis techniques, *SSRN*,
  [doi:10.2139/ssrn.3043076](https://doi.org/10.2139/ssrn.3043076).
- Smith, Heger & Sudbery 2017, UMI-tools: modeling sequencing errors in Unique Molecular
  Identifiers to improve quantification accuracy, *Genome Res* 27:491–499,
  [doi:10.1101/gr.209601.116](https://doi.org/10.1101/gr.209601.116).
- Lun et al. 2019, EmptyDrops, *Genome Biol* 20:63,
  [doi:10.1186/s13059-019-1662-y](https://doi.org/10.1186/s13059-019-1662-y): barcodeRanks'
  knee and inflection.
- [Same idea, different tools](../../overview/families.md) and
  [Reads, UMIs, ties and rounding](../../overview/arithmetic.md): the details that move a call
  within the family.

<script src="{{ '/assets/explained/widgets.js' | relative_url }}" defer></script>
