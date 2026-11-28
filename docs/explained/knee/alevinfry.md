---
title: alevin-fry's knee
parent: Knee geometry
grand_parent: Explained
nav_order: 3
explains: [alevinfry-permit-list]
---

# alevin-fry's knee
{: .no_toc }

`alevin-fry generate-permit-list` decides which barcodes alevin-fry quantifies as cells; simpleaf
runs it with `--knee`. Its knee is [UMI-tools'](umitools.md) rule with three changes:

- the window grows to **five** times the knee, not three;
- of two equally distant points, the **last** wins;
- the result is a count, and every barcode at or above it is kept.

The last two move a barcode or two. The first changes which fixed point the loop finds, and on a
gradual bend that decides everything: on the spinal cord nuclei, alevin-fry keeps 154,250
barcodes where UMI-tools keeps 3,554.

1. TOC
{:toc}

## The rule

1. Sort the barcodes with a count in decreasing order, $$f_0 \ge f_1 \ge \dots$$, and take the
   cumulative counts $$v_i$$.
2. For the first $$m$$ of them, scale both axes, with $$x_i = i/m$$ and $$y_i = v_i / v_{m-1}$$.
   Draw the line from $$(0, y_0)$$ to $$(1, 1)$$, and take the **last** point farthest from it.
3. Start with all of them. While the knee $$k$$ moves (at most 100 times), redraw over the first
   $$\min(n - 1, 5k)$$ and take the knee again.
4. The threshold is $$t = f_k$$, and the cells are every barcode with a count of at least $$t$$.

The package draws the nuclei's result this way: the rank curve with its cut at 128 UMIs, and the
last pass's cumulative curve and chord.

![alevin-fry's diagnostics on the spinal cord nuclei: the cut at 128 UMIs, and the last pass's cumulative curve over 771,245 barcodes]({{ '/assets/explained/alevinfry-diagnostics.png' | relative_url }})

## Why it lands where it does

Scaling the axes cannot move a knee (the [family page](index.md#the-chord-and-the-parallel-tangent)
proves it), so alevin-fry's unit square is the same drawing as UMI-tools'. Two details remain:

- **Its line ends at $$x = 1$$**, one step past the last point, which sits at $$(m-1)/m$$. In
  counts, the chord's slope is $$(v_{m-1} - v_0)/m$$: the mean of $$f_1, \dots, f_{m-1}$$ divided
  by $$m$$ rather than $$m - 1$$.
- **The last maximum wins**, so a barcode whose count equals that mean counts as a cell.

Each pass therefore keeps the barcodes with

$$f_i \ \ge\ \frac{f_1 + \dots + f_{m-1}}{m}$$

The script behind these pages checks this against alevin-fry's knee (ported, and validated against
programs built from alevin-fry's source) on every pass of both samples.

On PBMC 1k the passes run 1,565 → 1,168 → 1,145 → 1,144 → 1,144 cells, over windows of 234,600,
7,820, 5,835, 5,720 and 5,715 barcodes. UMI-tools ends at 1,089. alevin-fry's windows are larger,
so they hold more empty droplets, which pull the mean down and let a few more cells in.

### A second fixed point

On the [family page](index.md#moving-the-window-fixed-points), the refit loop is a map
$$k \mapsto f(ck)$$, and the answer is a fixed point of that map. Here is every fixed point, for
every multiplier $$c$$ from 2 to 8:

![The fixed points of the refit loop against the multiplier, on PBMC 1k and on the spinal cord nuclei, with where the loop ends when it starts from every barcode]({{ '/assets/explained/alevinfry-bifurcation.png' | relative_url }})

On **PBMC 1k** there is one stable branch, at the cliff. The loop reaches it for every
multiplier from 2.05 up, and from 2.5 to 8 it ends between 1,065 and 1,174 cells: the multiplier
hardly matters.

On **the nuclei** the picture has three regimes:

- **Below $$c = 2.70$$** the loop collapses: there is no fixed point for it to reach.
- **From 2.70 to about 4.1** there is one stable branch, inside the nuclei. It rises from 2,635
  cells to about 7,000, and UMI-tools' 3 sits on it, at 3,554.
- **Near $$c = 4.1$$** a second stable branch appears at the second cliff, together with an
  unstable branch between the two. This is a *saddle-node bifurcation*: the map, which used to
  stay below the diagonal there, now touches it and then crosses it twice.

The loop always starts from the knee of every barcode, 158,878, at the second cliff. That start
lies above the unstable branch, so from $$c = 4.1$$ on the loop never comes down: at alevin-fry's
5 it ends at 154,250. The lower branch is still there. Started at 50,000, below the unstable
branch, the same loop would settle at 10,390.

So UMI-tools and alevin-fry, with the same idea, sit on opposite sides of a bifurcation. Which
side a multiplier falls on depends on the curve: PBMC 1k has no bifurcation at all.

## Where it goes wrong

On the nuclei, alevin-fry's knee keeps 154,250 barcodes, with a cut at 128 UMIs. That is every
barcode above the error tail: the nuclei and the whole ambient plateau. The
[comparison on public data](../../overview/comparison.md) finds the same on this sample.

The cause is not the multiplier alone. It is the multiplier together with the start: the first
pass looks at every barcode, the mean of all of them lies below the ambient plateau, and that
puts the start in the upper branch's basin. When the shape of the curve cannot be trusted,
alevin-fry offers other rules: `--expect-cells`, `--force-cells` and an unfiltered permit list.
The first is a robust-maximum rule, with a twist: it takes the count at rank $$0.99E$$, not
$$0.01E$$ ([Same idea, different tools](../../overview/families.md#a-robust-maximum-divided-by-ten)).

## Details that matter

Each is on the [method page](../../methods/alevinfry/permit-list.md).

- **Reads, not UMIs.** alevin-fry counts reads per barcode, before barcode correction. This
  package counts UMIs by default, or reads with `--unit reads`.
- **The window** is $$\min(n - 1, 5k)$$, so the last barcode never enters a refit.
- **Errors.** alevin-fry stops when the knee is 0 or fewer than two counts are left.
- **Ties never matter**, because the cut is a count: every barcode at or above it is a cell.
- **Flex input, 0.12 – 0.17.** Multiplexed input (a sample barcode and a cell barcode) had rules
  of its own. Its knee cut one rank higher, $$t = f_{k-1}$$
  ([versions](../../methods/alevinfry/permit-list.md#versions)).

## Further reading

- He, Zakeri, Sarkar, Soneson, Srivastava & Patro 2022, Alevin-fry unlocks rapid, accurate and
  memory-frugal quantification of single-cell RNA-seq data, *Nature Methods* 19:316–322,
  [doi:10.1038/s41592-022-01408-3](https://doi.org/10.1038/s41592-022-01408-3).
- Source: [`knee_finding.rs`](https://github.com/COMBINE-lab/alevin-fry/blob/v0.19.0/src/knee_finding.rs)
  and [`cellfilter.rs`](https://github.com/COMBINE-lab/alevin-fry/blob/v0.19.0/src/cellfilter.rs)
  (0.19.0).
- Strogatz 2015, *Nonlinear Dynamics and Chaos*, 2nd edition, chapters 3 and 10: saddle-node
  bifurcations, and fixed points of one-dimensional maps and their cobwebs.
- The reference page: [alevin-fry generate-permit-list](../../methods/alevinfry/permit-list.md).
- [UMI-tools' distance knee](umitools.md), the same rule with a multiplier of 3.
