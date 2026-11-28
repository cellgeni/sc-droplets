---
title: Kneedle's difference curve
parent: Knee geometry
grand_parent: Explained
nav_order: 5
explains: [kneedle]
---

# Kneedle's difference curve
{: .no_toc }

Kneedle (Satopää et al. 2011) is a general-purpose knee finder, and kneed packages it for Python.
Pipelines that need a knee of the rank curve reach for it. Sanger's yascp runs it on the
barcodes with at least 1,000 UMIs, and lowers its sensitivity $$S$$ from 5,000 in steps of 100
until it finds a knee. In the family's terms:

- it draws counts against rank on linear axes, scaled to the unit square;
- it measures each point's height above the chord between the curve's ends (the *difference
  curve*);
- but it does not take the highest point. A walk along the curve picks a local maximum, and the
  sensitivity decides which.

That last step makes Kneedle the odd member of the family. The highest point is where the curve
runs parallel to its chord, as for every chord rule. The point Kneedle takes is set by something
else: how much neighbouring barcodes differ.

1. TOC
{:toc}

## The rule

1. Take the barcodes with at least 1,000 UMIs in rank order: $$x$$ = rank, $$y$$ = count.
2. Scale both to $$[0, 1]$$, and turn the falling, convex curve into a rising, concave one:
   $$y_n \to 1 - y_n$$.
3. The difference curve $$d_i = y_{n,i} - x_{n,i}$$ is the height above the diagonal, the chord
   between the scaled curve's ends.
4. At each local maximum of $$d$$, set a threshold $$T = d_{\max} - S \cdot \overline{\Delta x_n}$$,
   where $$\overline{\Delta x_n} = 1/(n-1)$$ is the mean step.
5. Walk along the curve. At a maximum the threshold becomes its $$T$$. At a minimum it resets
   (and from kneed 0.8.6 nothing is detected until the next maximum). When the next point falls
   below the threshold, the last maximum is the knee.
6. yascp's search: start at $$S = 5{,}000$$ and lower $$S$$ by 100 until a knee appears.

## Why it lands where it does

![Kneedle's difference curve with its maxima, the highest point and the search's knee, and the gaps between neighbouring barcodes against the chord's slope, on PBMC 1k and the spinal cord nuclei]({{ '/assets/explained/kneedle-gaps.png' | relative_url }})

**The highest point** of $$d$$ is the farthest point from the chord. On PBMC 1k the chord falls
47.6 UMIs per rank, and the curve falls that fast only among its first few dozen barcodes: the
highest point is rank 69. Kneedle does not take it.

**A maximum becomes a knee** when, after it, $$d$$ falls by more than $$S/(n-1)$$ before the next
minimum. A large $$S$$ asks for a long fall. yascp's search lowers $$S$$ by 100 at a time until
some maximum has a fall that long, and the first such maximum along the curve is the knee. When
one fall dwarfs the others, as on both samples here, that is **the maximum followed by the longest
fall**. (With few points a step of 100 in $$S$$ is coarse, several maxima can qualify at once, and
the first wins.)

**What makes a fall?** From one barcode to the next, $$d$$ changes by

$$d_{i+1} - d_i = \frac{c_i - c_{i+1}}{c_{\max} - c_{\min}} - \frac{1}{n - 1}$$

so it rises when the gap between neighbouring counts exceeds the chord's slope,
$$(c_{\max} - c_{\min})/(n-1)$$, and falls otherwise. A fall is a stretch in which no two
neighbouring barcodes differ by more than the chord's slope.

So on these samples yascp's Kneedle picks **the barcode that starts the longest stretch in which
no neighbour gap exceeds the chord's slope**. The script behind these pages checks this against
kneed's own walk (ported, and validated against kneed itself) on both samples:

| | PBMC 1k | spinal cord nuclei |
|---|---|---|
| barcodes with at least 1,000 UMIs | 1,174 | 12,666 |
| the chord's slope (UMIs per rank) | 47.6 | 7.1 |
| the highest point of $$d$$ | rank 69 | rank 680 |
| the longest fall | ranks 321 – 1,053, depth 0.535 | ranks 4,389 – 12,666, depth 0.615 |
| the next longest | depth 0.013 | depth 0.029 |
| the sensitivity found | 600 | 5,000 |
| Kneedle's cells | 321 | 4,389 |
| 10x's cells | 1,222 | 15,257 |

Among the top few hundred barcodes, neighbours differ by tens to thousands of UMIs, and the gaps
cross the chord's slope again and again: many small maxima. Further down, neighbours differ by a
few UMIs, and $$d$$ falls without a break. The knee is where that smooth descent begins: the last
barcode in the noisy top whose gap to the next exceeds the chord's slope. On PBMC 1k the one long
fall is 40 times deeper than any other, so the search's answer is robust. It is robust about the
wrong thing, though: where the sampling noise between the largest barcodes dies down, not where
the cells end.

## Where it goes wrong

- **High on the curve.** 321 cells on PBMC 1k and 4,389 on the nuclei, against 10x's 1,222 and
  15,257. The [comparison on public data](../../overview/comparison.md) finds Kneedle stopping at
  the first sharp bend at the top of every whole-cell sample.
- **The floor of 1,000 sets the chord.** The bottom of the curve is the chord's lower end, so the
  floor decides the slope that a neighbour gap must beat.

## Details that matter

Each is on the [method page](../../methods/kneedle/kneedle.md).

- **kneed's releases** changed how plateaus count as extrema (0.5.1), how the polynomial is fitted
  (0.7.0), the pause after a minimum (0.8.6), and the mirroring and interpolation (0.9.0). On
  these two samples every release gives the same knee: 321 and 4,389
  ([versions](../../methods/kneedle/kneedle.md#versions)).
- **No smoothing by default.** `interp1d` evaluates the line through the points at the points
  themselves. The optional degree-7 polynomial is ill-conditioned on ranks in the thousands.
- **`online`** keeps walking and takes the last knee instead of the first; cellsweep 1.0.2 runs
  it that way ([cellsweep's bend](cellsweep.md)).
- **The search can fail.** With $$S$$ below 0, every threshold lies above its maximum, so the
  answer can no longer change. The package stops there with an error; yascp's loop would never
  end.

## Further reading

- Satopää, Albrecht, Irwin & Raghavan 2011, Finding a "Kneedle" in a haystack: detecting knee
  points in system behavior, *ICDCS Workshops* 166–171,
  [doi:10.1109/ICDCSW.2011.20](https://doi.org/10.1109/ICDCSW.2011.20).
- Source: [`kneed/knee_locator.py`](https://github.com/arvkevi/kneed/blob/v0.9.0/kneed/knee_locator.py)
  (0.9.0), and yascp's
  [`get_estimates_from_umi_counts.py`](https://github.com/wtsi-hgi/yascp/blob/main/bin/get_estimates_from_umi_counts.py).
- The reference page: [Kneedle on the barcode-rank curve](../../methods/kneedle/kneedle.md).
