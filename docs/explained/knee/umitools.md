---
title: UMI-tools' distance knee
parent: Knee geometry
grand_parent: Explained
nav_order: 2
explains: [umitools-distance]
---

# UMI-tools' distance knee
{: .no_toc }

`umi_tools whitelist` chooses a run's cell barcodes from the reads, before any alignment. Drop-seq
and inDrop workflows built on UMI-tools take their cells from it, and since 1.0.0 its default
rule is the *distance* knee. In the family's terms:

- it draws the cumulative curve over every barcode with a count;
- it takes the point farthest from the chord between the curve's ends;
- it redraws the curve over three times as many barcodes as the knee, and repeats until the knee
  settles.

By the [family's identity](index.md#on-the-cumulative-curve-a-chord-is-a-mean), the answer is a
self-consistent cell set: the barcodes above the mean count of the top three times as many.

1. TOC
{:toc}

## The rule

1. Sort the barcodes with a count in decreasing order and take the cumulative counts
   $$v_0, v_1, \dots, v_{n-1}$$.
2. Find the knee $$k$$: the first point of the curve farthest from the chord between
   $$(0, v_0)$$ and $$(n - 1, v_{n-1})$$.
3. Redraw the curve over the first $$\min(3k, n)$$ barcodes, with the chord between its new ends,
   and find the knee again. Repeat until $$k$$ no longer changes, at most 100 times.
4. The cells are the first $$k + 1$$ barcodes.

The package draws the result this way. On the left are the rank curve and the cut; on the right,
the last pass's cumulative curve with its chord, and each earlier pass's knee in grey:

![UMI-tools' diagnostics on PBMC 1k: the rank curve with the cut at 2,868 UMIs, and the last pass's cumulative curve with its chord]({{ '/assets/explained/umitools-diagnostics.png' | relative_url }})

## Why it lands where it does

Each pass is a threshold. On the cumulative curve the farthest point from the chord is the last
barcode whose count exceeds the window's mean, here the mean of $$c_1, \dots, c_{m-1}$$. (The
chord starts at the first barcode, so that barcode is not in the mean.) So UMI-tools computes a
sequence of means, each over a window three times as long as the cells the last mean found:

![Each pass's mean drawn as a horizontal line over its window, on PBMC 1k and the spinal cord nuclei]({{ '/assets/explained/umitools-window.png' | relative_url }})

On PBMC 1k, the passes run:

| Pass | Window (barcodes) | Cells (k + 1) |
|---|---|---|
| 1 | 234,600, every barcode with a count | 1,565 |
| 2 | 4,692 | 1,126 |
| 3 | 3,375 | 1,097 |
| 4 | 3,288 | 1,090 |
| 5 | 3,267 | 1,089 |
| 6 | 3,264 | 1,089 |

The first mean is 43 UMIs. The error tail drags it down, so the first knee lands at the foot of
the cliff. From the second pass on, the window holds the cells and two to three times as many
barcodes below them. Its mean jumps to 1,985 UMIs, then 2,753, and the cut climbs to the top of
the cliff. The last window's mean is 2,846, and the last cell has 2,868 UMIs.

The loop stops at a fixed point of the map $$k \mapsto f(3k)$$, where $$f(m)$$ is the knee of
the first $$m$$ barcodes. In cobweb form:

![The knee after each pass on both samples, and the cobweb of the refit loop on PBMC 1k and on the spinal cord nuclei]({{ '/assets/explained/umitools-passes.png' | relative_url }})

On PBMC 1k the map is almost flat past the cliff: any window from 3,000 barcodes to all of them
gives a knee between 1,081 and 1,565. Wherever the loop starts, it lands at once near the cliff,
and it settles in six passes.

A clean step shows why this is robust. Suppose $$K$$ cells hold about $$C$$ UMIs each and the
barcodes below them almost none. Then the mean of the top $$3K$$ is about $$C/3$$. Every cell is
above it, nothing else is, and $$K$$ is a fixed point. The cut sits at a third of a typical cell,
deep inside the cliff, where no barcode lives.

## Where it goes wrong

On the spinal cord nuclei, the first pass's mean is 84 UMIs. That is below the ambient plateau,
about 200 UMIs high, so the first knee lands at the second cliff, 158,878 barcodes down, where
the empty droplets end.

The refits then slide back up. The next window, 476,631 barcodes, holds the nuclei, the ambient
plateau and part of the tail; its mean is 205, and the knee falls to 105,898. The next window holds more nuclei, its mean rises, and
the knee falls again. Once the window is mostly nuclei, the gentle bend never steepens enough to
hold it: the multiplier 3 needs a log-log slope of only 0.17 to stop, and the bend reaches that
far inside the nuclei. After 24 passes the loop settles at **3,554 cells, a quarter of 10x's
15,257**. Those are the nuclei above the mean count of the top 10,659 barcodes, 5,678 UMIs.

The [family page](index.md#moving-the-window-fixed-points) works through why: the map runs along
the diagonal on a gradual bend, so its single fixed point sits wherever the bend first steepens
past the critical slope. With a multiplier of 5, alevin-fry finds three fixed points on the same
curve, and lands on another; see [alevin-fry's knee](alevinfry.md).

## Details that matter

These move single barcodes or decide whether the rule runs at all. Each is on the
[method page](../../methods/umitools/distance.md).

- **What is counted.** UMI-tools counts reads, or UMIs with `--method umis`, from the first 100
  million reads. This package takes the counts from the matrix, and the same counts in the same
  order give the same barcodes
  ([differences](../../methods/umitools/distance.md#differences-from-the-original)).
- **Ties.** Barcodes are sorted with `Counter.most_common`, so equal counts keep the order in
  which they were first seen. A barcode tied with the last cell can fall outside the cells.
- **The first maximum wins** (numpy's `argmax`). alevin-fry takes the last.
- **A knee at the first barcode** is an error in UMI-tools, and here too.
- **`--set-cell-number N`** keeps the barcodes strictly above the $$(N+1)$$-th count. It still
  computes the knee, so it can still fail
  ([how it works](../../methods/umitools/distance.md#how-it-works)).
- **Before 1.0.0** the default rule was the density valley, which belongs to the density family
  ([UMI-tools density](../../methods/umitools/density.md)).

## Further reading

- Smith, Heger & Sudbery 2017, UMI-tools: modeling sequencing errors in Unique Molecular
  Identifiers to improve quantification accuracy, *Genome Res* 27:491–499,
  [doi:10.1101/gr.209601.116](https://doi.org/10.1101/gr.209601.116).
- Source: `getKneeEstimateDistance` in
  [`umi_tools/whitelist_methods.py`](https://github.com/CGATOxford/UMI-tools/blob/1.1.6/umi_tools/whitelist_methods.py)
  (1.1.6).
- The reference page: [UMI-tools whitelist, distance knee](../../methods/umitools/distance.md),
  with its releases, options and validation.
- [Same idea, different tools](../../overview/families.md#the-farthest-point-from-a-chord):
  UMI-tools beside the other chord rules.
