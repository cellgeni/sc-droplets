---
title: zUMIs' UIK
parent: Knee geometry
grand_parent: Explained
nav_order: 4
explains: [zumis-cell-barcodes]
---

# zUMIs' UIK
{: .no_toc }

zUMIs processes UMI-based single-cell data from many protocols. When it is given no list of
barcodes and no number of cells, it finds the cells itself. Since 2.5.6b it does so with the R
package inflection's *unit invariant knee* (UIK). In the family's terms:

- it draws the cumulative reads of the barcodes with at least 100 reads, against rank;
- it takes the point farthest above the chord between the curve's ends;
- every barcode with at least that point's count is a cell.

The name says what the [family page](index.md#the-chord-and-the-parallel-tangent) proves: the
farthest point from a chord does not move when either axis changes units. The rest follows from
the [family's identity](index.md#on-the-cumulative-curve-a-chord-is-a-mean): the cells are **the
barcodes above the mean count of the barcodes with at least 100 reads**.

1. TOC
{:toc}

## The rule

1. Keep the barcodes with at least `nReadsperCell` reads (100). Sort them in decreasing order,
   $$n_1 \ge \dots \ge n_N$$, and take the cumulative reads $$s_i$$, divided by 1,000.
2. inflection's `check_curve` decides whether the curve is concave. It looks at the signs of
   trapezoid integrals of the curve minus its chords, at rank 1, the three quartiles and $$N$$. A
   concave curve is negated.
3. The knee is the first minimum of the (negated) curve minus its chord: the point farthest
   above the chord.
4. The cells are the barcodes with at least the knee's count. If fewer than 10 pass, zUMIs takes
   the top 100 instead.

## Why it lands where it does

Two of the steps do less than they seem to.

- **Dividing by 1,000** is a change of units on one axis. The farthest point from a chord is the
  point of the largest triangle, and every triangle's area is scaled by the same factor, so
  nothing moves. That is the "unit invariance" of the name.
- **`check_curve`** always says "concave" here. The cumulative sum of sorted counts is concave by
  construction, since its slope is the count itself, which only falls. So the negation always
  happens, and the knee is always the point farthest *above* the chord.

What remains is the identity: the cells are the barcodes whose count exceeds the mean of the
window, the first barcode excluded. The window is set by the floor. The script behind these pages
checks the identity against zUMIs' rule (ported, and validated against zUMIs and inflection in R)
at floors of 30, 100 and 1,000, on both samples.

| | PBMC 1k | spinal cord nuclei |
|---|---|---|
| barcodes with at least 100 | 1,367 | 158,149 |
| their mean | 6,761 | 614 |
| cells | 594 | 15,092 |
| 10x's cells | 1,222 | 15,257 |

The floor of 100 decides what the mean is a mean of:

- **On PBMC 1k** the floor sits well above the empty droplets, which hold about 13 UMIs. The
  window is almost only cells, so the rule keeps the cells above the average cell: 594, about
  half.
- **On the nuclei** the empty droplets hold about 200 UMIs, above the floor. The window is the
  nuclei plus the whole ambient plateau, the mean of 614 falls inside the cliff, and the call
  lands 165 cells short of 10x's.

![zUMIs' number of cells as its floor moves, on PBMC 1k and the spinal cord nuclei]({{ '/assets/explained/zumis-floor.png' | relative_url }})

The figure moves the floor:

- **On the nuclei**, every floor from 4 to 200 gives about 15,100 cells. Those floors differ only
  by a few thousand barcodes on the second cliff, which barely move the mean. Below 4 the error
  tail comes in, and the answer jumps to the second cliff. Above the ambient plateau the window
  is nuclei only: 6,733 cells at a floor of 300, 3,989 at 1,000.
- **On PBMC 1k** the answer falls steadily as the floor rises: 1,258 cells at 10, 1,071 at 20,
  711 at 50, 594 at 100. Every floor above the ambient level leaves a window of cells, and their
  mean cuts them roughly in half.

The rule works when the floor admits the empty droplets and rejects the error tail. It fails when
the floor is above the empty droplets.

## Where it goes wrong

- **Whole cells with little ambient RNA.** Half of PBMC 1k's cells fall below the mean of the
  cells. The [comparison on public data](../../overview/comparison.md) finds zUMIs calling about
  half of 10x's cells on each of its whole-cell samples (Jaccard 0.39 to 0.57 with 10x's
  call), and a Jaccard of 0.93 on the nuclei.
- **Reads against UMIs.** zUMIs applies its floor to reads, and there are several reads per UMI.
  A floor of 100 reads therefore sits lower on the curve than a floor of 100 UMIs, and lets in
  more of the empty droplets. This package applies the floor to what it counts: UMIs by default,
  or reads with `--unit reads`. The PBMC 1k numbers here are for UMIs. On reads the window would
  hold more empty droplets, and the call would rise.
- **Fewer than 10 cells.** zUMIs falls back to the top 100 barcodes.

## Details that matter

Each is on the [method page](../../methods/zumis/cell-barcodes.md).

- **Before 2.5.6b** zUMIs fitted a Gaussian mixture with `mclust`, initialised on an unseeded
  random subset, so it cannot be reproduced. It belongs to the mixture family
  ([versions](../../methods/zumis/cell-barcodes.md#versions)).
- **`check_curve`'s quartiles** are R's type-7 quantiles of $$1, \dots, N$$. They matter only for
  curves that are not concave.
- **Ties** at the knee's count are all cells, because the cut is a count.
- **`barcode_file` and `barcode_num`**, a list or a number given in zUMIs' configuration,
  replace the detection. They are not rules, and are not offered.

## Further reading

- Parekh, Ziegenhain, Vieth, Enard & Hellmann 2018, zUMIs – A fast and flexible pipeline to
  process RNA sequencing data with UMIs, *GigaScience* 7(6),
  [doi:10.1093/gigascience/giy059](https://doi.org/10.1093/gigascience/giy059).
- Christopoulos 2016, Introducing Unit Invariant Knee (UIK) as an objective choice for elbow
  point in multivariate data analysis techniques, *SSRN*,
  [doi:10.2139/ssrn.3043076](https://doi.org/10.2139/ssrn.3043076).
- Source: [`barcodeIDFUN.R`](https://github.com/sdparekh/zUMIs/blob/2.9.7/barcodeIDFUN.R)
  (`cellBC`, `.FindBCcut`), and inflection's `uik` and `check_curve`.
- The reference page: [zUMIs cell barcodes](../../methods/zumis/cell-barcodes.md).
- [dropbead's diagonal](dropbead.md): the same mean, over a fixed number of barcodes instead of a
  floor.
