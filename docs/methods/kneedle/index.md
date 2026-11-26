---
title: Kneedle
parent: Methods
nav_order: 6
has_children: true
---

# Kneedle (kneed)
{: .no_toc }

Kneedle (Satopää et al. 2011) finds the knee of a curve as the point where it departs most from
the straight line between its ends, after scaling both axes to [0, 1]. It is not a single-cell
method; pipelines apply it to the barcode-rank curve through `kneed`, its Python implementation.
Sanger's yascp pipeline, for example, runs it on rank against UMIs to estimate the expected cells
it gives CellBender, and cellsweep 1.0.2 uses it for its knee (see its own page when it lands).

1. TOC
{:toc}

## Which release does what

`kneed` changed how it finds the knee several times. The releases that change the answer:

| kneed | What changed |
|---|---|
| 0.5.1 (Nov 2019) | local maxima and minima of the difference curve include plateaus (`>=`, `<=`) |
| 0.7.0 (Aug 2020) | the polynomial fit uses `numpy.polyfit` (scikit-learn before) and takes a degree; a convex curve is no longer mirrored in x before the search |
| 0.8.6 (Mar 2026) | after a local minimum, no knee is detected until the next local maximum |
| 0.9.0 (Oct 2026) | no interpolation step; a concave decreasing or convex increasing curve is mirrored in x |

Releases before 0.5.1 (2017–2019) resampled the curve on an even grid and compared strictly; they
are not offered.

## Algorithms

| Page | Id | Rule |
|---|---|---|
| [kneedle](kneedle.md) | `kneedle` | Kneedle on the barcode-rank curve, with yascp's search for the sensitivity by default |
