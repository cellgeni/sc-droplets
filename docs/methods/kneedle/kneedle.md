---
title: kneedle
parent: Kneedle
grand_parent: Methods
nav_order: 1
method_ids: [kneedle]
status: implemented
---

# Kneedle on the barcode-rank curve
{: .no_toc }

`kneed.KneeLocator` run on rank against count, as pipelines run it. By default this is what
Sanger's yascp does (`get_estimates_from_umi_counts.py`): the barcodes with at least 1,000 UMIs,
linear axes, a convex decreasing curve, and the sensitivity $$S$$ lowered from 5,000 in steps of
100 until a knee appears. The cells are the barcodes ranked up to the knee.

Why the rule lands where it does, drawn on real data:
[Kneedle's difference curve](../../explained/knee/kneedle.md), in *Explained*.

1. TOC
{:toc}

## At a glance

| Id | `kneedle` |
|---|---|
| Tool | kneed 0.5.1 and later (`KneeLocator`); yascp's use of it |
| Category | knee |
| Counts | UMIs, or reads with `--unit reads` |
| Needs expected cells | no |
| Deterministic | yes |
| Original licence | BSD-3-Clause (kneed); GPL-3 (yascp) |
| Validated against | `kneed` 0.6.0, 0.8.5, 0.8.6 and 0.9.0 (exact) |

## How it works

### The curve

Barcodes are sorted by decreasing count, ties in input order, and ranked $$r = 1, 2, \dots$$. The
barcodes with count $$\ge$$ `--lower` (1,000) make the curve: $$x = r$$ and $$y$$ = the count, or
their $$\log_{10}$$ with `--log-x` and `--log-y`. At least two points are needed.

### Kneedle

With $$S$$ the sensitivity, and `--curve` and `--direction` describing the curve's shape (convex
and decreasing for a linear rank curve):

1. **Smoothing.** `interp1d`: the linear interpolant of the points, evaluated at their own $$x$$
   (the points themselves, up to rounding); from 0.9.0 the points themselves. `polynomial`: the
   least-squares polynomial of degree 7 (`--polynomial-degree` from 0.7.0) evaluated at $$x$$;
   scikit-learn's `LinearRegression` on `PolynomialFeatures` before 0.7.0, `numpy.polyfit` from
   0.7.0.
2. **Scaling.** $$x_n = (x - \min x)/(\max x - \min x)$$, and $$y_n$$ likewise from the smoothed
   values.
3. **Turning the curve into a concave increasing one.**
   - From 0.7.0, $$y_n$$ only: reversed for concave decreasing; $$\max y_n - y_n$$ for convex
     decreasing; reversed $$\max y_n - y_n$$ for convex increasing. From 0.9.0, for concave
     decreasing and convex increasing, $$x_n$$ also becomes the reverse of $$1 - x_n$$.
   - Before 0.7.0: for a convex curve $$x_n \to \max x_n - x_n$$ and $$y_n \to \max y_n - y_n$$; for
     a decreasing curve $$y_n$$ is reversed; for a convex curve both are then reversed.
4. **Difference curve** $$d_i = y_{n,i} - x_{n,i}$$.
5. **Local maxima** of $$d$$: points at least as high as both neighbours; **local minima**: at
   most as high as both (scipy's `argrelextrema`, which compares an end with itself, so an end
   can be one).
6. **Thresholds**: at the $$k$$-th maximum, $$T_k = d_{\max,k} - S \cdot |\overline{\Delta x_n}|$$,
   with $$\overline{\Delta x_n}$$ the mean step of $$x_n$$.
7. **Walk** from the first maximum to the second-last point. At a maximum, the threshold becomes
   its $$T_k$$; at a minimum, 0 (and from 0.8.6, nothing is detected until the next maximum).
   When the next point's $$d$$ is below the threshold, the last maximum is a knee. Without
   `--online`, the first knee is the answer; with it, the last one found.
8. **Knee point.** For a convex decreasing or concave increasing curve, the point at the last
   maximum's index; for the other two shapes, the point at the mirrored index (counted from the end).

No maximum, or no knee: none is found.

### The search and the call

Without `--sensitivity`, $$S$$ starts at `--s-start` (5,000) and is lowered by `--s-step` (100)
until Kneedle finds a knee, as yascp does. Below 0 every threshold is above its maximum, so the
first maximum the walk reaches is a knee whatever $$S$$ is: if there is still none, there never
will be, and the search stops with an error (yascp's would not end). With
`--sensitivity S`, only that $$S$$ is tried. The cells are the barcodes ranked up to the knee's
rank: yascp's estimate of the number of cells, as a set of barcodes.

### Arithmetic

The interpolation, polynomial fits, local extrema and means use the same numpy, scipy and
scikit-learn calls as kneed, so given the same libraries the knee is the same. The degree-7
polynomial on ranks in the thousands is ill-conditioned: across numpy or LAPACK builds its fitted
values can differ beyond the last bits, and so, rarely, can the knee.

## Versions

| kneed | Variant | `--tool-version` | Behaviour |
|---|---|---|---|
| before 0.5.1 | — | — | not offered |
| 0.5.1 – 0.6.0 | `v0.5.1` | `0.5.1` | plateaus are extrema; polynomial by scikit-learn; convex curves mirrored in x |
| 0.7.0 – 0.8.5 | `v0.7` | `0.7.0` | polynomial by `numpy.polyfit`, degree settable; no mirroring in x |
| 0.8.6 | `v0.8.6` | `0.8.6` | no detection between a minimum and the next maximum |
| 0.9.0 | `v0.9` | `0.9.0` | no interpolation; concave decreasing and convex increasing curves mirrored in x |

## Parameters

| Option | Default | Original | Meaning |
|---|---|---|---|
| `--lower` | 1000 | yascp `--lower_bound_umis_estimate_ncells` | barcodes with fewer counts are not on the curve |
| `--sensitivity` | none (search) | `S` | one sensitivity, instead of yascp's search |
| `--s-start` | 5000 | yascp | the search's first sensitivity |
| `--s-step` | 100 | yascp | how much the search lowers it each time |
| `--curve` | `convex` | `curve` | `convex` or `concave` |
| `--direction` | `decreasing` | `direction` | `decreasing` or `increasing` |
| `--interp-method` | `interp1d` | `interp_method` | `interp1d` or `polynomial` |
| `--polynomial-degree` | 7 | `polynomial_degree` (0.7.0+) | the polynomial's degree |
| `--online` / `--no-online` | off | `online` | keep looking after the first knee; the last one wins |
| `--log-x` / `--no-log-x` | off | | rank on a log10 axis |
| `--log-y` / `--no-log-y` | off | | counts on a log10 axis |
| `--unit` | `umis` | | count `umis` or `reads` (needs `--reads`) |

## Run it

```bash
sc-droplets call kneedle raw_feature_bc_matrix.h5 --json
sc-droplets call kneedle raw_feature_bc_matrix.h5 --tool-version 0.8.5 --interp-method polynomial
sc-droplets call kneedle raw_feature_bc_matrix.h5 --log-x --log-y --curve concave --sensitivity 1 --lower 100
```

```python
import sc_droplets as scd

adata = scd.read("raw_feature_bc_matrix.h5")
result = scd.call(adata, "kneedle", tool_version="0.8.5")
result.thresholds["knee_rank"], result.thresholds["sensitivity"]
```

The thresholds hold `knee_rank`, the count at the knee in the unit used, and the sensitivity that
found it. `sc-droplets plot diagnostics` draws the rank curve with the cells, and the difference
curve with its maxima, their thresholds and the knee.

## Differences from the original

- yascp takes the number of cells; this calls the barcodes ranked up to it (ties at the knee's
  count in input order).
- When no sensitivity can give a knee, yascp's search never ends; this stops with an error once
  the sensitivity is below 0, where every value gives the same answer.
- Releases before 0.5.1 are not offered.

## Validation

Compared with `kneed.KneeLocator` itself in 0.6.0, 0.8.5, 0.8.6 and 0.9.0 (each in its own
environment from conda-forge, with the numpy, scipy and scikit-learn of its time), running yascp's
search on the per-barcode UMI totals of 10x's pbmc_1k_v3 and pbmc4k raw matrices: linear axes
with `interp1d` and with `polynomial`, and log-log axes with a fixed sensitivity, concave and
convex, online and not. Identical knees.

## References

- Satopää, Albrecht, Irwin & Raghavan 2011, Finding a "Kneedle" in a haystack: detecting knee
  points in system behavior, *31st International Conference on Distributed Computing Systems
  Workshops*, 166–171, [doi:10.1109/ICDCSW.2011.20](https://doi.org/10.1109/ICDCSW.2011.20).
- Source: [`kneed/knee_locator.py`](https://github.com/arvkevi/kneed/blob/v0.9.0/kneed/knee_locator.py)
  (0.9.0) and the releases on [PyPI](https://pypi.org/project/kneed/#history).
- yascp: [`bin/get_estimates_from_umi_counts.py`](https://github.com/wtsi-hgi/yascp/blob/main/bin/get_estimates_from_umi_counts.py)
  (`kneedle_cutoff`).
