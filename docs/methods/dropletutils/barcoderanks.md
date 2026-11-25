---
title: barcodeRanks
parent: DropletUtils
grand_parent: Methods
nav_order: 1
method_ids: [dropletutils-barcoderanks]
status: implemented
---

# DropletUtils barcodeRanks
{: .no_toc }

`barcodeRanks()` computes the barcode-rank curve and two points on it: the **knee**, where the
curve turns from the plateau of cells into the cliff, and the **inflection**, the steepest point
of the cliff. It calls no cells itself; `emptyDrops` keeps every barcode at or above the knee
without testing it, and many pipelines use the knee or the inflection as a cut. Here a barcode is
a cell when its total is at least the chosen point.

The knee has had four definitions: the least curvature of a smoothing spline (1.0), the same with
the first barcodes skipped and an indexing bug fixed (1.10), the farthest point from a chord
(1.28), and the end of a sliding window traced along the curve (1.30).

1. TOC
{:toc}

## At a glance

| Id | `dropletutils-barcoderanks` |
|---|---|
| Tool | DropletUtils 1.0 and later (`barcodeRanks()`) |
| Category | knee |
| Counts | column sums of the matrix (UMIs) |
| Needs expected cells | no |
| Deterministic | yes |
| Original licence | GPL-3 |
| Validated against | `DropletUtils::barcodeRanks` in R, every variant (exact) |

## How it works

Every variant starts the same way:

1. Totals $$t_b$$ are the column sums; from 1.12, `as.integer(round(colSums(m)))`.
2. Sort the totals in decreasing order and collapse runs of equal totals. A run of length
   $$\ell$$ ending at rank $$c$$ gets the mid-rank $$r = c - (\ell - 1)/2$$.
3. Keep the runs with total above `lower` (100): $$x = \log_{10} r$$, $$y = \log_{10} t$$, in
   decreasing order of total. At least 3 runs are needed (2 from 1.30).

### Knee by spline curvature (1.0 – 1.26)

4. Numerical slopes $$d_i = (y_{i+1} - y_i)/(x_{i+1} - x_i)$$. From 1.10, the first
   $$s = \min(n_d - 1,\ \#\{x \le \log_{10} 50\})$$ slopes are skipped (`exclude.from = 50`).
5. The **right edge** is the least slope (the cliff); the **left edge** is the largest slope up
   to it (the plateau). Both are indices into the kept runs (plus the skipped count).
6. The **inflection** is $$10^{y}$$ at the right edge.
7. Fit `smooth.spline(x, y, df = 20)` to the runs from the left edge to the right edge (or, with
   `fit.bounds`, to the runs with $$y$$ strictly between the logs of the bounds).
8. At those runs, from the spline's derivatives, the signed curvature
   $$\kappa = y'' / (1 + y'^2)^{3/2}$$. The knee is $$10^{y}$$ at the least $$\kappa$$ (the first,
   if tied).
   - In 1.0 – 1.8 the position of the least curvature, counted within the fitted runs, is used
     as a position among *all* kept runs, so the knee lands (left edge − 1) runs too early: on a
     higher total.
   - From 1.10, fewer than 4 fitted runs give the knee at the first of them, with no fit.

`smooth.spline(df = 20)` searches the smoothing parameter for the fit with 20 equivalent degrees
of freedom (when there are more than 20 distinct x; otherwise it uses generalised
cross-validation). sc-droplets carries a port of R's own code for it: the knot rule, the B-spline
basis, the banded Cholesky solve, the leverages and R's golden-section search, which returns the
fit of the last smoothing parameter it evaluated.

### Knee by chord distance (1.28)

Steps 4 – 6 as above. Over the runs from the left edge to the right edge, draw the chord between
the first and the last; among the runs on or above the chord, the knee is the one farthest from
it:

$$\text{knee} = 10^{y_k},\quad k = \arg\max_{i:\ y_i \ge g x_i + c}\ \frac{|g x_i - y_i + c|}{\sqrt{g^2 + 1}}$$

with $$g$$ and $$c$$ the chord's slope and intercept. Fewer than 4 runs give the first of them.
`df` is ignored.

### Knee by curve tracing (1.30 on)

Runs with mid-rank at most `exclude.from` (50) are dropped before anything else.

1. Arc length along the curve: $$s_0 = 0$$, $$s_{i} = s_{i-1} + \sqrt{\Delta x^2 + \Delta y^2}$$.
2. A window of length `window` (1, in arc length) starts at every run with
   $$s_i + 1 \le s_{\max}$$. Its right end and its midpoint ($$s_i + 0.5$$) are interpolated
   linearly along the curve (between the runs around them; a point exactly on a run uses the
   segment before it).
3. For each window: the gap between its ends (straight-line distance), the slope $$g$$ of the
   line joining them, and whether the midpoint lies above that line.
4. An **elbow** window has its midpoint on or below the line and $$g < $$ `gradient.threshold`
   (−1). Only the windows up to the first elbow window are knee candidates; the inflection is
   searched among the windows that start at or before the midpoint of that first elbow window.
   Without an elbow window, every window counts.
5. The **knee** is $$10^{y}$$ at the midpoint of the candidate window (midpoint above its line)
   with the smallest gap; the **inflection** is $$10^{y}$$ at the midpoint of the window with the
   least slope. If no window has its midpoint above its line, the knee is the inflection. If no
   window fits on the curve, both are the last run.

### The call

A barcode is a cell when $$t_b \ge$$ the knee (`--use knee`, as `emptyDrops`' `retain`) or the
inflection (`--use inflection`). The knee and inflection are $$10^{\log_{10} t}$$, computed in
floating point; when that is not exactly $$t$$, the barcodes with total $$t$$ fall on one side of
the cut as they do in R.

## Versions

| DropletUtils | Variant | `--tool-version` | Knee | Notes |
|---|---|---|---|---|
| 1.0 – 1.8 | `v1.0` | `1.0` | spline curvature, mis-indexed | totals are unrounded column sums |
| 1.10 | `v1.10` | `1.10` | spline curvature | `exclude.from = 50`; fallback below 4 runs |
| 1.12 – 1.26 | `v1.12` | `1.12` | spline curvature | totals rounded to integers |
| 1.28 | `v1.28` | `1.28` | chord distance | `df` ignored |
| 1.30 – 1.34 | `v1.30` | `1.30` | curve tracing | `window`, `gradient.threshold`; `fit.bounds` and `df` ignored |

## Parameters

| Option | Default | Original | Meaning |
|---|---|---|---|
| `--lower` | 100 | `lower` | runs with total at or below this are not on the curve |
| `--exclude-from` | 50 (1.10 on) | `exclude.from` | barcodes ranked this high or higher are skipped |
| `--df` | 20 (up to 1.26) | `df` | degrees of freedom of the smoothing spline |
| `--fit-lower` | none (up to 1.28) | `fit.bounds[1]` | fit the runs with total above this instead of between the edges |
| `--fit-upper` | none (up to 1.28) | `fit.bounds[2]` | ... and below this; give both or neither |
| `--window` | 1 (1.30 on) | `window` | window length along the curve |
| `--gradient-threshold` | -1 (1.30 on) | `gradient.threshold` | slope below which a window with its midpoint under the line is an elbow |
| `--use` | `knee` | | which point is the cut: `knee` or `inflection` |

An option a release does not have is refused with that release's `--tool-version`.

## Run it

```bash
sc-droplets call dropletutils-barcoderanks raw_feature_bc_matrix.h5 --json
sc-droplets call dropletutils-barcoderanks raw_feature_bc_matrix.h5 --tool-version 1.26 --use inflection --barcodes cells.tsv
```

```python
import sc_droplets as scd

adata = scd.read("raw_feature_bc_matrix.h5")
result = scd.call(adata, "dropletutils-barcoderanks", tool_version="1.26")
result.thresholds["knee"], result.thresholds["inflection"]
```

The thresholds hold `knee`, `inflection` and `umis` (the cut used). `sc-droplets plot
diagnostics` draws the curve with the fitted spline (1.0 – 1.26), the chord (1.28) or the
windows (1.30).

## Differences from the original

- R's `smooth.spline` and `bvalue` are ported operation for operation and agree with R to the
  last bit. R linked to an optimised BLAS gives the same fits here (the banded solve uses only
  dot products of length 3 or less, which OpenBLAS sums in the same order).
- $$10^y$$ is computed correctly rounded, as R's `powl` on x86-64 gives it. It is often not the
  total again: R's knee for a run of total 1,852 is 1851.9999999999993, so those barcodes are
  above it. R on Apple silicon uses `pow`, which can differ in the last bit.
- `fit.bounds` is two options here (`--fit-lower`, `--fit-upper`).

## Validation

Compared with `DropletUtils::barcodeRanks` from bioconda (1.0.3 to 1.30.0, twelve releases;
1.28's function evaluated from its release branch inside 1.30) on 10x's pbmc_1k_v3 raw matrix
and on a simulated matrix: identical knee and inflection, to the last bit. The spline port
alone is compared with R's `smooth.spline` and `predict` on curves of 30 to 2,500 points:
identical smoothing parameter, fitted values and derivatives.

## References

- Lun et al. 2019, EmptyDrops: distinguishing cells from empty droplets in droplet-based
  single-cell RNA sequencing data, *Genome Biol* 20:63,
  [doi:10.1186/s13059-019-1662-y](https://doi.org/10.1186/s13059-019-1662-y).
- Source: [`R/barcodeRanks.R`](https://code.bioconductor.org/browse/DropletUtils/blob/RELEASE_3_22/R/barcodeRanks.R)
  (1.30); the earlier variants at `RELEASE_3_7`, `RELEASE_3_12`, `RELEASE_3_13` and
  `RELEASE_3_21`.
- R's `smooth.spline`: `src/library/stats/R/smspline.R`, `src/library/stats/src/sbart.c` and
  the Fortran it calls (R Core Team, GPL-2 or later).
