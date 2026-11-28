"""alevin's initial whitelist (salmon 0.11.0 - 1.10.3): the knee with its density correction,
expected cells or forced cells, then the low-confidence band.

Page: docs/methods/alevin/initial-whitelist.md. Ported from salmon's `src/Alevin.cpp`
(`getLeftBoundary`, `gaussianKDE`, `sampleTrueBarcodes`, each release where they changed), which
salmon distributes under the GNU General Public License v3, as sc-droplets is.

The arithmetic follows alevin's: the cumulative sums in double, sequentially; the C library's
`log` and `log10`; the kernel density summed barcode by barcode in alevin's order, on a grid built
by repeated subtraction. Ties at the knee fall in the matrix's order (alevin's depend on its hash
map and the C++ library's sort; see the page).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional

import numpy as np
import pandas as pd

from sc_droplets.methods._common import result
from sc_droplets.primitives.rmath import clog, clog10
from sc_droplets.registry import Method, MethodInfo, Reference, Variant, param, register

TOOL = "alevin"
LATEST = "1.10.3"
CLASSES = ["cell", "low_confidence", "background"]
MAX_NUM_BARCODES = 100_000
LOW_REGION_MIN = 200
LOW_REGION_MAX = 1000
FREQ_THRESHOLD = 10
GRID = 10_000
KEEP_CAP_OLD = 737_000  # --keepCBFraction's cap before 0.14.1
ORDER = ("v0.11", "v0.12", "v0.13", "v0.14", "v0.14.1", "v1.2")


class AlevinError(ValueError):
    """Where alevin itself stops (an error, or reading past its arrays)."""


def _from(variant: str, first: str) -> bool:
    return ORDER.index(variant) >= ORDER.index(first)


def left_boundary(f: np.ndarray, top: int, shifted: bool = False) -> int:
    """`getLeftBoundary`: the tangent from the origin to the log cumulative curve of the top
    `top` counts (`f` sorted in decreasing order). Returns the number of cells, or 0.
    `shifted`: 0.11.x's sum, one barcode further down (`sortedIdx[topx - i]`)."""
    if shifted:
        if top >= f.size:
            raise AlevinError("alevin 0.11 reads one barcode past the top --maxNumBarcodes, past the end here")
        values = f[1:top + 1][::-1]
    else:
        values = f[:top][::-1]
    c = clog(np.cumsum(values.astype(np.float64)))  # double cumCount += f; freqs[i] = log(cumCount)
    m = c.size
    if m < 2:
        return 0
    i = np.arange(m, dtype=np.float64)
    ratio = np.full(m, np.inf)
    ratio[1:] = c[1:] / i[1:]
    # the smallest ratio after each j: a j whose ratio is clearly above it has a point below its line
    after = np.minimum.accumulate(ratio[::-1])[::-1]
    after = np.append(after[1:], np.inf)
    for j in np.flatnonzero((c != c[0]) & (ratio <= after * (1 + 1e-12))):
        slope = c[j] / float(j)
        if not np.any(i[j + 1:] * slope > c[j + 1:]):  # lineY > curveY: the curve dips below
            cutoff = top - int(j)
            above = np.flatnonzero(f[:cutoff] != f[cutoff])  # step left past the barcodes tied with it
            if above.size == 0:
                raise AlevinError("every barcode up to the knee has the same count (alevin reads before its array)")
            return int(above[-1])
    return 0


@dataclass
class Density:
    """What the density correction saw: its grid (log10 counts), density, minima and choice."""

    grid: np.ndarray
    density: np.ndarray
    minima: np.ndarray  # grid indices, increasing (decreasing count)
    boundaries: np.ndarray  # barcodes at or above each minimum's count, as far as tried
    cells: int  # the prediction, 0 if none
    accepted: bool


def gaussian_kde(f: np.ndarray, knee: int) -> Density:
    """`gaussianKDE`: a Gaussian kernel density of log10 counts (variance 0.01 x the sample
    variance) on 10,000 points; the first minimum whose barcodes number between 0.1 x knee and the
    knee replaces it."""
    threshold = 0.001 * float(f[0])
    kept = f[f.astype(np.float64) > threshold]  # sorted, so the loop's break keeps a prefix
    logs = clog10(kept.astype(np.float64))
    n = logs.size
    mean = np.cumsum(logs)[-1] / n
    deviation = logs - mean
    covariance = (np.cumsum(deviation * deviation)[-1] * 0.01) / (n - 1.0)
    if covariance == 0:
        raise AlevinError("0 Covariance error for Gaussian kde")
    inverse = 1.0 / covariance
    decrement = (logs[0] - logs[-1]) / float(GRID)
    grid = np.empty(GRID)
    pred = logs[0]
    for j in range(GRID):
        grid[j] = pred
        pred -= decrement
    density = np.zeros(GRID)
    for value in logs:  # alevin sums barcode by barcode
        diff = value - grid
        density += np.exp(-((diff * diff * inverse) / 2.0))
    inner = density[1:-1]
    minima = np.flatnonzero((density[:-2] > inner) & (inner < density[2:])) + 1
    boundaries = []
    for t in minima:
        tau = math.pow(10.0, logs[0] - (float(t) * decrement))
        boundary = int(np.sum(f >= tau))
        boundaries.append(boundary)
        if boundary > knee:
            return Density(grid, density, minima, np.array(boundaries), boundary, False)
        if knee * 0.1 > boundary:
            continue
        return Density(grid, density, minima, np.array(boundaries), boundary, True)
    return Density(grid, density, minima, np.array(boundaries), 0, False)


def expected_cells(f: np.ndarray, expect: int, plus_one: bool) -> int:
    """Cell Ranger 2's rule as alevin has it: the first rank in [r, min(n, 2E)) below a tenth of
    the count at rank r = max(1, 0.01 E)."""
    r = max(1, int(expect * 0.01))
    if r >= f.size:
        raise AlevinError(f"--expectCells {expect}: rank {r} is past the {f.size} barcodes")
    cut = max(1.0, float(f[r]) * 0.1)
    most = min(f.size, expect * 2)
    below = np.flatnonzero(f[r:most] < cut)
    cells = most if below.size == 0 else r + int(below[0]) + (1 if plus_one else 0)
    if cells == most:
        raise AlevinError("Can't find right Boundary (alevin's --expectCells rule)")
    return cells


def forced_cells(f: np.ndarray, force: int, variant: str, freq_threshold: int) -> int:
    """`--forceCells` as each release has it."""
    if not _from(variant, "v0.14"):
        return force
    k = force if variant == "v0.14" else force - 1
    if k >= f.size or k < 0:
        raise AlevinError(f"--forceCells {force} is past the {f.size} barcodes")
    if f[k] < freq_threshold:
        k = int(np.sum(f >= freq_threshold)) - 1  # lowered until f[k] reaches the floor
        if k < 0:
            raise AlevinError(f"no barcode has --freqThreshold {freq_threshold} reads")
    return k + 1 if _from(variant, "v1.2") else k


def whitelist(f: np.ndarray, k: int, variant: str, low_min: int):
    """The low-confidence band and the whitelist's size. Returns (whitelist size, alevin's
    reported low-confidence count)."""
    n = f.size
    low = int(0.5 * k)
    low = low_min if low < low_min else LOW_REGION_MAX if low > LOW_REGION_MAX else low
    if _from(variant, "v0.14.1") and low + k + 1 > n:
        low = n - k - 1
    t = k + low
    if t >= n or t < 0:
        raise AlevinError(f"the band ends past the {n} barcodes (alevin reads past its array)")
    cutoff = f[t]
    if _from(variant, "v0.14.1"):
        nearest = cutoff
        while nearest == cutoff and low != 0:
            t -= 1
            if t < 0:
                raise AlevinError("every barcode down to the band's end has the same count")
            nearest = f[t]
            if low > low_min:
                low -= 1
        if low != 0:
            low += 1
            t += 1
        return t, low
    above = int(np.sum(f > cutoff))
    if above == 0:
        raise AlevinError("every barcode down to the band's end has the same count")
    return above, above - k


@register
class InitialWhitelist(Method):
    info = MethodInfo(
        id="alevin-initial-whitelist",
        title="alevin initial whitelist",
        tool=TOOL,
        category="knee",
        summary=(
            "alevin's initial whitelist: a tangent knee on the log cumulative curve corrected by a "
            "kernel-density valley, Cell Ranger 2's expected cells, or forced cells; then a band of "
            "low-confidence barcodes."
        ),
        docs="methods/alevin/initial-whitelist.md",
        counts="either",
        variants=(
            Variant("v0.11", since="0.11.0", until="0.12.0", summary="the knee only, off by one barcode"),
            Variant("v0.12", since="0.12.0", until="0.13.0", summary="knee fixed; expected and forced cells"),
            Variant("v0.13", since="0.13.0", until="0.14.0", summary="adds --keepCBFraction"),
            Variant("v0.14", since="0.14.0", until="0.14.1", summary="forced cells below the read floor dropped"),
            Variant("v0.14.1", since="0.14.1", until="1.2.0", summary="forced cells from F - 1; the band stops at the end"),
            Variant("v1.2", since="1.2.0", summary="forced cells F; expected cells one more"),
        ),
        latest_known=LATEST,
        references=(
            Reference(
                "Srivastava, Malik, Smith, Sudbery & Patro 2019, Alevin efficiently estimates accurate gene "
                "abundances from dscRNA-seq data, Genome Biology 20:65",
                doi="10.1186/s13059-019-1670-y",
            ),
        ),
        source="https://github.com/COMBINE-lab/salmon/blob/v1.10.3/src/Alevin.cpp",
        original_license="GPL-3.0",
    )

    @dataclass
    class Params:
        expect_cells: Optional[int] = param(None, "Cell Ranger 2's rule from expected cells", original="--expectCells", minimum=1)
        force_cells: Optional[int] = param(None, "Force this many cells", original="--forceCells", minimum=1)
        keep_cb_fraction: Optional[float] = param(None, "Force this fraction of the barcodes seen", original="--keepCBFraction",
                                                  minimum=0, maximum=1)
        freq_threshold: int = param(FREQ_THRESHOLD, "Forced cells need at least this count (0.14.0 on)", original="--freqThreshold", minimum=0)
        max_num_barcodes: int = param(MAX_NUM_BARCODES, "The knee looks at this many barcodes at most", original="--maxNumBarcodes",
                                      minimum=2)
        low_region_min_num_barcodes: int = param(LOW_REGION_MIN, "The low-confidence band's least width",
                                                 original="--lowRegionMinNumBarcodes", minimum=0)
        unit: str = param("umis", "Count UMIs (the matrix) or reads (--reads)", original="reads per barcode", choices=("umis", "reads"))

    def run(self, droplets, params, variant, random_state):
        name = variant.name
        if name == "v0.11" and (params.expect_cells or params.force_cells or params.keep_cb_fraction):
            raise AlevinError("alevin before 0.12.0 has only the knee (use --tool-version 0.12.0 or later)")
        if params.keep_cb_fraction and not _from(name, "v0.13"):
            raise AlevinError("alevin before 0.13.0 has no --keepCBFraction (use --tool-version 0.13.0 or later)")
        values = droplets.values(params.unit)
        index = np.flatnonzero(values > 0)
        counts = values[index]
        if not np.all(counts == np.floor(counts)):
            raise AlevinError("alevin counts whole reads; these counts have fractions")
        counts = counts.astype(np.int64)
        order = np.argsort(-counts, kind="stable")  # ties in matrix order (see the page)
        f = counts[order]
        n = f.size
        force = params.force_cells
        if params.keep_cb_fraction:  # alevin turns the fraction into forced cells
            cap = params.max_num_barcodes if _from(name, "v0.14.1") else KEEP_CAP_OLD
            force = min(int(params.keep_cb_fraction * n), cap)
        thresholds = {}
        diagnostics = {"unit": params.unit}
        if force:
            k = forced_cells(f, force, name, params.freq_threshold)
        elif params.expect_cells:
            k = expected_cells(f, params.expect_cells, plus_one=_from(name, "v1.2"))
        else:
            top = min(params.max_num_barcodes, n)
            k = left_boundary(f, top, shifted=name == "v0.11")
            if k == 0:
                raise AlevinError("Can't find left Boundary")
            thresholds["knee_rank"] = float(k)
            found = gaussian_kde(f, k)
            if found.cells:
                thresholds["gauss_rank"] = float(found.cells)
            if found.accepted:
                k = found.cells
            chosen = found.minima[len(found.boundaries) - 1] if found.accepted else None
            diagnostics.update({
                "grid": found.grid, "density": found.density,
                "valleys": found.grid[found.minima], "valley_passing": np.asarray(found.boundaries, dtype=np.int64),
                "density_title": "Density of log10 counts (alevin's correction)",
                **({"chosen": float(found.grid[chosen])} if chosen is not None else {}),
            })
        size, reported_low = whitelist(f, k, name, params.low_region_min_num_barcodes)
        cells = min(k, size)
        classes = np.full(values.shape[0], CLASSES.index("background"), dtype=np.int8)
        classes[index[order[cells:size]]] = CLASSES.index("low_confidence")
        classes[index[order[:cells]]] = CLASSES.index("cell")
        thresholds.update({
            params.unit: float(f[cells - 1]) if cells else float("nan"),
            "cells": float(k),
            "whitelist": float(size),
            "low_confidence": float(reported_low),
        })
        return result(
            type(self), variant, params, droplets, classes == CLASSES.index("cell"),
            classes=pd.Categorical.from_codes(classes, categories=CLASSES),
            thresholds=thresholds, diagnostics=diagnostics,
        )
