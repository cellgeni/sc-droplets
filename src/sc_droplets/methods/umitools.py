"""UMI-tools `whitelist`: the distance knee (1.0.0+) and the density valley (0.5.0+).

Pages: docs/methods/umitools/distance.md, docs/methods/umitools/density.md.

Ported from UMI-tools' `getKneeEstimateDistance` and `getKneeEstimateDensity`
(`umi_tools/whitelist_methods.py`, 1.1.6; `umi_methods.getKneeEstimate`, 0.5),
under the MIT licence:

    Copyright (c) 2015 CGAT

    Permission is hereby granted, free of charge, to any person obtaining a copy of
    this software and associated documentation files (the "Software"), to deal in
    the Software without restriction, including without limitation the rights to
    use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
    the Software, and to permit persons to whom the Software is furnished to do so,
    subject to the following conditions: The above copyright notice and this
    permission notice shall be included in all copies or substantial portions of
    the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND.

The arithmetic follows UMI-tools' numpy and scipy calls one for one, so the same
counts in the same order give the same barcodes.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Tuple

import numpy as np

from sc_droplets.methods._common import result
from sc_droplets.registry import Method, MethodInfo, Reference, Variant, param, register

TOOL = "UMI-tools"
LATEST = "1.1.6"
_SOURCE = "https://github.com/CGATOxford/UMI-tools/blob/1.1.6/umi_tools/whitelist_methods.py"
_PAPER = Reference(
    "Smith, Heger & Sudbery 2017, UMI-tools: modeling sequencing errors in Unique Molecular "
    "Identifiers to improve quantification accuracy, Genome Res 27:491-499",
    doi="10.1101/gr.209601.116",
)
#: UMI-tools gives up on a knee that keeps moving after this many passes.
MAX_PASSES = 100
#: Points on the density's grid.
GRID = 10000


class UmiToolsError(ValueError):
    """Where UMI-tools itself stops with an error."""


def _observed(values: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """The barcodes UMI-tools would see (count > 0), as indices, and their counts.

    Counts are 64-bit integers when they are whole numbers, as UMI-tools' are.
    """
    index = np.flatnonzero(values > 0)
    counts = values[index]
    if np.all(counts == np.floor(counts)) and (counts.size == 0 or counts.max() < 2.0**63):
        counts = counts.astype(np.int64)
    return index, counts


def _cut_above_nth(counts: np.ndarray, n_cells: int) -> float:
    """`--set-cell-number N`: the (N+1)-th largest count; cells are strictly above it."""
    if n_cells >= counts.shape[0]:
        raise UmiToolsError(
            f"--n-cells {n_cells} is not below the number of barcodes with counts ({counts.shape[0]}); "
            "UMI-tools takes the (N+1)-th largest count as the cut"
        )
    return np.sort(counts)[::-1][n_cells]


# ---------------------------------------------------------------------------
# Distance
# ---------------------------------------------------------------------------


def knee_distance(values: np.ndarray) -> Tuple[np.ndarray, int]:
    """Distances of (i, values[i]) from the chord of its ends, and the first farthest i.

    UMI-tools' `getKneeDistance`, operation for operation: `np.sum(x, axis=1)` over
    two columns is one addition, so the result does not depend on memory order.
    """
    n = values.shape[0]
    coords = np.vstack((np.arange(n), values)).T
    first = coords[0]
    line = coords[-1] - coords[0]
    with np.errstate(invalid="ignore", divide="ignore"):
        unit = line / np.sqrt(np.sum(line * line))
        from_first = coords - first
        along = np.sum(from_first * unit, axis=1)
        to_line = from_first - np.outer(along, unit)
        distance = np.sqrt(np.sum(to_line * to_line, axis=1))
    return distance, int(np.argmax(distance))


def distance_knee(counts_sorted: np.ndarray) -> Tuple[int, List[Tuple[int, int]]]:
    """The knee index k (cells are the first k + 1) and each pass as (points, knee)."""
    n = counts_sorted.shape[0]
    if n == 0:
        raise UmiToolsError("no barcode has a count; UMI-tools needs at least one")
    cumulative = np.cumsum(counts_sorted)
    _, knee = knee_distance(cumulative)
    passes = [(n, knee)]
    if knee == 0:
        raise UmiToolsError("the distance knee is the first barcode (UMI-tools: \"Something's gone wrong here!!\")")
    previous = 0
    iterations = 0
    while knee - previous != 0:
        previous = knee
        iterations += 1
        if iterations > MAX_PASSES:
            break
        points = min(knee * 3, n)
        if points == 0:
            raise UmiToolsError("the distance knee fell to the first barcode; UMI-tools fails here too")
        _, knee = knee_distance(cumulative[:points])
        passes.append((points, knee))
    return knee, passes


@register
class Distance(Method):
    info = MethodInfo(
        id="umitools-distance",
        title="UMI-tools whitelist, distance knee",
        tool=TOOL,
        category="knee",
        summary=(
            "The point of the cumulative count curve farthest from the chord joining its ends, "
            "searched again on three times the knee's rank until it settles."
        ),
        docs="methods/umitools/distance.md",
        counts="either",
        variants=(
            Variant(
                "v1.0",
                since="1.0.0",
                summary="whitelist --knee-method=distance (the default from 1.0.0)",
            ),
        ),
        latest_known=LATEST,
        references=(_PAPER,),
        source=_SOURCE,
        original_license="MIT",
    )

    @dataclass
    class Params:
        n_cells: Optional[int] = param(
            None,
            "Call the barcodes above the (N+1)-th largest count instead of the knee",
            original="--set-cell-number",
            minimum=1,
        )
        unit: str = param("umis", "Count UMIs (the matrix) or reads (--reads)", original="--method", choices=("umis", "reads"))

    def run(self, droplets, params, variant, random_state):
        values = droplets.values(params.unit)
        index, counts = _observed(values)
        order = np.argsort(-counts, kind="stable")
        knee, passes = distance_knee(counts[order])
        is_cell = np.zeros(values.shape[0], dtype=bool)
        thresholds = {"knee_rank": float(knee + 1)}
        if params.n_cells is not None:
            cut = _cut_above_nth(counts, params.n_cells)
            is_cell[index[counts > cut]] = True
            thresholds[params.unit] = float(cut)
        else:
            is_cell[index[order[: knee + 1]]] = True
            thresholds[params.unit] = float(counts[order[knee]])
        return result(
            type(self),
            variant,
            params,
            droplets,
            is_cell,
            thresholds=thresholds,
            diagnostics={"passes": np.asarray(passes, dtype=np.int64), "unit": params.unit,
                         **_last_pass_curve(np.cumsum(counts[order]), passes)},
        )


def _last_pass_curve(cumulative: np.ndarray, passes: List[Tuple[int, int]], most: int = 2000) -> dict:
    """The last pass's cumulative curve, thinned to about `most` points (its ends and every
    knee kept), so the plot needs no barcodes: a call saved with only its cells still draws."""
    points = passes[-1][0]
    step = max(1, points // most)
    knees = [k for _, k in passes if k < points]
    x = np.unique(np.concatenate([np.arange(0, points, step), [points - 1], knees])).astype(np.int64)
    return {"cumulative_x": x, "cumulative_y": cumulative[x]}


# ---------------------------------------------------------------------------
# Density
# ---------------------------------------------------------------------------


@dataclass
class DensityValleys:
    """What the density rule saw: its grid, density, valleys and choice."""

    grid: np.ndarray
    density: np.ndarray
    valleys: np.ndarray  # grid indices of the local minima, increasing
    passing: np.ndarray  # barcodes above each valley, in the order tried (highest first)
    chosen: Optional[int]  # grid index of the chosen valley


def density_valleys(counts: np.ndarray, expect_cells: Optional[int] = None) -> DensityValleys:
    """UMI-tools' density rule up to the choice of a valley (`counts` > 0, any order)."""
    from scipy.signal import argrelextrema
    from scipy.stats import gaussian_kde

    if counts.shape[0] == 0:
        raise UmiToolsError("no barcode has a count; UMI-tools needs at least one")
    # UMI-tools keeps the counts sorted in decreasing order: the kernel sums follow it.
    ordered = np.sort(counts)[::-1]
    floor = 0.001 * ordered[0]
    log_counts = np.log10(ordered[ordered > floor])
    try:
        kde = gaussian_kde(log_counts, bw_method=0.1)
    except (ValueError, np.linalg.LinAlgError) as error:
        raise UmiToolsError(f"the kernel density cannot be fitted to these counts: {error}") from None
    top = log_counts.max()
    grid = np.linspace(log_counts.min(), top, GRID)
    density = kde(grid)
    valleys = argrelextrema(density, np.less)[0]
    passing = []
    chosen = None
    for valley in valleys[::-1]:
        n_pass = int(np.count_nonzero(counts > np.power(10, grid[valley])))
        passing.append(n_pass)
        if chosen is not None:
            continue
        if expect_cells:
            if n_pass > expect_cells * 0.1 and n_pass <= expect_cells:
                chosen = int(valley)
        elif valley >= 0.2 * GRID and (top - grid[valley] > 0.5 or grid[valley] < top / 2):
            chosen = int(valley)
    return DensityValleys(grid, density, valleys, np.asarray(passing, dtype=np.int64), chosen)


@register
class Density(Method):
    info = MethodInfo(
        id="umitools-density",
        title="UMI-tools whitelist, density valley",
        tool=TOOL,
        category="density",
        summary=(
            "A valley of a Gaussian kernel density of log10 counts, chosen by fixed rules or by "
            "an expected number of cells."
        ),
        docs="methods/umitools/density.md",
        counts="either",
        variants=(
            Variant(
                "v0.5",
                since="0.5.0",
                summary="whitelist's only rule in 0.5; --knee-method=density from 1.0.0",
            ),
        ),
        latest_known=LATEST,
        references=(_PAPER,),
        source=_SOURCE,
        original_license="MIT",
    )

    @dataclass
    class Params:
        expect_cells: Optional[int] = param(
            None,
            "Choose the first valley with more than a tenth of this and at most this many barcodes above it",
            original="--expect-cells",
            minimum=1,
        )
        n_cells: Optional[int] = param(
            None,
            "Call the barcodes above the (N+1)-th largest count instead of a valley",
            original="--set-cell-number",
            minimum=1,
        )
        allow_threshold_error: bool = param(
            False, "Call no cells, instead of failing, when no valley passes", original="--allow-threshold-error"
        )
        unit: str = param("umis", "Count UMIs (the matrix) or reads (--reads)", original="--method", choices=("umis", "reads"))

    def run(self, droplets, params, variant, random_state):
        if params.expect_cells is not None and params.n_cells is not None:
            raise UmiToolsError("--expect-cells and --n-cells cannot be given together (UMI-tools refuses both)")
        values = droplets.values(params.unit)
        index, counts = _observed(values)
        found = density_valleys(counts, params.expect_cells)
        is_cell = np.zeros(values.shape[0], dtype=bool)
        thresholds = {}
        notes = []
        if params.n_cells is not None:
            cut = _cut_above_nth(counts, params.n_cells)
        elif found.chosen is not None:
            cut = np.power(10, found.grid[found.chosen])
        else:
            if not params.allow_threshold_error:
                raise UmiToolsError(
                    "no valley of the density passes UMI-tools' rules (\"No local minima was accepted\"); "
                    "give --expect-cells or --n-cells, or --allow-threshold-error to call no cells"
                )
            cut = None
            notes.append("no valley passed; no cells called (--allow-threshold-error)")
        if cut is not None:
            is_cell[index[counts > cut]] = True
            thresholds[params.unit] = float(cut)
        step = max(1, GRID // 1000)
        diagnostics = {
            "unit": params.unit,
            "grid": found.grid[::step],
            "density": found.density[::step],
            "valleys": found.grid[found.valleys],
            "valley_density": found.density[found.valleys],
            "valley_passing": found.passing[::-1],
        }
        if found.chosen is not None and params.n_cells is None:
            diagnostics["chosen"] = float(found.grid[found.chosen])
        return result(
            type(self),
            variant,
            params,
            droplets,
            is_cell,
            thresholds=thresholds,
            diagnostics=diagnostics,
            notes=notes,
        )
