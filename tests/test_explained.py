"""The identities the Explained pages derive, checked against the methods' own code.

The knee-geometry pages reduce several rules to simpler statements: a chord rule on the cumulative
curve keeps the barcodes above a mean, a refit loop ends at a fixed point of a map, Kneedle's
searched sensitivity picks the maximum followed by the longest fall. The pages teach from those
statements, so they must stay true of the code. Each identity is recomputed here independently of
the method and compared with it. docs/_scripts/explained.py checks the same on the two public
matrices the pages draw.
"""

from __future__ import annotations

import numpy as np
import pytest

import sc_droplets as scd
from sc_droplets.counts import Droplets
from sc_droplets.methods import alevinfry, cellsweep, dropseq, umitools
from sc_droplets.primitives import kneed


@pytest.fixture(scope="module")
def totals(sim):
    return Droplets.from_anndata(sim).values("umis")


@pytest.fixture(scope="module")
def counts(totals):
    """The totals above 0 in decreasing order, ties in matrix order, as the knee rules sort them."""
    _, observed = umitools._observed(totals)
    return observed[np.argsort(-observed, kind="stable")]


def _farthest(x: np.ndarray, y: np.ndarray) -> int:
    """The first point farthest from the chord between the ends: the largest triangle."""
    area = np.abs((x[-1] - x[0]) * (y - y[0]) - (y[-1] - y[0]) * (x - x[0]))
    return int(np.argmax(area))


def test_the_farthest_point_from_a_chord_ignores_the_scale_of_either_axis():
    rng = np.random.default_rng(0)
    for _ in range(50):
        steps = np.sort(rng.lognormal(0, 1.5, size=400))[::-1]
        x, y = np.arange(400.0), np.cumsum(steps)
        found = _farthest(x, y)
        for sx, sy in ((1e-3, 1.0), (1.0, 1e4), (7.0, 0.01), (1 / 399, 1 / y[-1])):
            assert _farthest(sx * x, sy * y) == found


@pytest.mark.parametrize("window", [500, 2000, 7000])
def test_dropbeads_cells_are_the_barcodes_above_the_mean_of_its_window(counts, window):
    top = counts[:window].astype(np.float64)
    assert dropseq.estimate_cell_number(counts.astype(np.float64), window) == np.sum(top > top.mean())


def test_each_umi_tools_pass_keeps_the_barcodes_above_the_mean_of_its_window(counts):
    cumulative = np.cumsum(counts)
    _, passes = umitools.distance_knee(counts)
    windows = [m for m, _ in passes] + [10, 100, 1000, 5000]
    for m in windows:
        _, knee = umitools.knee_distance(cumulative[:m])
        mean = (cumulative[m - 1] - cumulative[0]) / (m - 1)  # the chord starts at the first barcode
        assert knee == np.sum(counts[1:m] > mean), m


def test_each_alevin_fry_pass_keeps_the_barcodes_at_or_above_its_mean(counts):
    cumulative = np.cumsum(counts.astype(np.uint64))
    for m in (10, 100, 1000, 5000, counts.size):
        knee = alevinfry.max_distance_index(cumulative[:m])
        mean = float(cumulative[m - 1] - cumulative[0]) / m  # its line ends one step past the last point
        assert knee == np.sum(counts[1:m] >= mean), m


@pytest.mark.parametrize("floor", [10, 100, 1000])
def test_zumis_cells_are_the_barcodes_above_the_mean_of_its_window(sim, counts, floor):
    window = counts[counts >= floor].astype(np.float64)
    result = scd.call(sim, f"zumis-cell-barcodes:n_reads_per_cell={floor}", inplace=False)
    assert result.n_cells == np.sum(counts > window[1:].mean())


def test_the_refit_loops_end_at_fixed_points_of_their_maps(counts):
    cumulative = np.cumsum(counts)
    knee, _ = umitools.distance_knee(counts)
    assert umitools.knee_distance(cumulative[: min(3 * knee, counts.size)])[1] == knee
    knee, _ = alevinfry.get_knee(counts)
    wide = np.cumsum(counts.astype(np.uint64))
    assert alevinfry.max_distance_index(wide[: min(counts.size - 1, 5 * knee)]) == knee


def test_kneedles_search_takes_the_first_maximum_whose_fall_beats_the_sensitivity(sim, totals):
    """yascp lowers S by 100 until some maximum of the difference curve falls by more than
    S / (n - 1) before the next minimum; the first such maximum is the knee. (The page's "longest
    fall" is the case where one fall dwarfs the rest; with few points a step of S is coarse, and on
    this matrix several maxima qualify at once.)"""
    result = scd.call(sim, "kneedle", inplace=False)
    ranked = np.sort(totals)[::-1]
    y = ranked[ranked >= 1000]
    found = kneed.prepare(np.arange(1, y.size + 1), y, curve="convex", direction="decreasing", behaviour="0.9")
    d, minima = found.y_difference, found.minima
    step = 1 / (y.size - 1)
    falls = []
    for top in found.maxima:
        after = minima[minima > top]
        falls.append((int(top), d[top] - d[after[0] if after.size else d.size - 1]))
    S = result.thresholds["sensitivity"]
    assert max(fall for _, fall in falls) <= (S + 100) * step  # no knee one step of S earlier
    first = next(top for top, fall in falls if fall > S * step)
    assert first + 1 == result.thresholds["knee_rank"]


def test_cellsweeps_online_knee_is_the_last_maximum_of_its_difference_curve(totals):
    for version in ("0.8.5", "0.8.6", "0.9.0"):
        _, found, index = cellsweep.kneed_cutoff(totals, version)
        assert index == found.maxima[-1], version
