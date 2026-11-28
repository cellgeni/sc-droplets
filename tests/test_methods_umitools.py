"""UMI-tools' whitelist rules on small and simulated matrices."""

from __future__ import annotations

import re

import anndata as ad
import numpy as np
import pytest
import scipy.sparse as sp
from typer.testing import CliRunner

import sc_droplets as scd
from sc_droplets.cli import app
from sc_droplets.methods.umitools import (
    UmiToolsError,
    _observed,
    density_valleys,
    distance_knee,
    knee_distance,
)
from sc_droplets.registry import UnsupportedVersion

runner = CliRunner()


def _out(result) -> str:
    return re.sub(r"\s+", " ", re.sub(r"\x1b\[[0-9;]*m", "", result.output))


def _matrix(totals) -> ad.AnnData:
    """One feature: each barcode's total is its count."""
    adata = ad.AnnData(sp.csr_matrix(np.asarray(totals, dtype=np.float32).reshape(-1, 1)))
    adata.obs_names = [f"B{i}" for i in range(len(totals))]
    adata.var_names = ["G0"]
    adata.var["feature_types"] = "Gene Expression"
    return adata


def _cells(adata, result) -> list:
    return list(adata.obs_names[result.is_cell])


# ---------------------------------------------------------------------------
# Distance
# ---------------------------------------------------------------------------


def test_knee_distance_is_the_perpendicular_distance_from_the_chord():
    values = np.cumsum([50, 40, 30, 5, 4, 3, 2, 1, 1, 1])
    distance, knee = knee_distance(values)
    x = np.arange(values.size, dtype=float)
    dx, dy = x[-1] - x[0], float(values[-1] - values[0])
    expected = np.abs(dy * x - dx * (values - values[0])) / np.hypot(dx, dy)
    np.testing.assert_allclose(distance, expected, rtol=1e-12, atol=1e-12)
    assert knee == int(np.argmax(expected)) == 2


def test_knee_distance_takes_the_first_of_equal_distances():
    # A symmetric bump: two points at the same distance from the chord.
    values = np.array([0, 2, 2, 0], dtype=np.int64)
    distance, knee = knee_distance(values)
    assert distance[1] == distance[2]
    assert knee == 1


def test_distance_knee_searches_again_on_three_times_the_knee():
    counts = np.array([1000] * 20 + [100] * 30 + [3] * 3000, dtype=np.int64)
    knee, passes = distance_knee(counts)
    assert passes[0][0] == counts.size
    for (points, _), (_, previous) in zip(passes[1:], passes[:-1]):
        assert points == min(3 * previous, counts.size)
    assert passes[-1][1] == passes[-2][1] == knee


def test_distance_knee_at_the_first_barcode_is_an_error():
    with pytest.raises(UmiToolsError, match="first barcode"):
        distance_knee(np.array([10], dtype=np.int64))
    with pytest.raises(UmiToolsError, match="no barcode"):
        distance_knee(np.array([], dtype=np.int64))


def test_distance_calls_the_barcodes_up_to_the_knee(sim):
    result = scd.call(sim, "umitools-distance", inplace=False)
    truth = sim.obs["truth"].to_numpy()
    assert result.n_cells == result.thresholds["knee_rank"]
    assert set(truth[result.is_cell]) == {"cell"}
    assert result.n_cells > 0.8 * (truth == "cell").sum()


def test_barcodes_without_counts_are_not_on_umitools_curve():
    """A raw matrix lists barcodes with no counts of the chosen features; UMI-tools never sees
    them, and adding them to the cumulative curve would move its chord and the knee."""
    totals = [900, 800, 700, 600, 50, 40, 30, 20, 10, 5, 3, 2, 1, 1, 1]
    alone = scd.call(_matrix(totals), "umitools-distance", inplace=False)
    padded = scd.call(_matrix(totals + [0] * 200), "umitools-distance", inplace=False)
    assert padded.n_cells == alone.n_cells
    assert _cells(_matrix(totals), alone) == list(np.asarray(padded.barcodes)[padded.is_cell])


def test_distance_breaks_ties_at_the_knee_by_input_order():
    # The chord's slope is exactly 13, so along the run of 13s the distances are equal but for
    # rounding, and the knee lands inside the run: which 13s are cells follows input order.
    counts = [59, 44, 34, 21] + [13] * 8 + [1] * 5
    order = np.random.default_rng(3).permutation(len(counts))
    adata = _matrix([counts[i] for i in order])
    result = scd.call(adata, "umitools-distance", inplace=False)
    assert result.thresholds["knee_rank"] == 7
    totals = np.asarray(result.totals)
    thirteens = [f"B{i}" for i in np.flatnonzero(totals == 13)]
    assert _cells(adata, result) == sorted(
        [f"B{i}" for i in np.flatnonzero(totals > 13)] + thirteens[:3], key=lambda b: int(b[1:])
    )


def test_n_cells_calls_the_barcodes_above_the_next_count():
    adata = _matrix([90, 80, 70, 70, 60, 5, 4, 3, 2, 1])
    result = scd.call(adata, "umitools-distance", n_cells=3, inplace=False)
    # the 4th largest count is 70, tied with the 3rd: only 2 barcodes are above it
    assert _cells(adata, result) == ["B0", "B1"]
    assert result.thresholds["umis"] == 70
    result = scd.call(adata, "umitools-distance", n_cells=4, inplace=False)
    assert _cells(adata, result) == ["B0", "B1", "B2", "B3"]


def test_n_cells_must_leave_a_count_below_it():
    with pytest.raises(UmiToolsError, match="not below the number of barcodes"):
        scd.call(_matrix([90, 80, 70, 0]), "umitools-distance", n_cells=3, inplace=False)


def test_distance_predates_umitools_1_0():
    with pytest.raises(UnsupportedVersion, match="predates"):
        scd.call(_matrix([90, 80, 70]), "umitools-distance@0.5.5", inplace=False)


def test_whole_counts_are_integers_as_in_umitools():
    _, counts = _observed(np.array([0.0, 3.0, 2.0]))
    assert counts.dtype == np.int64 and counts.tolist() == [3, 2]
    _, counts = _observed(np.array([0.5, 3.0]))
    assert counts.dtype == np.float64


# ---------------------------------------------------------------------------
# Density
# ---------------------------------------------------------------------------


def test_density_cuts_at_a_valley_between_empty_droplets_and_cells(sim):
    result = scd.call(sim, "umitools-density", inplace=False)
    truth = sim.obs["truth"].to_numpy()
    totals = np.asarray(result.totals)
    cut = result.thresholds["umis"]
    assert np.array_equal(result.is_cell, totals > cut)
    assert set(truth[result.is_cell]) == {"cell"}
    assert 10 ** result.diagnostics["chosen"] == cut


def test_density_tries_the_valleys_from_the_highest_down(sim):
    totals = np.asarray(scd.call(sim, "umitools-density", inplace=False).totals)
    _, counts = _observed(totals)
    found = density_valleys(counts)
    assert len(found.valleys) >= 1
    for valley, n_pass in zip(found.valleys[::-1], found.passing):
        assert n_pass == np.count_nonzero(counts > 10 ** found.grid[valley])
    assert found.passing.tolist() == sorted(found.passing.tolist())


def test_expect_cells_picks_the_first_valley_in_range(sim):
    totals = np.asarray(scd.call(sim, "umitools-density", inplace=False).totals)
    _, counts = _observed(totals)
    found = density_valleys(counts)
    target = int(found.passing[0])  # barcodes above the highest valley
    chosen = density_valleys(counts, expect_cells=max(target, 1))
    assert chosen.chosen == found.valleys[-1]


def test_no_valley_in_range_is_umitools_error_unless_allowed(sim):
    with pytest.raises(UmiToolsError, match="No local minima was accepted"):
        scd.call(sim, "umitools-density", expect_cells=1, inplace=False)
    result = scd.call(sim, "umitools-density", expect_cells=1, allow_threshold_error=True, inplace=False)
    assert result.n_cells == 0
    assert result.thresholds == {}
    assert any("no valley" in note for note in result.notes)


def test_expect_cells_and_n_cells_together_are_refused(sim):
    with pytest.raises(UmiToolsError, match="cannot be given together"):
        scd.call(sim, "umitools-density", expect_cells=300, n_cells=300, inplace=False)


def test_density_with_n_cells_uses_the_next_count(sim):
    result = scd.call(sim, "umitools-density", n_cells=100, inplace=False)
    totals = np.asarray(result.totals)
    nth = np.sort(totals[totals > 0])[::-1][100]
    assert np.array_equal(result.is_cell, totals > nth)
    assert "chosen" not in result.diagnostics


def test_density_of_equal_counts_cannot_be_fitted():
    with pytest.raises(UmiToolsError, match="kernel density"):
        scd.call(_matrix([5, 5, 5, 5]), "umitools-density", inplace=False)


# ---------------------------------------------------------------------------
# Reads, storage, plots, CLI
# ---------------------------------------------------------------------------


def test_reads_are_counted_when_asked(sim):
    adata = sim.copy()
    adata.obs["reads"] = np.asarray(adata.X.sum(axis=1)).ravel() * 3
    by_umis = scd.call(adata, "umitools-distance", inplace=False)
    by_reads = scd.call(adata, "umitools-distance", unit="reads", reads="reads", inplace=False)
    assert np.array_equal(by_umis.is_cell, by_reads.is_cell)
    assert by_reads.thresholds["reads"] == 3 * by_umis.thresholds["umis"]


def test_calls_and_their_diagnostics_survive_h5ad(adata, tmp_path):
    scd.call(adata, "umitools-distance")
    scd.call(adata, "umitools-density")
    adata.write_h5ad(tmp_path / "calls.h5ad")
    back = ad.read_h5ad(tmp_path / "calls.h5ad")
    record = back.uns["sc_droplets"]["umitools_distance"]
    assert np.asarray(record["diagnostics"]["passes"]).shape[1] == 2
    assert "grid" in back.uns["sc_droplets"]["umitools_density"]["diagnostics"]


def test_diagnostics_draw_the_cumulative_curve_and_the_density(adata):
    import matplotlib

    matplotlib.use("Agg")
    scd.call(adata, "umitools-distance")
    scd.call(adata, "umitools-density")
    fig = scd.pl.diagnostics(adata, "umitools_distance")
    assert len(fig.axes) == 2
    assert "Cumulative" in fig.axes[1].get_title()
    assert any("knee" in (line.get_label() or "") for line in fig.axes[1].get_lines())
    fig = scd.pl.diagnostics(adata, "umitools_density")
    assert "Density" in fig.axes[1].get_title()
    chosen = adata.uns["sc_droplets"]["umitools_density"]["diagnostics"]["chosen"]
    assert any(line.get_xdata()[0] == chosen for line in fig.axes[1].get_lines()[1:])


def test_distance_diagnostics_draw_after_keeping_only_the_cells(adata):
    """The cumulative panel used to rebuild the curve from the AnnData's barcodes; once a call
    was saved with --cells-only (or subset to its cells) the last pass had more points than
    were left, and the plot failed on arrays of unequal length."""
    import matplotlib

    matplotlib.use("Agg")
    scd.call(adata, "umitools-distance")
    cells = adata[adata.obs["umitools_distance"].to_numpy(dtype=bool)].copy()
    passes = np.asarray(cells.uns["sc_droplets"]["umitools_distance"]["diagnostics"]["passes"])
    assert passes[-1][0] > cells.n_obs
    fig = scd.pl.diagnostics(cells, "umitools_distance")
    curve = fig.axes[1].get_lines()[0]
    assert curve.get_xdata()[-1] == passes[-1][0] - 1


def test_cli_calls_umitools_distance(sim, tmp_path):
    path = tmp_path / "sim.h5ad"
    sim.write_h5ad(path)
    result = runner.invoke(app, ["call", "umitools-distance", str(path), "--n-cells", "100", "--json"])
    assert result.exit_code == 0, _out(result)
    expected = scd.call(sim, "umitools-distance", n_cells=100, inplace=False).n_cells
    assert f'"n_cells": {expected}' in _out(result)


def test_cli_reports_umitools_density_failure_in_one_line(sim, tmp_path):
    path = tmp_path / "sim.h5ad"
    sim.write_h5ad(path)
    result = runner.invoke(app, ["call", "umitools-density", str(path), "--expect-cells", "1"])
    assert result.exit_code == 1
    assert "Error:" in _out(result) and "No local minima was accepted" in _out(result)
