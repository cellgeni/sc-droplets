"""Kneedle (kneed's KneeLocator) and its use on the barcode-rank curve."""

from __future__ import annotations

import re

import anndata as ad
import numpy as np
import pytest
import scipy.sparse as sp
from typer.testing import CliRunner

import sc_droplets as scd
from sc_droplets.cli import app
from sc_droplets.methods.kneedle import KneedleError, search
from sc_droplets.primitives import kneed

runner = CliRunner()
BEHAVIOURS = kneed.BEHAVIOURS


def _out(result) -> str:
    return re.sub(r"\s+", " ", re.sub(r"\x1b\[[0-9;]*m", "", result.output))


def _matrix(totals) -> ad.AnnData:
    adata = ad.AnnData(sp.csr_matrix(np.asarray(totals, dtype=np.float32).reshape(-1, 1)))
    adata.obs_names = [f"B{i}" for i in range(len(totals))]
    adata.var_names = ["G0"]
    adata.var["feature_types"] = "Gene Expression"
    return adata


# kneed's own sample curves (kneed/data_generator.py) and the knees its tests expect.
X10 = np.arange(0, 10)
SHAPES = [
    ("concave", "increasing", [0, 60, 80, 85, 90, 95, 96, 97, 98, 99], 2),
    ("concave", "decreasing", [99, 98, 97, 96, 95, 90, 85, 80, 60, 0], 7),
    ("convex", "increasing", [1, 2, 3, 4, 5, 10, 15, 20, 40, 100], 7),
    ("convex", "decreasing", [100, 40, 20, 15, 10, 5, 4, 3, 2, 1], 2),
]


@pytest.mark.parametrize("behaviour", BEHAVIOURS)
@pytest.mark.parametrize("curve, direction, y, expected", SHAPES)
def test_kneed_sample_curves_have_kneeds_knees(behaviour, curve, direction, y, expected):
    knee = kneed.knee_locator(X10, np.array(y), curve=curve, direction=direction, behaviour=behaviour)
    assert knee == expected


@pytest.mark.parametrize("behaviour", ["0.7", "0.8.6", "0.9"])
def test_flat_maxima_follow_the_sensitivity(behaviour):
    # kneed's test_flat_maxima: S = 0 stops at the first maximum, S = 1 at the global one.
    y = [1, 0.787701317715959, 0.7437774524158126, 0.6559297218155198, 0.5065885797950219,
         0.36749633967789164, 0.2547584187408492, 0.16251830161054173, 0.10395314787701318,
         0.06734992679355783, 0.043923865300146414, 0.027818448023426062, 0.01903367496339678,
         0.013177159590043924, 0.010248901903367497, 0.007320644216691069, 0.005856515373352855,
         0.004392386530014641]
    x = np.arange(18, dtype=float)
    assert kneed.knee_locator(x, y, S=0.0, curve="convex", direction="decreasing", behaviour=behaviour) == 1.0
    assert kneed.knee_locator(x, y, S=1.0, curve="convex", direction="decreasing", behaviour=behaviour) == 8.0


@pytest.mark.parametrize("behaviour", ["0.7", "0.8.6", "0.9"])
def test_online_search_keeps_the_last_knee(behaviour):
    # kneed's test_gamma_online_offline (numpy's legacy generator is frozen).
    rng = np.random.RandomState(23)
    y = sorted(rng.gamma(0.5, 1.0, 1000), reverse=True)
    x = np.arange(1, 1001)
    assert kneed.knee_locator(x, y, curve="convex", direction="decreasing", online=False, behaviour=behaviour) == 22
    assert kneed.knee_locator(x, y, curve="convex", direction="decreasing", online=True, behaviour=behaviour) == 482


def test_from_0_9_uneven_x_gives_the_same_knee_in_every_orientation():
    # kneed's test_uneven_x_mirrored (0.9.0): mirroring x by value, not by index.
    x = np.linspace(0, 1, 30) ** 2
    y = -1 / (x + 0.1) + 5
    expected = kneed.knee_locator(x, y, curve="concave", direction="increasing", behaviour="0.9")
    for curve, direction, flip_x, flip_y in [("concave", "decreasing", True, False), ("convex", "increasing", True, True)]:
        xs = 1 - x if flip_x else x
        ys = -y if flip_y else y
        order = np.argsort(xs)
        knee = kneed.knee_locator(xs[order], ys[order], curve=curve, direction=direction, behaviour="0.9")
        assert np.isclose(1 - knee if flip_x else knee, expected)


def test_from_0_8_6_no_knee_is_detected_between_a_minimum_and_the_next_maximum():
    # After the minimum at 2 the threshold is 0 and the next point is below it: before 0.8.6
    # that makes the maximum at 1 a knee; from 0.8.6 detection waits for the next maximum (4),
    # whose threshold nothing after it falls below.
    d = np.array([0.0, 0.5, -0.2, -0.1, 0.3, 0.0])
    fields = dict(x=np.arange(6), x_normalized=np.linspace(0, 1, 6), y_difference=d,
                  maxima=np.array([1, 4]), minima=np.array([2]), curve="convex", direction="decreasing")
    assert kneed.Kneedle(**fields, behaviour="0.7").find(S=10.0) == 1
    assert kneed.Kneedle(**fields, behaviour="0.8.6").find(S=10.0) is None


def test_before_0_7_polynomials_are_of_degree_7_only():
    with pytest.raises(ValueError, match="degree 7 only"):
        kneed.prepare(X10, np.arange(10.0), interp_method="polynomial", polynomial_degree=5, behaviour="0.5.1")


def test_bad_shape_or_smoothing_is_refused():
    with pytest.raises(ValueError, match="curve and direction"):
        kneed.prepare(X10, X10, curve="wavy")
    with pytest.raises(ValueError, match="invalid interp_method"):
        kneed.prepare(X10, X10, interp_method="spline")
    with pytest.raises(ValueError, match="same length"):
        kneed.prepare(X10, X10[:5])


# ---------------------------------------------------------------------------
# The method
# ---------------------------------------------------------------------------


def test_yascp_search_lowers_the_sensitivity_until_a_knee(sim):
    result = scd.call(sim, "kneedle", lower=10, inplace=False)
    S = result.thresholds["sensitivity"]
    assert S <= 5000 and (5000 - S) % 100 == 0
    totals = np.asarray(result.totals)
    found = kneed.prepare(np.arange(1, np.count_nonzero(totals >= 10) + 1),
                          np.sort(totals)[::-1][: np.count_nonzero(totals >= 10)],
                          curve="convex", direction="decreasing")
    if S < 5000:
        assert found.find(S + 100) is None
    assert found.find(S) == result.thresholds["knee_rank"] - 1


def test_cells_are_the_barcodes_ranked_up_to_the_knee():
    totals = [5, 300, 280, 0, 250, 20, 260, 15, 10, 8, 270, 6]
    adata = _matrix(totals)
    result = scd.call(adata, "kneedle", lower=1, sensitivity=1.0, inplace=False)
    rank = int(result.thresholds["knee_rank"])
    order = np.argsort(-np.asarray(totals), kind="stable")
    assert sorted(np.flatnonzero(result.is_cell)) == sorted(order[:rank])
    assert result.thresholds["umis"] == totals[order[rank - 1]]


def test_lower_keeps_the_barcodes_at_or_above_it():
    totals = [100, 90, 50, 40, 40, 3, 2]
    adata = _matrix(totals)
    result = scd.call(adata, "kneedle", lower=40, sensitivity=1.0, inplace=False)
    assert result.diagnostics["points"] == 5


def test_a_curve_needs_two_points_and_log_counts_need_a_positive_lower():
    with pytest.raises(KneedleError, match="needs two points"):
        scd.call(_matrix([5000, 10, 5]), "kneedle", inplace=False)
    with pytest.raises(KneedleError, match="--log-y"):
        scd.call(_matrix([50, 40, 0]), "kneedle", lower=0, log_y=True, inplace=False)


def test_no_knee_at_any_sensitivity_is_an_error():
    # Two points: the difference curve is exactly 0, both points are a maximum and a minimum,
    # the threshold falls to 0 at each, and nothing drops below it whatever S is.
    with pytest.raises(KneedleError, match="at any sensitivity"):
        scd.call(_matrix([100, 50]), "kneedle", lower=1, inplace=False)


def test_search_stops_below_zero():
    class Never:
        tried = []

        def find(self, S, online):
            self.tried.append(S)
            return None

    never = Never()
    assert search(never, None, 300.0, 100.0, False) == (None, -100.0)
    assert never.tried == [300.0, 200.0, 100.0, 0.0, -100.0]


def test_releases_differ_where_kneed_changed(sim):
    # Convex curve on log axes: kneed before 0.7.0 mirrored x as well, which moves the knee.
    params = {"lower": 10, "log_x": True, "log_y": True, "sensitivity": 1.0}
    old = scd.call(sim, "kneedle@0.6.0", inplace=False, **params)
    new = scd.call(sim, "kneedle@0.8.5", inplace=False, **params)
    assert old.variant == "v0.5.1" and new.variant == "v0.7"
    assert old.n_cells != new.n_cells


def test_kneedle_predates_0_5_1():
    with pytest.raises(scd.registry.UnsupportedVersion, match="predates"):
        scd.call(_matrix([90, 80, 70]), "kneedle@0.4.2", inplace=False)


def test_diagnostics_survive_h5ad_and_draw_the_difference_curve(adata, tmp_path):
    import matplotlib

    matplotlib.use("Agg")
    scd.call(adata, "kneedle", lower=10)
    adata.write_h5ad(tmp_path / "calls.h5ad")
    back = ad.read_h5ad(tmp_path / "calls.h5ad")
    cells = back[back.obs["kneedle"].to_numpy(dtype=bool)].copy()
    fig = scd.pl.diagnostics(cells, "kneedle")
    assert len(fig.axes) == 2
    assert "difference" in fig.axes[1].get_title()
    knee = [line for line in fig.axes[1].get_lines() if (line.get_label() or "").startswith("knee")]
    assert knee and knee[0].get_xdata()[0] == back.uns["sc_droplets"]["kneedle"]["diagnostics"]["knee_x"]


def test_cli_calls_kneedle_with_a_release(sim, tmp_path):
    path = tmp_path / "sim.h5ad"
    sim.write_h5ad(path)
    result = runner.invoke(app, ["call", "kneedle", str(path), "--tool-version", "0.8.5", "--lower", "10", "--json"])
    assert result.exit_code == 0, _out(result)
    expected = scd.call(sim, "kneedle@0.8.5", lower=10, inplace=False).n_cells
    assert f'"n_cells": {expected}' in _out(result)
