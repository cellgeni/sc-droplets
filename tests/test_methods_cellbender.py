"""CellBender's priors on small and simulated matrices."""

from __future__ import annotations

import re

import anndata as ad
import numpy as np
import pytest
import scipy.sparse as sp
from typer.testing import CliRunner

import sc_droplets as scd
from sc_droplets.cli import app
from sc_droplets.methods.cellbender import CellBenderError, threshold_otsu
from sc_droplets.primitives.legacysort import legacy_argsort

runner = CliRunner()
RELEASES = ("0.1.0", "0.2.2", "0.4.0")


def _out(result) -> str:
    return re.sub(r"\s+", " ", re.sub(r"\x1b\[[0-9;]*m", "", result.output))


def _matrix(totals) -> ad.AnnData:
    adata = ad.AnnData(sp.csr_matrix(np.asarray(totals, dtype=np.int32).reshape(-1, 1)))
    adata.obs_names = [f"B{i}" for i in range(len(totals))]
    adata.var_names = ["G0"]
    return adata


def _classes(result) -> np.ndarray:
    return np.asarray(result.classes.astype(str))


@pytest.mark.parametrize("release", RELEASES)
def test_classes_follow_the_ranking_of_the_totals(sim, release):
    result = scd.call(sim, f"cellbender-prior@{release}", inplace=False)
    # 0.3+ ranks on the totals over the features it keeps; on `sim` it keeps them all.
    totals, classes = np.asarray(result.totals), _classes(result)
    cell, transition, empty = (totals[classes == c] for c in ("cell", "transition", "empty"))
    assert cell.size == result.thresholds["expected_cells"] == result.n_cells
    assert cell.min() >= transition.max()
    if empty.size:
        assert transition.min() >= empty.max()
        assert empty.min() > result.thresholds["low_count_cutoff"]


def test_ties_at_a_class_boundary_fall_as_numpys_classic_quicksort_leaves_them():
    totals = [3000] * 60 + [2000] * 40 + [21] * 100 + [20] * 400 + [0] * 100
    adata = _matrix(totals)
    result = scd.call(adata, "cellbender-prior@0.1.0", total_droplets_included=300, inplace=False)
    order = legacy_argsort(np.asarray(totals, dtype=np.int64))[::-1]
    classes = _classes(result)
    assert set(np.flatnonzero(classes == "transition")) == set(order[100:300])
    assert (classes == "empty").sum() == 300  # the other 20s


def test_before_0_3_the_empty_prior_is_the_mode_of_rounded_log_counts():
    # log1p(20) rounds to 3.0 (400 barcodes), log1p(21) to 3.1 (100): the empty prior is
    # int(expm1(3.0)) = 19; the cell prior is the log-median of the barcodes above 5 * 19.
    totals = [3000] * 60 + [2000] * 40 + [21] * 100 + [20] * 400 + [0] * 100
    result = scd.call(_matrix(totals), "cellbender-prior@0.1.0", inplace=False)
    big = np.array([3000] * 60 + [2000] * 40, dtype=np.int64)
    assert result.thresholds["empty_counts"] == 19
    assert result.thresholds["cell_counts"] == int(np.expm1(np.median(np.log1p(big))))
    assert result.n_cells == 100
    assert result.thresholds["low_count_cutoff"] == 15  # max(15, int(0.8 * 19))


def test_the_low_cut_is_0_8_of_the_empty_prior_in_0_1_and_0_5_from_0_2(sim):
    old = scd.call(sim, "cellbender-prior@0.1.0", inplace=False).thresholds
    new = scd.call(sim, "cellbender-prior@0.2.0", inplace=False).thresholds
    assert old["low_count_cutoff"] == max(15, int(0.8 * old["empty_counts"]))
    assert new["low_count_cutoff"] == max(15, int(0.5 * new["empty_counts"]))


def test_given_expected_cells_set_the_cell_prior_by_release(sim):
    totals = np.sort(np.asarray(sim.X.sum(axis=1)).ravel().astype(np.int64))[::-1]
    v1 = scd.call(sim, "cellbender-prior@0.1.0", expected_cells=100, inplace=False).thresholds
    v2 = scd.call(sim, "cellbender-prior@0.2.0", expected_cells=100, inplace=False).thresholds
    v3 = scd.call(sim, "cellbender-prior@0.3.0", expected_cells=100, total_droplets_included=1000, inplace=False).thresholds
    auto = scd.call(sim, "cellbender-prior@0.1.0", inplace=False).thresholds
    assert v1["cell_counts"] == auto["cell_counts"]  # 0.1.0 keeps its log-median
    assert v2["cell_counts"] == int(np.median(totals[:100]))
    assert v3["cell_counts"] == np.exp(np.mean(np.log(totals[:100])))
    assert v1["expected_cells"] == v2["expected_cells"] == v3["expected_cells"] == 100
    assert v3["total_droplets"] == 1000


def test_from_0_3_the_expected_cells_are_the_barcodes_at_or_above_the_cell_prior(sim):
    result = scd.call(sim, "cellbender-prior", inplace=False)
    totals = np.asarray(result.totals)
    assert result.thresholds["expected_cells"] == (totals >= result.thresholds["cell_counts"]).sum()
    assert result.thresholds["surely_empty_counts"] == np.sort(totals)[::-1][int(result.thresholds["total_droplets"])]


def test_forced_priors_belong_to_0_3_and_set_the_upper_limit(sim):
    result = scd.call(sim, "cellbender-prior", force_cell_umi_prior=3000, force_empty_umi_prior=40, inplace=False)
    assert result.thresholds["cell_counts"] == 3000 and result.thresholds["empty_counts"] == 40
    assert result.thresholds["empty_count_upper_limit"] == min(np.sqrt(3000 * 40), 80)
    with pytest.raises(CellBenderError, match="before 0.3.0 has no --force-cell-umi-prior"):
        scd.call(sim, "cellbender-prior@0.2.2", force_cell_umi_prior=3000, inplace=False)


def test_features_with_almost_no_ambient_counts_are_left_out_from_0_3(sim):
    base = scd.call(sim, "cellbender-prior", inplace=False)
    trimmed = scd.call(sim, "cellbender-prior", projected_ambient_count_threshold=150, inplace=False)
    totals = np.asarray(sim.X.sum(axis=1)).ravel()
    upper = trimmed.thresholds["empty_count_upper_limit"]
    rows = (totals < upper) & (totals > 5)
    mean = np.asarray(sim.X[rows].astype(np.float64).mean(axis=0)).ravel()
    keep = max(trimmed.thresholds["expected_cells"], 5000) * mean > 150
    assert base.thresholds["n_features"] == sim.n_vars
    assert trimmed.thresholds["n_features"] == keep.sum() < sim.n_vars
    trimmed_totals = np.asarray(sim.X[:, keep].sum(axis=1)).ravel()
    classes = _classes(trimmed)
    assert trimmed_totals[classes == "cell"].min() >= trimmed_totals[classes == "transition"].max()


def test_total_droplets_must_exceed_the_cells(sim):
    # From 0.3 a total at or below the cells first empties the window of the empty prior.
    with pytest.raises(CellBenderError, match=r"expected cells \(500\) must be fewer than the total droplets \(400\)"):
        scd.call(sim, "cellbender-prior", expected_cells=500, total_droplets_included=400, inplace=False)
    with pytest.raises(CellBenderError, match="Increase --total-droplets-included above 300"):
        scd.call(sim, "cellbender-prior@0.2.2", total_droplets_included=200, inplace=False)


@pytest.mark.parametrize("release, totals, message", [
    ("0.1.0", [10] * 20, "no barcode has more than 15 counts"),
    ("0.1.0", [20] * 50, "more than five times the empty prior"),
    ("0.4.0", [1000] * 10 + [3] * 50, "fewer than 100 barcodes above the low-count threshold"),
])
def test_matrices_cellbender_cannot_take_stop_with_its_reason(release, totals, message):
    with pytest.raises(CellBenderError, match=message):
        scd.call(_matrix(totals), f"cellbender-prior@{release}", inplace=False)


def test_otsu_maximises_the_between_class_variance():
    rng = np.random.default_rng(1)
    x = np.concatenate([rng.normal(2, 0.5, 500), rng.normal(7, 0.7, 200)])
    counts, edges = np.histogram(x, 256)
    centers = (edges[:-1] + edges[1:]) / 2
    variance = []
    for i in range(255):
        w1, w2 = counts[: i + 1].sum(), counts[i + 1 :].sum()
        m1 = (counts[: i + 1] * centers[: i + 1]).sum() / w1
        m2 = (counts[i + 1 :] * centers[i + 1 :]).sum() / w2
        variance.append(w1 * w2 * (m1 - m2) ** 2)
    # Empty bins between the modes give a plateau of equal variance; any bin on it is a maximum.
    chosen = int(np.flatnonzero(centers == threshold_otsu(x))[0])
    assert variance[chosen] == pytest.approx(max(variance), rel=1e-9)
    assert 3 < centers[chosen] < 6


def test_cli_calls_cellbender_prior(sim, tmp_path):
    path = tmp_path / "sim.h5ad"
    sim.write_h5ad(path)
    result = runner.invoke(app, ["call", "cellbender-prior", str(path), "--tool-version", "0.2.2",
                                 "--expected-cells", "120", "--json"])
    assert result.exit_code == 0, _out(result)
    assert '"n_cells": 120' in _out(result)
