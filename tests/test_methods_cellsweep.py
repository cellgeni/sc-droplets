"""cellsweep's empty-droplet step on small and simulated matrices."""

from __future__ import annotations

import re

import anndata as ad
import numpy as np
import pytest
import scipy.sparse as sp
from typer.testing import CliRunner

import sc_droplets as scd
from sc_droplets.cli import app
from sc_droplets.primitives import kneed
from sc_droplets.methods.cellsweep import (
    CellsweepError,
    kneed_cutoff,
    mx_filter_cutoff,
    threshold_cutoff,
)

runner = CliRunner()


def _out(result) -> str:
    return re.sub(r"\s+", " ", re.sub(r"\x1b\[[0-9;]*m", "", result.output))


def _matrix(totals) -> ad.AnnData:
    adata = ad.AnnData(sp.csr_matrix(np.asarray(totals, dtype=np.float32).reshape(-1, 1)))
    adata.obs_names = [f"B{i}" for i in range(len(totals))]
    adata.var_names = ["G0"]
    adata.var["feature_types"] = "Gene Expression"
    return adata


def _totals(adata) -> np.ndarray:
    return np.asarray(adata.X.sum(axis=1), dtype=np.float32).ravel()


def test_cells_are_the_barcodes_at_or_above_the_cutoff(sim):
    result = scd.call(sim, "cellsweep-empty", inplace=False)
    totals = np.asarray(result.totals)
    assert np.array_equal(result.is_cell, totals >= result.thresholds["umis"])


def test_threshold_rule_is_the_most_negative_second_derivative(sim):
    totals = _totals(sim)
    cutoff, counts, smooth, d2 = threshold_cutoff(totals)
    assert counts.dtype == np.float32 and smooth.dtype == np.float32
    assert (counts > 10).all() and np.all(np.diff(counts) <= 0)
    assert cutoff == counts[np.argmin(d2)]
    result = scd.call(sim, "cellsweep-empty@1.0.1", inplace=False)
    assert result.variant == "v0.1" and result.diagnostics["rule"] == "threshold"
    assert result.thresholds["umis"] == cutoff


def test_kneed_rule_places_each_distinct_total_at_its_last_rank():
    totals = np.array([500, 500, 400, 300, 300, 300, 20, 15, 15, 12, 11, 11, 5, 0], dtype=np.float32)
    cutoff, found, index = kneed_cutoff(totals, "0.9.0")
    assert found.x.tolist() == [2, 3, 6, 7, 9, 10, 12]  # 500, 400, 300, 20, 15, 12, 11
    assert cutoff == [500, 400, 300, 20, 15, 12, 11][index]


@pytest.mark.parametrize("version, behaviour", [("0.8.5", "0.7"), ("0.8.6", "0.8.6"), ("0.9.0", "0.9")])
def test_kneed_rule_runs_the_installed_kneeds_online_search(sim, version, behaviour):
    totals = _totals(sim)
    t = totals[totals > 10]
    umis, n = np.unique(t, return_counts=True)
    found = kneed.prepare(np.cumsum(n[::-1]), umis[::-1], curve="convex", direction="decreasing", behaviour=behaviour)
    expected = umis[::-1][found.find(1.0, online=True)]
    result = scd.call(sim, "cellsweep-empty", kneed_version=version, inplace=False)
    assert result.thresholds["umis"] == expected


def test_a_single_distinct_total_has_no_knee():
    with pytest.raises(CellsweepError, match="Could not detect a knee"):
        scd.call(_matrix([50, 50, 50, 3]), "cellsweep-empty", inplace=False)


def test_mx_filter_cuts_where_the_mixture_is_least_sure(sim):
    totals = _totals(sim)
    totals = totals[totals >= 100]
    cutoff, x, ent, gm = mx_filter_cutoff(totals)
    assert x.dtype == np.float32
    assert cutoff == x[np.argmax(ent)]
    assert gm.means_.shape == (2, 1)
    assert min(gm.means_.ravel()) < np.log1p(cutoff) < max(gm.means_.ravel())


def test_manual_cutoffs_win_and_expected_cells_ranks_the_totals():
    adata = _matrix([90, 80, 70, 70, 60, 5, 0])
    result = scd.call(adata, "cellsweep-empty", umi_cutoff=65, inplace=False)
    assert result.n_cells == 4
    result = scd.call(adata, "cellsweep-empty", expected_cells=3, inplace=False)
    assert result.thresholds["umis"] == 70 and result.n_cells == 4  # ties with the 3rd are cells
    with pytest.raises(CellsweepError, match="more than"):
        scd.call(adata, "cellsweep-empty", expected_cells=8, inplace=False)


def test_before_1_0_2_mx_filter_ignores_expected_cells(sim):
    """cellsweep 0.1.0 - 1.0.1 checks only umi_cutoff in its mx_filter branch; from 1.0.2 any
    expected_cells selects the manual rule."""
    old = scd.call(sim, "cellsweep-empty@1.0.1", empty_method="mx_filter", expected_cells=100, inplace=False)
    assert old.diagnostics["rule"] == "mx_filter"
    new = scd.call(sim, "cellsweep-empty@1.0.2", empty_method="mx_filter", expected_cells=100, inplace=False)
    assert new.diagnostics["rule"] == "manual"


def test_rules_belong_to_their_releases(sim):
    with pytest.raises(CellsweepError, match="no 'kneed' rule"):
        scd.call(sim, "cellsweep-empty@1.0.1", empty_method="kneed", inplace=False)
    # 1.0.2 renamed threshold to manual, which needs a cutoff or expected cells
    with pytest.raises(CellsweepError, match="requires either umi_cutoff or expected_cells"):
        scd.call(sim, "cellsweep-empty@1.0.2", empty_method="threshold", inplace=False)


def test_float64_matrices_compute_in_float64(sim):
    result = scd.call(sim, "cellsweep-empty", matrix_dtype="float64", inplace=False)
    assert np.array_equal(result.is_cell, np.asarray(result.totals) >= result.thresholds["umis"])


@pytest.mark.parametrize("params", [{"tool_version": "1.0.1"}, {}, {"empty_method": "mx_filter"}])
def test_diagnostics_survive_h5ad_and_draw_each_rule(adata, tmp_path, params):
    import matplotlib

    matplotlib.use("Agg")
    scd.call(adata, "cellsweep-empty", **params)
    adata.write_h5ad(tmp_path / "calls.h5ad")
    back = ad.read_h5ad(tmp_path / "calls.h5ad")
    key = next(iter(back.uns["sc_droplets"]))
    fig = scd.pl.diagnostics(back, key)
    rule = back.uns["sc_droplets"][key]["diagnostics"]["rule"]
    assert fig.axes[1].get_title() == f"cellsweep: {rule}"


def test_cli_calls_cellsweep_empty(sim, tmp_path):
    path = tmp_path / "sim.h5ad"
    sim.write_h5ad(path)
    result = runner.invoke(app, ["call", "cellsweep-empty", str(path), "--kneed-version", "0.8.5", "--json"])
    assert result.exit_code == 0, _out(result)
    expected = scd.call(sim, "cellsweep-empty", kneed_version="0.8.5", inplace=False).n_cells
    assert f'"n_cells": {expected}' in _out(result)
