"""alevin-fry's generate-permit-list rules on small and simulated matrices."""

from __future__ import annotations

import re

import anndata as ad
import numpy as np
import pytest
import scipy.sparse as sp
from typer.testing import CliRunner

import sc_droplets as scd
from sc_droplets.cli import app
from sc_droplets.methods.alevinfry import AlevinFryError, get_knee, threshold

runner = CliRunner()


def _out(result) -> str:
    return re.sub(r"\s+", " ", re.sub(r"\x1b\[[0-9;]*m", "", result.output))


def _matrix(totals) -> ad.AnnData:
    adata = ad.AnnData(sp.csr_matrix(np.asarray(totals, dtype=np.float64).reshape(-1, 1)))
    adata.obs_names = [f"B{i}" for i in range(len(totals))]
    adata.var_names = ["G0"]
    adata.var["feature_types"] = "Gene Expression"
    return adata


def test_the_knee_takes_the_last_of_equally_distant_points():
    # alevin-fry's own test (knee_finding.rs, 0.19.0): the index moves down, 4 -> 2
    knee, _ = get_knee(np.array([100, 100, 100, 100]))
    assert knee == 2


def test_cells_are_the_barcodes_at_or_above_the_knee_count(sim):
    result = scd.call(sim, "alevinfry-permit-list", inplace=False)
    totals = np.asarray(result.totals)
    freq = np.sort(totals[totals > 0])[::-1]
    knee = int(result.thresholds["knee_rank"]) - 1
    assert result.thresholds["umis"] == freq[knee]
    assert np.array_equal(result.is_cell, totals >= freq[knee])


def test_expected_cells_round_halves_away_from_zero():
    # rank round(0.99 * 10) = 10 holds 25: a tenth is 2.5, which Rust rounds to 3
    totals = [1000] * 10 + [25] + [2] * 5 + [3] * 4
    result = scd.call(_matrix(totals), "alevinfry-permit-list", expect_cells=10, inplace=False)
    assert result.thresholds["umis"] == 3
    assert result.n_cells == 15  # the 3s count too


def test_forced_cells_keep_ties_with_the_last():
    result = scd.call(_matrix([50, 40, 30, 30, 30, 5]), "alevinfry-permit-list", force_cells=3, inplace=False)
    assert result.thresholds["umis"] == 30 and result.n_cells == 5


def test_flex_rules_belong_to_0_12_to_0_17():
    freq = np.array([1000, 900, 800, 50, 40, 30, 20, 10, 5, 2, 1])
    single, *_ = threshold(freq, expect_cells=3)
    flex, *_ = threshold(freq, expect_cells=3, flex_rules=True)
    assert single == max(1, round(freq[3] / 10))  # rank round(2.97) = 3
    # floor(1000 / 3) = 333; the first count below it is at rank 3; 10 * 3 = 30 cells, capped at 11
    assert flex == freq[10]
    adata = _matrix(freq)
    assert scd.call(adata, "alevinfry-permit-list@0.17.1", flex=True, expect_cells=3, inplace=False).n_cells == 11
    assert scd.call(adata, "alevinfry-permit-list@0.19.0", flex=True, expect_cells=3, inplace=False).n_cells == 9  # >= 5
    with pytest.raises(AlevinFryError, match="before 0.12.0 takes no multi-barcode"):
        scd.call(adata, "alevinfry-permit-list@0.11.2", flex=True, inplace=False)


def test_the_flex_knee_is_one_rank_higher(sim):
    single = scd.call(sim, "alevinfry-permit-list@0.17.1", inplace=False)
    flex = scd.call(sim, "alevinfry-permit-list@0.17.1", flex=True, inplace=False)
    freq = np.sort(np.asarray(single.totals))[::-1]
    knee = int(single.thresholds["knee_rank"]) - 1
    assert flex.thresholds["umis"] == freq[knee - 1]


def test_the_unfiltered_floor_came_in_0_4_2(sim):
    result = scd.call(sim, "alevinfry-permit-list", unfiltered=True, min_reads=20, inplace=False)
    assert np.array_equal(result.is_cell, np.asarray(result.totals) >= 20)
    with pytest.raises(AlevinFryError, match="before 0.4.2 has no unfiltered"):
        scd.call(sim, "alevinfry-permit-list@0.4.1", unfiltered=True, inplace=False)


def test_one_rule_at_a_time(sim):
    with pytest.raises(AlevinFryError, match="give one of --expect-cells, --force-cells"):
        scd.call(sim, "alevinfry-permit-list", expect_cells=100, force_cells=100, inplace=False)


def test_fractional_counts_are_refused():
    with pytest.raises(AlevinFryError, match="whole reads"):
        scd.call(_matrix([10.5, 3.0, 2.0]), "alevinfry-permit-list", inplace=False)


def test_diagnostics_draw_the_knees_cumulative_curve(adata):
    import matplotlib

    matplotlib.use("Agg")
    scd.call(adata, "alevinfry-permit-list", key_added="af")
    fig = scd.pl.diagnostics(adata[adata.obs["af"].to_numpy()].copy(), "af")
    assert fig.axes[1].get_title() == "Cumulative counts (alevin-fry knee)"


def test_cli_calls_alevinfry_permit_list(sim, tmp_path):
    path = tmp_path / "sim.h5ad"
    sim.write_h5ad(path)
    result = runner.invoke(app, ["call", "alevinfry-permit-list", str(path), "--force-cells", "120", "--json"])
    assert result.exit_code == 0, _out(result)
    expected = scd.call(sim, "alevinfry-permit-list", force_cells=120, inplace=False).n_cells
    assert f'"n_cells": {expected}' in _out(result)
