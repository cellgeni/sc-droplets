"""zUMIs' cell barcodes on small and simulated matrices."""

from __future__ import annotations

import re

import anndata as ad
import numpy as np
import pytest
import scipy.sparse as sp
from typer.testing import CliRunner

import sc_droplets as scd
from sc_droplets.cli import app
from sc_droplets.methods.zumis import check_curve, uik
from sc_droplets.registry import UnsupportedVersion

runner = CliRunner()


def _out(result) -> str:
    return re.sub(r"\s+", " ", re.sub(r"\x1b\[[0-9;]*m", "", result.output))


def _matrix(totals) -> ad.AnnData:
    adata = ad.AnnData(sp.csr_matrix(np.asarray(totals, dtype=np.float64).reshape(-1, 1)))
    adata.obs_names = [f"B{i}" for i in range(len(totals))]
    adata.var_names = ["G0"]
    adata.var["feature_types"] = "Gene Expression"
    return adata


def test_cumulative_reads_are_a_concave_curve_and_the_knee_is_farthest_above_the_chord():
    x = np.arange(1, 51, dtype=np.float64)
    y = np.cumsum(np.concatenate([np.full(10, 100.0), np.full(40, 5.0)]))
    assert check_curve(x, y) == 1
    chord = y[0] + (y[-1] - y[0]) * (x - x[0]) / (x[-1] - x[0])
    assert uik(x, y) == x[int(np.argmax(y - chord))]


def test_zumis_cells_are_the_barcodes_with_the_knees_count(sim):
    result = scd.call(sim, "zumis-cell-barcodes", inplace=False)
    totals = np.asarray(result.totals)
    assert np.array_equal(result.is_cell, totals >= result.thresholds["umis"])


def test_zumis_falls_back_to_the_top_100_below_10_cells():
    totals = [10000, 9000] + [150] * 200
    result = scd.call(_matrix(totals), "zumis-cell-barcodes", inplace=False)
    assert result.n_cells == 100 and result.notes


def test_zumis_before_2_5_6b_is_refused():
    with pytest.raises(UnsupportedVersion, match="random subset"):
        scd.call(_matrix([500] * 20), "zumis-cell-barcodes@2.5.4", inplace=False)


def test_cli_calls_zumis_cell_barcodes(sim, tmp_path):
    path = tmp_path / "sim.h5ad"
    sim.write_h5ad(path)
    result = runner.invoke(app, ["call", "zumis-cell-barcodes", str(path), "--json", "--n-reads-per-cell", "50"])
    assert result.exit_code == 0, _out(result)
    assert '"n_cells":' in _out(result)
