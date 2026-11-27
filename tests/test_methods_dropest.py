"""dropEst's cell-number bounds on small and simulated matrices."""

from __future__ import annotations

import re

import anndata as ad
import numpy as np
import pytest
import scipy.sparse as sp
from typer.testing import CliRunner

import sc_droplets as scd
from sc_droplets.cli import app
from sc_droplets.methods.dropest import DropEstError, longest_true_start

runner = CliRunner()


def _out(result) -> str:
    return re.sub(r"\s+", " ", re.sub(r"\x1b\[[0-9;]*m", "", result.output))


def _matrix(totals) -> ad.AnnData:
    adata = ad.AnnData(sp.csr_matrix(np.asarray(totals, dtype=np.float64).reshape(-1, 1)))
    adata.obs_names = [f"B{i}" for i in range(len(totals))]
    adata.var_names = ["G0"]
    adata.var["feature_types"] = "Gene Expression"
    return adata


@pytest.mark.parametrize("runs, start", [
    ([False, True, True, False, True, True], 2),  # two longest runs: the first
    ([True, False, True, True, True], 3),
    ([False, False, False], 1),  # R's which.max of all FALSE is 1: the first run
])
def test_the_longest_run_of_true_is_found_as_r_finds_it(runs, start):
    assert longest_true_start(np.array(runs)) == start


def test_dropest_bands_follow_the_bounds(sim):
    result = scd.call(sim, "dropest-bounds", inplace=False)
    t = result.thresholds
    classes = np.asarray(result.classes.astype(str))
    totals = np.asarray(result.totals)
    order = np.argsort(-totals, kind="stable")
    assert (classes[order[: int(t["min"])]] == "high").all()
    assert (classes[order[int(t["max"]) - 1:]] == "low").all()
    assert (classes == "unknown").sum() == int(t["max"]) - int(t["min"]) - 1
    assert result.n_cells == int(t["expected"])


def test_dropest_needs_a_lag_of_a_barcode_at_least():
    with pytest.raises(DropEstError, match="lag"):
        scd.call(_matrix([50, 40, 30, 20, 10, 5]), "dropest-bounds", lag=0.05, inplace=False)


def test_cli_calls_dropest_bounds(sim, tmp_path):
    path = tmp_path / "sim.h5ad"
    sim.write_h5ad(path)
    result = runner.invoke(app, ["call", "dropest-bounds", str(path), "--json"])
    assert result.exit_code == 0, _out(result)
    assert '"n_cells":' in _out(result)
