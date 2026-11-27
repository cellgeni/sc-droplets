"""MAD outliers."""

from __future__ import annotations

import math
import re

import anndata as ad
import numpy as np
import pytest
import scipy.sparse as sp
from typer.testing import CliRunner

import sc_droplets as scd
from sc_droplets.cli import app
from sc_droplets.methods.generic import RuleError
from sc_droplets.primitives.rmath import clog2

runner = CliRunner()


def _out(result) -> str:
    return re.sub(r"\s+", " ", re.sub(r"\x1b\[[0-9;]*m", "", result.output))


def _matrix(totals) -> ad.AnnData:
    adata = ad.AnnData(sp.csr_matrix(np.asarray(totals, dtype=np.float64).reshape(-1, 1)))
    adata.obs_names = [f"B{i}" for i in range(len(totals))]
    adata.var_names = ["G0"]
    adata.var["feature_types"] = "Gene Expression"
    return adata


# ---------------------------------------------------------------------------
# Every rule
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("spec", ['mad'])
def test_generic_rules_need_a_barcode_at_the_floor(spec):
    with pytest.raises(RuleError, match="--min-counts"):
        scd.call(_matrix([1, 2, 3]), spec, inplace=False)


@pytest.mark.parametrize("name", ['mad'])
def test_cli_calls_the_generic_rule(name, sim, tmp_path):
    path = tmp_path / "sim.h5ad"
    sim.write_h5ad(path)
    result = runner.invoke(app, ["call", name, str(path), "--json"])
    assert result.exit_code == 0, _out(result)
    assert '"n_cells": 30' in _out(result)


# ---------------------------------------------------------------------------
# mad
# ---------------------------------------------------------------------------


def test_mad_limit_is_the_log2_median_less_k_mads_of_the_barcodes_above_the_floor():
    totals = [5, 50, 150, 300, 1000, 1100, 1200, 1300, 1500, 2000, 4000]
    result = scd.call(_matrix(totals), "mad", min_counts=100, nmads=2, inplace=False)
    candidates = np.array([v for v in totals if v >= 100], dtype=float)
    metric = clog2(candidates)
    median = np.median(metric)
    mad = 1.4826 * np.median(np.abs(metric - median))
    lower = median - 2 * mad
    assert result.thresholds["umis"] == 2.0 ** lower
    assert result.thresholds["median"] == 2.0 ** median and result.thresholds["mad"] == mad
    expected = [v >= 100 and math.log2(v) >= lower for v in totals]
    assert result.is_cell.tolist() == expected
    assert not result.is_cell[:2].any()  # below the floor: empty, whatever the limit


def test_mad_min_diff_widens_a_narrow_limit_and_no_log_works_on_counts():
    totals = [100, 990, 1000, 1000, 1010, 2000]
    narrow = scd.call(_matrix(totals), "mad", inplace=False)
    wide = scd.call(_matrix(totals), "mad", min_diff=4.0, inplace=False)
    assert narrow.n_cells < wide.n_cells == 6
    assert wide.thresholds["umis"] == 2.0 ** (math.log2(1000) - 4.0)
    plain = scd.call(_matrix(totals), "mad", log=False, inplace=False)
    assert plain.thresholds["median"] == 1000.0
