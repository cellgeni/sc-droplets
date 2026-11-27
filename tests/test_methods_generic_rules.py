"""MAD outliers, Otsu and multi-Otsu thresholds."""

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


def _two_populations(low=(2.3, 0.15, 3000), high=(3.8, 0.2, 500)):
    """Totals of two log-normal populations, at evenly spaced quantiles (no sampling noise), and
    which barcodes are in the upper one."""
    from scipy.stats import norm

    def population(mean, sd, n):
        return 10 ** norm.ppf((np.arange(n) + 0.5) / n, mean, sd)

    values = np.concatenate([population(*low), population(*high)])
    return np.round(values), np.r_[np.zeros(low[2], bool), np.ones(high[2], bool)]


# ---------------------------------------------------------------------------
# Every rule
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("spec", ['mad', 'otsu', 'multiotsu:min_counts=10'])
def test_generic_rules_need_a_barcode_at_the_floor(spec):
    with pytest.raises(RuleError, match="--min-counts"):
        scd.call(_matrix([1, 2, 3]), spec, inplace=False)


@pytest.mark.parametrize("name", ['mad', 'otsu', 'multiotsu'])
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


# ---------------------------------------------------------------------------
# Otsu and multi-Otsu
# ---------------------------------------------------------------------------


def test_otsu_splits_two_clusters_between_them_and_calls_the_upper_one():
    values, upper = _two_populations()
    result = scd.call(_matrix(values), "otsu", inplace=False)
    cut = result.thresholds["transformed"]
    assert 2.6 < cut < 3.4 and result.thresholds["umis"] == 10.0 ** cut
    assert np.array_equal(result.is_cell, np.log10(values) > cut)
    assert (result.is_cell == upper).mean() > 0.995


def test_otsu_of_identical_values_calls_nothing():
    assert scd.call(_matrix([500] * 10), "otsu", inplace=False).n_cells == 0


def test_otsu_arcsinh_transform_reports_the_cut_back_in_counts():
    values, _ = _two_populations()
    result = scd.call(_matrix(values), "otsu", transform="arcsinh", inplace=False)
    assert result.thresholds["umis"] == pytest.approx(math.sinh(result.thresholds["transformed"]))


def test_multiotsu_stores_each_barcodes_class_and_calls_the_top_one(sim):
    result = scd.call(sim, "multiotsu", inplace=False)
    classes = np.asarray(result.classes.astype(str))
    assert set(classes) <= {"below floor", "class 1", "class 2", "class 3"}
    assert np.array_equal(result.is_cell, classes == "class 3")
    assert (classes == "below floor").sum() == (np.asarray(result.totals) < 1).sum()
    top = 10.0 ** result.thresholds["threshold_2"]
    assert result.thresholds["umis"] == top


def test_multiotsu_refuses_more_than_four_classes(sim):
    with pytest.raises(ValueError, match="classes"):
        scd.call(sim, "multiotsu", classes=5, inplace=False)
