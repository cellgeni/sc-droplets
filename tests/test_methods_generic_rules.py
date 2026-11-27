"""MAD outliers, Otsu and multi-Otsu thresholds, and the kernel-density valley."""

from __future__ import annotations

import math
import re

import anndata as ad
import matplotlib
import numpy as np
import pytest
import scipy.sparse as sp
from typer.testing import CliRunner

import sc_droplets as scd
from sc_droplets.cli import app
from sc_droplets.methods.generic import GRID_POINTS, RuleError, density_valley
from sc_droplets.primitives.rmath import clog2

matplotlib.use("Agg")
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


@pytest.mark.parametrize("spec", ['mad', 'otsu', 'multiotsu:min_counts=10', 'kde-valley'])
def test_generic_rules_need_a_barcode_at_the_floor(spec):
    with pytest.raises(RuleError, match="--min-counts"):
        scd.call(_matrix([1, 2, 3]), spec, inplace=False)


@pytest.mark.parametrize("spec", ["mad", "otsu", "multiotsu", "kde-valley"])
def test_rules_on_a_log_scale_leave_out_zero_counts_whatever_the_floor(spec):
    """A floor of 0 used to make zero-count barcodes candidates: their log, -inf, broke the
    histogram and the density."""
    values, _ = _two_populations()
    with_zeros = _matrix(np.concatenate([np.zeros(50), values]))
    at_zero = scd.call(with_zeros, spec, min_counts=0, inplace=False)
    at_one = scd.call(with_zeros, spec, min_counts=1, inplace=False)
    assert at_zero.thresholds == at_one.thresholds
    assert not at_zero.is_cell[:50].any()


def test_rules_off_a_log_scale_keep_zero_counts_as_candidates():
    values, _ = _two_populations()
    with_zeros = _matrix(np.concatenate([np.zeros(50), values]))
    plain = scd.call(with_zeros, "otsu", min_counts=0, transform="none", inplace=False)
    assert plain.diagnostics["histogram_counts"].sum() == values.size + 50
    unlogged = scd.call(with_zeros, "mad", min_counts=0, log=False, inplace=False)
    assert unlogged.diagnostics["histogram_counts"].sum() == values.size + 50


@pytest.mark.parametrize("name", ['mad', 'otsu', 'multiotsu', 'kde-valley'])
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


# ---------------------------------------------------------------------------
# kde-valley
# ---------------------------------------------------------------------------


def test_kde_valley_lies_between_two_populations_and_calls_the_upper_one():
    values, upper = _two_populations()
    result = scd.call(_matrix(values), "kde-valley", inplace=False)
    assert 2.6 < result.thresholds["log10"] < 3.4
    assert np.array_equal(result.is_cell, values > result.thresholds["umis"])
    assert (result.is_cell == upper).mean() > 0.995


def test_kde_valley_density_is_scipys_on_a_10000_point_grid():
    from scipy.stats import gaussian_kde

    values, _ = _two_populations()
    x = np.log10(values[values >= 100])
    found = density_valley(x, 0.1)
    grid = np.linspace(x.min(), x.max(), GRID_POINTS)
    assert np.array_equal(found["grid"], grid)
    assert np.array_equal(found["density"], gaussian_kde(x, bw_method=0.1)(grid))
    low, high = found["bounds"]
    assert found["density"][found["valley"]] == found["density"][low:high + 1].min()


def test_kde_valley_needs_two_peaks_and_a_valid_bandwidth():
    values, _ = _two_populations(high=(3.8, 0.2, 0))
    with pytest.raises(RuleError, match="1 peak"):
        scd.call(_matrix(values), "kde-valley", inplace=False)
    with pytest.raises(RuleError, match="scott, silverman or a number"):
        scd.call(_matrix(values), "kde-valley", bandwidth="wide", inplace=False)
    with pytest.raises(RuleError, match="positive"):
        scd.call(_matrix(values), "kde-valley", bandwidth="-1", inplace=False)


def test_cli_reports_a_bad_bandwidth_as_an_error(sim, tmp_path):
    path = tmp_path / "sim.h5ad"
    sim.write_h5ad(path)
    bad = runner.invoke(app, ["call", "kde-valley", str(path), "--json", "--bandwidth", "wide"])
    assert bad.exit_code == 1 and "Error:" in _out(bad)


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------


def test_diagnostics_of_the_histogram_rules_draw_their_cuts(adata):
    for spec, n_cuts in (("mad", 2), ("otsu", 1), ("multiotsu", 2)):
        figure = scd.pl.diagnostics(adata, spec)
        panel = figure.axes[1]
        dashed = [line for line in panel.get_lines() if line.get_linestyle() == "--"]
        assert len(dashed) == n_cuts and panel.get_yscale() == "log"


def test_diagnostics_of_kde_valley_marks_the_peaks_and_the_valley(adata):
    figure = scd.pl.diagnostics(adata, "kde-valley")
    panel = figure.axes[1]
    labels = [line.get_label() for line in panel.get_lines()]
    assert "peaks" in labels and any(label.startswith("chosen valley") for label in labels)


def test_density_plot_has_a_line_at_each_calls_smallest_cell(adata):
    scd.call(adata, "otsu")
    scd.call(adata, "top-n", n_cells=100)
    ax = scd.pl.density(adata, ["otsu", "top_n"])
    totals = np.asarray(adata.X.sum(axis=1)).ravel()
    lines = {line.get_label().split(":")[0]: line.get_xdata()[0] for line in ax.get_lines()}
    assert lines["otsu"] == np.log10(totals[adata.obs["otsu"].to_numpy()].min())
    assert lines["top_n"] == np.log10(totals[adata.obs["top_n"].to_numpy()].min())
    assert ax.get_yscale() == "log" and ax.get_ylabel() == "barcodes"
    weighted = scd.pl.density(adata, "otsu", weighted=True)
    assert weighted.get_ylabel() == "share of UMIs"


def test_cli_plots_the_density_with_each_calls_cut(sim, tmp_path):
    path = tmp_path / "sim.h5ad"
    sim.write_h5ad(path)
    figure = tmp_path / "density.png"
    result = runner.invoke(app, ["plot", "density", str(path), "-m", "otsu", "-m", "kde-valley", "-o", str(figure)])
    assert result.exit_code == 0, _out(result)
    assert figure.stat().st_size > 0
