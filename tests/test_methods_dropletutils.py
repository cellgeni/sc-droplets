"""DropletUtils' methods on small and simulated matrices: the rules each release follows."""

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
from sc_droplets.methods.dropletutils import (
    DropletUtilsError,
    _BR_VARIANTS,
    _runs,
    barcode_ranks,
)

runner = CliRunner()


def _out(result) -> str:
    return re.sub(r"\s+", " ", re.sub(r"\x1b\[[0-9;]*m", "", result.output))


def _matrix(totals, n_features=20, seed=0) -> ad.AnnData:
    """Barcodes with the given totals, counts spread over features at random."""
    rng = np.random.default_rng(seed)
    rows = []
    for t in totals:
        rows.append(np.bincount(rng.integers(0, n_features, int(t)), minlength=n_features))
    adata = ad.AnnData(sp.csr_matrix(np.array(rows, dtype=np.float32)))
    adata.obs_names = [f"B{i}" for i in range(len(totals))]
    adata.var_names = [f"G{i}" for i in range(n_features)]
    adata.var["feature_types"] = "Gene Expression"
    return adata


def _variant(name):
    return next(v for v in _BR_VARIANTS if v.name == name)


# ---------------------------------------------------------------------------
# barcodeRanks
# ---------------------------------------------------------------------------


def test_runs_get_their_mid_rank():
    totals, rank = _runs(np.array([5, 9, 5, 5, 2, 9]))
    assert totals.tolist() == [9, 5, 2]
    assert rank.tolist() == [1.5, 4.0, 6.0]


def test_barcode_ranks_needs_three_runs_above_lower():
    with pytest.raises(DropletUtilsError, match="insufficient unique points"):
        barcode_ranks(np.array([500, 400, 50, 10]), _variant("v1.12").settings, lower=100, exclude_from=50, df=20)


def test_version_1_0_lands_the_knee_too_early(adata):
    """1.0 - 1.8 index the full curve with a position counted inside the fitted section."""
    old = scd.call(adata, "dropletutils-barcoderanks", tool_version="1.8", inplace=False)
    v110 = barcode_ranks(adata.X.sum(axis=1).A1, _variant("v1.10").settings, lower=100, exclude_from=None, df=20)
    # the same edges without exclude.from, but the knee from the fitted section
    assert old.thresholds["knee"] >= v110.knee


def test_each_knee_definition_gives_a_cut_on_the_curve(adata):
    knees = {}
    for version in ("1.0", "1.10", "1.12", "1.28", "1.30"):
        result = scd.call(adata, "dropletutils-barcoderanks", tool_version=version, inplace=False)
        knee, inflection = result.thresholds["knee"], result.thresholds["inflection"]
        assert 100 < inflection < knee < adata.X.sum(axis=1).max()
        assert result.n_cells == int((np.rint(adata.X.sum(axis=1).A1) >= knee).sum())
        knees[version] = knee
    assert len(set(knees.values())) >= 3  # the definitions do differ


def test_inflection_can_be_the_cut(adata):
    result = scd.call(adata, "dropletutils-barcoderanks", use="inflection", inplace=False)
    assert result.thresholds["umis"] == result.thresholds["inflection"]


@pytest.mark.parametrize(
    "version, option",
    [("1.26", {"window": 2.0}), ("1.30", {"df": 10.0}), ("1.8", {"exclude_from": 20}), ("1.30", {"fit_lower": 200.0, "fit_upper": 5000.0})],
)
def test_options_a_release_lacks_are_refused(adata, version, option):
    with pytest.raises(ValueError, match="has no"):
        scd.call(adata, "dropletutils-barcoderanks", tool_version=version, inplace=False, **option)


def test_fit_bounds_need_both_ends(adata):
    with pytest.raises(ValueError, match="both"):
        scd.call(adata, "dropletutils-barcoderanks", tool_version="1.26", fit_lower=200.0, inplace=False)


def test_fit_bounds_choose_the_fitted_section(adata):
    result = scd.call(adata, "dropletutils-barcoderanks", tool_version="1.26", fit_lower=300.0, fit_upper=20000.0, inplace=False)
    assert result.thresholds["knee"] > 300


def test_curve_tracing_falls_back_to_the_last_point_when_no_window_fits():
    totals = np.array([1000, 999, 998, 997] + [0] * 3)
    ranks = barcode_ranks(totals, _variant("v1.30").settings, lower=100, exclude_from=0, window=50.0, gradient_threshold=-1.0)
    assert ranks.knee == ranks.inflection == pytest.approx(997)


# ---------------------------------------------------------------------------
# emptyDrops
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def emptydrops_130(sim):
    adata = sim.copy()
    return adata, scd.call(adata, "dropletutils-emptydrops", niters=500, inplace=False)


def test_emptydrops_classes_follow_the_test(emptydrops_130):
    adata, result = emptydrops_130
    totals = np.rint(adata.X.sum(axis=1).A1)
    classes = np.asarray(result.classes)
    retain = result.thresholds["retain"]
    assert set(classes[totals >= retain]) == {"retained"}
    assert set(classes[totals <= 100]) == {"ambient"}
    fdr = result.columns["fdr"]
    assert np.all(fdr[classes == "rescued"] <= 0.001)
    assert np.all(fdr[classes == "tested"] > 0.001)
    assert np.all(np.isnan(result.columns["pvalue"][totals <= 100]))


def test_emptydrops_p_values_are_monte_carlo_fractions(emptydrops_130):
    _, result = emptydrops_130
    p = result.columns["pvalue"]
    p = p[~np.isnan(p)]
    n_above = p * 501 - 1
    assert np.allclose(n_above, np.round(n_above))
    assert np.array_equal(result.columns["limited"][~np.isnan(result.columns["pvalue"])], np.round(n_above) == 0)


def test_the_seed_changes_the_p_values_and_the_default_is_100(sim):
    a = scd.call(sim.copy(), "dropletutils-emptydrops", niters=200, inplace=False)
    b = scd.call(sim.copy(), "dropletutils-emptydrops", niters=200, random_state=100, inplace=False)
    c = scd.call(sim.copy(), "dropletutils-emptydrops", niters=200, random_state=1, inplace=False)
    assert np.array_equal(a.columns["pvalue"], b.columns["pvalue"], equal_nan=True)
    assert not np.array_equal(a.columns["pvalue"], c.columns["pvalue"], equal_nan=True)


def test_ignore_leaves_small_droplets_untested(sim):
    result = scd.call(sim.copy(), "dropletutils-emptydrops", niters=100, ignore=300.0, inplace=False)
    totals = np.rint(sim.X.sum(axis=1).A1)
    assert np.all(np.isnan(result.columns["pvalue"][totals <= 300]))
    assert "ignored" in set(np.asarray(result.classes))


def test_test_ambient_na_needs_1_12(sim):
    with pytest.raises(ValueError, match="1.12"):
        scd.call(sim.copy(), "dropletutils-emptydrops", tool_version="1.10", test_ambient="na", niters=50, inplace=False)


def test_ambient_droplets_tested_lose_their_fdr_from_1_12(sim):
    adata = sim.copy()
    totals = np.rint(adata.X.sum(axis=1).A1)
    new = scd.call(adata, "dropletutils-emptydrops", tool_version="1.30", test_ambient="true", niters=50, inplace=False)
    small = (totals > 0) & (totals <= 100)
    assert np.all(~np.isnan(new.columns["pvalue"][small]))
    assert np.all(np.isnan(new.columns["fdr"][small]))
    keep = scd.call(adata, "dropletutils-emptydrops", tool_version="1.30", test_ambient="na", niters=50, inplace=False)
    assert np.all(~np.isnan(keep.columns["fdr"][small]))


@pytest.mark.parametrize("version, option", [("1.8", {"by_rank": 100}), ("1.0", {"alpha": "50"}), ("1.6", {"round_counts": False})])
def test_emptydrops_options_a_release_lacks_are_refused(sim, version, option):
    with pytest.raises(ValueError, match="has no"):
        scd.call(sim.copy(), "dropletutils-emptydrops", tool_version=version, niters=10, inplace=False, **option)


def test_version_1_2s_dirichlet_draws_are_refused_but_its_multinomial_runs(sim):
    with pytest.raises(NotImplementedError, match="rgamma"):
        scd.call(sim.copy(), "dropletutils-emptydrops", tool_version="1.2", niters=10, inplace=False)
    result = scd.call(sim.copy(), "dropletutils-emptydrops", tool_version="1.2", alpha="inf", niters=10, inplace=False)
    assert result.diagnostics["alpha"] == "Inf"


def test_alpha_must_be_a_word_or_a_positive_number(sim):
    with pytest.raises(ValueError, match="positive"):
        scd.call(sim.copy(), "dropletutils-emptydrops", alpha="-3", niters=10, inplace=False)
    with pytest.raises(ValueError, match="estimate, inf"):
        scd.call(sim.copy(), "dropletutils-emptydrops", alpha="lots", niters=10, inplace=False)


def test_a_given_alpha_uses_the_dirichlet_multinomial(sim):
    result = scd.call(sim.copy(), "dropletutils-emptydrops", alpha="200", niters=20, inplace=False)
    assert result.diagnostics["alpha"] == 200.0


def test_by_rank_counts_ties_differently_from_1_20():
    """1.10 - 1.18 take every droplet tied with the (k+1)-th as ambient; 1.20 takes exactly those ranked below k."""
    totals = [3000, 2500, 2000] + [300] * 6 + [40] * 400 + [5] * 600
    adata = _matrix(totals, n_features=30, seed=4)
    old = scd.call(adata, "dropletutils-emptydrops", tool_version="1.18", by_rank=5, niters=20, retain=10000.0, inplace=False)
    new = scd.call(adata, "dropletutils-emptydrops", tool_version="1.20", by_rank=5, niters=20, retain=10000.0, inplace=False)
    assert old.diagnostics["n_ambient"] == 1006  # every 300 is <= the 6th largest, 300
    assert new.diagnostics["n_ambient"] == 1004  # the first two 300s are ranked 4 and 5
    assert new.thresholds["lower"] == 300.0


def test_by_rank_beyond_the_matrix_is_an_error(sim):
    with pytest.raises(DropletUtilsError, match="by.rank"):
        scd.call(sim.copy(), "dropletutils-emptydrops", by_rank=10**7, niters=10, inplace=False)


def test_an_empty_ambient_profile_stops_emptydrops():
    adata = _matrix([3000, 2000, 1500, 1200] + [0] * 50, n_features=10)
    with pytest.raises(DropletUtilsError, match="no counts"):
        scd.call(adata, "dropletutils-emptydrops", retain=10.0, niters=10, inplace=False)


def test_non_integer_counts_are_rounded_from_1_6_only_when_needed():
    adata = _matrix([3000, 2000] + [150] * 20 + [30] * 300, n_features=25, seed=2)
    adata.X = adata.X.astype(np.float64) * 1.0
    adata.X.data[::7] += 0.4  # rounds back
    a = scd.call(adata, "dropletutils-emptydrops", tool_version="1.8", niters=20, retain=10000.0, inplace=False)
    b = scd.call(adata, "dropletutils-emptydrops", tool_version="1.8", niters=20, retain=10000.0, round_counts=False, inplace=False)
    assert a.diagnostics["ambient_total"] == float(np.rint(adata.X[np.rint(adata.X.sum(axis=1).A1) <= 100].toarray()).sum())
    assert b.diagnostics["ambient_total"] != a.diagnostics["ambient_total"]


# ---------------------------------------------------------------------------
# emptyDropsCellRanger and defaultDrops
# ---------------------------------------------------------------------------


def test_emptydropscellranger_needs_1_14(sim):
    with pytest.raises(Exception, match="predates"):
        scd.call(sim.copy(), "dropletutils-emptydropscellranger", tool_version="1.12", inplace=False)


def test_emptydropscellranger_tests_the_droplets_below_the_cut():
    rng = np.random.default_rng(1)
    totals = list(rng.integers(5000, 20000, 60)) + list(rng.integers(600, 900, 40)) + list(rng.integers(1, 60, 3000))
    adata = _matrix(totals, n_features=40, seed=1)
    result = scd.call(adata, "dropletutils-emptydropscellranger", ind_min=200, ind_max=3000, niters=100, inplace=False)
    t = np.array(totals)
    retain = result.thresholds["umis"]
    order = np.argsort(-t, kind="stable")
    assert retain == max(round(t[order[29]] / 10), 1)
    tested = ~np.isnan(result.columns["pvalue"])
    assert tested.sum() == int((t > 500).sum()) - int((t >= retain).sum())
    assert set(np.asarray(result.classes)[t >= retain]) == {"initial"}


def test_the_1_14_bug_reads_the_median_from_the_unsorted_totals():
    """1.14 - 1.18.0 took totals[n/2] in matrix order: here a small droplet, so a lower minimum."""
    totals = [20] * 25 + [20000] * 40 + [800] * 50 + [20] * 2000
    adata = _matrix(totals, n_features=30, seed=3)
    common = dict(umi_min=50.0, ind_min=100, ind_max=2000, niters=10, inplace=False)
    old = scd.call(adata, "dropletutils-emptydropscellranger", tool_version="1.14", **common)
    new = scd.call(adata, "dropletutils-emptydropscellranger", tool_version="1.18.1", **common)
    # 40 initial cells; position 20 is 20 UMIs unsorted, 20000 sorted
    assert old.thresholds["candidate_umis"] == 50.0  # max(50, round(0.01 * 20))
    assert new.thresholds["candidate_umis"] == 200.0  # max(50, round(0.01 * 20000))


def test_defaultdrops_is_ten_percent_of_the_99th_percentile():
    totals = list(range(1, 4001))
    adata = _matrix(totals, n_features=5)
    result = scd.call(adata, "dropletutils-defaultdrops", inplace=False)
    top = np.arange(1001, 4001, dtype=float)
    h = 1 + 2999 * 0.99
    q = (1 - (h - math.floor(h))) * top[math.floor(h) - 1] + (h - math.floor(h)) * top[math.ceil(h) - 1]
    assert result.thresholds["umis"] == q * 0.1
    assert result.n_cells == int((np.array(totals) > q * 0.1).sum())


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------


def test_each_dropletutils_method_has_a_call_command():
    result = runner.invoke(app, ["call", "--help"])
    text = _out(result)
    for name in ("dropletutils-barcoderanks", "dropletutils-emptydrops", "dropletutils-emptydropscellranger", "dropletutils-defaultdrops"):
        assert name in text


def test_emptydrops_runs_from_the_command_line(tmp_path, sim):
    path = tmp_path / "sim.h5ad"
    sim.write_h5ad(path)
    result = runner.invoke(app, ["call", "dropletutils-emptydrops", str(path), "--niters", "50", "--tool-version", "1.28", "--json"])
    assert result.exit_code == 0, result.output
    assert '"method": "dropletutils-emptydrops"' in result.output


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("version, label", [("1.26", "smoothing spline (df 20)"), ("1.28", "chord"), ("1.30", "knee's window")])
def test_diagnostics_draw_what_found_the_knee(adata, version, label):
    import matplotlib

    matplotlib.use("Agg")
    fig = scd.pl.diagnostics(adata, f"dropletutils-barcoderanks@{version}")
    labels = fig.axes[0].get_legend_handles_labels()[1]
    assert label in labels
    assert any(text.startswith("cut:") for text in labels)


def test_emptydrops_diagnostics_show_retain_and_the_log_probabilities(adata):
    import matplotlib

    matplotlib.use("Agg")
    fig = scd.pl.diagnostics(adata, "dropletutils-emptydrops:niters=100")
    assert len(fig.axes) == 3
    first = fig.axes[0].get_legend_handles_labels()[1]
    assert any(text.startswith("retain (knee)") for text in first)
    assert any(text.startswith("called") for text in fig.axes[1].get_legend_handles_labels()[1])


def test_a_zero_ambient_probability_gives_minus_infinity_as_in_r():
    """Before 1.8, a profile without singletons left unseen features at probability 0; R's
    log(0) is -Inf, so droplets with those features get the smallest p-value. Python's math.log(0)
    raised instead."""
    rng = np.random.default_rng(5)
    ambient = [2 * np.bincount(rng.integers(0, 10, 40), minlength=20) for _ in range(300)]  # even counts only
    cells = [np.bincount(rng.integers(0, 20, 2000), minlength=20) for _ in range(5)]
    adata = ad.AnnData(sp.csr_matrix(np.array(ambient + cells, dtype=np.float32)))
    adata.var["feature_types"] = "Gene Expression"
    result = scd.call(adata, "dropletutils-emptydrops", tool_version="1.6", alpha="inf", retain=10**6, niters=50, inplace=False)
    tested = ~np.isnan(result.columns["pvalue"])
    assert tested.sum() == 5
    assert np.all(result.columns["logprob"][tested] == -np.inf)
    assert np.all(result.columns["pvalue"][tested] == 1 / 51)
