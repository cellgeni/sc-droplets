"""Plots: what they draw, checked through the Axes, not pixels."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

import sc_droplets as scd
from sc_droplets.pl._data import curve_from_totals


@pytest.fixture(autouse=True)
def _close_figures():
    yield
    plt.close("all")


def test_the_curve_steps_through_each_run_of_equal_totals():
    curve = curve_from_totals(np.array([1, 5, 0, 5, 2]))
    x, y = curve.steps()
    assert x.tolist() == [1, 2, 3, 3, 4, 4]
    assert y.tolist() == [5, 5, 2, 2, 1, 1]
    assert curve.rank_of().tolist() == [4, 1, 0, 2, 3]


def test_one_call_is_marked_on_the_curve_with_its_cut_off(adata):
    result = scd.call(adata, "threshold", min_counts=150)
    ax = scd.pl.knee(adata, "threshold")
    assert ax.get_xscale() == "log" and ax.get_yscale() == "log"
    cut_lines = [line for line in ax.get_lines() if line.get_linestyle() == "--"]
    assert cut_lines and cut_lines[0].get_ydata()[0] == pytest.approx(result.min_cell_total)
    assert len(ax.collections) == 1, "the called barcodes are drawn"
    assert "threshold: " in ax.get_legend().get_texts()[1].get_text()


def test_several_calls_each_get_a_line_and_a_legend_entry(adata):
    scd.compare(adata, ["top-n:n_cells=300", "threshold:min_counts=500"])
    ax = scd.pl.knee(adata)
    labels = [t.get_text() for t in ax.get_legend().get_texts()]
    assert labels[0] == "all barcodes"
    assert any(l.startswith("top_n_n_cells_300: 300 cells") for l in labels)
    assert len(ax.collections) == 0, "no per-barcode overlay with several calls"


def test_specs_not_yet_stored_are_run_first(adata):
    scd.pl.knee(adata, ["top-n:n_cells=50"])
    assert "top_n_n_cells_50" in adata.uns["sc_droplets"]


def test_rotated_swaps_the_axes(adata):
    scd.call(adata, "top-n", n_cells=100)
    ax = scd.pl.knee(adata, rotated=True)
    assert ax.get_xlabel() == "UMIs per barcode" and ax.get_ylabel() == "Barcode rank"


def test_cumulative_ends_at_one_and_marks_where_cells_end(adata):
    scd.call(adata, "top-n", n_cells=100)
    ax = scd.pl.knee(adata, cumulative=True)
    curve = ax.get_lines()[0]
    assert curve.get_ydata()[-1] == pytest.approx(1.0)
    assert ax.get_lines()[1].get_xdata()[0] == 100


def test_a_call_with_no_cells_is_still_listed(adata):
    scd.call(adata, "threshold", min_counts=1e9)
    ax = scd.pl.knee(adata)
    assert any("(none)" in t.get_text() for t in ax.get_legend().get_texts())


def test_an_unknown_key_with_running_off_is_an_error(adata):
    from sc_droplets.pl._data import resolve_calls

    with pytest.raises(KeyError, match="no stored call"):
        resolve_calls(adata, ["nope"], run_missing=False)


def test_the_terminal_knee_has_the_asked_size_and_names_its_marks():
    totals = np.concatenate([np.full(100, 5000), np.full(2000, 20), np.ones(10000)])
    lines = scd.pl.knee_text(totals, width=40, height=10, marks={"cut": 500})
    plot_rows = [l for l in lines if "┤" in l]
    assert len(plot_rows) == 10
    assert all(len(l.split("┤", 1)[1]) == 40 for l in plot_rows)
    assert lines[-1].strip() == "···· cut: 500"
    assert scd.pl.knee_text(np.zeros(3)) == ["(no barcode has any counts)"]


def test_diagnostics_of_a_threshold_method_shows_its_cut_and_robust_maximum(adata):
    scd.call(adata, "starsolo-cellranger2.2")
    figure = scd.pl.diagnostics(adata, "starsolo_cellranger2_2")
    (ax,) = figure.axes
    labels = [t.get_text() for t in ax.get_legend().get_texts()]
    assert any(l.startswith("cut: ") for l in labels) and any(l.startswith("robust maximum") for l in labels)


def test_diagnostics_of_ordmag_shades_its_bootstraps(adata):
    figure = scd.pl.diagnostics(adata, "cellranger-ordmag@2.1")
    labels = [t.get_text() for t in figure.axes[0].get_legend().get_texts()]
    assert any("OrdMag bootstraps" in l for l in labels)


def test_diagnostics_of_a_test_has_classes_likelihoods_and_p_values():
    adata = scd.simulate_droplets(n_cells=300, n_empty=3000, n_low=3000, n_genes=200, cell_umis=1500, cell_sigma=1.0, seed=7)
    scd.call(adata, "starsolo-emptydrops-cr", ind_min=2000, ind_max=5000, umi_min=150, sim_n=200)
    figure = scd.pl.diagnostics(adata, "starsolo_emptydrops_cr")
    assert len(figure.axes) == 3
    classes, loglik, pvalues = figure.axes
    assert any("ambient ranks" in t.get_text() for t in classes.get_legend().get_texts())
    assert loglik.get_xlabel() == "UMIs per candidate"
    assert "adjusted p ≤ 0.01" in pvalues.get_title()


def test_diagnostics_of_a_test_that_tested_nothing_says_so(adata):
    """Cell Ranger on a matrix with fewer barcodes than its ambient range tests no barcode; its
    candidates' panel used to put an empty scatter on a log axis, which matplotlib cannot draw,
    and the diagnostics and the report failed."""
    scd.call(adata, "cellranger")
    assert not np.isfinite(adata.obs["cellranger_loglik" if "cellranger_loglik" in adata.obs else "cellranger_logprob"]).any()
    figure = scd.pl.diagnostics(adata, "cellranger")
    figure.canvas.draw()
    assert any(t.get_text() == "no barcode was tested" for t in figure.axes[1].texts)
    assert 'id="call-cellranger"' in scd.report(adata)
