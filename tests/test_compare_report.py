"""Several calls side by side (`pl.compare`) and quality measures against the calls (`pl.qc`)."""

from __future__ import annotations

import matplotlib
import numpy as np
import pytest

import sc_droplets as scd
from sc_droplets.pl._data import curve_from_totals
from sc_droplets.pl.compare import strip_shares
from sc_droplets.pl.qc import disagreement_band

matplotlib.use("Agg")


@pytest.fixture
def called(spliced):
    """Three calls on the simulated matrix with splicing: two exact, one that takes in empties."""
    adata = spliced.copy()
    names = list(adata.var_names)
    names[:2] = ["MT-CO1", "MT-ND1"]
    adata.var_names = names
    scd.add_metrics(adata)
    for spec in ("otsu", "kde-valley", "top-n:n_cells=600"):
        scd.call(adata, spec)
    return adata


# ---------------------------------------------------------------------------
# pl.compare
# ---------------------------------------------------------------------------


def test_strip_shares_are_the_called_share_of_each_rank_bin():
    totals = np.array([9.0, 8, 7, 6, 5, 4, 3, 2, 1, 0])
    curve = curve_from_totals(totals)
    cells = totals >= 5  # the top five ranks
    edges, shares = strip_shares(cells, curve, bins=4)
    assert edges[0] == 1 and edges[-1] == curve.n
    for lo, hi, share in zip(edges[:-1], edges[1:], shares):
        ranks = np.arange(lo, hi + 1)
        assert share == np.mean(ranks <= 5)


def test_compare_plot_shows_strips_counts_and_the_jaccard_index(called):
    keys = ["otsu", "kde_valley", "top_n_n_cells_600"]
    figure = scd.pl.compare(called, keys)
    knee_ax, strip_ax, bar_ax, jaccard_ax = figure.axes[:4]
    assert [t.get_text() for t in strip_ax.get_yticklabels()] == keys
    assert [p.get_width() for p in bar_ax.patches] == [int(called.obs[k].sum()) for k in keys]
    image = jaccard_ax.get_images()[0].get_array()
    assert np.allclose(image, scd.agreement(called, keys).to_numpy())
    assert len([line for line in knee_ax.get_lines() if line.get_linestyle() == "--"]) == 3


def test_compare_plot_needs_a_call(sim):
    with pytest.raises(ValueError, match="no calls to compare"):
        scd.pl.compare(sim.copy())


# ---------------------------------------------------------------------------
# pl.qc
# ---------------------------------------------------------------------------


def test_disagreement_band_runs_between_the_calls_smallest_cells(called):
    totals = np.asarray(called.X.sum(axis=1)).ravel()
    keys = ["otsu", "top_n_n_cells_600"]
    smallest = sorted(float(totals[called.obs[k].to_numpy()].min()) for k in keys)
    assert disagreement_band(called, keys, totals) == tuple(smallest)
    assert disagreement_band(called, ["otsu"], totals) is None


def test_qc_plot_has_genes_mito_and_the_disputed_barcodes(called):
    figure = scd.pl.qc(called, ["otsu", "top_n_n_cells_600"])
    assert [ax.get_ylabel() for ax in figure.axes] == ["log10 genes detected", "mitochondrial fraction", "density"]
    legend = [t.get_text() for t in figure.axes[2].get_legend().get_texts()]
    assert any(text.startswith("no call in the band") for text in legend)
    one = scd.pl.qc(called, "otsu")  # one call: no disagreement panel
    assert len(one.axes) == 2


def test_qc_plot_without_mito_has_no_mito_panel(sim):
    adata = sim.copy()
    scd.call(adata, "otsu")
    assert [ax.get_ylabel() for ax in scd.pl.qc(adata).axes] == ["log10 genes detected"]
