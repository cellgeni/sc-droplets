"""Several calls side by side: `pl.compare`."""

from __future__ import annotations

import matplotlib
import numpy as np
import pytest

import sc_droplets as scd
from sc_droplets.pl._data import curve_from_totals
from sc_droplets.pl.compare import strip_shares

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
