"""Several calls side by side (`pl.compare`), quality measures against the calls (`pl.qc`), and the
HTML report."""

from __future__ import annotations

import re

import anndata as ad
import matplotlib
import numpy as np
import pytest
from typer.testing import CliRunner

import sc_droplets as scd
from sc_droplets.cli import app
from sc_droplets.pl._data import curve_from_totals
from sc_droplets.pl.compare import strip_shares
from sc_droplets.pl.qc import disagreement_band

matplotlib.use("Agg")
runner = CliRunner()


def _out(result) -> str:
    return re.sub(r"\s+", " ", re.sub(r"\x1b\[[0-9;]*m", "", result.output))


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


# ---------------------------------------------------------------------------
# The report
# ---------------------------------------------------------------------------


def test_report_holds_every_section_inline(called, tmp_path):
    path = tmp_path / "report.html"
    page = scd.report(called, path=path, title="simulated")
    assert path.read_text() == page
    for anchor in ("summary", "knee", "density", "compare", "qc", "splicing", "calls"):
        assert f'id="{anchor}"' in page
    for key in scd.stored(called):
        assert f'id="call-{key}"' in page and f'id="call-{key}-figure"' in page
    assert page.count("<svg") == 5 + len(scd.stored(called))
    assert "unspliced: AUROC" in page and "mito: AUROC" in page
    assert "<script" not in page and 'src="http' not in page and "href=\"http" in page  # links only to method pages


def test_report_lists_parameters_thresholds_and_notes(sim):
    adata = sim.copy()
    scd.call(adata, "gmm", max_iter=1, tol=1e-12)
    page = scd.report(adata)
    assert "<td>max_iter</td>" in page and "<td>cell_component</td>" not in page
    assert 'class="note"' in page and "--max-iter 1" in page
    assert 'id="compare"' not in page and 'id="splicing"' not in page  # one call, no splicing


def test_report_names_the_input_file_not_its_path(sim, tmp_path):
    path = tmp_path / "deep" / "sample.h5ad"
    path.parent.mkdir()
    sim.write_h5ad(path)
    adata = scd.read(path)
    adata.uns["sc_droplets_input"] = {"path": str(path), "format": "h5ad"}
    scd.call(adata, "otsu")
    page = scd.report(adata)
    assert "sample.h5ad" in page and str(path.parent) not in page


def test_report_escapes_call_keys_in_its_anchors(sim):
    """A key with quotes and tags used to go into a figure's id as it was: a shared report could run
    a script."""
    adata = sim.copy()
    key = 'x"><img src=x onerror=alert(1)><figure id="x'
    scd.call(adata, "otsu", key_added=key)
    page = scd.report(adata)
    assert "<img" not in page and "onerror=alert(1)>" not in page
    assert 'id="call-x&quot;&gt;&lt;img src=x onerror=alert(1)&gt;&lt;figure id=&quot;x-figure"' in page


def test_report_needs_a_call(sim):
    with pytest.raises(ValueError, match="nothing to report"):
        scd.report(sim.copy())


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------


def test_cli_writes_the_compare_and_qc_plots_and_the_report(called, tmp_path):
    path = tmp_path / "calls.h5ad"
    called.write_h5ad(path)
    for args in (["plot", "compare", str(path)], ["plot", "qc", str(path), "--key", "otsu", "--key", "kde_valley"]):
        figure = tmp_path / f"{args[1]}.png"
        result = runner.invoke(app, args + ["-o", str(figure)])
        assert result.exit_code == 0, _out(result)
        assert figure.stat().st_size > 0
    page = tmp_path / "report.html"
    result = runner.invoke(app, ["report", str(path), "--key", "otsu", "-m", "gmm", "-o", str(page), "--title", "Simulated"])
    assert result.exit_code == 0, _out(result)
    text = page.read_text()
    assert "<title>Simulated</title>" in text and 'id="call-gmm"' in text and 'id="call-otsu"' in text
    assert 'id="call-kde_valley"' not in text  # stored, but not asked for


def test_cli_report_refuses_to_overwrite_its_input(called, tmp_path):
    path = tmp_path / "calls.h5ad"
    called.write_h5ad(path)
    result = runner.invoke(app, ["report", str(path), "-o", str(path)])
    assert result.exit_code == 1 and "also the input" in _out(result)
    assert ad.read_h5ad(path).n_obs == called.n_obs
