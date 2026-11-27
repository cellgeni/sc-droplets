"""The unspliced and mitochondrial fractions: their sources, their separation scores, and the plots
that use them."""
from __future__ import annotations


import anndata as ad
import matplotlib
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from conftest import write_alevin_fry, write_kb, write_loom, write_starsolo

import sc_droplets as scd
from sc_droplets.io.loom import is_loom
from sc_droplets.metrics import MetricError, auroc, barcode_key

matplotlib.use("Agg")

def _expected(adata) -> np.ndarray:
    s = np.asarray(adata.layers["spliced"].sum(axis=1)).ravel().astype(float)
    u = np.asarray(adata.layers["unspliced"].sum(axis=1)).ravel().astype(float)
    out = np.full(s.shape, np.nan)
    np.divide(u, s + u, out=out, where=(s + u) > 0)
    return out


# ---------------------------------------------------------------------------
# Sources
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name, key", [
    ("ACGTACGT-1", "ACGTACGT"),
    ("sample1:ACGTACGTx", "ACGTACGT"),
    ("ACGTACGT", "ACGTACGT"),
    ("run_2:ACGTACGT-3", "ACGTACGT"),
])
def test_barcode_key_drops_velocyto_and_gem_group_decorations(name, key):
    assert barcode_key(name) == key


def test_unspliced_fraction_from_layers_counts_gene_expression_only(spliced):
    adata = spliced.copy()
    values = scd.unspliced_fraction(adata)
    assert np.allclose(values.to_numpy(), _expected(adata), equal_nan=True)
    assert values.isna().sum() == (_expected(adata) != _expected(adata)).sum()
    # An antibody column must not count.
    extra = ad.AnnData(sp.csr_matrix(np.full((adata.n_obs, 1), 50)), obs=adata.obs[[]])
    extra.var_names = ["CD3"]
    extra.var["feature_types"] = "Antibody Capture"
    for name in ("spliced", "unspliced"):
        extra.layers[name] = sp.csr_matrix(np.full((adata.n_obs, 1), 50))
    both = ad.concat([adata, extra], axis=1, merge="first")
    assert np.allclose(scd.unspliced_fraction(both).to_numpy(), _expected(adata), equal_nan=True)


def test_a_loom_reads_as_velocyto_writes_it(spliced, tmp_path):
    path = write_loom(spliced, tmp_path / "sample.loom")
    assert is_loom(path)
    adata = scd.read(path)
    assert adata.obs_names[0] == "sample1:" + spliced.obs_names[0].split("-")[0] + "x"
    assert (adata.layers["spliced"] != spliced.layers["spliced"]).nnz == 0
    assert (adata.layers["unspliced"] != spliced.layers["unspliced"]).nnz == 0
    assert (adata.X != spliced.layers["spliced"] + spliced.layers["unspliced"]).nnz == 0


def test_a_loom_needs_spliced_and_unspliced_layers(spliced, tmp_path):
    import h5py

    path = write_loom(spliced, tmp_path / "sample.loom")
    with h5py.File(path, "a") as handle:
        del handle["layers/unspliced"]
    with pytest.raises(ValueError, match="layers/unspliced"):
        scd.read(path)


@pytest.mark.parametrize("layout", ["loom", "velocyto", "kb", "alevin-fry"])
def test_every_splicing_layout_gives_the_same_fraction(spliced, tmp_path, layout):
    if layout == "loom":
        source = write_loom(spliced, tmp_path / "sample.loom")
    elif layout == "velocyto":
        source = write_starsolo(spliced, tmp_path / "star", velocyto=True)
    elif layout == "kb":
        source = write_kb(spliced, tmp_path / "counts_unfiltered", nac=True)
    else:
        source = write_alevin_fry(spliced, tmp_path / "alevin", usa=True)
    values = scd.unspliced_fraction(source)
    assert np.allclose(values.to_numpy(), _expected(spliced), equal_nan=True)


def test_starsolo_without_velocyto_estimates_from_genefull_less_gene(spliced, tmp_path):
    exonic = spliced.copy()
    exonic.X = spliced.layers["spliced"]
    write_starsolo(exonic, tmp_path / "star", feature="Gene")
    write_starsolo(spliced, tmp_path / "star_full", feature="GeneFull")
    (tmp_path / "star_full" / "Solo.out" / "GeneFull").rename(tmp_path / "star" / "Solo.out" / "GeneFull")
    values = scd.unspliced_fraction(tmp_path / "star")
    assert np.allclose(values.to_numpy(), _expected(spliced), equal_nan=True)


def test_a_source_without_splicing_is_refused(sim, tmp_path):
    with pytest.raises(MetricError, match="no 'spliced' and 'unspliced' layers"):
        scd.unspliced_fraction(sim)
    source = write_kb(sim, tmp_path / "counts_unfiltered")
    with pytest.raises(MetricError, match="no spliced and unspliced counts"):
        scd.unspliced_fraction(source)


# ---------------------------------------------------------------------------
# add_metrics
# ---------------------------------------------------------------------------


def test_add_metrics_matches_a_separate_source_by_barcode(spliced, sim, tmp_path):
    source = write_loom(spliced, tmp_path / "sample.loom")  # CellIDs like sample1:ACGTx
    adata = sim.copy()
    assert scd.add_metrics(adata, splicing=source) == ["unspliced_fraction"]
    assert np.allclose(adata.obs["unspliced_fraction"].to_numpy(), _expected(spliced), equal_nan=True)


def test_add_metrics_leaves_unmatched_barcodes_unmeasured(spliced, sim):
    values = scd.unspliced_fraction(spliced).iloc[:100]
    adata = sim.copy()
    scd.add_metrics(adata, splicing=values)
    assert adata.obs["unspliced_fraction"].iloc[100:].isna().all()
    unrelated = pd.Series([0.5], index=["TTTTTTTTTTTTTTTT-1"])
    with pytest.raises(MetricError, match="no barcode of the splicing source matches"):
        scd.add_metrics(sim.copy(), splicing=unrelated)


def test_add_metrics_takes_the_mitochondrial_fraction_from_mt_genes(spliced):
    adata = spliced.copy()
    names = list(adata.var_names)
    names[:3] = ["MT-CO1", "MT-ND1", "MT-ATP6"]
    adata.var_names = names
    assert scd.add_metrics(adata) == ["unspliced_fraction", "mito_fraction"]
    x = sp.csr_matrix(adata.X)
    mito = np.asarray(x[:, :3].sum(axis=1)).ravel()
    total = np.asarray(x.sum(axis=1)).ravel()
    expected = np.where(total > 0, mito / np.where(total > 0, total, 1), np.nan)
    assert np.allclose(adata.obs["mito_fraction"].to_numpy(), expected, equal_nan=True)
    assert scd.add_metrics(scd.simulate_droplets(n_cells=10, n_empty=10, n_low=10, n_genes=20, seed=1)) == []


# ---------------------------------------------------------------------------
# Separation scores
# ---------------------------------------------------------------------------


def test_auroc_is_the_share_of_pairs_ordered_ties_counting_half():
    assert auroc(np.array([3.0, 4.0]), np.array([1.0, 2.0])) == 1.0
    assert auroc(np.array([1.0, 2.0]), np.array([3.0, 4.0])) == 0.0
    assert auroc(np.array([2.0]), np.array([2.0, 1.0])) == 0.75
    assert np.isnan(auroc(np.array([]), np.array([1.0])))


def test_separation_compares_a_calls_cells_with_the_barcodes_it_left_above_the_floor():
    totals = [500, 400, 300, 50, 20, 5]
    adata = ad.AnnData(sp.csr_matrix(np.asarray(totals, float).reshape(-1, 1)))
    adata.obs_names = [f"B{i}" for i in range(6)]
    adata.var_names = ["G0"]
    adata.obs["unspliced_fraction"] = [0.3, 0.25, 0.1, 0.05, 0.2, 0.9]
    scd.call(adata, "top-n", n_cells=2)
    row = scd.separation(adata, min_counts=10).loc["top_n"]
    assert row["unspliced_fraction_cells"] == 0.275
    assert row["unspliced_fraction_others"] == 0.1  # 0.1, 0.05, 0.2: the barcode with 5 is below the floor
    assert row["unspliced_fraction_auroc"] == 1.0


def test_separation_needs_the_metric(sim):
    adata = sim.copy()
    scd.call(adata, "otsu")
    with pytest.raises(MetricError, match="add_metrics"):
        scd.separation(adata)


def test_compare_scores_each_call_when_the_fractions_are_there(spliced):
    adata = spliced.copy()
    scd.add_metrics(adata)
    summary = scd.compare(adata, ["otsu", "top-n:n_cells=600"])
    assert summary.loc["otsu", "unspliced_fraction_auroc"] > summary.loc["top_n_n_cells_600", "unspliced_fraction_auroc"]
    assert "mito_fraction_auroc" not in summary.columns
    plain = scd.compare(scd.simulate_droplets(n_cells=50, n_empty=200, n_low=500, n_genes=50, seed=2), ["otsu"])
    assert "unspliced_fraction_auroc" not in plain.columns


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------


@pytest.fixture
def measured(spliced):
    adata = spliced.copy()
    scd.add_metrics(adata)
    return adata


def test_splicing_plot_has_a_disagreement_panel_only_for_several_calls(measured):
    one = scd.pl.splicing(measured, "otsu")
    assert len(one.axes) == 2
    several = scd.pl.splicing(measured, ["otsu", "top-n:n_cells=600"])
    assert len(several.axes) == 3
    labels = [t.get_text() for t in several.axes[1].get_xticklabels()]
    assert all("AUROC" in label for label in labels)
    legend = [t.get_text() for t in several.axes[0].get_legend().get_texts()]
    assert any(text.startswith("some calls") for text in legend)


def test_splicing_plot_needs_the_metric(sim):
    with pytest.raises(KeyError, match="add_metrics"):
        scd.pl.splicing(sim.copy(), "otsu")
