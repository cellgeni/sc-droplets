"""The unspliced and mitochondrial fractions and their sources."""
from __future__ import annotations


import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
from conftest import write_alevin_fry, write_kb, write_loom, write_starsolo

import sc_droplets as scd
from sc_droplets.io.loom import is_loom
from sc_droplets.metrics import MetricError, barcode_key


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
