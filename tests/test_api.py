"""Counting, the two generic methods, and how calls are stored on AnnData."""

from __future__ import annotations

import json

import anndata as ad
import numpy as np
import pytest
import scipy.sparse as sp

import sc_droplets as scd
from sc_droplets.counts import Droplets, looks_filtered


def _tiny(totals, feature_types=None):
    """An AnnData whose row sums are `totals` (one count per gene column)."""
    totals = np.asarray(totals, dtype=np.int64)
    n_genes = 3
    rows = np.repeat(np.arange(len(totals)), totals)
    cols = np.arange(rows.size) % n_genes
    matrix = sp.csr_matrix((np.ones(rows.size, dtype=np.int64), (rows, cols)), shape=(len(totals), n_genes))
    adata = ad.AnnData(X=matrix)
    adata.obs_names = [f"B{i}" for i in range(len(totals))]
    if feature_types is not None:
        adata.var["feature_types"] = feature_types
    return adata


# --- counting ---------------------------------------------------------------


def test_only_the_requested_feature_types_are_counted():
    adata = _tiny([6, 3], feature_types=["Gene Expression", "Gene Expression", "Antibody Capture"])
    gex = Droplets.from_anndata(adata, feature_types="Gene Expression")
    assert gex.totals.tolist() == [4, 2]
    assert gex.feature_types == ("Gene Expression",)
    assert Droplets.from_anndata(adata, feature_types="all").totals.tolist() == [6, 3]
    with pytest.raises(ValueError, match="no features of type 'CRISPR'.*'Gene Expression'"):
        Droplets.from_anndata(adata, feature_types=["CRISPR"])


def test_a_matrix_without_feature_types_counts_every_feature():
    adata = _tiny([6, 3])
    assert Droplets.from_anndata(adata, feature_types="Gene Expression").totals.tolist() == [6, 3]


def test_genome_selects_one_genomes_features(sim):
    adata = sim.copy()
    adata.var["genome"] = ["hg"] * 100 + ["mm"] * 100
    human = Droplets.from_anndata(adata, genome="hg")
    assert human.n_features == 100
    with pytest.raises(ValueError, match="no genome 'rn'"):
        Droplets.from_anndata(adata, genome="rn")
    with pytest.raises(ValueError, match="no 'genome' column"):
        Droplets.from_anndata(_tiny([1]), genome="hg")


def test_layers_reads_and_genes_are_available_to_methods():
    adata = _tiny([6, 3])
    adata.layers["double"] = adata.X * 2
    adata.obs["n_reads"] = [60.0, 30.0]
    droplets = Droplets.from_anndata(adata, layer="double", reads="n_reads")
    assert droplets.totals.tolist() == [12, 6]
    assert droplets.values("reads").tolist() == [60, 30]
    assert droplets.values("genes").tolist() == [3, 3]
    with pytest.raises(ValueError, match="needs read counts"):
        Droplets.from_anndata(adata).values("reads")
    with pytest.raises(ValueError, match="no obs column 'x'"):
        Droplets.from_anndata(adata, reads="x")


def test_explicit_zeros_are_not_genes_detected():
    matrix = sp.csr_matrix((np.array([0, 2]), np.array([0, 1]), np.array([0, 2])), shape=(1, 2))
    adata = ad.AnnData(X=matrix)
    assert Droplets.from_anndata(adata).values("genes").tolist() == [1]


def test_a_filtered_matrix_is_recognised(sim):
    assert looks_filtered(Droplets.from_anndata(sim).totals) is None
    cells = sim[sim.obs["truth"] == "cell"].copy()
    assert "filtered" in looks_filtered(Droplets.from_anndata(cells).totals)
    assert "zero counts" in looks_filtered(np.zeros(5))


# --- generic methods --------------------------------------------------------


def test_top_n_breaks_ties_by_input_order_or_keeps_them_all():
    adata = _tiny([5, 3, 3, 3, 1])
    first = scd.call(adata, "top-n", n_cells=2, key_added="first")
    assert adata.obs["first"].tolist() == [True, True, False, False, False]
    assert first.thresholds == {"umis": 3.0}
    scd.call(adata, "top-n", n_cells=2, ties="include", key_added="all")
    assert adata.obs["all"].tolist() == [True, True, True, True, False]


def test_top_n_of_zero_or_more_than_there_are():
    adata = _tiny([5, 3])
    assert scd.call(adata, "top-n", n_cells=0, inplace=False).n_cells == 0
    assert scd.call(adata, "top-n", n_cells=10, inplace=False).n_cells == 2


def test_threshold_is_inclusive_unless_strict_and_counts_any_unit():
    adata = _tiny([5, 3, 2])
    assert scd.call(adata, "threshold", min_counts=3, inplace=False).is_cell.tolist() == [True, True, False]
    assert scd.call(adata, "threshold", min_counts=3, strict=True, inplace=False).is_cell.tolist() == [True, False, False]
    genes = scd.call(adata, "threshold", min_counts=3, unit="genes", inplace=False)
    assert genes.is_cell.tolist() == [True, True, False]


def test_simple_methods_separate_simulated_cells_from_empties(sim):
    truth = (sim.obs["truth"] == "cell").to_numpy()
    result = scd.call(sim.copy(), "top-n", n_cells=int(truth.sum()))
    assert (result.is_cell == truth).mean() > 0.999


# --- storing calls ----------------------------------------------------------


def test_a_call_is_stored_in_obs_and_uns(adata):
    result = scd.call(adata, "threshold", min_counts=100)
    assert result.key == "threshold"
    assert adata.obs["threshold"].dtype == bool
    record = adata.uns["sc_droplets"]["threshold"]
    assert record["method"] == "threshold"
    assert record["n_cells"] == result.n_cells
    assert record["params"]["min_counts"] == 100.0
    assert record["counting"]["feature_types"] == ["Gene Expression"]


def test_stored_calls_survive_writing_and_reading_h5ad(adata, tmp_path):
    scd.call(adata, "top-n", n_cells=10)
    scd.call(adata, "threshold", min_counts=50)
    adata.write_h5ad(tmp_path / "calls.h5ad")
    back = ad.read_h5ad(tmp_path / "calls.h5ad")
    assert set(scd.stored(back)) == {"top_n", "threshold"}
    assert back.obs["top_n"].sum() == 10
    assert back.uns["sc_droplets"]["top_n"]["params"]["n_cells"] == 10


def test_storing_again_under_a_key_replaces_the_old_call(adata):
    scd.call(adata, "top-n", n_cells=10, key_added="k")
    scd.call(adata, "threshold", min_counts=50, key_added="k")
    assert adata.uns["sc_droplets"]["k"]["method"] == "threshold"


def test_a_result_cannot_be_stored_on_different_barcodes(adata, sim):
    result = scd.call(adata, "top-n", n_cells=5, inplace=False)
    other = sim[:10].copy()
    with pytest.raises(ValueError, match="different barcodes"):
        result.store(other, "x")


def test_the_table_and_json_summary_describe_the_call(adata):
    result = scd.call(adata, "top-n", n_cells=5, inplace=False)
    table = result.table
    assert list(table.columns) == ["total", "is_cell"]
    assert table["is_cell"].sum() == 5
    summary = json.loads(result.to_json())
    assert summary["n_cells"] == 5 and summary["params"]["n_cells"] == 5
    assert result.key is None, "not stored with inplace=False"


def test_spec_parameters_and_keyword_parameters_combine(adata):
    result = scd.call(adata, "top-n:n_cells=7,ties=include")
    assert result.key == "top_n_n_cells_7_ties_include"
    assert result.params["ties"] == "include"
    override = scd.call(adata, "top-n:n_cells=7", n_cells=3, inplace=False)
    assert override.n_cells == 3


def test_a_tool_version_on_a_generic_method_is_refused(adata):
    with pytest.raises(scd.UnsupportedVersion, match="takes no tool version"):
        scd.call(adata, "top-n@1.0", n_cells=3)
    with pytest.raises(ValueError, match="names version"):
        scd.call(adata, "top-n@1.0", tool_version="2.0", n_cells=3)


def test_compare_runs_each_spec_and_agreement_compares_them(adata):
    summary = scd.compare(adata, ["top-n:n_cells=300", "threshold:min_counts=150", "threshold:min_counts=1000"])
    assert list(summary.index) == ["top_n_n_cells_300", "threshold_min_counts_150", "threshold_min_counts_1000"]
    jaccard = scd.agreement(adata)
    assert np.allclose(np.diag(jaccard), 1.0)
    assert jaccard.loc["top_n_n_cells_300", "threshold_min_counts_1000"] < 1.0
    with pytest.raises(KeyError, match="no stored call nope"):
        scd.agreement(adata, ["nope"])


def test_list_methods_describes_every_registered_method():
    table = scd.list_methods()
    assert {"top-n", "threshold"} <= set(table.index)
    assert table.loc["top-n", "docs"].endswith("/methods/generic/top-n/")


def test_the_simulator_is_deterministic_and_labelled():
    a = scd.simulate_droplets(n_cells=20, n_empty=50, n_low=80, n_genes=30, seed=1)
    b = scd.simulate_droplets(n_cells=20, n_empty=50, n_low=80, n_genes=30, seed=1)
    assert (a.X != b.X).nnz == 0 and list(a.obs_names) == list(b.obs_names)
    assert a.obs["truth"].value_counts().to_dict() == {"low": 80, "empty": 50, "cell": 20}
    spliced = scd.simulate_droplets(n_cells=20, n_empty=50, n_low=80, n_genes=30, seed=1, splicing=True)
    assert ((spliced.layers["spliced"] + spliced.layers["unspliced"]) != spliced.X).nnz == 0
    big = scd.simulate_droplets(n_cells=20, n_empty=50, n_low=80, n_genes=30, seed=1, ensure_ambient_range=500)
    assert big.n_obs == 500
