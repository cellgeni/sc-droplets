"""Readers, checked against scanpy where scanpy reads the same format, and writers."""

from __future__ import annotations

import gzip

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

import sc_droplets as scd
from sc_droplets.io import find_matrix
from sc_droplets.io.write import check_outputs, write_barcodes, write_table

from conftest import (
    write_10x_h5_v2,
    write_10x_h5_v3,
    write_10x_mtx,
    write_alevin_fry,
    write_kb,
    write_starsolo,
)

sc = pytest.importorskip("scanpy")


def _same_counts(a, b) -> None:
    assert list(a.obs_names) == list(b.obs_names)
    assert (sp.csr_matrix(a.X) != sp.csr_matrix(b.X)).nnz == 0


def test_10x_h5_v3_reads_as_scanpy_reads_it(sim, tmp_path):
    path = write_10x_h5_v3(sim, tmp_path / "raw_feature_bc_matrix.h5")
    ours = scd.read(path)
    theirs = sc.read_10x_h5(path, gex_only=False)
    _same_counts(ours, theirs)
    assert list(ours.var_names) == list(theirs.var_names)
    assert list(ours.var["gene_ids"]) == list(theirs.var["gene_ids"])
    assert list(ours.var["feature_types"]) == list(theirs.var["feature_types"])
    assert ours.X.dtype.kind == "i", "counts stay integers"


def test_10x_h5_v2_reads_as_scanpy_reads_it(sim, tmp_path):
    path = write_10x_h5_v2(sim, tmp_path / "raw_gene_bc_matrices_h5.h5")
    ours = scd.read(path)
    theirs = sc.read_10x_h5(path, genome="GRCh38")
    _same_counts(ours, theirs)
    assert list(ours.var["gene_ids"]) == list(theirs.var["gene_ids"])


def test_a_barnyard_v2_file_is_read_whole_with_each_genome_marked(sim, tmp_path):
    path = write_10x_h5_v2(sim, tmp_path / "barnyard.h5", genomes=("hg19", "mm10"))
    both = scd.read(path)
    assert both.n_vars == sim.n_vars
    assert set(both.var["genome"]) == {"hg19", "mm10"}
    _same_counts(both, sim)
    mouse = scd.read(path, genome="mm10")
    assert set(mouse.var["genome"]) == {"mm10"}
    with pytest.raises(ValueError, match="no genome 'rn6'"):
        scd.read(path, genome="rn6")


@pytest.mark.parametrize("version", [2, 3])
def test_10x_mtx_folders_read_as_scanpy_reads_them(sim, tmp_path, version):
    folder = write_10x_mtx(sim, tmp_path / "raw_feature_bc_matrix", version=version)
    ours = scd.read(folder)
    theirs = sc.read_10x_mtx(folder, var_names="gene_symbols", gex_only=False)
    _same_counts(ours, theirs)
    assert list(ours.var_names) == list(theirs.var_names)


def test_a_starsolo_folder_is_found_from_its_parent_solo_out_or_feature(sim, tmp_path):
    solo = write_starsolo(sim, tmp_path / "sample", feature="GeneFull")
    for start in (tmp_path / "sample", solo, solo / "GeneFull"):
        assert find_matrix(start, solo_feature="GeneFull") == solo / "GeneFull" / "raw"
    _same_counts(scd.read(solo, solo_feature="GeneFull"), sim)


def test_starsolo_velocyto_reads_spliced_and_unspliced(tmp_path):
    adata = scd.simulate_droplets(n_cells=50, n_empty=200, n_low=300, n_genes=50, splicing=True, seed=3)
    solo = write_starsolo(adata, tmp_path / "s", velocyto=True)
    velo = scd.read(solo, solo_feature="Velocyto")
    assert (velo.layers["spliced"] != adata.layers["spliced"]).nnz == 0
    assert (velo.layers["unspliced"] != adata.layers["unspliced"]).nnz == 0
    assert (sp.csr_matrix(velo.X) != sp.csr_matrix(adata.X)).nnz == 0


@pytest.mark.parametrize("nac", [False, True])
def test_kb_counts_unfiltered_is_read_cells_by_genes(tmp_path, nac):
    adata = scd.simulate_droplets(n_cells=50, n_empty=200, n_low=300, n_genes=50, splicing=True, seed=4)
    write_kb(adata, tmp_path / "out" / "counts_unfiltered", nac=nac)
    read = scd.read(tmp_path / "out")
    _same_counts(read, adata)
    assert list(read.var_names) == list(adata.var_names)
    if nac:
        assert (read.layers["unspliced"] != adata.layers["unspliced"]).nnz == 0


@pytest.mark.parametrize("usa", [False, True])
def test_alevin_fry_quants_are_read_and_usa_blocks_summed(tmp_path, usa):
    adata = scd.simulate_droplets(n_cells=50, n_empty=200, n_low=300, n_genes=50, splicing=True, seed=5)
    write_alevin_fry(adata, tmp_path / "af_quant" / "alevin", usa=usa)
    read = scd.read(tmp_path)
    _same_counts(read, adata)
    if usa:
        # anndata 0.13 lists X as the layer None
        assert {k for k in read.layers.keys() if k is not None} == {"spliced", "unspliced"}


def _split_ambiguous(adata):
    """Move about half of each spliced count to an `ambiguous` layer, as kb nac and alevin-fry USA
    split the reads within exons from those across a junction; X stays the same."""
    spliced = sp.csr_matrix(adata.layers["spliced"])
    ambiguous = spliced.copy()
    ambiguous.data = ambiguous.data // 2
    ambiguous.eliminate_zeros()
    adata.layers["ambiguous"] = ambiguous
    adata.layers["spliced"] = (spliced - ambiguous).tocsr()
    return spliced


@pytest.mark.parametrize("layout", ["kb", "alevin-fry"])
def test_kb_nac_and_alevin_fry_count_ambiguous_reads_as_spliced(tmp_path, layout):
    """kb nac's ambiguous counts and alevin-fry's A block (reads within exons, mostly) used to be
    left out of the unspliced fraction: on 10x's PBMC 1k, kb's counts put the cells at 0.63 and the
    barcodes left out at 0.50; with those reads spliced, as velocyto counts them, 0.32 and 0.14."""
    adata = scd.simulate_droplets(n_cells=50, n_empty=200, n_low=300, n_genes=50, splicing=True, seed=6)
    whole_spliced = _split_ambiguous(adata)
    if layout == "kb":
        read = scd.read(write_kb(adata, tmp_path / "counts_unfiltered", nac=True))
    else:
        read = scd.read(write_alevin_fry(adata, tmp_path / "alevin", usa=True))
    _same_counts(read, adata)
    assert {k for k in read.layers.keys() if k is not None} == {"spliced", "unspliced"}
    assert (read.layers["spliced"] != whole_spliced).nnz == 0
    unspliced = np.asarray(adata.layers["unspliced"].sum(axis=1)).ravel()
    whole = unspliced + np.asarray(whole_spliced.sum(axis=1)).ravel()
    expected = np.divide(unspliced, whole, out=np.full(whole.shape, np.nan), where=whole > 0)
    assert np.allclose(scd.unspliced_fraction(read).to_numpy(), expected, equal_nan=True)


def test_a_cell_ranger_outs_folder_gives_its_raw_matrix(sim, tmp_path):
    outs = tmp_path / "run" / "outs"
    outs.mkdir(parents=True)
    write_10x_h5_v3(sim, outs / "raw_feature_bc_matrix.h5")
    write_10x_h5_v3(sim[:10].copy(), outs / "filtered_feature_bc_matrix.h5")
    assert find_matrix(tmp_path / "run") == outs / "raw_feature_bc_matrix.h5"
    assert scd.read(tmp_path / "run").n_obs == sim.n_obs


def test_a_folder_with_only_a_filtered_matrix_says_so(sim, tmp_path):
    write_10x_h5_v3(sim[:10].copy(), tmp_path / "filtered_feature_bc_matrix.h5")
    with pytest.raises(FileNotFoundError, match="only a filtered matrix"):
        scd.read(tmp_path)


def test_h5ad_is_read_with_anndata(sim, tmp_path):
    sim.write_h5ad(tmp_path / "a.h5ad")
    _same_counts(scd.read(tmp_path / "a.h5ad"), sim)


def test_zarr_is_read_with_anndata(sim, tmp_path):
    # anndata installs zarr only from 0.12; at the 0.11 floor it is optional.
    pytest.importorskip("zarr")
    sim.write_zarr(tmp_path / "a.zarr")
    _same_counts(scd.read(tmp_path / "a.zarr"), sim)


def test_an_hdf5_file_that_is_not_10x_is_refused(tmp_path):
    import h5py

    with h5py.File(tmp_path / "other.h5", "w") as f:
        f.create_dataset("x", data=[1, 2])
    with pytest.raises(ValueError, match="not a 10x"):
        scd.read(tmp_path / "other.h5")


def test_unknown_inputs_are_refused_with_what_was_expected(tmp_path):
    (tmp_path / "notes.txt").write_text("hello")
    with pytest.raises(ValueError, match="cannot tell"):
        scd.read(tmp_path / "notes.txt")
    with pytest.raises(FileNotFoundError, match="does not exist"):
        scd.read(tmp_path / "missing")
    with pytest.raises(ValueError, match="solo_feature"):
        scd.read(tmp_path, solo_feature="Genes")


def test_barcodes_are_written_one_per_line_gzipped_or_not(tmp_path):
    write_barcodes(["A-1", "C-1"], tmp_path / "cells.tsv")
    assert (tmp_path / "cells.tsv").read_text() == "A-1\nC-1\n"
    write_barcodes(["A-1"], tmp_path / "cells.tsv.gz")
    assert gzip.open(tmp_path / "cells.tsv.gz").read() == b"A-1\n"


def test_tables_are_csv_or_tsv_by_name(tmp_path):
    frame = pd.DataFrame({"total": [3, 1]}, index=pd.Index(["A", "B"], name="barcode"))
    write_table(frame, tmp_path / "t.csv")
    write_table(frame, tmp_path / "t.tsv.gz")
    assert (tmp_path / "t.csv").read_text().splitlines()[0] == "barcode,total"
    assert gzip.open(tmp_path / "t.tsv.gz").read().decode().splitlines()[0] == "barcode\ttotal"


def test_a_failed_write_leaves_no_partial_file(tmp_path):
    class Broken:
        def write_h5ad(self, path):
            path.write_text("half")
            raise OSError("disk full")

    from sc_droplets.io.write import write_h5ad

    with pytest.raises(OSError):
        write_h5ad(Broken(), tmp_path / "out.h5ad")
    assert list(tmp_path.iterdir()) == []


def test_outputs_may_not_be_the_input_or_each_other(tmp_path):
    source = tmp_path / "in.h5ad"
    source.write_text("x")
    with pytest.raises(ValueError, match="also the input"):
        check_outputs([source], [source])
    with pytest.raises(ValueError, match="two outputs"):
        check_outputs([tmp_path / "a.csv", tmp_path / "a.csv"], [source])
    check_outputs([None, tmp_path / "a.csv"], [source])
