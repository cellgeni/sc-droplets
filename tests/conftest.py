"""Fixtures: simulated raw matrices, and the files each pipeline would write for them."""

from __future__ import annotations

import gzip
from pathlib import Path

import h5py
import numpy as np
import pytest
import scipy.io
import scipy.sparse as sp

import sc_droplets as scd


def pytest_addoption(parser):
    parser.addoption(
        "--shard",
        default=None,
        metavar="I/N",
        help="run every N-th selected test, starting at the I-th (1-based), so N CI jobs share a suite",
    )


@pytest.hookimpl(trylast=True)  # after -m and -k have deselected
def pytest_collection_modifyitems(config, items):
    shard = config.getoption("--shard")
    if not shard:
        return
    index, count = (int(part) for part in shard.split("/"))
    if not 1 <= index <= count:
        raise pytest.UsageError(f"--shard {shard}: I must be between 1 and N")
    # Round robin in collection order, so each shard gets a share of every slow test function.
    keep = [item for i, item in enumerate(items) if i % count == index - 1]
    config.hook.pytest_deselected(items=[item for i, item in enumerate(items) if i % count != index - 1])
    items[:] = keep


@pytest.fixture(scope="session")
def sim():
    """A small raw matrix: 300 cells, 2,000 empty droplets, 5,000 low barcodes."""
    return scd.simulate_droplets(n_cells=300, n_empty=2000, n_low=5000, n_genes=200, seed=7)


@pytest.fixture(scope="session")
def spliced():
    """`sim` with `spliced` and `unspliced` layers: cells about 25 % unspliced, ambient RNA 10 %."""
    return scd.simulate_droplets(n_cells=300, n_empty=2000, n_low=5000, n_genes=200, seed=7, splicing=True)


@pytest.fixture
def adata(sim):
    """A fresh copy of `sim` for tests that store calls on it."""
    return sim.copy()


def write_10x_h5_v3(adata, path: Path, feature_types=None, genomes=None) -> Path:
    """Write `adata.X` as Cell Ranger 3+ writes raw_feature_bc_matrix.h5."""
    matrix = sp.csr_matrix(adata.X)
    n = adata.n_vars
    types = feature_types if feature_types is not None else ["Gene Expression"] * n
    genome = genomes if genomes is not None else ["GRCh38"] * n
    with h5py.File(path, "w") as f:
        group = f.create_group("matrix")
        group.create_dataset("data", data=matrix.data.astype(np.int32))
        group.create_dataset("indices", data=matrix.indices.astype(np.int64))
        group.create_dataset("indptr", data=matrix.indptr.astype(np.int64))
        group.create_dataset("shape", data=np.array([n, adata.n_obs], dtype=np.int32))
        group.create_dataset("barcodes", data=np.array(adata.obs_names, dtype="S"))
        feats = group.create_group("features")
        feats.create_dataset("id", data=np.array(adata.var["gene_ids"], dtype="S"))
        feats.create_dataset("name", data=np.array(adata.var_names, dtype="S"))
        feats.create_dataset("feature_type", data=np.array(types, dtype="S"))
        feats.create_dataset("genome", data=np.array(genome, dtype="S"))
        feats.create_dataset("_all_tag_keys", data=np.array([b"genome"]))
    return path


def write_10x_h5_v2(adata, path: Path, genomes=("GRCh38",)) -> Path:
    """Write as Cell Ranger 2 does: one group per genome, features split between them."""
    matrix = sp.csr_matrix(adata.X)
    parts = np.array_split(np.arange(adata.n_vars), len(genomes))
    with h5py.File(path, "w") as f:
        for genome, features in zip(genomes, parts):
            sub = sp.csr_matrix(matrix[:, features])
            group = f.create_group(genome)
            group.create_dataset("data", data=sub.data.astype(np.int32))
            group.create_dataset("indices", data=sub.indices.astype(np.int64))
            group.create_dataset("indptr", data=sub.indptr.astype(np.int64))
            group.create_dataset("shape", data=np.array([len(features), adata.n_obs], dtype=np.int32))
            group.create_dataset("barcodes", data=np.array(adata.obs_names, dtype="S"))
            group.create_dataset("genes", data=np.array(adata.var["gene_ids"].iloc[features], dtype="S"))
            group.create_dataset("gene_names", data=np.array(adata.var_names[features], dtype="S"))
    return path


def _write_text(path: Path, lines, gz: bool) -> None:
    text = "".join(f"{line}\n" for line in lines).encode()
    if gz:
        with gzip.open(path, "wb") as handle:
            handle.write(text)
    else:
        path.write_bytes(text)


def _write_mtx(path: Path, matrix, gz: bool) -> None:
    if gz:
        with gzip.open(path, "wb") as handle:
            scipy.io.mmwrite(handle, sp.coo_matrix(matrix), field="integer")
    else:
        scipy.io.mmwrite(str(path), sp.coo_matrix(matrix), field="integer")


def write_10x_mtx(adata, folder: Path, version: int = 3) -> Path:
    """Cell Ranger's MTX folder: v3 gzipped with features.tsv.gz, v2 plain with genes.tsv."""
    folder.mkdir(parents=True, exist_ok=True)
    gz = version == 3
    suffix = ".gz" if gz else ""
    _write_mtx(folder / f"matrix.mtx{suffix}", sp.csr_matrix(adata.X).T, gz)
    _write_text(folder / f"barcodes.tsv{suffix}", adata.obs_names, gz)
    if version == 3:
        lines = [f"{i}\t{n}\tGene Expression" for i, n in zip(adata.var["gene_ids"], adata.var_names)]
        _write_text(folder / "features.tsv.gz", lines, True)
    else:
        lines = [f"{i}\t{n}" for i, n in zip(adata.var["gene_ids"], adata.var_names)]
        _write_text(folder / "genes.tsv", lines, False)
    return folder


def write_starsolo(adata, root: Path, feature: str = "Gene", velocyto: bool = False) -> Path:
    """A STARsolo Solo.out/ with <feature>/raw (and Velocyto/raw when asked)."""
    raw = root / "Solo.out" / feature / "raw"
    raw.mkdir(parents=True)
    _write_mtx(raw / "matrix.mtx", sp.csr_matrix(adata.X).T, False)
    _write_text(raw / "barcodes.tsv", adata.obs_names, False)
    features = [f"{i}\t{n}\tGene Expression" for i, n in zip(adata.var["gene_ids"], adata.var_names)]
    _write_text(raw / "features.tsv", features, False)
    if velocyto:
        vraw = root / "Solo.out" / "Velocyto" / "raw"
        vraw.mkdir(parents=True)
        for name in ("spliced", "unspliced", "ambiguous"):
            layer = adata.layers[name] if name in adata.layers else sp.csr_matrix(adata.shape, dtype=np.int64)
            _write_mtx(vraw / f"{name}.mtx", sp.csr_matrix(layer).T, False)
        _write_text(vraw / "barcodes.tsv", adata.obs_names, False)
        _write_text(vraw / "features.tsv", features, False)
    return root / "Solo.out"


def write_kb(adata, folder: Path, nac: bool = False) -> Path:
    """kb count's counts_unfiltered/ (cells × genes); nac writes mature/nascent/ambiguous."""
    folder.mkdir(parents=True)
    _write_text(folder / "cells_x_genes.barcodes.txt", adata.obs_names, False)
    _write_text(folder / "cells_x_genes.genes.txt", adata.var["gene_ids"], False)
    _write_text(folder / "cells_x_genes.genes.names.txt", adata.var_names, False)
    if nac:
        _write_mtx(folder / "cells_x_genes.mature.mtx", adata.layers["spliced"], False)
        _write_mtx(folder / "cells_x_genes.nascent.mtx", adata.layers["unspliced"], False)
        ambiguous = adata.layers["ambiguous"] if "ambiguous" in adata.layers else sp.csr_matrix(adata.shape, dtype=np.int64)
        _write_mtx(folder / "cells_x_genes.ambiguous.mtx", ambiguous, False)
    else:
        _write_mtx(folder / "cells_x_genes.mtx", adata.X, False)
    return folder


def write_loom(adata, path: Path, sample: str = "sample1") -> Path:
    """A velocyto loom: dense genes × cells layers, CellID as `sample:BARCODEx` (no `-1`)."""
    import h5py

    def dense(name):
        layer = adata.layers[name] if name in adata.layers else sp.csr_matrix(adata.shape, dtype=np.int64)
        return sp.csr_matrix(layer).T.toarray()

    cells = [f"{sample}:{b.split('-')[0]}x" for b in adata.obs_names]
    with h5py.File(path, "w") as handle:
        handle["matrix"] = dense("spliced")
        for name in ("spliced", "unspliced", "ambiguous"):
            handle[f"layers/{name}"] = dense(name)
        handle["row_attrs/Gene"] = np.array(adata.var_names, dtype="S")
        handle["row_attrs/Accession"] = np.array(adata.var["gene_ids"], dtype="S")
        handle["col_attrs/CellID"] = np.array(cells, dtype="S")
    return path


def write_alevin_fry(adata, folder: Path, usa: bool = False) -> Path:
    """alevin-fry's alevin/ quant folder; USA mode has S, U, A column blocks."""
    folder.mkdir(parents=True)
    if usa:
        ambiguous = adata.layers["ambiguous"] if "ambiguous" in adata.layers else sp.csr_matrix(adata.shape, dtype=np.int64)
        matrix = sp.hstack([adata.layers["spliced"], adata.layers["unspliced"], ambiguous], format="csr")
    else:
        matrix = adata.X
    _write_mtx(folder / "quants_mat.mtx", matrix, False)
    _write_text(folder / "quants_mat_rows.txt", adata.obs_names, False)
    _write_text(folder / "quants_mat_cols.txt", adata.var["gene_ids"], False)
    return folder
