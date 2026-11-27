"""velocyto's `.loom` files: spliced, unspliced and ambiguous counts per gene and barcode.

A loom is HDF5 with `/matrix` (genes × cells, dense and chunked), `/layers/<name>` of the same
shape, `/row_attrs` (velocyto writes `Gene` and `Accession`) and `/col_attrs` (`CellID`, as
`sample:BARCODEx`). velocyto's `/matrix` is the spliced layer; the total here is spliced +
unspliced + ambiguous, as STARsolo's Velocyto/raw is read.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

LAYERS = ("spliced", "unspliced", "ambiguous")
#: Genes read at a time: a loom's layers are dense.
CHUNK = 2048


def is_loom(path: Path) -> bool:
    """Whether `path` is an HDF5 file laid out as a loom (`/matrix`, `/row_attrs`, `/col_attrs`)."""
    import h5py

    try:
        with h5py.File(path, "r") as handle:
            return all(key in handle for key in ("matrix", "row_attrs", "col_attrs"))
    except OSError:
        return False


def _text(values) -> np.ndarray:
    return np.array([v.decode() if isinstance(v, bytes) else str(v) for v in values], dtype=object)


def _sparse_layer(dataset) -> sp.csr_matrix:
    """A dense genes × cells dataset as cells × genes CSR, a block of genes at a time."""
    n_genes = dataset.shape[0]
    blocks = [sp.csr_matrix(dataset[start:start + CHUNK, :]) for start in range(0, n_genes, CHUNK)]
    return sp.vstack(blocks, format="csr").T.tocsr() if blocks else sp.csr_matrix(dataset.shape[::-1])


def read_loom(path: Path):
    """A velocyto loom as an AnnData (barcodes × genes) with `spliced`, `unspliced` and, when
    present, `ambiguous` layers; X is their sum."""
    import anndata as ad
    import h5py

    path = Path(path)
    with h5py.File(path, "r") as handle:
        rows, cols = handle["row_attrs"], handle["col_attrs"]
        if "CellID" not in cols:
            raise ValueError(f"{path}: a loom without col_attrs/CellID has no barcodes")
        barcodes = _text(cols["CellID"][:])
        genes = _text(rows["Gene"][:]) if "Gene" in rows else _text(rows["Accession"][:])
        ids = _text(rows["Accession"][:]) if "Accession" in rows else genes
        present = [name for name in LAYERS if name in handle.get("layers", {})]
        if not {"spliced", "unspliced"} <= set(present):
            raise ValueError(f"{path}: a loom needs layers/spliced and layers/unspliced (velocyto's), it has {present}")
        layers = {name: _sparse_layer(handle["layers"][name]) for name in present}
    total = layers["spliced"] + layers["unspliced"]
    if "ambiguous" in layers:
        total = total + layers["ambiguous"]
    var = pd.DataFrame({"gene_ids": ids, "feature_types": "Gene Expression"}, index=pd.Index(genes.astype(str)))
    adata = ad.AnnData(X=total.tocsr(), obs=pd.DataFrame(index=pd.Index(barcodes.astype(str))), var=var, layers=layers)
    adata.var_names_make_unique()
    adata.uns["sc_droplets_input"] = {"path": str(path), "format": "loom"}
    return adata
