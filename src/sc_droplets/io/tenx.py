"""10x Genomics feature-barcode matrices in HDF5 (`*_feature_bc_matrix.h5`).

Two layouts:

- v3 (Cell Ranger 3 onwards, also CellBender's output):
  `/matrix/{data,indices,indptr,shape,barcodes}` and
  `/matrix/features/{id,name,feature_type,genome,...}`;
- v2 (Cell Ranger 2): one group per genome,
  `/<genome>/{data,indices,indptr,shape,barcodes,genes,gene_names}`.

The matrix is stored features × barcodes in CSC, which is barcodes × features
in CSR, AnnData's layout, so the arrays are taken as they are. Counts stay
integers. `var` has scanpy's columns (`gene_ids`, `feature_types`, `genome`),
indexed by gene symbol made unique, so code written against
`scanpy.read_10x_h5(gex_only=False)` works on the result.

A v2 file with several genomes (a barnyard run) is read whole: the genomes'
features side by side, `var["genome"]` saying which is which, as Cell Ranger 2
counted each genome separately.
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp


def _strings(dataset: h5py.Dataset) -> np.ndarray:
    values = dataset[()]
    if values.dtype.kind in ("S", "O"):
        return np.array([v.decode("utf-8") if isinstance(v, bytes) else str(v) for v in values], dtype=object)
    return values.astype(str).astype(object)


def is_10x_h5(path: Path) -> bool:
    """Whether an HDF5 file holds a 10x matrix (either layout)."""
    try:
        with h5py.File(path, "r") as f:
            if "matrix" in f and "barcodes" in f["matrix"]:
                return True
            return any(
                isinstance(f[name], h5py.Group) and {"barcodes", "data", "indptr"} <= set(f[name].keys())
                for name in f.keys()
            )
    except OSError:
        return False


def _csr(group: h5py.Group) -> sp.csr_matrix:
    n_features, n_barcodes = (int(x) for x in group["shape"][()])
    return sp.csr_matrix(
        (group["data"][()], group["indices"][()], group["indptr"][()]),
        shape=(n_barcodes, n_features),
    )


def read_10x_h5(path: Path, genome: Optional[str] = None):
    """Read a 10x HDF5 matrix into an AnnData (barcodes × features)."""
    import anndata as ad

    path = Path(path)
    with h5py.File(path, "r") as f:
        if "matrix" in f:
            group = f["matrix"]
            matrix = _csr(group)
            barcodes = _strings(group["barcodes"])
            feats = group["features"]
            var = pd.DataFrame(
                {
                    "gene_ids": _strings(feats["id"]),
                    "feature_types": _strings(feats["feature_type"]),
                    "genome": _strings(feats["genome"]) if "genome" in feats else "",
                },
            )
            names = _strings(feats["name"])
            layout = "10x-h5-v3"
        else:
            genomes = [n for n in f.keys() if isinstance(f[n], h5py.Group) and "barcodes" in f[n]]
            if genome is not None:
                if genome not in genomes:
                    raise ValueError(f"no genome {genome!r} in {path}; it has {', '.join(genomes)}")
                genomes = [genome]
            matrices: List[sp.csr_matrix] = []
            frames: List[pd.DataFrame] = []
            barcodes = None
            for name in genomes:
                group = f[name]
                these = _strings(group["barcodes"])
                if barcodes is None:
                    barcodes = these
                elif not np.array_equal(barcodes, these):
                    raise ValueError(
                        f"the genomes of {path} list different barcodes; read one with genome="
                    )
                matrices.append(_csr(group))
                frames.append(
                    pd.DataFrame(
                        {
                            "gene_ids": _strings(group["genes"]),
                            "gene_names": _strings(group["gene_names"]),
                            "feature_types": "Gene Expression",
                            "genome": name,
                        }
                    )
                )
            matrix = sp.hstack(matrices, format="csr") if len(matrices) > 1 else matrices[0]
            var = pd.concat(frames, ignore_index=True)
            names = var.pop("gene_names").to_numpy()
            layout = "10x-h5-v2"

    var.index = pd.Index(names.astype(str), name=None)
    adata = ad.AnnData(X=matrix, obs=pd.DataFrame(index=pd.Index(barcodes.astype(str))), var=var)
    adata.var_names_make_unique()
    adata.uns["sc_droplets_input"] = {"path": str(path), "format": layout}
    return adata
