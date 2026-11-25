"""Matrix Market folders, as each pipeline writes them.

| Layout | Files | Orientation |
|---|---|---|
| Cell Ranger 3+ | `matrix.mtx.gz`, `barcodes.tsv.gz`, `features.tsv.gz` | features × barcodes |
| Cell Ranger 2 | `matrix.mtx`, `barcodes.tsv`, `genes.tsv` | features × barcodes |
| STARsolo | `matrix.mtx`, `barcodes.tsv`, `features.tsv` (uncompressed); `Velocyto/raw/` has `spliced.mtx`, `unspliced.mtx`, `ambiguous.mtx` | features × barcodes |
| kb (kallisto bustools) | `cells_x_genes.mtx`, `cells_x_genes.barcodes.txt`, `cells_x_genes.genes.txt`; `--workflow nac` adds `.mature`, `.nascent`, `.ambiguous` matrices | barcodes × features |
| alevin-fry | `quants_mat.mtx`, `quants_mat_rows.txt`, `quants_mat_cols.txt`; USA mode has three column blocks (S, U, A) | barcodes × features |

The counts used for calling are the matrix the folder's own tool calls its
count matrix: `matrix.mtx` for 10x and STARsolo; for kb the plain
`cells_x_genes.mtx`, or mature + nascent + ambiguous for nac; for alevin-fry
USA the sum of the three blocks. Spliced and unspliced counts, where a layout
has them, go to `layers["spliced"]` and `layers["unspliced"]`.
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp


def _open(path: Path):
    return gzip.open(path, "rb") if path.suffix == ".gz" else open(path, "rb")


def _first(folder: Path, *names: str) -> Optional[Path]:
    for name in names:
        for candidate in (folder / name, folder / (name + ".gz")):
            if candidate.exists():
                return candidate
    return None


def _read_mtx(path: Path) -> sp.csr_matrix:
    with _open(path) as handle:
        return sp.csr_matrix(scipy.io.mmread(handle))


def _lines(path: Path) -> List[str]:
    with _open(path) as handle:
        return [line.decode("utf-8").rstrip("\r\n") for line in handle if line.strip()]


def _table(path: Path) -> pd.DataFrame:
    rows = [line.split("\t") for line in _lines(path)]
    width = max(len(r) for r in rows)
    return pd.DataFrame([r + [""] * (width - len(r)) for r in rows])


@dataclass
class MtxLayout:
    name: str
    matrix: Path
    barcodes: Path
    features: Path
    #: True when rows are barcodes (kb, alevin-fry), False for features × barcodes.
    barcodes_are_rows: bool
    extra: Dict[str, Path]


def detect_mtx(folder: Path) -> Optional[MtxLayout]:
    """The layout of a Matrix Market folder, or None if it is not one."""
    folder = Path(folder)
    if not folder.is_dir():
        return None
    matrix = _first(folder, "matrix.mtx")
    if matrix is not None:
        barcodes = _first(folder, "barcodes.tsv")
        features = _first(folder, "features.tsv", "genes.tsv")
        if barcodes is None or features is None:
            return None
        name = "10x-mtx-v2" if features.name.startswith("genes") else "10x-mtx"
        return MtxLayout(name, matrix, barcodes, features, False, {})
    spliced = _first(folder, "spliced.mtx")
    if spliced is not None and _first(folder, "barcodes.tsv") is not None:
        # STARsolo Velocyto/raw: no matrix.mtx, so the total is S + U + A.
        extra = {k: p for k in ("spliced", "unspliced", "ambiguous") if (p := _first(folder, k + ".mtx"))}
        return MtxLayout(
            "starsolo-velocyto", spliced, _first(folder, "barcodes.tsv"), _first(folder, "features.tsv"), False, extra
        )
    kb = _first(folder, "cells_x_genes.mtx")
    kb_nac = _first(folder, "cells_x_genes.mature.mtx")
    if kb is not None or kb_nac is not None:
        barcodes = _first(folder, "cells_x_genes.barcodes.txt")
        features = _first(folder, "cells_x_genes.genes.txt")
        if barcodes is None or features is None:
            return None
        extra = {
            k: p
            for k in ("mature", "nascent", "ambiguous")
            if (p := _first(folder, f"cells_x_genes.{k}.mtx"))
        }
        names = _first(folder, "cells_x_genes.genes.names.txt")
        if names is not None:
            extra["names"] = names
        return MtxLayout("kb", kb or kb_nac, barcodes, features, True, extra)
    fry = _first(folder, "quants_mat.mtx")
    if fry is not None:
        barcodes = _first(folder, "quants_mat_rows.txt")
        features = _first(folder, "quants_mat_cols.txt")
        if barcodes is None or features is None:
            return None
        return MtxLayout("alevin-fry", fry, barcodes, features, True, {})
    return None


def read_mtx_folder(folder: Path):
    """Read any layout `detect_mtx` knows into an AnnData (barcodes × features)."""
    import anndata as ad

    folder = Path(folder)
    layout = detect_mtx(folder)
    if layout is None:
        raise ValueError(f"{folder} is not a Matrix Market folder this tool can read")

    barcodes = np.array(_lines(layout.barcodes), dtype=object)
    barcodes = np.array([b.split("\t")[0] for b in barcodes], dtype=object)
    layers: Dict[str, sp.csr_matrix] = {}

    def oriented(path: Path) -> sp.csr_matrix:
        m = _read_mtx(path)
        return m if layout.barcodes_are_rows else m.T.tocsr()

    if layout.name == "starsolo-velocyto":
        for key in ("spliced", "unspliced", "ambiguous"):
            if key in layout.extra:
                layers[key] = oriented(layout.extra[key])
        matrix = sum(layers.values(), sp.csr_matrix(layers["spliced"].shape, dtype=layers["spliced"].dtype))
    elif layout.name == "kb" and {"mature", "nascent"} <= set(layout.extra):
        mature = oriented(layout.extra["mature"])
        nascent = oriented(layout.extra["nascent"])
        ambiguous = oriented(layout.extra["ambiguous"]) if "ambiguous" in layout.extra else None
        layers["spliced"], layers["unspliced"] = mature, nascent
        matrix = mature + nascent
        if ambiguous is not None:
            layers["ambiguous"] = ambiguous
            matrix = matrix + ambiguous
    else:
        matrix = oriented(layout.matrix)

    if layout.name in ("kb", "alevin-fry"):
        ids = np.array(_lines(layout.features), dtype=object)
        if layout.name == "alevin-fry" and matrix.shape[1] == 3 * len(ids):
            n = len(ids)
            blocks = [matrix[:, i * n : (i + 1) * n] for i in range(3)]
            layers["spliced"], layers["unspliced"], layers["ambiguous"] = blocks
            matrix = (blocks[0] + blocks[1] + blocks[2]).tocsr()
        names = np.array(_lines(layout.extra["names"]), dtype=object) if "names" in layout.extra else ids
        var = pd.DataFrame({"gene_ids": ids, "feature_types": "Gene Expression"})
    else:
        table = _table(layout.features)
        ids = table[0].to_numpy(dtype=object)
        names = table[1].to_numpy(dtype=object) if table.shape[1] > 1 else ids
        types = table[2].to_numpy(dtype=object) if table.shape[1] > 2 else np.array(["Gene Expression"] * len(ids), dtype=object)
        var = pd.DataFrame({"gene_ids": ids, "feature_types": types})
        if table.shape[1] > 3:
            var["genome"] = table[3].to_numpy(dtype=object)

    if matrix.shape != (len(barcodes), len(var)):
        raise ValueError(
            f"{folder}: the matrix is {matrix.shape[0]} × {matrix.shape[1]} but there are "
            f"{len(barcodes)} barcodes and {len(var)} features"
        )
    var.index = pd.Index(np.asarray(names).astype(str))
    adata = ad.AnnData(
        X=matrix.tocsr(),
        obs=pd.DataFrame(index=pd.Index(barcodes.astype(str))),
        var=var,
        layers={k: v.tocsr() for k, v in layers.items()},
    )
    adata.var_names_make_unique()
    adata.uns["sc_droplets_input"] = {"path": str(folder), "format": layout.name}
    return adata
