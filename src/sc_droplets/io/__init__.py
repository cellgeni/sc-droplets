"""Reading any count matrix a droplet pipeline writes, and writing results.

`read(path)` takes a file or a folder and finds the raw (unfiltered) matrix in
it. See `describe_inputs` for every layout it recognises.
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional

from sc_droplets.io.mtx import detect_mtx, read_mtx_folder
from sc_droplets.io.tenx import is_10x_h5, read_10x_h5
from sc_droplets.io.write import write_barcodes, write_h5ad, write_table

__all__ = [
    "read",
    "find_matrix",
    "write_barcodes",
    "write_h5ad",
    "write_table",
    "SOLO_FEATURES",
]

#: STARsolo feature folders, in the order `--solo-feature` documents them.
SOLO_FEATURES = ("Gene", "GeneFull", "GeneFull_Ex50pAS", "GeneFull_ExonOverIntron", "Velocyto")

#: Where a raw matrix sits inside a pipeline's output folder, most specific first.
_RAW_CANDIDATES = (
    "raw_feature_bc_matrix.h5",  # Cell Ranger 3+ outs/
    "raw_gene_bc_matrices_h5.h5",  # Cell Ranger 2 outs/
    "raw_feature_bc_matrix",  # Cell Ranger 3+ MTX
    "outs/raw_feature_bc_matrix.h5",
    "outs/raw_gene_bc_matrices_h5.h5",
    "outs/multi/count/raw_feature_bc_matrix.h5",  # Cell Ranger multi
    "multi/count/raw_feature_bc_matrix.h5",
    "count/raw_feature_bc_matrix.h5",
    "counts_unfiltered",  # kb count
    "alevin",  # alevin-fry quant folder
    "af_quant/alevin",  # simpleaf
)


def find_matrix(path: Path, solo_feature: str = "Gene") -> Path:
    """The file or folder holding the raw matrix inside `path`.

    `path` itself when it already is one. Raises with what was looked for.
    """
    path = Path(path)
    if path.is_file() or path.suffix == ".zarr":
        return path
    if not path.is_dir():
        raise FileNotFoundError(f"{path} does not exist")
    if detect_mtx(path) is not None:
        return path
    # STARsolo: Solo.out/<feature>/raw, from Solo.out, its parent or the feature folder.
    for base in (path, path / "Solo.out"):
        feature = base / solo_feature
        if (feature / "raw").is_dir() and detect_mtx(feature / "raw") is not None:
            return feature / "raw"
    if path.name == solo_feature and (path / "raw").is_dir():
        return path / "raw"
    for relative in _RAW_CANDIDATES:
        candidate = path / relative
        if candidate.is_file() or (candidate.is_dir() and detect_mtx(candidate) is not None):
            return candidate
    filtered = [p for p in ("filtered_feature_bc_matrix.h5", "outs/filtered_feature_bc_matrix.h5") if (path / p).exists()]
    hint = f" (only a filtered matrix, {filtered[0]}, which has no empty droplets)" if filtered else ""
    raise FileNotFoundError(
        f"no raw count matrix found in {path}{hint}. Give a .h5ad, a 10x .h5, an MTX folder, "
        "a Cell Ranger outs/, a STARsolo Solo.out/, a kb counts_unfiltered/ or an alevin-fry quant folder"
    )


def read(path, *, solo_feature: str = "Gene", genome: Optional[str] = None):
    """Read a raw count matrix into an AnnData (barcodes × features).

    `path` may be an `.h5ad`, a `.zarr` store, a 10x `.h5` (v2 or v3, also
    CellBender's output), any MTX folder `detect_mtx` knows, or a pipeline
    output folder holding one (Cell Ranger `outs/`, STARsolo `Solo.out/`, kb,
    alevin-fry). `solo_feature` picks the STARsolo feature folder.
    """
    import anndata as ad

    if solo_feature not in SOLO_FEATURES:
        raise ValueError(f"solo_feature must be one of {', '.join(SOLO_FEATURES)}")
    target = find_matrix(Path(path), solo_feature=solo_feature)
    suffix = target.suffix.lower()
    if target.is_dir() and suffix != ".zarr":
        adata = read_mtx_folder(target)
    elif suffix == ".h5ad":
        adata = ad.read_h5ad(target)
    elif suffix == ".zarr":
        adata = ad.read_zarr(target)
    elif suffix in (".h5", ".hdf5"):
        if not is_10x_h5(target):
            raise ValueError(f"{target} is HDF5 but not a 10x feature-barcode matrix")
        return read_10x_h5(target, genome=genome)
    else:
        raise ValueError(f"cannot tell what {target} is from its name (.h5ad, .zarr, .h5 or a folder)")
    return adata


def describe_inputs() -> List[str]:
    """One line per input this module reads, for help text and docs."""
    return [
        ".h5ad or .zarr (AnnData, barcodes × features, raw counts in X or a layer)",
        "10x .h5, v2 (Cell Ranger 2) or v3 (Cell Ranger 3+, CellBender)",
        "MTX folder: Cell Ranger 2 or 3+, STARsolo raw/ or Velocyto/raw/, kb, alevin-fry",
        "Cell Ranger outs/ (count or multi), STARsolo Solo.out/, kb output, simpleaf output",
    ]
