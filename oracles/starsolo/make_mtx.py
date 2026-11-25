"""Write a dataset's raw matrix as the uncompressed MTX folder STAR's soloCellFiltering reads.

    python oracles/starsolo/make_mtx.py DATASET OUT_DIR

STAR wants features × barcodes in `matrix.mtx`, a one-column `barcodes.tsv` and
`features.tsv`, none gzipped.
"""

from __future__ import annotations

import sys
from pathlib import Path

import scipy.io
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tests"))

import sc_droplets as scd  # noqa: E402
from datasets import raw_path  # noqa: E402


def main(name: str, out: Path) -> None:
    adata = scd.read(raw_path(name))
    out.mkdir(parents=True, exist_ok=True)
    scipy.io.mmwrite(str(out / "matrix.mtx"), sp.coo_matrix(sp.csr_matrix(adata.X).T), field="integer")
    (out / "barcodes.tsv").write_text("".join(f"{b}\n" for b in adata.obs_names))
    types = adata.var["feature_types"] if "feature_types" in adata.var else ["Gene Expression"] * adata.n_vars
    (out / "features.tsv").write_text(
        "".join(f"{i}\t{n}\t{t}\n" for i, n, t in zip(adata.var["gene_ids"], adata.var_names, types))
    )
    print(f"{name}: {adata.n_obs} barcodes × {adata.n_vars} features -> {out}")


if __name__ == "__main__":
    main(sys.argv[1], Path(sys.argv[2]))
