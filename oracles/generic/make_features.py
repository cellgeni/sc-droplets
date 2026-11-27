"""Write a dataset's per-barcode UMIs and genes detected, as the mixture oracles read them.

    python oracles/generic/make_features.py DATASET OUT.tsv

One line per barcode with a count, `barcode<TAB>umis<TAB>genes`, in the matrix's order. Gene
Expression features only, as sc-droplets counts them.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tests"))

import sc_droplets as scd  # noqa: E402
from datasets import raw_path  # noqa: E402
from sc_droplets.counts import Droplets  # noqa: E402


def main(name: str, out: Path) -> None:
    adata = scd.read(raw_path(name))
    droplets = Droplets.from_anndata(adata)
    umis, genes = droplets.values("umis"), droplets.values("genes")
    keep = np.flatnonzero(umis > 0)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as handle:
        for i in keep:
            handle.write(f"{adata.obs_names[i]}\t{int(umis[i])}\t{int(genes[i])}\n")
    print(f"{name}: {keep.size} barcodes -> {out}")


if __name__ == "__main__":
    main(sys.argv[1], Path(sys.argv[2]))
