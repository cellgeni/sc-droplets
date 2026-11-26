"""Write a dataset's raw matrix as a sorted BUS file, one record per UMI, for the bustools oracle.

    python oracles/bustools/make_bus.py DATASET OUT.bus

Each barcode with UMIs gets that many records (UMI 0, 1, ..., equivalence class 0, count 1), in
bustools' order: barcodes by their 2-bit code (A=0, C=1, G=2, T=3, first base highest), which is
their order as text. Gene Expression features only, as sc-droplets counts them; the barcode's
`-1` suffix is dropped.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tests"))

import sc_droplets as scd  # noqa: E402
from datasets import raw_path  # noqa: E402
from sc_droplets.counts import Droplets  # noqa: E402

RECORD = np.dtype([("barcode", "<u8"), ("umi", "<u8"), ("ec", "<i4"), ("count", "<u4"), ("flags", "<u4"), ("pad", "<u4")])
CODE = {"A": 0, "C": 1, "G": 2, "T": 3}


def encode(sequence: str) -> int:
    value = 0
    for base in sequence:
        value = (value << 2) | CODE[base]
    return value


def main(name: str, out: Path) -> None:
    adata = scd.read(raw_path(name))
    totals = Droplets.from_anndata(adata).totals.astype(np.int64)
    keep = np.flatnonzero(totals > 0)
    seqs = np.array([str(b).split("-", 1)[0] for b in np.asarray(adata.obs_names)[keep]])
    order = np.argsort(seqs, kind="stable")
    seqs, umis = seqs[order], totals[keep][order]
    length = len(seqs[0])
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("wb") as handle:
        handle.write(b"BUS\0")
        handle.write(np.array([1, length, 12, 0], dtype="<u4").tobytes())  # version, bclen, umilen, tlen
        for start in range(0, seqs.size, 50_000):
            codes = np.array([encode(s) for s in seqs[start:start + 50_000]], dtype=np.uint64)
            n = umis[start:start + 50_000]
            records = np.zeros(int(n.sum()), dtype=RECORD)
            records["barcode"] = np.repeat(codes, n)
            first = np.repeat(np.cumsum(n) - n, n)
            records["umi"] = np.arange(records.size, dtype=np.uint64) - first.astype(np.uint64)
            records["count"] = 1
            handle.write(records.tobytes())
    print(f"{name}: {seqs.size} barcodes, {int(umis.sum())} records, barcode length {length}")


if __name__ == "__main__":
    main(sys.argv[1], Path(sys.argv[2]))
