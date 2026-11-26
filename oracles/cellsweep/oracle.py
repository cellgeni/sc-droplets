"""Run one cellsweep release's own `infer_empty_droplets` on a dataset's totals.

    python oracle.py COUNTS.tsv OUT_DIR EXPECTED_CELLS

COUNTS.tsv is `barcode<TAB>count` for every barcode of the matrix, zeros included, in the
matrix's order (oracles/make_counts.py --all). cellsweep sums `adata.X` over features; it is
given here a float32 AnnData with one feature holding each barcode's total (the sum, in float32,
of whole counts below 2**24 is exact either way). Writes OUT_DIR/cutoffs.tsv
(`case<TAB>cutoff<TAB>cells`, or `error<TAB>message`) and OUT_DIR/versions.txt.

Cases: default (the release's automatic rule), mx_filter, and expected_cells, on every barcode;
mx_filter_detected, mx_filter on the barcodes with a count (as kb's matrices hold them); and
mx_filter_over100, on the barcodes with at least 100 (where the mixture separates cells from
background rather than zeros from the rest).
"""

import importlib.metadata as metadata
import os
import sys

import anndata as ad
import numpy as np
import scipy
import scipy.sparse as sp
import sklearn
from cellsweep.utils.data_utils import infer_empty_droplets

VERSION = metadata.version("cellsweep")


def read_counts(path):
    barcodes, counts = [], []
    with open(path) as handle:
        for line in handle:
            barcode, count = line.rstrip("\n").split("\t")
            barcodes.append(barcode)
            counts.append(int(count))
    return barcodes, np.asarray(counts, dtype=np.float32)


def adata_of(barcodes, counts):
    adata = ad.AnnData(sp.csr_matrix(counts.reshape(-1, 1)))
    adata.obs_names = barcodes
    adata.var_names = ["total"]
    return adata


def cases(expected):
    default = "kneed" if VERSION >= "1.0.2" else "threshold"
    # name, keyword arguments, the smallest total of the barcodes given
    return [
        ("default", dict(method=default), 0),
        ("mx_filter", dict(method="mx_filter"), 0),
        ("expected", dict(method=default, expected_cells=expected), 0),
        ("mx_filter_detected", dict(method="mx_filter"), 1),
        ("mx_filter_over100", dict(method="mx_filter"), 100),
    ]


def main(counts_path, out, expected):
    os.makedirs(out, exist_ok=True)
    barcodes, counts = read_counts(counts_path)
    try:
        kneed_version = metadata.version("kneed")
    except metadata.PackageNotFoundError:
        kneed_version = "-"
    with open(os.path.join(out, "versions.txt"), "w") as handle:
        handle.write("cellsweep %s\nkneed %s\npython %s\nnumpy %s\nscipy %s\nscikit-learn %s\nplatform %s\n" % (
            VERSION, kneed_version, sys.version.split()[0], np.__version__, scipy.__version__, sklearn.__version__, sys.platform))
    with open(os.path.join(out, "cutoffs.tsv"), "w") as handle:
        for name, kwargs, lowest in cases(expected):
            keep = counts >= lowest
            adata = adata_of([b for b, k in zip(barcodes, keep) if k], counts[keep])
            try:
                adata = infer_empty_droplets(adata, quiet=True, **kwargs)
            except ValueError as error:
                handle.write("%s\terror\t%s\n" % (name, str(error).splitlines()[0]))
                print("%s %s: error %s" % (VERSION, name, error))
                continue
            cells = int((~adata.obs["is_empty"].to_numpy(dtype=bool)).sum())
            # the cutoff: the smallest total among the cells (cells are total >= cutoff)
            cutoff = float(counts[keep][~adata.obs["is_empty"].to_numpy(dtype=bool)].min()) if cells else float("nan")
            handle.write("%s\t%r\t%d\n" % (name, cutoff, cells))
            print("%s (kneed %s) %s: cutoff %r, %d cells" % (VERSION, kneed_version, name, cutoff, cells))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]))
