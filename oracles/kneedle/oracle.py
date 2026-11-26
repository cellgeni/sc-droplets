"""Run one kneed release's KneeLocator on a dataset's rank curve, as yascp and others use it.

    python oracle.py COUNTS.tsv OUT_DIR

Runs inside a pixi environment with that release installed (Python 3.7 or later; no newer
syntax). COUNTS.tsv is `barcode<TAB>count` in the matrix's order (oracles/make_counts.py).
Writes OUT_DIR/knees.tsv (`case<TAB>knee rank<TAB>sensitivity`, or `none`) and OUT_DIR/versions.txt.

The curve: counts sorted in decreasing order, ties in input order; ranks 1, 2, ...; the
barcodes with count >= lower; x = rank and y = count, or their log10 (C library, `math.log10`).
yascp's search (bin/get_estimates_from_umi_counts.py, kneedle_cutoff): S = 5000, then 100 less
each time, until KneeLocator finds a knee.
"""

import math
import os
import sys

import numpy
import scipy
from kneed import KneeLocator
import kneed

# name: (lower, log axes, curve, direction, interp_method, S (None: yascp's search), online)
CASES = [
    ("yascp_interp1d", 1000, False, "convex", "decreasing", "interp1d", None, False),
    ("yascp_polynomial", 1000, False, "convex", "decreasing", "polynomial", None, False),
    ("loglog_concave", 10, True, "concave", "decreasing", "interp1d", 1.0, False),
    ("loglog_concave_online", 100, True, "concave", "decreasing", "interp1d", 1.0, True),
    ("loglog_convex", 100, True, "convex", "decreasing", "interp1d", 1.0, False),
]


def read_counts(path):
    counts = []
    with open(path) as handle:
        for line in handle:
            counts.append(int(line.rstrip("\n").split("\t")[1]))
    return counts


def curve(counts, lower, log):
    ranked = sorted(counts, key=lambda c: -c)  # stable: ties keep input order
    kept = [c for c in ranked if c >= lower]
    ranks = list(range(1, len(kept) + 1))
    if log:
        return [math.log10(r) for r in ranks], [math.log10(c) for c in kept]
    return ranks, kept


def knee(x, y, curve_, direction, interp, S, online):
    kwargs = dict(curve=curve_, direction=direction, interp_method=interp, online=online)
    if S is not None:
        return KneeLocator(numpy.array(x), numpy.array(y), S=S, **kwargs).knee, S
    sensitivity = 5000
    while True:
        found = KneeLocator(numpy.array(x), numpy.array(y), S=sensitivity, **kwargs).knee
        if found is not None:
            return found, sensitivity
        if sensitivity < 0:
            return None, sensitivity
        sensitivity -= 100


def main(counts_path, out):
    if not os.path.isdir(out):
        os.makedirs(out)
    counts = read_counts(counts_path)
    try:
        import sklearn
        sk = sklearn.__version__
    except ImportError:
        sk = "-"
    with open(os.path.join(out, "versions.txt"), "w") as handle:
        handle.write("kneed %s\npython %s\nnumpy %s\nscipy %s\nscikit-learn %s\nplatform %s\n" % (
            kneed.__version__, sys.version.split()[0], numpy.__version__, scipy.__version__, sk, sys.platform))
    with open(os.path.join(out, "knees.tsv"), "w") as handle:
        for name, lower, log, curve_, direction, interp, S, online in CASES:
            x, y = curve(counts, lower, log)
            value, used = knee(x, y, curve_, direction, interp, S, online)
            rank = "none" if value is None else str(x.index(value) + 1)
            handle.write("%s\t%s\t%r\n" % (name, rank, float(used)))
            print("%s %s: rank %s, S %r" % (kneed.__version__, name, rank, used))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
