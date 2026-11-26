"""Run one UMI-tools release's own whitelist functions on per-barcode counts.

    python oracle.py COUNTS.tsv OUT_DIR EXPECT_CELLS[,EXPECT_CELLS...]

Runs inside a pixi environment with that release installed (Python 3.6 or later; no newer
syntax). Writes OUT_DIR/<case>.txt, the called barcodes sorted, or the single line `none` when
UMI-tools accepts no threshold, and OUT_DIR/versions.txt.

Cases: distance and distance_n1000 (1.0.0+); density, density_expect<E> for each E, and
density_n1000.
"""

import collections
import os
import sys

import numpy
import scipy
import umi_tools.version

VERSION = umi_tools.version.__version__
N_CELLS = 1000


def read_counts(path):
    counts = collections.Counter()
    with open(path) as handle:
        for line in handle:
            barcode, count = line.rstrip("\n").split("\t")
            counts[barcode] = int(count)
    return counts


def cases(expects):
    if VERSION.startswith("0."):
        from umi_tools import umi_methods

        def density(counts, expect_cells=False, cell_number=False):
            return umi_methods.getKneeEstimate(counts, expect_cells, cell_number)

        return (
            [("density", lambda c: density(c))]
            + [("density_expect%d" % e, lambda c, e=e: density(c, expect_cells=e)) for e in expects]
            + [("density_n1000", lambda c: density(c, cell_number=N_CELLS))]
        )
    from umi_tools import whitelist_methods

    def call(counts, knee_method, expect_cells=False, cell_number=False):
        whitelist, _ = whitelist_methods.getCellWhitelist(
            counts, knee_method, expect_cells, cell_number, error_correct_threshold=0
        )
        return whitelist

    return (
        [
            ("distance", lambda c: call(c, "distance")),
            ("distance_n1000", lambda c: call(c, "distance", cell_number=N_CELLS)),
            ("density", lambda c: call(c, "density")),
        ]
        + [("density_expect%d" % e, lambda c, e=e: call(c, "density", expect_cells=e)) for e in expects]
        + [("density_n1000", lambda c: call(c, "density", cell_number=N_CELLS))]
    )


def main(counts_path, out, expects):
    if not os.path.isdir(out):
        os.makedirs(out)
    counts = read_counts(counts_path)
    with open(os.path.join(out, "versions.txt"), "w") as handle:
        handle.write("umi_tools %s\npython %s\nnumpy %s\nscipy %s\nplatform %s\n" % (
            VERSION, sys.version.split()[0], numpy.__version__, scipy.__version__, sys.platform))
    for name, run in cases(expects):
        called = run(counts)
        with open(os.path.join(out, name + ".txt"), "w") as handle:
            if not called:
                handle.write("none\n")
            else:
                handle.write("".join(b + "\n" for b in sorted(called)))
        print("%s %s %s: %s" % (VERSION, os.path.basename(out), name, len(called) if called else "none"))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], [int(e) for e in sys.argv[3].split(",")])
