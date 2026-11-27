"""Sum kb's nac matrices per barcode, as the splicing tests and the comparison page read them.

    python summarise.py COUNTS_UNFILTERED OUT.tsv.gz

COUNTS_UNFILTERED is kb count's `counts_unfiltered/` (`--workflow nac`): `cells_x_genes.mature.mtx`,
`.nascent.mtx`, `.ambiguous.mtx` and `.barcodes.txt`, barcodes as rows. Writes one line per
barcode with any count, `barcode<TAB>mature<TAB>nascent<TAB>ambiguous`, gzipped, in kb's order.
"""

import gzip
import sys
from pathlib import Path

import numpy as np
import scipy.io

folder, out = Path(sys.argv[1]), Path(sys.argv[2])
barcodes = folder.joinpath("cells_x_genes.barcodes.txt").read_text().split()
sums = []
for name in ("mature", "nascent", "ambiguous"):
    matrix = scipy.io.mmread(folder / f"cells_x_genes.{name}.mtx").tocsr()
    if matrix.shape[0] != len(barcodes):
        sys.exit(f"{name}: {matrix.shape[0]} rows, {len(barcodes)} barcodes")
    sums.append(np.asarray(matrix.sum(axis=1)).ravel().astype(np.int64))
mature, nascent, ambiguous = sums
keep = np.flatnonzero(mature + nascent + ambiguous > 0)
out.parent.mkdir(parents=True, exist_ok=True)
with gzip.open(out, "wt") as handle:
    for i in keep:
        handle.write(f"{barcodes[i]}\t{mature[i]}\t{nascent[i]}\t{ambiguous[i]}\n")
print(f"{keep.size} barcodes with a count -> {out}")
