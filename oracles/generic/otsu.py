"""Call scikit-image's threshold_otsu and threshold_multiotsu on a dataset's counts, as the `otsu`
and `multiotsu` rules do.

    python otsu.py COUNTS.tsv OUT_DIR

COUNTS.tsv is `barcode<TAB>count` (oracles/make_counts.py). The values are the counts at or above
each floor, through the C library's log10 or asinh (Python's math module), as sc-droplets
transforms them. Writes OUT_DIR/otsu.tsv (`floor<TAB>transform<TAB>threshold`),
OUT_DIR/multiotsu.tsv (`floor<TAB>transform<TAB>classes<TAB>thresholds`, comma-separated), each
threshold as Python's repr, and OUT_DIR/versions.txt.
"""

import math
import platform
import sys
from pathlib import Path

import numpy as np
import skimage
from skimage.filters import threshold_multiotsu, threshold_otsu

TRANSFORMS = {"log10": math.log10, "arcsinh": math.asinh}

counts = np.array([float(line.split("\t")[1]) for line in Path(sys.argv[1]).read_text().splitlines()])
out = Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)
otsu, multi = [], []
for floor in (1, 10, 100, 200):
    kept = counts[counts >= floor]
    for name, function in TRANSFORMS.items():
        x = np.array([function(v) for v in kept.tolist()], dtype=np.float64)
        otsu.append(f"{floor}\t{name}\t{float(threshold_otsu(x))!r}")
        for classes in (3, 4):
            cuts = threshold_multiotsu(x, classes=classes)
            multi.append(f"{floor}\t{name}\t{classes}\t{','.join(repr(float(c)) for c in cuts)}")
        print(otsu[-1], "|", multi[-2].split("\t")[-1])
(out / "otsu.tsv").write_text("\n".join(otsu) + "\n")
(out / "multiotsu.tsv").write_text("\n".join(multi) + "\n")
(out / "versions.txt").write_text(
    f"scikit-image {skimage.__version__}\nnumpy {np.__version__}\npython {platform.python_version()}\n"
    f"platform {platform.system()}-{platform.machine()}\n"
)
