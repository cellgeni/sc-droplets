"""Draw the figures of the *Visualising calls* guide, with the package's own plots.

    uv run python docs/_scripts/figures.py

One sample, 10x's pbmc_1k_v3 raw matrix, with five calls and the unspliced fraction from kb's nac
counts (`oracles/splicing/run.sh`; without them the splicing figures are left out). Writes
docs/assets/guide/*.png. Needs the dataset of tests/datasets.py (downloaded on first use).
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tests"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import sc_droplets as scd  # noqa: E402
from comparison import splicing_fractions  # noqa: E402
from datasets import raw_path  # noqa: E402

OUT = ROOT / "docs" / "assets" / "guide"
NAME = "pbmc_1k_v3_3_0"
CALLS = {
    "cellranger": "cellranger@3.0.0:expect_cells=1000",
    "emptydrops": "dropletutils-emptydrops",
    "otsu": "otsu",
    "gmm": "gmm",
    "zumis": "zumis-cell-barcodes",
}


def save(figure, name: str) -> None:
    figure.savefig(OUT / f"{name}.png", dpi=90, bbox_inches="tight")
    plt.close(figure)
    print(f"  {name}.png")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    adata = scd.read(raw_path(NAME))
    splicing = splicing_fractions(NAME)
    measured = scd.add_metrics(adata, splicing=splicing)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for key, spec in CALLS.items():
            scd.call(adata, spec, key_added=key)
    keys = list(CALLS)
    save(scd.pl.knee(adata, keys).figure, "knee")
    save(scd.pl.knee(adata, "cellranger", rotated=True).figure, "knee-rotated")
    save(scd.pl.density(adata, keys).figure, "density")
    save(scd.pl.compare(adata, keys), "compare")
    save(scd.pl.qc(adata, keys), "qc")
    save(scd.pl.diagnostics(adata, "emptydrops"), "diagnostics-emptydrops")
    save(scd.pl.diagnostics(adata, "gmm"), "diagnostics-gmm")
    if "unspliced_fraction" in measured:
        tracks = [m for m in ("unspliced_fraction", "mito_fraction") if m in measured]
        save(scd.pl.knee(adata, "cellranger", color="unspliced_fraction", track=tracks).figure, "knee-tracks")
        save(scd.pl.splicing(adata, keys), "splicing")


if __name__ == "__main__":
    main()
