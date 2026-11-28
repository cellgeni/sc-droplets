"""Run the members of each family of the *Same idea, different tools* page on the same matrices.

    uv run python docs/_scripts/families.py

Each family's methods run with the same expected cells (3,000, the common default) where they take
it, and otherwise with their defaults, on 10x's PBMC 1k (whole cells) and mouse spinal cord nuclei
raw matrices. Writes each family's table between its markers in docs/overview/families.md. Needs
the datasets of tests/datasets.py (downloaded on first use).
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tests"))

import sc_droplets as scd  # noqa: E402
from datasets import raw_path  # noqa: E402

PAGE = ROOT / "docs" / "overview" / "families.md"
DATASETS = {"pbmc_1k_v3_3_0": "PBMC 1k", "spinal_cord_10_0": "spinal cord nuclei"}
EXPECTED = 3000

#: family: [(what the row is, the spec)]
FAMILIES = {
    "robust-maximum": [
        ("Cell Ranger 2.x – 5.x OrdMag", f"cellranger-ordmag@2.1.0:expect_cells={EXPECTED}"),
        ("Cell Ranger 7.2+ OrdMag", f"cellranger-ordmag@10.1.0:expect_cells={EXPECTED}"),
        ("STARsolo CellRanger2.2", f"starsolo-cellranger2.2:expected_cells={EXPECTED}"),
        ("DropletUtils defaultDrops", f"dropletutils-defaultdrops:expected_cells={EXPECTED}"),
        ("alevin `--expectCells`", f"alevin-initial-whitelist:expect_cells={EXPECTED}"),
        ("alevin-fry `--expect-cells`", f"alevinfry-permit-list:expect_cells={EXPECTED}"),
    ],
    "chord": [
        ("UMI-tools distance", "umitools-distance"),
        ("alevin-fry `--knee-distance`", "alevinfry-permit-list"),
        ("zUMIs (inflection's `uik`)", "zumis-cell-barcodes"),
        ("dropbead", "dropseq-dropbead"),
        ("DropletUtils barcodeRanks 1.28", "dropletutils-barcoderanks@1.28"),
        ("Kneedle with yascp's search", "kneedle"),
    ],
    "density": [
        ("UMI-tools density", "umitools-density"),
        ("alevin's knee and density correction", "alevin-initial-whitelist"),
        ("kernel-density valley, Scott's bandwidth", "kde-valley"),
        ("kernel-density valley, UMI-tools' bandwidth", "kde-valley:bandwidth=0.1"),
    ],
    "ambient-test": [
        ("Cell Ranger 3.x", f"cellranger@3.0.0:expect_cells={EXPECTED}"),
        ("STARsolo EmptyDrops_CR", f"starsolo-emptydrops-cr:expected_cells={EXPECTED}"),
        ("DropletUtils emptyDropsCellRanger", f"dropletutils-emptydropscellranger:expected_cells={EXPECTED}"),
        ("DropletUtils emptyDrops", "dropletutils-emptydrops"),
    ],
}


def run() -> dict:
    results = {}
    for name in DATASETS:
        adata = scd.read(raw_path(name))
        print(name)
        for family, members in FAMILIES.items():
            for label, spec in members:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    result = scd.call(adata, spec, inplace=False)
                cells = set(np.asarray(result.barcodes)[np.asarray(result.is_cell)])
                results[name, spec] = (result.n_cells, result.min_cell_total, cells)
                print(f"  {spec}: {result.n_cells:,}")
    return results


def table(members, results) -> str:
    header = "| Method | Spec |" + "".join(f" {DATASETS[n]}: cells | smallest | Jaccard |" for n in DATASETS)
    lines = [header, "|---|---|" + "---|---|---|" * len(DATASETS)]
    for label, spec in members:
        row = f"| {label} | `{spec}` |"
        for name in DATASETS:
            n, smallest, cells = results[name, spec]
            first = results[name, members[0][1]][2]
            union = len(cells | first)
            jaccard = len(cells & first) / union if union else 1.0
            row += f" {n:,} | {smallest:g} | {jaccard:.2f} |"
        lines.append(row)
    return "\n".join(lines)


def write_page(results) -> None:
    text = PAGE.read_text()
    for family, members in FAMILIES.items():
        start, end = f"<!-- families:{family}:start -->", f"<!-- families:{family}:end -->"
        head, rest = text.split(start)
        tail = rest.split(end)[1]
        text = head + start + "\n\n" + table(members, results) + "\n\n" + end + tail
    PAGE.write_text(text)


if __name__ == "__main__":
    write_page(run())
