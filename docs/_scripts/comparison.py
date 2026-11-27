"""Run every method on public raw matrices and write the *Comparison on public data* page.

    uv run python docs/_scripts/comparison.py [DATASET ...]

For each dataset (default: all of DATASETS) every method in SPECS runs with its defaults on the raw
matrix, as `sc-droplets compare` would. Each call is set against 10x's own filtered matrix (the
call Cell Ranger made when 10x processed the sample), and scored by the mitochondrial fraction and,
for pbmc_1k_v3, the unspliced fraction from kb's nac counts (`oracles/splicing/run.sh`).

Writes docs/overview/comparison/<dataset>.csv, docs/assets/comparison/<dataset>.png (pl.compare
of a representative subset), and the tables between the markers of docs/overview/comparison.md.
Needs the datasets of tests/datasets.py (downloaded on first use).
"""

from __future__ import annotations

import gzip
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tests"))

import sc_droplets as scd  # noqa: E402
from datasets import DATASETS as KNOWN, filtered_barcodes, raw_path  # noqa: E402

PAGE = ROOT / "docs" / "overview" / "comparison.md"
TABLES = ROOT / "docs" / "overview" / "comparison"
FIGURES = ROOT / "docs" / "assets" / "comparison"
SPLICING = ROOT / "tests" / "oracles" / "splicing"
START, END = "<!-- comparison:start -->", "<!-- comparison:end -->"

#: dataset: (label, genome of a Cell Ranger 2 file, the Cell Ranger spec 10x's run matches)
DATASETS = {
    "pbmc4k_2_1": ("PBMC 4k, v2 chemistry, Cell Ranger 2.1", "GRCh38", "cellranger@2.1.0:expect_cells=5000"),
    "pbmc_1k_v3_3_0": ("PBMC 1k, v3 chemistry, Cell Ranger 3.0", None, "cellranger@3.0.0:expect_cells=1000"),
    "pbmc5k_gemx_9_0": ("PBMC 5k, GEM-X, Cell Ranger 9.0", None, "cellranger@9.0.0"),
    "spinal_cord_10_0": ("mouse spinal cord nuclei, Cell Ranger 10.0", None, "cellranger@10.0.0"),
    "flex_lymphnode1_7_2": ("mouse lymph node, Flex (1 of a 4-plex), Cell Ranger 7.2", None, "cellranger@7.2.0"),
}

#: Every method with its defaults (the three that need a number are left out: top-n, threshold,
#: starsolo-topcells), by tool page.
SPECS = [
    "cellranger", "cellranger-ordmag",
    "starsolo-cellranger2.2", "starsolo-emptydrops-cr", "starsolo-basic-awk",
    "dropletutils-barcoderanks", "dropletutils-emptydrops", "dropletutils-emptydropscellranger",
    "dropletutils-defaultdrops",
    "umitools-distance", "umitools-density", "kneedle", "cellsweep-empty", "cellbender-prior",
    "alevinfry-permit-list", "bustools-allowlist", "alevin-initial-whitelist",
    "dropseq-dropbead", "dropest-bounds", "zumis-cell-barcodes",
    "mad", "otsu", "multiotsu", "kde-valley", "gmm", "gmm-2d",
]
#: The calls drawn in each dataset's figure: one per family.
FIGURE_SPECS = ["dropletutils-emptydrops", "starsolo-emptydrops-cr", "cellbender-prior", "umitools-distance",
                "alevinfry-permit-list", "kneedle", "zumis-cell-barcodes", "otsu", "gmm"]


def splicing_fractions(name: str):
    """nascent / (mature + ambiguous + nascent) per barcode from kb's counts, or None when not made
    yet. Ambiguous reads (mostly within exons) count as spliced, as `scd.read` counts them."""
    found = sorted(SPLICING.glob(f"{name}/kb-*/counts.tsv.gz"))
    if not found:
        return None
    with gzip.open(found[-1], "rt") as handle:
        table = pd.read_csv(handle, sep="\t", header=None, names=["barcode", "mature", "nascent", "ambiguous"],
                            index_col=0)
    whole = (table["mature"] + table["ambiguous"] + table["nascent"]).to_numpy(dtype=float)
    fraction = np.divide(table["nascent"].to_numpy(dtype=float), whole, out=np.full(whole.shape, np.nan), where=whole > 0)
    return pd.Series(fraction, index=table.index)


def key_of(spec: str) -> str:
    return spec.split(":")[0].replace("@", "_").replace("-", "_").replace(".", "_")


def run(name: str) -> pd.DataFrame:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    label, genome, reference_spec = DATASETS[name]
    adata = scd.read(raw_path(name), genome=genome)
    splicing = splicing_fractions(name)
    measured = scd.add_metrics(adata, splicing=splicing)
    reference = set(filtered_barcodes(name, genome))
    barcodes = np.asarray(adata.obs_names)
    rows = []
    for spec in [reference_spec] + SPECS:
        key = key_of(spec)
        started = time.time()
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                result = scd.call(adata, spec, key_added=key)
        except Exception as error:  # a method that cannot run here says why on the page
            rows.append({"spec": spec, "key": key, "note": str(error).split("\n")[0][:160]})
            print(f"  {spec}: {error}")
            continue
        called = set(barcodes[result.is_cell])
        union = len(called | reference)
        rows.append({
            "spec": spec, "key": key, "cells": result.n_cells,
            "jaccard_10x": len(called & reference) / union if union else 1.0,
            "smallest": result.min_cell_total, "seconds": round(time.time() - started, 1),
            "note": "; ".join(str(n) for n in result.notes)[:160],
        })
        print(f"  {spec}: {result.n_cells:,} cells, Jaccard {rows[-1]['jaccard_10x']:.3f}, {rows[-1]['seconds']}s")
    frame = pd.DataFrame(rows).set_index("key")
    stored = [k for k in frame.index if k in adata.obs.columns]
    for metric in measured:
        frame = frame.join(scd.separation(adata, stored, metric=metric))
    frame.insert(0, "dataset", name)
    TABLES.mkdir(parents=True, exist_ok=True)
    frame.to_csv(TABLES / f"{name}.csv")
    FIGURES.mkdir(parents=True, exist_ok=True)
    shown = [key_of(reference_spec)] + [key_of(s) for s in FIGURE_SPECS if key_of(s) in stored]
    figure = scd.pl.compare(adata, shown)
    figure.suptitle(f"{label}: 10x called {len(reference):,} cells", fontsize="medium")
    figure.savefig(FIGURES / f"{name}.png", dpi=90)
    plt.close(figure)
    return frame


def _number(value, spec: str) -> str:
    return "—" if pd.isna(value) else format(value, spec)


def _cells(row) -> str:
    if pd.isna(row.get("cells")):
        return "fails"
    return f"{int(row['cells']):,} ({row['jaccard_10x']:.2f})"


def write_page(frames) -> None:
    lines = [START, "", "Cells called, and in brackets the Jaccard index of the cell set with 10x's own call:", ""]
    names = list(frames)
    lines.append("| Method | " + " | ".join(DATASETS[n][0] for n in names) + " |")
    lines.append("|---|" + "---|" * len(names))
    lines.append("| 10x's filtered matrix | " + " | ".join(f"{len(set(filtered_barcodes(n, DATASETS[n][1]))):,}" for n in names) + " |")
    specs = [None] + SPECS
    for spec in specs:
        label = "Cell Ranger, the release 10x ran" if spec is None else f"`{spec}`"
        cells = []
        for n in names:
            frame = frames[n]
            key = key_of(DATASETS[n][2]) if spec is None else key_of(spec)
            cells.append(_cells(frame.loc[key]) if key in frame.index else "—")
        lines.append(f"| {label} | " + " | ".join(cells) + " |")
    lines += ["", "Where a method fails on a dataset, its row on that dataset's table says why.", ""]
    for n in names:
        frame = frames[n]
        lines += [f"### {DATASETS[n][0]}", "", f"![Calls on {n}]({{{{ '/assets/comparison/{n}.png' | relative_url }}}})", ""]
        scored = [m for m in ("unspliced_fraction", "mito_fraction") if f"{m}_auroc" in frame.columns]
        header = ["Method", "Cells", "Jaccard with 10x", "Smallest cell total"]
        for m in scored:
            short = "unspliced" if m.startswith("unspliced") else "mito"
            header += [f"{short}: cells", f"{short}: others", f"{short}: AUROC"]
        header.append("Seconds")
        lines += ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
        for key, row in frame.iterrows():
            if pd.isna(row.get("cells")):
                lines.append(f"| `{row['spec']}` | fails: {row['note']} |" + " |" * (len(header) - 2))
                continue
            values = [f"`{row['spec']}`", f"{int(row['cells']):,}", f"{row['jaccard_10x']:.3f}", f"{row['smallest']:g}"]
            for m in scored:  # a call that leaves out no barcode above the floor has no outsiders: a dash
                values += [_number(row[f"{m}_cells"], ".3f"), _number(row[f"{m}_others"], ".3f"),
                           _number(row[f"{m}_auroc"], ".2f")]
            values.append(f"{row['seconds']:g}")
            lines.append("| " + " | ".join(values) + " |")
        lines.append("")
    lines.append(END)
    text = PAGE.read_text()
    head, tail = text.split(START)[0], text.split(END)[1]
    PAGE.write_text(head + "\n".join(lines) + tail)


def main(argv) -> None:
    names = argv or list(DATASETS)
    unknown = [n for n in names if n not in DATASETS or n not in KNOWN]
    if unknown:
        sys.exit(f"unknown dataset {', '.join(unknown)}")
    frames = {}
    for name in DATASETS:
        if name in names:
            print(name)
            frames[name] = run(name)
        elif (TABLES / f"{name}.csv").exists():
            frames[name] = pd.read_csv(TABLES / f"{name}.csv", index_col="key")
    write_page(frames)


if __name__ == "__main__":
    main(sys.argv[1:])
