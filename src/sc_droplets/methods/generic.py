"""Generic rules that belong to no single tool: top N, a fixed threshold and MAD outliers.

Pages: docs/methods/generic/{top-n,threshold,mad}.md. MAD follows scuttle's `isOutlier` (R's
`log2`, `median` and `mad`).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional

import numpy as np

from sc_droplets.methods._common import result, top_by_rank
from sc_droplets.primitives.otsu import histogram
from sc_droplets.primitives.rmath import clog2
from sc_droplets.primitives.rstats import r_mad, r_median
from sc_droplets.registry import Method, MethodInfo, Reference, param, register

_GENERIC = "Generic"
_UNITS = ("umis", "reads", "genes")
_UNIT_LABEL = {"umis": "UMIs", "reads": "reads", "genes": "genes"}


class RuleError(ValueError):
    """Where a generic rule cannot place its cut."""


def candidates(values: np.ndarray, min_counts: float) -> np.ndarray:
    """Indices of the barcodes with at least `min_counts`; the others are empty."""
    index = np.flatnonzero(values >= min_counts)
    if index.size == 0:
        raise RuleError(f"no barcode has at least --min-counts {min_counts:g}")
    return index


def _histogram_diagnostics(x: np.ndarray, n_bins: int, label: str, cuts) -> dict:
    counts, centers = histogram(x, n_bins)
    half = (centers[1] - centers[0]) / 2 if centers.size > 1 else 0.5
    return {
        "histogram_counts": counts.astype(np.int64),
        "histogram_edges": np.concatenate([centers - half, centers[-1:] + half]),
        "histogram_label": label,
        "cuts": np.asarray([c for c, _ in cuts], dtype=np.float64),
        "cut_labels": [str(t) for _, t in cuts],
    }


@register
class TopN(Method):
    info = MethodInfo(
        id="top-n",
        title="Top N barcodes",
        tool=_GENERIC,
        category="rank-threshold",
        summary="The N barcodes with the most counts are cells.",
        docs="methods/generic/top-n.md",
        counts="either",
        references=(
            Reference(
                "The rule behind Cell Ranger --force-cells, STARsolo TopCells, UMI-tools "
                "--set-cell-number and alevin --forceCells; each breaks ties its own way (see their pages)."
            ),
        ),
    )

    @dataclass
    class Params:
        n_cells: int = param(help="How many barcodes to call", original="--force-cells", minimum=0)
        ties: str = param(
            "first",
            "first: exactly N, ties broken by input order; include: every barcode with at least the N-th total",
            choices=("first", "include"),
        )
        unit: str = param("umis", "What to rank by", choices=("umis", "reads", "genes"))

    def run(self, droplets, params, variant, random_state):
        values = droplets.values(params.unit)
        n = min(params.n_cells, values.shape[0])
        if params.ties == "first":
            is_cell = top_by_rank(values, n)
        else:
            nth = np.sort(values)[::-1][n - 1] if n > 0 else np.inf
            is_cell = values >= nth if n > 0 else np.zeros(values.shape[0], dtype=bool)
        nth_value = float(np.sort(values)[::-1][n - 1]) if n > 0 else float("nan")
        return result(
            type(self),
            variant,
            params,
            droplets,
            is_cell,
            thresholds={params.unit: nth_value},
            diagnostics={"rank": n},
        )


@register
class Threshold(Method):
    info = MethodInfo(
        id="threshold",
        title="Fixed threshold",
        tool=_GENERIC,
        category="rank-threshold",
        summary="Barcodes with at least a fixed number of UMIs, reads or genes are cells.",
        docs="methods/generic/threshold.md",
        counts="either",
        references=(
            Reference(
                "Used as Drop-seq SelectCellsByNumTranscripts MIN_TRANSCRIPTS_PER_CELL, "
                "alevin-fry --min-reads, and the 'more than 200 genes' rule of the scanpy and Seurat tutorials."
            ),
        ),
    )

    @dataclass
    class Params:
        min_counts: float = param(help="The threshold", original="MIN_TRANSCRIPTS_PER_CELL", minimum=0)
        unit: str = param("umis", "What the threshold applies to", choices=("umis", "reads", "genes"))
        strict: bool = param(False, "Require more than the threshold, not at least it")

    def run(self, droplets, params, variant, random_state):
        values = droplets.values(params.unit)
        is_cell = values > params.min_counts if params.strict else values >= params.min_counts
        return result(
            type(self),
            variant,
            params,
            droplets,
            is_cell,
            thresholds={params.unit: float(params.min_counts)},
        )


@register
class Mad(Method):
    info = MethodInfo(
        id="mad",
        title="MAD outliers",
        tool=_GENERIC,
        category="density",
        summary="Among the barcodes above a floor, those whose log total is more than k MADs below the median "
                "are empty (scuttle's isOutlier, lower tail); the rest are cells.",
        docs="methods/generic/mad.md",
        counts="either",
        references=(
            Reference(
                "McCarthy, Campbell, Lun & Wills 2017, Scater: pre-processing, quality control, normalization "
                "and visualization of single-cell RNA-seq data in R, Bioinformatics 33:1179",
                doi="10.1093/bioinformatics/btw777",
            ),
        ),
        source="https://github.com/LTLA/scuttle/blob/master/R/isOutlier.R",
        original_license="GPL-3.0",
    )

    @dataclass
    class Params:
        nmads: float = param(3.0, "How many MADs below the median a barcode may be", original="nmads", minimum=0)
        min_counts: float = param(100.0, "Candidates have at least this value; the others are empty", minimum=0)
        log: bool = param(True, "Work on log2 values", original="log")
        min_diff: Optional[float] = param(None, "The least distance from the median to the limit", original="min.diff",
                                          minimum=0)
        unit: str = param("umis", "What the rule applies to", choices=_UNITS)

    def run(self, droplets, params, variant, random_state):
        values = droplets.values(params.unit)
        index = candidates(values, params.min_counts)
        metric = clog2(values[index]) if params.log else values[index].astype(np.float64)
        median = r_median(metric)
        mad = r_mad(metric, median)
        spread = params.nmads * mad
        if params.min_diff is not None and not math.isnan(spread):
            spread = max(params.min_diff, spread)
        elif params.min_diff is not None:
            spread = params.min_diff
        lower = median - spread
        is_cell = np.zeros(values.size, dtype=bool)
        is_cell[index[~(metric < lower)]] = True
        back = (lambda v: 2.0 ** v) if params.log else float
        label = ("log2 " if params.log else "") + _UNIT_LABEL[params.unit]
        return result(
            type(self), variant, params, droplets, is_cell,
            thresholds={params.unit: back(lower), "median": back(median), "mad": float(mad)},
            diagnostics=_histogram_diagnostics(metric, 100, label, [(median, "median"), (lower, f"limit: {params.nmads:g} MADs")]),
        )
