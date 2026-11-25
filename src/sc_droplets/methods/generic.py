"""Generic rules that belong to no single tool: top-N and a fixed threshold.

Pages: docs/methods/generic/top-n.md, docs/methods/generic/threshold.md.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from sc_droplets.methods._common import result, top_by_rank
from sc_droplets.registry import Method, MethodInfo, Reference, param, register

_GENERIC = "Generic"


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
