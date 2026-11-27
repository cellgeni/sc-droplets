"""Drop-seq's cell number by dropbead's `estimateCellNumber`.

Page: docs/methods/dropseq/dropbead.md. Ported from dropbead's `R/Tools.R` (commit 73fd121,
version 0.3.1), distributed under the GNU GPL v3, as sc-droplets is. R's `sum` and `cumsum`
accumulate in long double; `rmath.ld_sum` and `ld_cumsum` do the same.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from sc_droplets.methods._common import result
from sc_droplets.primitives.rmath import ld_cumsum, ld_sum
from sc_droplets.registry import Method, MethodInfo, Reference, Variant, param, register

MAX_CELLS = 15000


class DropbeadError(ValueError):
    """Where dropbead gives no number."""


def estimate_cell_number(counts: np.ndarray, max_cells: int = MAX_CELLS) -> int:
    """`estimateCellNumber`: the first barcode farthest from the diagonal of the cumulative
    fraction of the top `max_cells` counts (`counts` sorted in decreasing order)."""
    if counts.size < max_cells:
        raise DropbeadError(f"dropbead needs at least max.cells = {max_cells} barcodes (it has {counts.size}); "
                            "with fewer its sums are missing and it returns no number")
    top = counts[:max_cells].astype(np.float64)
    x = np.arange(1, max_cells + 1, dtype=np.float64) / max_cells
    cs = ld_cumsum(top / ld_sum(top))
    m = (cs[-1] - cs[0]) / (max_cells - x[0])
    dists = math.sin(math.atan(m)) * np.abs(cs - x)
    return int(np.argmax(dists)) + 1  # which.max: the first maximum, 1-based


@register
class Dropbead(Method):
    info = MethodInfo(
        id="dropseq-dropbead",
        title="Drop-seq cell number (dropbead)",
        tool="Drop-seq (dropbead)",
        category="knee",
        summary="dropbead's estimateCellNumber: the barcode farthest from the diagonal of the cumulative "
                "fraction of reads of the top 15,000.",
        docs="methods/dropseq/dropbead.md",
        counts="either",
        variants=(Variant("default", since="0.3.1", summary="estimateCellNumber"),),
        latest_known="0.3.1",
        references=(
            Reference(
                "Alles, Karaiskos, Praktiknjo et al. 2017, Cell fixation and preservation for droplet-based "
                "single-cell transcriptomics, BMC Biology 15:44",
                doi="10.1186/s12915-017-0383-5",
            ),
        ),
        source="https://github.com/rajewsky-lab/dropbead/blob/73fd121/R/Tools.R",
        original_license="GPL-3.0",
    )

    @dataclass
    class Params:
        max_cells: int = param(MAX_CELLS, "The top barcodes the curve is drawn over", original="max.cells", minimum=2)
        unit: str = param("umis", "Count UMIs (the matrix) or reads (--reads)", original="reads per barcode", choices=("umis", "reads"))

    def run(self, droplets, params, variant, random_state):
        values = droplets.values(params.unit)
        order = np.argsort(-values, kind="stable")  # ties keep the matrix's order
        counts = values[order]
        cells = estimate_cell_number(counts, params.max_cells)
        is_cell = np.zeros(values.size, dtype=bool)
        is_cell[order[:cells]] = True
        return result(type(self), variant, params, droplets, is_cell,
                      thresholds={"cells": float(cells), params.unit: float(counts[cells - 1])})
