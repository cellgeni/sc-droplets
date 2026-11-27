"""dropEst's cell-number bounds: dropestr's `EstimateCellsNumber`.

Page: docs/methods/dropest/bounds.md. Ported from dropestr's `R/cell_number.R`
(`EstimateCellsNumber`, `ArrayDerivative`, `GetLongestTrue`; unchanged from 0.6.0 to 0.7.9),
distributed under the GNU GPL v3, as sc-droplets is. `log` and `exp` are the C library's, as R
calls them; `round` halves to even, as R's does.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import pandas as pd

from sc_droplets.methods._common import result
from sc_droplets.primitives.rmath import clog
from sc_droplets.primitives.rounding import round_half_even
from sc_droplets.registry import Method, MethodInfo, Reference, Variant, param, register

CLASSES = ["high", "unknown", "low"]


class DropEstError(ValueError):
    """Where dropestr stops with an error."""


def longest_true_start(arr: np.ndarray) -> int:
    """`GetLongestTrue(arr)$start`: where the first of the longest runs of TRUE starts (1-based).
    With no TRUE at all, R's `which.max` of an all-FALSE vector gives the first run, at 1."""
    if arr.size == 0:
        raise DropEstError("no second derivative to look at: too few barcodes for the lag")
    change = np.flatnonzero(np.diff(arr.astype(np.int8))) + 1
    starts = np.concatenate([[0], change])
    lengths = np.diff(np.concatenate([starts, [arr.size]]))
    true_runs = arr[starts]
    if not true_runs.any():
        return 1
    longest = lengths[true_runs].max()
    return int(starts[np.flatnonzero(true_runs & (lengths == longest))[0]]) + 1


def estimate_cells_number(counts: np.ndarray, lag: float = 0.05) -> dict:
    """`EstimateCellsNumber(umi.counts, lag)` on counts sorted in decreasing order (all > 0).
    Returns expected, min, max and the lag in barcodes."""
    n = counts.size
    lag_n = int(round_half_even(n * lag))
    if lag_n < 1 or 2 * lag_n >= n:
        raise DropEstError(f"lag {lag} of {n} barcodes gives {lag_n}: R's diff needs a lag of at least 1 and more barcodes")
    log_u = clog(counts.astype(np.float64))
    log_rank = clog(np.arange(1, n + 1, dtype=np.float64))
    x = log_rank[lag_n:]
    y = (log_u[lag_n:] - log_u[:-lag_n]) / (log_rank[lag_n:] - log_rank[:-lag_n])
    x2 = x[lag_n:]
    y2 = (y[lag_n:] - y[:-lag_n]) / (x[lag_n:] - x[:-lag_n])
    start = longest_true_start(y2 > 0)
    max_num = int(round_half_even(math.exp(x2[start - 1])))
    count = max_num - lag_n
    if count < 0:
        raise DropEstError("the upper bound falls within the lag (R cannot index y[1:(max.num - lag)])")
    window = y[: max(count, 1)]  # y[1:0] is y[1] in R
    expected = int(round_half_even(math.exp(x[int(np.argmin(window))]) - lag_n / 2))
    return {"expected": expected, "max": max_num, "min": int(round_half_even(expected * 0.75)), "lag": lag_n}


@register
class Bounds(Method):
    info = MethodInfo(
        id="dropest-bounds",
        title="dropEst cell-number bounds",
        tool="dropEst",
        category="knee",
        summary="dropestr's EstimateCellsNumber: lagged derivatives of log count against log rank give the "
                "expected number of cells and the High / Unknown / Low bands.",
        docs="methods/dropest/bounds.md",
        variants=(Variant("default", since="0.6.0", summary="EstimateCellsNumber"),),
        latest_known="0.7.9",
        references=(
            Reference(
                "Petukhov, Guo, Baryawno, Severe, Scadden, Samsonova & Kharchenko 2018, dropEst: pipeline for "
                "accurate estimation of molecular counts in droplet-based single-cell RNA-seq experiments, "
                "Genome Biology 19:78",
                doi="10.1186/s13059-018-1449-6",
            ),
        ),
        source="https://github.com/kharchenkolab/dropEst/blob/v0.8.6/dropestr/R/cell_number.R",
        original_license="GPL-3.0",
    )

    @dataclass
    class Params:
        lag: float = param(0.05, "The derivatives' lag, as a fraction of the barcodes", original="lag", minimum=0, maximum=0.5)

    def run(self, droplets, params, variant, random_state):
        values = droplets.totals
        index = np.flatnonzero(values > 0)
        order = index[np.argsort(-values[index], kind="stable")]  # ties keep the matrix's order
        found = estimate_cells_number(values[order], params.lag)
        classes = np.full(values.size, CLASSES.index("low"), dtype=np.int8)
        classes[order[found["min"]:found["max"] - 1]] = CLASSES.index("unknown")
        classes[order[: found["min"]]] = CLASSES.index("high")
        is_cell = np.zeros(values.size, dtype=bool)
        is_cell[order[: found["expected"]]] = True
        return result(
            type(self), variant, params, droplets, is_cell,
            classes=pd.Categorical.from_codes(classes, categories=CLASSES),
            thresholds={**{k: float(v) for k, v in found.items()}, "umis": float(values[order[found["expected"] - 1]])},
        )
