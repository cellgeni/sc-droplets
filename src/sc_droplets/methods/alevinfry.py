"""alevin-fry `generate-permit-list`: the knee, expected cells, forced cells and the read floor.

Page: docs/methods/alevinfry/permit-list.md. Ported from alevin-fry's `knee_finding.rs`
(`distance_to_line`, `get_max_distance_index`, `get_knee`) and `cellfilter.rs`
(`select_retained_barcodes`, 0.19.0; the multi-barcode rules of 0.12.0 - 0.17.1), under the
BSD 3-Clause licence:

    Copyright (c) 2020-, University of Maryland. All rights reserved. Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met: 1.
    Redistributions of source code must retain the above copyright notice, this list of
    conditions and the following disclaimer. 2. Redistributions in binary form must reproduce the
    above copyright notice, this list of conditions and the following disclaimer in the
    documentation and/or other materials provided with the distribution. 3. Neither the name of
    the copyright holder nor the names of its contributors may be used to endorse or promote
    products derived from this software without specific prior written permission. THIS
    SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
    IMPLIED WARRANTIES ARE DISCLAIMED.

The distances are computed in alevin-fry's order of operations, in float64, on 64-bit integer
counts, so the same counts give the same knee.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Tuple

import numpy as np

from sc_droplets.methods._common import result
from sc_droplets.methods.umitools import _last_pass_curve
from sc_droplets.primitives.rounding import round_half_away
from sc_droplets.registry import Method, MethodInfo, Reference, Variant, param, register

TOOL = "alevin-fry"
LATEST = "0.19.0"
MAX_ITERATIONS = 100
ITER_SLACK = 5


class AlevinFryError(ValueError):
    """Where alevin-fry itself stops (a failed assertion)."""


def max_distance_index(cumulative: np.ndarray) -> int:
    """`get_max_distance_index(cfreq, true)`: the last point farthest from the line through
    (0, first / last) and (1, 1), with x = index / length and y = count / last."""
    n = cumulative.size
    if n < 2:
        raise AlevinFryError(f"when attempting to find a knee-distance threshold, the list of putative cells is only of length {n}")
    max_x = float(n)
    max_y = float(cumulative[-1])
    x1, y1 = 0.0, float(cumulative[0]) / max_y
    x2, y2 = 1.0, float(cumulative[-1]) / max_y
    x0 = np.arange(n, dtype=np.float64) / max_x
    y0 = cumulative.astype(np.float64) / max_y
    numer = np.abs((y2 - y1) * x0 - (x2 - x1) * y0 + x2 * y1 - y2 * x1)
    denom = np.sqrt((y2 - y1) * (y2 - y1) + (x2 - x1) * (x2 - x1))
    distance = numer / denom
    return int(n - 1 - np.argmax(distance[::-1]))  # `d >= max_d`: the last maximum wins


def get_knee(freq: np.ndarray) -> Tuple[int, List[Tuple[int, int]]]:
    """`get_knee(freq, 100)` on counts sorted in decreasing order. Returns the knee's index and
    each pass as (points searched, knee)."""
    cumulative = np.cumsum(freq.astype(np.uint64))
    knee = max_distance_index(cumulative)
    passes = [(cumulative.size, knee)]
    if knee == 0:
        raise AlevinFryError("get_knee determined a knee index of 0")
    previous, iterations = 0, 0
    while knee != previous:
        previous = knee
        iterations += 1
        if iterations > MAX_ITERATIONS:
            break
        last = min(cumulative.size - 1, knee * ITER_SLACK)
        knee = max_distance_index(cumulative[:last])
        passes.append((last, knee))
        if knee == 0:
            raise AlevinFryError("get_knee determined a knee index of 0")
    return knee, passes


def threshold(freq: np.ndarray, *, expect_cells: Optional[int] = None, force_cells: Optional[int] = None,
              min_reads: Optional[int] = None, flex_rules: bool = False):
    """The count a barcode needs, by the rule given (the knee when none is). Returns
    (threshold, knee or None, passes)."""
    n = freq.size
    if min_reads is not None:
        return int(min_reads), None, []
    if n == 0:
        raise AlevinFryError("no barcode has a count")
    if force_cells is not None:
        return int(freq[min(force_cells, n) - 1]), None, []
    if expect_cells is not None:
        if flex_rules:
            start = int(freq[0]) // expect_cells
            below = np.flatnonzero(freq < start)
            position = int(below[0]) if below.size else n
            cells = min(int(position * 10.0), n)
            return (int(freq[cells - 1]) if cells > 0 else 0), None, []
        robust = min(n - 1, int(round_half_away(expect_cells * 0.99)))
        return max(1, int(round_half_away(freq[robust] / 10.0))), None, []
    knee, passes = get_knee(freq)
    if flex_rules:
        return int(freq[knee - 1]), knee, passes
    return int(freq[min(knee, n - 1)]), knee, passes


@register
class PermitList(Method):
    info = MethodInfo(
        id="alevinfry-permit-list",
        title="alevin-fry generate-permit-list",
        tool=TOOL,
        category="knee",
        summary=(
            "alevin-fry's permit list: a distance knee (repeated on five times the knee's rank), "
            "a tenth of the count at rank 0.99 x expected cells, forced cells, or a read floor."
        ),
        docs="methods/alevinfry/permit-list.md",
        counts="either",
        variants=(
            Variant("v0.1", since="0.1.0", until="0.4.2", summary="knee, expected cells, forced cells"),
            Variant("v0.4.2", since="0.4.2", until="0.12.0", summary="adds the unfiltered permit list"),
            Variant("v0.12", since="0.12.0", until="0.18.0",
                    summary="multi-barcode (Flex) input has rules of its own"),
            Variant("v0.18", since="0.18.0", summary="one rule set for single- and multi-barcode input"),
        ),
        latest_known=LATEST,
        references=(
            Reference(
                "He, Zakeri, Sarkar, Soneson, Srivastava & Patro 2022, Alevin-fry unlocks rapid, "
                "accurate and memory-frugal quantification of single-cell RNA-seq data, Nature Methods 19:316-322",
                doi="10.1038/s41592-022-01408-3",
            ),
        ),
        source="https://github.com/COMBINE-lab/alevin-fry/blob/v0.19.0/src/cellfilter.rs",
        original_license="BSD-3-Clause",
    )

    @dataclass
    class Params:
        expect_cells: Optional[int] = param(None, "Threshold from expected cells: a tenth of the count at rank 0.99 E",
                                            original="--expect-cells", minimum=1)
        force_cells: Optional[int] = param(None, "Threshold at the count of the F-th barcode",
                                           original="--force-cells", minimum=1)
        unfiltered: bool = param(False, "Keep every barcode with at least --min-reads", original="--unfiltered-pl")
        min_reads: int = param(10, "The floor of --unfiltered", original="--min-reads", minimum=1)
        flex: bool = param(False, "Multi-barcode (10x Flex) input: 0.12-0.17 use rules of their own",
                           original="multi-barcode RAD input")
        unit: str = param("umis", "Count UMIs (the matrix) or reads (--reads)", original="reads per barcode",
                          choices=("umis", "reads"))

    def run(self, droplets, params, variant, random_state):
        chosen = [name for name, given in (("--expect-cells", params.expect_cells is not None),
                                           ("--force-cells", params.force_cells is not None),
                                           ("--unfiltered", params.unfiltered)) if given]
        if len(chosen) > 1:
            raise AlevinFryError(f"give one of {', '.join(chosen)} (alevin-fry takes one rule)")
        if params.flex and variant.name in ("v0.1", "v0.4.2"):
            raise AlevinFryError("alevin-fry before 0.12.0 takes no multi-barcode input (use --tool-version 0.12.0 or later)")
        if params.unfiltered and variant.name == "v0.1":
            raise AlevinFryError("alevin-fry before 0.4.2 has no unfiltered permit list (use --tool-version 0.4.2 or later)")
        values = droplets.values(params.unit)
        index = np.flatnonzero(values > 0)
        counts = values[index]
        if not np.all(counts == np.floor(counts)):
            raise AlevinFryError("alevin-fry counts whole reads; these counts have fractions")
        counts = counts.astype(np.int64)
        freq = np.sort(counts)[::-1]
        cut, knee, passes = threshold(
            freq,
            expect_cells=params.expect_cells,
            force_cells=params.force_cells,
            min_reads=params.min_reads if params.unfiltered else None,
            flex_rules=params.flex and variant.name == "v0.12",
        )
        is_cell = np.zeros(values.shape[0], dtype=bool)
        is_cell[index[counts >= cut]] = True
        thresholds = {params.unit: float(cut)}
        diagnostics = {"unit": params.unit}
        if knee is not None:
            thresholds["knee_rank"] = float(knee + 1)
            diagnostics.update({"passes": np.asarray(passes, dtype=np.int64),
                                "curve_title": "Cumulative counts (alevin-fry knee)",
                                **_last_pass_curve(np.cumsum(freq), passes)})
        return result(type(self), variant, params, droplets, is_cell, thresholds=thresholds, diagnostics=diagnostics)
