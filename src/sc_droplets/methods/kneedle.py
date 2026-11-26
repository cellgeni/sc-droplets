"""Kneedle (kneed's KneeLocator) on the barcode-rank curve, by default as yascp runs it.

Page: docs/methods/kneedle/kneedle.md. Kneedle itself is `primitives/kneed.py` (BSD-3, from
kneed); yascp's search for the sensitivity (GPL-3, `bin/get_estimates_from_umi_counts.py`,
`kneedle_cutoff`) is reimplemented here: S from 5000 down by 100 until a knee is found.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from sc_droplets.methods._common import result
from sc_droplets.primitives import kneed
from sc_droplets.primitives.rmath import clog10
from sc_droplets.registry import Method, MethodInfo, Reference, Variant, param, register

TOOL = "Kneedle (kneed)"
LATEST = "0.9.0"


class KneedleError(ValueError):
    """No knee, or a curve Kneedle cannot be run on."""


_VARIANTS = (
    Variant(
        "v0.5.1",
        since="0.5.1",
        until="0.7.0",
        summary="extrema include plateaus; polynomial by scikit-learn; convex curves mirrored in x",
        settings={"behaviour": "0.5.1"},
    ),
    Variant(
        "v0.7",
        since="0.7.0",
        until="0.8.6",
        summary="polynomial by numpy.polyfit with a degree; only y transformed",
        settings={"behaviour": "0.7"},
    ),
    Variant(
        "v0.8.6",
        since="0.8.6",
        until="0.9.0",
        summary="no detection between a local minimum and the next maximum",
        settings={"behaviour": "0.8.6"},
    ),
    Variant(
        "v0.9",
        since="0.9.0",
        summary="no interpolation; concave decreasing and convex increasing curves mirrored in x",
        settings={"behaviour": "0.9"},
    ),
)


def rank_curve(values: np.ndarray, lower: float, log_x: bool, log_y: bool):
    """The barcodes on the curve (in rank order), and its x and y."""
    order = np.argsort(-values, kind="stable")
    ranked = values[order]
    points = int(np.count_nonzero(ranked >= lower))
    if log_y and points and ranked[points - 1] <= 0:
        raise KneedleError("--log-y needs every barcode on the curve to have a count: raise --lower above 0")
    if points < 2:
        raise KneedleError(f"only {points} barcode(s) have at least {lower:g} counts; Kneedle needs two points")
    ranks = np.arange(1, points + 1, dtype=np.int64)
    x = clog10(ranks.astype(float)) if log_x else ranks
    y = clog10(ranked[:points]) if log_y else ranked[:points]
    return order[:points], x, y


def search(found: kneed.Kneedle, sensitivity: Optional[float], s_start: float, s_step: float, online: bool):
    """yascp's search: S from `s_start` down by `s_step` until a knee; or one given S."""
    if sensitivity is not None:
        return found.find(sensitivity, online), sensitivity
    S = s_start
    while True:
        index = found.find(S, online)
        if index is not None:
            return index, S
        if S < 0:
            # Below 0 every threshold is above its maximum: no S will find a knee.
            return None, S
        S -= s_step


@register
class Kneedle(Method):
    info = MethodInfo(
        id="kneedle",
        title="Kneedle on the barcode-rank curve",
        tool=TOOL,
        category="knee",
        summary=(
            "kneed's KneeLocator on rank against count; by default yascp's use of it (convex, "
            "decreasing, sensitivity lowered from 5000 until a knee appears)."
        ),
        docs="methods/kneedle/kneedle.md",
        counts="either",
        variants=_VARIANTS,
        latest_known=LATEST,
        references=(
            Reference(
                "Satopää, Albrecht, Irwin & Raghavan 2011, Finding a \"Kneedle\" in a haystack: "
                "detecting knee points in system behavior, ICDCS Workshops 166-171",
                doi="10.1109/ICDCSW.2011.20",
            ),
            Reference(
                "yascp, bin/get_estimates_from_umi_counts.py (kneedle_cutoff)",
                url="https://github.com/wtsi-hgi/yascp/blob/main/bin/get_estimates_from_umi_counts.py",
            ),
        ),
        source="https://github.com/arvkevi/kneed/blob/v0.9.0/kneed/knee_locator.py",
        original_license="BSD-3-Clause",
    )

    @dataclass
    class Params:
        lower: float = param(1000.0, "Barcodes with fewer counts are not on the curve",
                             original="yascp --lower_bound_umis_estimate_ncells", minimum=0)
        sensitivity: Optional[float] = param(None, "One sensitivity S, instead of yascp's search", original="S")
        s_start: float = param(5000.0, "The search's first sensitivity", original="yascp")
        s_step: float = param(100.0, "How much the search lowers the sensitivity each time", original="yascp", minimum=1e-9)
        curve: str = param("convex", "Shape of the curve", original="curve", choices=kneed.CURVES)
        direction: str = param("decreasing", "Direction of the curve", original="direction", choices=kneed.DIRECTIONS)
        interp_method: str = param("interp1d", "Smoothing before the search", original="interp_method",
                                   choices=("interp1d", "polynomial"))
        polynomial_degree: int = param(7, "Degree of the polynomial (kneed 0.7.0+)", original="polynomial_degree", minimum=1)
        online: bool = param(False, "Keep looking after the first knee; the last one wins", original="online")
        log_x: bool = param(False, "Rank on a log10 axis")
        log_y: bool = param(False, "Counts on a log10 axis")
        unit: str = param("umis", "Count UMIs (the matrix) or reads (--reads)", choices=("umis", "reads"))

    def run(self, droplets, params, variant, random_state):
        behaviour = variant.settings["behaviour"]
        values = droplets.values(params.unit)
        barcodes, x, y = rank_curve(values, params.lower, params.log_x, params.log_y)
        found = kneed.prepare(
            x,
            y,
            curve=params.curve,
            direction=params.direction,
            interp_method=params.interp_method,
            polynomial_degree=params.polynomial_degree,
            behaviour=behaviour,
        )
        index, S = search(found, params.sensitivity, params.s_start, params.s_step, params.online)
        if index is None:
            raise KneedleError(
                f"Kneedle finds no knee on this curve ({barcodes.size} barcodes with at least "
                f"{params.lower:g} counts" + (f", S = {S:g}" if params.sensitivity is not None else ", at any sensitivity") + ")"
            )
        rank = index + 1
        is_cell = np.zeros(values.shape[0], dtype=bool)
        is_cell[barcodes[:rank]] = True
        return result(
            type(self),
            variant,
            params,
            droplets,
            is_cell,
            thresholds={"knee_rank": float(rank), params.unit: float(values[barcodes[index]]), "sensitivity": float(S)},
            diagnostics=_diagnostics(found, index, S),
        )


def _diagnostics(found: kneed.Kneedle, index: int, S: float, most: int = 2000, most_maxima: int = 500) -> dict:
    """The difference curve, thinned to about `most` points (the knee's maximum kept), its
    first maxima and their thresholds at the sensitivity used, and the knee on it."""
    n = found.x_normalized.shape[0]
    at = found.maximum_of(index)
    step = max(1, n // most)
    keep = np.unique(np.concatenate([np.arange(0, n, step), [n - 1, at]])).astype(np.int64)
    maxima = found.maxima[:most_maxima]
    return {
        "difference_x": found.x_normalized[keep],
        "difference_y": found.y_difference[keep],
        "maxima_x": found.x_normalized[maxima],
        "maxima_y": found.y_difference[maxima],
        "maxima_threshold": found.thresholds(S)[:most_maxima],
        "knee_x": float(found.x_normalized[at]),
        "knee_y": float(found.y_difference[at]),
        "points": n,
    }
