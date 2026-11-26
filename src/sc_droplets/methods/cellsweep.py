"""cellsweep's empty-droplet step (`infer_empty_droplets`), release by release.

Page: docs/methods/cellsweep/empty.md. Ported from cellsweep's `cellsweep/utils/data_utils.py`
(`automatic_umi_cutoff_detection` to 1.0.1, `get_umi_cutoff_from_kneed` in 1.0.2, `gmm` and
`get_umi_cutoff_from_adata` throughout), under the BSD 2-Clause licence:

    Copyright (c) 2024, Pachter Lab
    Redistribution and use in source and binary forms, with or without modification, are
    permitted provided that the following conditions are met: 1. Redistributions of source code
    must retain the above copyright notice, this list of conditions and the following
    disclaimer. 2. Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation and/or other
    materials provided with the distribution. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS
    AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES ARE DISCLAIMED.

Every step calls the numpy, scipy and scikit-learn functions cellsweep calls, in the dtype of
cellsweep's matrix (float32 by default), so with the same libraries the cutoff is the same.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from sc_droplets.methods._common import result
from sc_droplets.primitives import kneed
from sc_droplets.registry import Method, MethodInfo, Reference, Variant, param, register

TOOL = "cellsweep"
LATEST = "1.0.2"
#: kneed release -> the Kneedle behaviour it has (primitives/kneed.py).
KNEED_BEHAVIOUR = {"0.8.5": "0.7", "0.8.6": "0.8.6", "0.9.0": "0.9"}
#: cellsweep's `min_counts` for its automatic rules.
MIN_COUNTS = 10
#: GaussianMixture's seed in `gmm`.
GMM_SEED = 42


class CellsweepError(ValueError):
    """Where cellsweep itself stops with an error."""


def threshold_cutoff(totals: np.ndarray):
    """`automatic_umi_cutoff_detection` (0.1.0 – 1.0.1): the total at the most negative second
    derivative of the smoothed log10 rank curve. Returns (cutoff, sorted totals, smoothed, d2)."""
    from scipy.ndimage import gaussian_filter1d

    counts = np.sort(totals)[::-1]
    counts = counts[counts > MIN_COUNTS]
    if counts.size < 2:
        raise CellsweepError(f"fewer than two barcodes have more than {MIN_COUNTS} counts")
    ranks = np.arange(1, len(counts) + 1)
    y = np.log10(counts)
    smooth = gaussian_filter1d(y, sigma=2)
    d2 = np.gradient(np.gradient(smooth, ranks), ranks)
    return counts[np.argmin(d2)], counts, smooth, d2


def kneed_cutoff(totals: np.ndarray, kneed_version: str):
    """`get_umi_cutoff_from_kneed` (1.0.2): kneed's online Kneedle on (last rank of each distinct
    total, total). Returns (cutoff, the Kneedle curves, knee index)."""
    t = totals[totals > MIN_COUNTS]
    umis, n_barcodes = np.unique(t, return_counts=True)
    umis = umis[::-1]
    ranks = np.cumsum(n_barcodes[::-1])
    if umis.size <= 1:
        raise CellsweepError("Could not detect a knee in the barcode rank curve (fewer than two distinct totals)")
    found = kneed.prepare(ranks, umis, curve="convex", direction="decreasing", behaviour=KNEED_BEHAVIOUR[kneed_version])
    index = found.find(1.0, online=True)
    if index is None:
        raise CellsweepError("Could not detect a knee in the barcode rank curve. Please provide either umi_cutoff or expected_cells.")
    # knee_y = y[x == knee][0]: ranks are distinct, so that is the total at the knee's index
    return umis[index], found, index


def mx_filter_cutoff(totals: np.ndarray):
    """`gmm(x, log1p(x), [2])`: the total whose two-component posterior has the highest entropy.
    Returns (cutoff, sorted totals, entropy, means, standard deviations, weights)."""
    from scipy.stats import entropy
    from sklearn.mixture import GaussianMixture

    x = np.sort(totals)[::-1]
    v = np.log1p(x).reshape(-1, 1)
    gm = GaussianMixture(n_components=2, random_state=GMM_SEED)
    gm.fit_predict(v)
    ent = entropy(gm.predict_proba(v), axis=1)
    return x[np.argmax(ent)], x, ent, gm


@register
class Empty(Method):
    info = MethodInfo(
        id="cellsweep-empty",
        title="cellsweep empty droplets",
        tool=TOOL,
        category="knee",
        summary=(
            "cellsweep's infer_empty_droplets: a cutoff on the totals (second derivative of the "
            "smoothed log curve, kneed's online knee, or a two-component mixture); below it, empty."
        ),
        docs="methods/cellsweep/empty.md",
        variants=(
            Variant("v0.1", since="0.1.0", until="1.0.2", summary="default rule: threshold (second derivative)",
                    defaults={"empty_method": "threshold"}),
            Variant("v1.0.2", since="1.0.2", summary="default rule: kneed (online Kneedle)",
                    defaults={"empty_method": "kneed"}),
        ),
        latest_known=LATEST,
        references=(
            Reference(
                "Caskey, Rich, Weber, Mortazavi, Pachter & Hallgrimsdottir 2026, Single-cell genomics "
                "decontamination with CellSweep, bioRxiv",
                doi="10.64898/2026.03.04.709349",
            ),
        ),
        source="https://github.com/pachterlab/cellsweep/blob/c81ab360/cellsweep/utils/data_utils.py",
        original_license="BSD-2-Clause",
    )

    @dataclass
    class Params:
        empty_method: Optional[str] = param(
            None,
            "The automatic rule: threshold (to 1.0.1), kneed (1.0.2) or mx_filter",
            original="empty_droplet_method",
            choices=("threshold", "kneed", "mx_filter"),
        )
        umi_cutoff: Optional[float] = param(None, "The cutoff itself (wins over everything)", original="umi_cutoff", minimum=0)
        expected_cells: Optional[int] = param(None, "The cutoff is the total of the barcode ranked this",
                                              original="expected_cells", minimum=1)
        kneed_version: str = param("0.9.0", "The installed kneed, whose behaviour the kneed rule follows",
                                   choices=tuple(KNEED_BEHAVIOUR))
        matrix_dtype: str = param("float32", "The dtype of cellsweep's matrix, in which it sums and computes",
                                  choices=("float32", "float64"))

    def run(self, droplets, params, variant, random_state):
        totals = droplets.totals.astype(params.matrix_dtype)
        method = params.empty_method
        if variant.name == "v0.1" and method == "kneed":
            raise CellsweepError("cellsweep before 1.0.2 has no 'kneed' rule (use --tool-version 1.0.2)")
        if variant.name == "v1.0.2" and method == "threshold":
            method = "manual"  # 1.0.2 renamed it; it then needs a cutoff or expected cells
        diagnostics = {}
        if params.umi_cutoff is not None:
            method = "manual"
        elif params.expected_cells is not None and not (variant.name == "v0.1" and method == "mx_filter"):
            method = "manual"  # before 1.0.2, mx_filter ignores expected_cells
        if method == "manual":
            if params.umi_cutoff is not None:
                cutoff = np.asarray(params.umi_cutoff, dtype=totals.dtype)[()]
            elif params.expected_cells is not None:
                ranked = np.sort(totals)[::-1]
                if params.expected_cells > ranked.size:
                    raise CellsweepError(f"expected_cells {params.expected_cells} is more than the {ranked.size} barcodes")
                cutoff = ranked[params.expected_cells - 1]
            else:
                raise CellsweepError("Empty droplet method 'manual' requires either umi_cutoff or expected_cells.")
        elif method == "threshold":
            cutoff, counts, smooth, d2 = threshold_cutoff(totals)
            diagnostics = _thin({"rank_curve": counts, "smoothed_log10": smooth, "second_derivative": d2})
        elif method == "kneed":
            cutoff, found, index = kneed_cutoff(totals, params.kneed_version)
            at = found.maximum_of(index)
            diagnostics = {
                **_thin({"difference_x": found.x_normalized, "difference_y": found.y_difference}),
                "knee_x": float(found.x_normalized[at]),
                "knee_y": float(found.y_difference[at]),
            }
        else:
            cutoff, x, ent, gm = mx_filter_cutoff(totals)
            diagnostics = {
                **_thin({"mixture_totals": x, "mixture_entropy": ent}),
                "mixture_means": gm.means_.ravel().astype(float),
                "mixture_sds": np.sqrt(gm.covariances_.ravel()).astype(float),
                "mixture_weights": gm.weights_.astype(float),
            }
        diagnostics["rule"] = method
        is_cell = ~(totals < cutoff)
        return result(
            type(self),
            variant,
            params,
            droplets,
            is_cell,
            thresholds={"umis": float(cutoff)},
            diagnostics=diagnostics,
        )


def _thin(curves: dict, most: int = 2000) -> dict:
    """Every curve at about `most` evenly spaced points, with their positions (`index`)."""
    n = len(next(iter(curves.values())))
    index = np.unique(np.concatenate([np.arange(0, n, max(1, n // most)), [n - 1]])).astype(np.int64)
    return {"curve_index": index, **{k: np.asarray(v, dtype=float)[index] for k, v in curves.items()}}
