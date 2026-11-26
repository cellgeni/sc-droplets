"""CellBender's priors: what `remove-background` decides before it trains, release by release.

Page: docs/methods/cellbender/prior.md. Ported from CellBender's
`cellbender/remove_background/data/dataset.py` (0.1.0, 0.2.2; 0.3.0 for `__init__`,
`_trim_noiseless_features` and `_trim_droplets`) and `cellbender/remove_background/data/priors.py`
(0.3.0, unchanged in 0.3.2 and 0.4.0), under the BSD 3-Clause licence:

    Copyright (c) 2019, Broad Institute, Inc. All rights reserved.
    Redistribution and use in source and binary forms, with or without modification, are
    permitted provided that the following conditions are met: * Redistributions of source code
    must retain the above copyright notice, this list of conditions and the following
    disclaimer. * Redistributions in binary form must reproduce the above copyright notice, this
    list of conditions and the following disclaimer in the documentation and/or other materials
    provided with the distribution. * Neither the name Broad Institute, Inc. nor the names of its
    contributors may be used to endorse or promote products derived from this software without
    specific prior written permission. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
    CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES ARE DISCLAIMED.

CellBender's totals are integers (int64 row sums of its int32 matrix); every step runs the numpy
and scipy calls CellBender makes, on int64 totals where it has them (np.mean of int64 sums in
buffered blocks, which can differ in the last bit from a float64 array). Ties are ranked as
numpy's classic quicksort leaves them (`legacy_argsort`).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional

import numpy as np
import pandas as pd

from sc_droplets.methods._common import result
from sc_droplets.primitives.legacysort import legacy_argsort
from sc_droplets.registry import Method, MethodInfo, Reference, Variant, param, register

TOOL = "CellBender"
LATEST = "0.4.0"
CLASSES = ["cell", "transition", "empty", "excluded"]

#: 0.1.0 and 0.2.x (consts.py).
TOTAL_DROPLET_DEFAULT = 25000
LOW_UMI_CUTOFF_V1 = 15
FRACTION_EMPTIES_V1 = 0.5
#: EMPIRICAL_LOW_UMI_TO_EMPTY_DROPLET_THRESHOLD: 0.8 in 0.1.0, 0.5 from 0.2.0.
LOW_CUT_FACTOR = {"v0.1": 0.8, "v0.2": 0.5, "v0.3": 0.5}
#: 0.3.0 and later (consts.py).
LOW_UMI_CUTOFF = 5
FRACTION_EMPTIES = 0.2
AMBIENT_COUNTS_IN_CELLS_LOW_LIMIT = 0.1
MAX_TOTAL_DROPLETS_GUESSED = 70000
MAX_EMPTIES_TO_INCLUDE = 20000
NUM_EMPTIES_INCREMENT = 20000


class CellBenderError(ValueError):
    """Where CellBender itself stops (an assertion, or no data for a step)."""


# --- 0.3.0 and later: priors.py ---------------------------------------------------------------


def threshold_otsu(values: np.ndarray, n_bins: int = 256) -> float:
    """`_threshold_otsu`: scikit-image's fast Otsu on a 256-bin histogram (float32 bin counts)."""
    counts, edges = np.histogram(values.reshape(-1), n_bins)
    centers = (edges[:-1] + edges[1:]) / 2
    counts = counts.astype("float32", copy=False)
    weight1 = np.cumsum(counts)
    weight2 = np.cumsum(counts[::-1])[::-1]
    mean1 = np.cumsum(counts * centers) / weight1
    mean2 = (np.cumsum((counts * centers)[::-1]) / weight2[::-1])[::-1]
    gap = mean1[:-1] - mean2[1:]
    variance12 = weight1[:-1] * weight2[1:] * (gap * gap)
    return centers[np.argmax(variance12)]


def peak_density(counts: np.ndarray, cutoff: float, cell_count_low_limit: float):
    """`_peak_density_given_cutoff`: the empty prior (the density peak of 1000 log-spaced
    resampled log counts at or below `cutoff`) and its upper limit. Returns
    (empty prior, upper limit, grid, density)."""
    from scipy.stats import gaussian_kde

    noncell = counts[counts <= cutoff]
    n_cells = (counts > cell_count_low_limit).sum()
    n_empties = len(noncell)
    if n_empties == 0:
        raise CellBenderError(f"no barcode at or below the empty-droplet cutoff {cutoff:g}")
    with np.errstate(divide="ignore", invalid="ignore"):
        inds = np.logspace(np.log10(n_cells), np.log10(n_cells + n_empties), num=1000, base=10)
    inds = [max(0, min(int(ind - n_cells), n_empties - 1)) for ind in inds]
    logs = np.log(np.sort(noncell)[::-1][inds])
    grid = np.arange(np.floor(logs.min()) - 0.01, np.ceil(logs.max()) + 0.01, 0.1)
    try:
        density = gaussian_kde(logs).evaluate(grid)
    except np.linalg.LinAlgError as error:
        raise CellBenderError(f"the counts below {cutoff:g} have no spread for a kernel density ({error})") from None
    peak = np.argmax(density)
    empty_count_prior = np.exp(grid[peak])
    below = np.where(density[peak:] < 0.6 * np.max(density))[0]
    step = below[0] if len(below) > 0 else len(density[peak:]) - 1
    return empty_count_prior, np.exp(grid[peak:][step]), grid, density


def cell_count_empty_count(u: np.ndarray, low_count_threshold: float) -> Dict[str, object]:
    """`get_cell_count_empty_count`: Otsu's lower limit for cells, then the kernel-density peak
    of the empty droplets, refined up to five times; then the cell prior as a median."""
    descending = np.sort(u)[::-1]
    for_otsu = descending[: (u > low_count_threshold).sum() // 4]
    if for_otsu.size == 0:
        raise CellBenderError(f"too few barcodes above the low-count threshold {low_count_threshold:g} for Otsu's threshold")
    cell_count_low_limit = np.exp(threshold_otsu(np.log(for_otsu)))
    cell_count_prior = np.mean(u[u > cell_count_low_limit])
    for_kde = descending[descending > low_count_threshold]
    if for_kde.size < 100:
        raise CellBenderError(f"fewer than 100 barcodes above the low-count threshold {low_count_threshold:g}")
    cutoff = max(0.1 * cell_count_low_limit, for_kde[-100])
    empty_count_prior = -100.0
    delta = np.inf
    rounds = 0
    while delta > 10:
        new_empty, upper, grid, density = peak_density(for_kde, cutoff, cell_count_low_limit)
        cutoff = 0.75 * np.sqrt(cell_count_prior * new_empty)
        delta = np.abs(new_empty - empty_count_prior)
        empty_count_prior = new_empty
        rounds += 1
        if rounds >= 5:
            break
    crossover = np.sqrt(cell_count_prior * empty_count_prior)
    cell_count_prior = np.median(u[u > crossover])
    return {
        "cell_counts": cell_count_prior,
        "empty_counts": empty_count_prior,
        "empty_count_upper_limit": upper,
        "_cell_count_low_limit": cell_count_low_limit,
        "_kde": (grid, density, rounds),
    }


def expected_cells_and_total_droplets(u, cell_counts, empty_counts, empty_count_upper_limit,
                                      max_empties: int = MAX_EMPTIES_TO_INCLUDE) -> Dict[str, object]:
    """`get_expected_cells_and_total_droplets`."""
    expected_cells = (u >= cell_counts).sum()
    total_droplets = (u >= np.sqrt(empty_counts * empty_count_upper_limit)).sum()
    transition_point = (u >= np.sqrt(cell_counts * empty_counts)).sum()
    total_droplets = min(total_droplets, transition_point + max_empties)
    return {"expected_cells": expected_cells, "total_droplets": total_droplets, "transition_point": transition_point}


def heuristic_priors(u: np.ndarray, low_count_threshold: float) -> Dict[str, object]:
    """`get_priors` without its last step: rerun with a higher low-count threshold (at most six
    times) while the total droplets exceed 70,000."""
    def once(threshold):
        priors = cell_count_empty_count(u, threshold)
        priors.update(expected_cells_and_total_droplets(
            u, priors["cell_counts"], priors["empty_counts"], priors["empty_count_upper_limit"]))
        return priors

    priors = once(low_count_threshold)
    retries = 0
    while priors["total_droplets"] > MAX_TOTAL_DROPLETS_GUESSED:
        priors = once(priors["empty_count_upper_limit"])
        retries += 1
        if retries > 5:
            break
    priors["_retries"] = retries
    return priors


def cell_count_given_expected_cells(u: np.ndarray, expected_cells: int) -> float:
    """`get_cell_count_given_expected_cells`: the geometric mean of the top `expected_cells`."""
    top = np.sort(u)[::-1][:expected_cells]
    return np.exp(np.mean(np.log(top))).item()


def empty_count_given_expected_cells_and_total_droplets(u, expected_cells, total_droplets):
    """`get_empty_count_given_expected_cells_and_total_droplets`: (empty prior, upper limit)."""
    descending = np.sort(u)[::-1]
    start = max(expected_cells, total_droplets - 500)
    window = descending[int(start):int(total_droplets)]
    if window.size == 0:
        # CellBender takes the median of the empty slice (NaN) and then divides by zero rows.
        raise CellBenderError(f"expected cells ({int(expected_cells)}) must be fewer than the total droplets "
                              f"({int(total_droplets)}); CellBender stops here with a division by zero")
    empty_counts = np.median(window).item()
    cell_counts = np.exp(np.mean(np.log(descending[:expected_cells]))).item()
    return empty_counts, min(np.sqrt(cell_counts * empty_counts), 1.5 * empty_counts)


def crossover_surely_empty_and_stds(u: np.ndarray, priors: dict) -> None:
    """`compute_crossover_surely_empty_and_stds`, in place."""
    descending = np.sort(u)[::-1]
    if priors["total_droplets"] >= descending.size:
        raise CellBenderError(f"total droplets {int(priors['total_droplets'])} is not below the {descending.size} barcodes")
    surely_empty_counts = descending[priors["total_droplets"]]
    log_counts_crossover = (np.log(surely_empty_counts) + np.log(priors["cell_counts"])) / 2
    logs = np.log(u[u > 0])
    priors["log_counts_crossover"] = log_counts_crossover
    priors["surely_empty_counts"] = surely_empty_counts
    priors["d_std"] = np.std(logs[logs > log_counts_crossover]).item() / 5.0


def kept_features(matrix, u, priors, low_count_threshold, threshold) -> np.ndarray:
    """`_trim_noiseless_features`: features whose mean count in the empty droplets
    (L < u < upper limit), times max(expected cells, 5000), is above `threshold`."""
    rows = (u < priors["empty_count_upper_limit"]) & (u > low_count_threshold)
    # scipy's mean sums x * (1 / n) in float64 for integer and float64 matrices alike; a float32
    # matrix would sum in float32, so it is widened first.
    empties = matrix[rows, :].astype(np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        mean_per_feature = np.asarray(empties.mean(axis=0)).ravel()
    return max(priors["expected_cells"], 5000) * mean_per_feature > threshold


# --- the method ------------------------------------------------------------------------------


@register
class Prior(Method):
    info = MethodInfo(
        id="cellbender-prior",
        title="CellBender priors",
        tool=TOOL,
        category="density",
        summary=(
            "CellBender remove-background's priors: probable cells, the other analysed droplets "
            "and the surely-empty droplets, from mode and median (0.1, 0.2) or Otsu and a kernel "
            "density peak (0.3+)."
        ),
        docs="methods/cellbender/prior.md",
        variants=(
            Variant("v0.1", since="0.1.0", until="0.2.0",
                    summary="mode and log-median priors; low cut 0.8 x empty; 25,000 droplets",
                    defaults={"total_droplets_included": TOTAL_DROPLET_DEFAULT,
                              "low_count_threshold": LOW_UMI_CUTOFF_V1,
                              "empty_drop_training_fraction": FRACTION_EMPTIES_V1}),
            Variant("v0.2", since="0.2.0", until="0.3.0",
                    summary="low cut 0.5 x empty; median of the top E when E is given",
                    defaults={"total_droplets_included": TOTAL_DROPLET_DEFAULT,
                              "low_count_threshold": LOW_UMI_CUTOFF_V1,
                              "empty_drop_training_fraction": FRACTION_EMPTIES_V1}),
            Variant("v0.3", since="0.3.0",
                    summary="Otsu and kernel-density priors; expected cells and total droplets estimated",
                    defaults={"low_count_threshold": LOW_UMI_CUTOFF,
                              "empty_drop_training_fraction": FRACTION_EMPTIES,
                              "projected_ambient_count_threshold": AMBIENT_COUNTS_IN_CELLS_LOW_LIMIT}),
        ),
        latest_known=LATEST,
        references=(
            Reference(
                "Fleming, Chaffin, Arduini, Akkad, Banks, Marioni, Philippakis, Ellinor & Babadi 2023, "
                "Unsupervised removal of systematic background noise from droplet-based single-cell "
                "experiments using CellBender, Nature Methods 20:1323-1335",
                doi="10.1038/s41592-023-01943-7",
            ),
        ),
        source="https://github.com/broadinstitute/CellBender/blob/v0.3.0/cellbender/remove_background/data/priors.py",
        original_license="BSD-3-Clause",
        feature_types=None,
    )

    @dataclass
    class Params:
        expected_cells: Optional[int] = param(None, "The number of probable cells (estimated if not given)",
                                              original="--expected-cells", minimum=1)
        total_droplets_included: Optional[int] = param(
            None, "Droplets analysed, cells included (25,000 to 0.2; estimated from 0.3)",
            original="--total-droplets-included", minimum=1)
        low_count_threshold: Optional[int] = param(None, "Barcodes at or below this count are left out (15 to 0.2; 5 from 0.3)",
                                                   original="--low-count-threshold", minimum=0)
        force_cell_umi_prior: Optional[float] = param(None, "The cell prior itself (0.3+)",
                                                      original="--force-cell-umi-prior", minimum=0)
        force_empty_umi_prior: Optional[float] = param(None, "The empty prior itself (0.3+)",
                                                       original="--force-empty-umi-prior", minimum=0)
        empty_drop_training_fraction: Optional[float] = param(
            None, "Enters the prior cell logit only (0.5 to 0.2; 0.2 from 0.3)",
            original="--empty-drop-training-fraction", minimum=0, maximum=1)
        projected_ambient_count_threshold: Optional[float] = param(
            None, "Keep features projected to have more ambient counts in cells than this (0.3+; 0.1)",
            original="--projected-ambient-count-threshold", minimum=0)

    def run(self, droplets, params, variant, random_state):
        if variant.name != "v0.3":
            given = [name for name in ("force_cell_umi_prior", "force_empty_umi_prior", "projected_ambient_count_threshold")
                     if getattr(params, name) is not None]
            if given:
                raise CellBenderError(f"CellBender before 0.3.0 has no --{given[0].replace('_', '-')} (use --tool-version 0.3.0)")
        u = _integer(droplets.totals)
        if variant.name == "v0.3":
            out = _priors_v03(droplets.counts, u, params)
        else:
            out = _priors_v01(u, params, variant.name)
        classes = np.full(u.size, CLASSES.index("excluded"), dtype=np.int8)
        classes[out.pop("_empty")] = CLASSES.index("empty")
        classes[out.pop("_transition")] = CLASSES.index("transition")
        classes[out.pop("_cell")] = CLASSES.index("cell")
        diagnostics = out.pop("_diagnostics")
        notes = out.pop("_notes")
        return result(
            type(self),
            variant,
            params,
            droplets,
            classes == CLASSES.index("cell"),
            classes=pd.Categorical.from_codes(classes, categories=CLASSES),
            thresholds={k: float(v) for k, v in out.items()},
            diagnostics=diagnostics,
            notes=notes,
        )


def _integer(totals: np.ndarray) -> np.ndarray:
    """CellBender's totals: int64 when every total is whole (as from a 10x matrix)."""
    if np.all(np.isfinite(totals)) and np.all(totals == np.floor(totals)) and (totals.size == 0 or totals.max() < 2**53):
        return totals.astype(np.int64)
    return np.asarray(totals, dtype=np.float64)


def _descending_order(values: np.ndarray, top: int) -> np.ndarray:
    """`np.argsort(values)[::-1][:top]` with numpy's classic quicksort's ties."""
    size = values.size
    top = int(min(max(top, 0), size))
    if top == 0:
        return np.zeros(0, dtype=np.int64)
    return legacy_argsort(values, window=(size - top, size))[::-1][:top]


def _classes(values: np.ndarray, n_cells: int, n_analysed: int, low_cut, with_empties: bool) -> dict:
    """Probable cells, the other analysed droplets, and the surely-empty droplets: ranks
    [0, n_cells), [n_cells, n_analysed), and the rest above `low_cut` (when `with_empties`)."""
    order = _descending_order(values, n_analysed)
    analysed = np.zeros(values.size, dtype=bool)
    analysed[order] = True
    empty = np.flatnonzero((values > low_cut) & ~analysed) if with_empties else np.zeros(0, dtype=np.int64)
    return {"_cell": order[:n_cells], "_transition": order[n_cells:n_analysed], "_empty": empty,
            "_class_ranks": np.array([n_cells, order.size, order.size + empty.size], dtype=np.int64)}


def _priors_v03(matrix, u: np.ndarray, params) -> dict:
    """`SingleCellRNACountsDataset.__init__` of 0.3.0: priors, user values, features, droplets."""
    low = params.low_count_threshold
    priors = heuristic_priors(u, low)
    cell_count_low_limit = priors.pop("_cell_count_low_limit")
    grid, density, rounds = priors.pop("_kde")
    retries = priors.pop("_retries")
    if params.expected_cells is not None:
        priors["expected_cells"] = params.expected_cells
        priors["cell_counts"] = cell_count_given_expected_cells(u, params.expected_cells)
        if params.expected_cells + NUM_EMPTIES_INCREMENT > priors["total_droplets"]:
            priors["total_droplets"] = params.expected_cells + NUM_EMPTIES_INCREMENT
            if params.total_droplets_included is None:
                priors["empty_counts"], priors["empty_count_upper_limit"] = \
                    empty_count_given_expected_cells_and_total_droplets(u, priors["expected_cells"], priors["total_droplets"])
    if params.total_droplets_included is not None:
        priors["total_droplets"] = params.total_droplets_included
        priors["empty_counts"], priors["empty_count_upper_limit"] = \
            empty_count_given_expected_cells_and_total_droplets(u, int(priors["expected_cells"]), params.total_droplets_included)
    if params.force_cell_umi_prior is not None:
        priors["cell_counts"] = params.force_cell_umi_prior
    if params.force_empty_umi_prior is not None:
        priors["empty_counts"] = params.force_empty_umi_prior
        middle = np.sqrt(priors["cell_counts"] * params.force_empty_umi_prior)
        priors["empty_count_upper_limit"] = min(middle, 2 * params.force_empty_umi_prior)
    crossover_surely_empty_and_stds(u, priors)

    keep = kept_features(matrix, u, priors, low, params.projected_ambient_count_threshold)
    trimmed = _integer(np.asarray(matrix @ keep.astype(np.float64)).ravel())

    n_cells = int(priors["expected_cells"])
    total = int(priors["total_droplets"])
    if total - n_cells <= 0:
        raise CellBenderError(f"The number of cells is {n_cells}, but the number of total droplets included is {total}. "
                              f"Increase --total-droplets-included above {n_cells}, or specify a different number of "
                              "expected cells using --expected-cells")
    low_cut = max(low, int(priors["empty_counts"] * LOW_CUT_FACTOR["v0.3"]))
    n_above = int(np.sum(trimmed > low_cut))
    if n_above == 0:
        raise CellBenderError(f"There are no barcodes with UMI counts over the lower cutoff of {low_cut}")
    if n_above <= n_cells:
        raise CellBenderError(f"There are no empty droplets with UMI counts over the lower cutoff of {low_cut}. "
                              "Reduce the --low-count-threshold parameter.")
    if total >= trimmed.size:
        raise CellBenderError(f"total droplets {total} is not below the {trimmed.size} barcodes")
    out = _classes(trimmed, n_cells, total, low_cut, with_empties=True)
    empty_umi_threshold = np.sort(trimmed)[::-1][total]
    cell_prob = (n_cells / total) * (1.0 - params.empty_drop_training_fraction)
    notes = []
    if low_cut == low and out["_empty"].size == 0:
        notes.append("few empty droplets identified: CellBender suggests a lower --low-count-threshold")
    above = u[u > low]
    hist_y, hist_edges = np.histogram(np.log(above), bins=100) if above.size else (np.zeros(0), np.zeros(1))
    out.update({
        "cell_counts": priors["cell_counts"],
        "empty_counts": priors["empty_counts"],
        "expected_cells": n_cells,
        "total_droplets": total,
        "low_count_cutoff": low_cut,
        "empty_umi_threshold": empty_umi_threshold,
        "log_counts_crossover": priors["log_counts_crossover"],
        "d_std": priors["d_std"],
        "cell_logit": np.log(cell_prob) - np.log(1.0 - cell_prob),
        "empty_count_upper_limit": priors["empty_count_upper_limit"],
        "surely_empty_counts": priors["surely_empty_counts"],
        "cell_count_low_limit": cell_count_low_limit,
        "n_features": int(keep.sum()),
        "_notes": notes,
        "_diagnostics": {
            "prior_rule": "otsu-kde",
            "class_ranks": out.pop("_class_ranks"),
            "log_histogram_x": (hist_edges[:-1] + hist_edges[1:]) / 2,
            "log_histogram_y": hist_y.astype(np.int64),
            "kde_x": np.asarray(grid, dtype=float),
            "kde_density": np.asarray(density, dtype=float),
            "kde_rounds": int(rounds),
            "prior_retries": int(retries),
        },
    })
    return out


def _priors_v01(u: np.ndarray, params, variant: str) -> dict:
    """0.1.0 and 0.2.x: `get_d_priors_from_dataset`, `estimate_cell_count_from_dataset`,
    `_trim_dataset_for_analysis` and `_estimate_priors`."""
    low = params.low_count_threshold
    expected = params.expected_cells
    total = params.total_droplets_included
    above = u[u > low]
    if above.size == 0:
        raise CellBenderError(f"no barcode has more than {low} counts (--low-count-threshold)")
    # scipy.stats.mode: the most frequent rounded value, the smallest if several
    values, frequency = np.unique(np.round(np.log1p(above), decimals=1), return_counts=True)
    empty_counts = int(np.expm1(values[np.argmax(frequency)]).item())
    if variant == "v0.2" and expected is not None:
        cell_counts = int(np.median(np.sort(u)[::-1][:expected]).item())
    else:
        big = u[u > 5 * empty_counts]
        if big.size == 0:
            raise CellBenderError(f"no barcode has more than five times the empty prior ({5 * empty_counts} counts)")
        cell_counts = int(np.expm1(np.median(np.log1p(big))).item())
    if expected is None:
        midway = np.mean([np.log1p(cell_counts), np.log1p(empty_counts)])
        n_prior = int(np.sum(u > np.expm1(midway)).item())
    else:
        n_prior = expected
    n_cells = min(n_prior, int(np.sum(u > 0).item()))
    if n_cells <= 0:
        raise CellBenderError("No cells identified. Try to use the option --expected-cells to specify a prior on cell count.")
    low_cut = max(low, int(empty_counts * LOW_CUT_FACTOR[variant]))
    n_above = int(np.sum(u > low_cut).item())
    if n_above == 0:
        raise CellBenderError(f"There are no barcodes with UMI counts over the lower cutoff of {low_cut}")
    if n_above <= n_cells:
        raise CellBenderError(f"There are no empty droplets with UMI counts over the lower cutoff of {low_cut}. "
                              "Reduce the --low-count-threshold parameter.")
    wanted = total - n_cells
    if wanted <= 0:
        raise CellBenderError(f"The number of cells is {n_cells}, but the number of total droplets included is {total}. "
                              f"Increase --total-droplets-included above {n_cells}, or specify a different number of "
                              "expected cells using --expected-cells")
    n_transition = max(0, min(wanted, n_above - n_cells))
    if n_transition == 0:
        raise CellBenderError(f"There are no barcodes identified from the transition region between cell and empty. "
                              f"The low UMI count cutoff, {low_cut}, was likely too high. Try to reduce --low-count-threshold")
    analysed = n_cells + n_transition
    out = _classes(u, n_cells, analysed, low_cut, with_empties=n_transition == wanted)
    crossover = np.mean(np.log1p([cell_counts, empty_counts])).item()
    cell_prob = (1 - params.empty_drop_training_fraction) * (n_prior / analysed)
    if not 0 < cell_prob <= 1:
        raise CellBenderError(f"Fraction of trimmed dataset containing cells should be in (0, 1], but is {cell_prob}.")
    with np.errstate(divide="ignore"):
        cell_logit = np.log(cell_prob / (1 - cell_prob)).item()
    notes = []
    if low_cut == low and out["_empty"].size == 0:
        notes.append("few empty droplets identified: CellBender suggests a lower --low-count-threshold")
    out.update({
        "cell_counts": cell_counts,
        "empty_counts": empty_counts,
        "expected_cells": n_cells,
        "total_droplets": total,
        "low_count_cutoff": low_cut,
        "log_counts_crossover": crossover,
        "d_std": (np.log1p(cell_counts) - crossover) / 5,
        "cell_logit": cell_logit,
        "_notes": notes,
        "_diagnostics": {
            "prior_rule": "mode-median",
            "class_ranks": out.pop("_class_ranks"),
            "mode_values": values.astype(float),
            "mode_counts": frequency.astype(np.int64),
        },
    })
    if variant == "v0.2":
        if analysed >= u.size:
            raise CellBenderError(f"every barcode is analysed ({analysed}): no surely-empty droplet is left")
        out["empty_umi_threshold"] = np.sort(u)[::-1][analysed]
    return out
