"""Otsu's threshold and its multi-class form, as scikit-image 0.20 and later compute them.

Ported from scikit-image's `skimage/filters/thresholding.py` (`threshold_otsu`,
`threshold_multiotsu`, `_validate_image_histogram`) and `skimage/filters/_multiotsu.pyx`
(`_get_multiotsu_thresh_indices_lut`), BSD-3-Clause:

    Copyright (C) 2019, the scikit-image team. All rights reserved.
    Redistribution and use in source and binary forms, with or without modification, are permitted
    provided that the conditions of the BSD 3-Clause licence are met (see scikit-image's
    LICENSE.txt).

The histogram is numpy's (equal bins from the minimum to the maximum); every sum after it is in
32-bit floats, in scikit-image's order, so that ties between two splits resolve as there.
"""

from __future__ import annotations

from typing import Tuple

import numpy as np


class OtsuError(ValueError):
    """Where scikit-image raises."""


def histogram(values: np.ndarray, n_bins: int) -> Tuple[np.ndarray, np.ndarray]:
    """`skimage.exposure.histogram(values, n_bins, source_range="image")` for floats: numpy's
    histogram over the values' range, and the bin centres."""
    counts, edges = np.histogram(np.asarray(values).reshape(-1), bins=n_bins)
    return counts, (edges[:-1] + edges[1:]) / 2.0


def threshold_otsu(values: np.ndarray, n_bins: int = 256) -> float:
    """`skimage.filters.threshold_otsu(values, nbins)`: the bin centre that maximises the
    variance between the two classes (the first one, when several do). When every value is the
    same, that value."""
    values = np.asarray(values).reshape(-1)
    if values.size == 0:
        raise OtsuError("no values to threshold")
    if np.all(values == values[0]):
        return float(values[0])
    counts, centers = histogram(values, n_bins)
    counts = counts.astype("float32", copy=False)
    weight1 = np.cumsum(counts)
    weight2 = np.cumsum(counts[::-1])[::-1]
    mean1 = np.cumsum(counts * centers) / weight1
    mean2 = (np.cumsum((counts * centers)[::-1]) / weight2[::-1])[::-1]
    gap = mean1[:-1] - mean2[1:]
    variance12 = weight1[:-1] * weight2[1:] * (gap * gap)
    return float(centers[np.argmax(variance12)])


def _between_class_table(prob: np.ndarray) -> np.ndarray:
    """`_set_var_btwcls_lut`, as a full table: H[i, j] for the bins i..j (upper triangle; H[0, 0]
    is left at 0 as there). First moments start at prob[0], not 0 * prob[0], as there."""
    n = prob.size
    zeroth = np.empty(n, dtype=np.float32)
    first = np.empty(n, dtype=np.float32)
    zeroth[0] = prob[0]
    first[0] = prob[0]
    for i in range(1, n):  # sequential float32 sums, as the C loop
        zeroth[i] = zeroth[i - 1] + prob[i]
        first[i] = first[i - 1] + np.float32(i) * prob[i]
    table = np.zeros((n, n), dtype=np.float32)
    with np.errstate(divide="ignore", invalid="ignore"):
        row = np.where(zeroth[1:] > 0, (first[1:] * first[1:]) / zeroth[1:], np.float32(0))
    table[0, 1:] = row
    for i in range(1, n):
        z = zeroth[i:] - zeroth[i - 1]
        f = first[i:] - first[i - 1]
        with np.errstate(divide="ignore", invalid="ignore"):
            table[i, i:] = np.where(z > 0, (f * f) / z, np.float32(0))
    return table


def multiotsu_indices(prob: np.ndarray, n_thresholds: int) -> np.ndarray:
    """`_get_multiotsu_thresh_indices_lut(prob, n_thresholds)`: the split bins whose classes have
    the largest variance between them. The search runs over every increasing choice in
    lexicographic order and keeps the first strict maximum, as the recursive C search does."""
    prob = np.asarray(prob, dtype=np.float32)
    n = prob.size
    table = _between_class_table(prob)
    last_col = table[:, n - 1]
    best = np.float32(0)
    best_indices = None
    k = n_thresholds

    def search(prefix: list, start: int) -> None:
        nonlocal best, best_indices
        depth = len(prefix)
        if depth == k - 1:
            # the last threshold, vectorised: t from start to n - 2
            t = np.arange(start, n - 1)
            if t.size == 0:
                return
            if depth == 0:
                sigma = table[0, t] + last_col[t + 1]
            else:
                sigma = table[0, prefix[0]] + last_col[t + 1]
                for a, b in zip(prefix[:-1], prefix[1:]):
                    sigma = sigma + table[a + 1, b]
                sigma = sigma + table[prefix[-1] + 1, t]
            sigma = sigma.astype(np.float32, copy=False)
            top = sigma.max()
            if top > best:
                best = top
                best_indices = prefix + [int(t[int(np.argmax(sigma))])]
            return
        for idx in range(start, n - k + depth):
            search(prefix + [idx], idx + 1)

    search([], 0)
    if best_indices is None:
        raise OtsuError("no split separates the values: the classes have no variance between them")
    return np.asarray(best_indices, dtype=np.intp)


def threshold_multiotsu(values: np.ndarray, classes: int = 3, n_bins: int = 256) -> np.ndarray:
    """`skimage.filters.threshold_multiotsu(values, classes, nbins)`: classes - 1 thresholds."""
    values = np.asarray(values).reshape(-1)
    if values.size == 0:
        raise OtsuError("no values to threshold")
    counts, centers = histogram(values, n_bins)
    prob = (counts / np.sum(counts)).astype("float32", copy=False)
    nvalues = int(np.count_nonzero(prob))
    if nvalues < classes:
        raise OtsuError(
            f"after binning, the values take only {nvalues} different values: they cannot be split into "
            f"{classes} classes (more bins may help)"
        )
    if nvalues == classes:
        indices = np.flatnonzero(prob)[:-1]
    else:
        indices = multiotsu_indices(prob, classes - 1)
    return centers[indices]
