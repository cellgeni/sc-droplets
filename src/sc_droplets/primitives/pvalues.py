"""Monte Carlo p-values and Benjamini-Hochberg adjustment, as each tool computes them."""

from __future__ import annotations

import numpy as np


def monte_carlo_pvalues(n_lower: np.ndarray, num_sims: int) -> np.ndarray:
    """(1 + number of simulations with a lower log-likelihood) / (1 + simulations).

    The form Cell Ranger and STAR use, with a strict `<`.
    """
    return (1.0 + np.asarray(n_lower, dtype=np.float64)) / (1.0 + num_sims)


def bh_cellranger(p: np.ndarray) -> np.ndarray:
    """Cell Ranger's `adjust_pvalue_bh` (3.0.2 to 10.1), capped at 1.

    Port of `lib/python/cellranger/analysis/diffexp.py` (Cell Ranger 3.0.2,
    MIT licence, Copyright (c) 2018 10X Genomics).
    """
    p = np.asarray(p, dtype=np.float64)
    if p.size == 0:
        return p.copy()
    descending = np.argsort(p)[::-1]
    scale = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(scale * p[descending]))
    return q[np.argsort(descending)]


def bh_starsolo(p: np.ndarray) -> np.ndarray:
    """STAR's EmptyDrops_CR adjustment: p·n/rank, made non-decreasing, not capped at 1.

    Tied p-values get the same adjusted value whatever order the sort leaves
    them in, so the result does not depend on the sort's tie order.
    """
    p = np.asarray(p, dtype=np.float64)
    n = p.size
    if n == 0:
        return p.copy()
    order = np.argsort(p, kind="stable")
    ranked = p[order] * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = adjusted
    return out
