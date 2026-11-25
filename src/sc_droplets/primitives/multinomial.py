"""Multinomial log-likelihoods of droplets under an ambient profile.

For counts x (summing to n) and profile p,

    log P(x | n, p) = log n! - sum_g log x_g! + sum_g x_g log p_g

computed over each row's stored entries only, so a sparse matrix of many
barcodes costs its number of nonzeros. The tools compute the same value with
different summation orders, which changes the last bits only; a p-value
compares a simulated value with an observed one, so those bits never decide a
call.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from scipy.special import gammaln


def multinomial_loglik_rows(counts: sp.csr_matrix, log_p: np.ndarray) -> np.ndarray:
    """log P(row | its total, p) for every row of a CSR matrix.

    `log_p[g]` is log p_g for each column. Columns where a row has counts must
    have p_g > 0.
    """
    counts = sp.csr_matrix(counts)
    data = counts.data.astype(np.float64)
    rows = np.repeat(np.arange(counts.shape[0]), np.diff(counts.indptr))
    n = np.bincount(rows, weights=data, minlength=counts.shape[0])
    per_entry = data * log_p[counts.indices] - gammaln(data + 1.0)
    return gammaln(n + 1.0) + np.bincount(rows, weights=per_entry, minlength=counts.shape[0])


def scipy_adjusted_profile(p: np.ndarray) -> np.ndarray:
    """The profile as scipy's `multinomial` (up to 1.4) used it: last entry 1 - sum(rest).

    Cell Ranger 3.x called `scipy.stats.multinomial.logpmf` and `.rvs`, which
    replaced the last probability this way before using it.
    """
    adjusted = np.array(p, dtype=np.float64, copy=True)
    adjusted[-1] = 1.0 - adjusted[:-1].sum()
    return adjusted


def multinomial_logpmf_dense(x: np.ndarray, n: int, p: np.ndarray) -> float:
    """scipy's `multinomial.logpmf(x, n, p)` for one dense count vector."""
    x = np.asarray(x, dtype=np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        terms = np.where(x > 0, x * np.log(p), 0.0) - gammaln(x + 1.0)
    return float(gammaln(n + 1.0) + np.sum(terms))
