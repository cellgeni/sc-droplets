"""Monte Carlo engines: log-likelihoods of droplets simulated from the ambient profile.

Each engine reproduces one tool's random-number stream, so that with the
tool's seed the simulated values, and so the p-values, are the tool's own.
Each returns, for every candidate, how many simulations gave a lower
log-likelihood at the candidate's own total (`n_lower`).

- `cellranger_v1`: Cell Ranger 3.0 to 8.x (`simulate_multinomial_loglikelihoods`),
  numpy's legacy global stream: a block of 10**6 pre-drawn features,
  multinomial draws for the first total and for jumps of 1,000 or more.
  Port of Cell Ranger 3.0.2 `lib/python/cellranger/stats.py` (MIT licence,
  Copyright (c) 2015-2018 10X Genomics).
- `cellranger_v9`: Cell Ranger 9.0 onwards: a PCG64 stream drawing every
  feature, written from the specification on the docs page.
- `starsolo`: STAR 2.7.8a onwards (`SoloFeature_emptyDrops_CR.cpp`, MIT
  licence, Copyright (c) 2019 Alexander Dobin): one `std::mt19937` per
  simulation, seeded 19760110·(i+1), drawing from libstdc++'s
  `discrete_distribution`.
"""

from __future__ import annotations

from typing import Callable, Optional

import numpy as np

from sc_droplets.primitives.multinomial import multinomial_logpmf_dense, scipy_adjusted_profile
from sc_droplets.primitives.stdrandom import LibstdcxxDiscrete, mt19937_raw


def occurrence_rank(values: np.ndarray, offset: Optional[np.ndarray] = None) -> np.ndarray:
    """For each position, 1 + how many earlier positions hold the same value.

    With `offset[v]`, that many earlier occurrences are assumed for value v.
    """
    values = np.asarray(values)
    n = values.size
    if n == 0:
        return np.zeros(0, dtype=np.int64)
    order = np.argsort(values, kind="stable")
    ordered = values[order]
    starts = np.flatnonzero(np.concatenate(([True], ordered[1:] != ordered[:-1])))
    run_start = np.repeat(starts, np.diff(np.concatenate((starts, [n]))))
    rank = np.empty(n, dtype=np.int64)
    rank[order] = np.arange(n) - run_start + 1
    if offset is not None:
        rank += offset[values]
    return rank


def _count_lower(simulated_at: np.ndarray, observed: np.ndarray) -> np.ndarray:
    """Per observed value, how many simulated values are strictly lower."""
    ordered = np.sort(simulated_at)
    return np.searchsorted(ordered, observed, side="left")


# ---------------------------------------------------------------------------
# Cell Ranger 3.0 - 8.x
# ---------------------------------------------------------------------------


def cellranger_v1(
    profile: np.ndarray,
    totals: np.ndarray,
    observed: np.ndarray,
    num_sims: int,
    state: np.random.RandomState,
    *,
    jump: int = 1000,
    block_size: int = 1_000_000,
    batch_small: Optional[int] = None,
) -> np.ndarray:
    """Cell Ranger's `simulate_multinomial_loglikelihoods`, then the p-value counts.

    `state` is the legacy numpy stream (Cell Ranger seeds the global one with
    0). For each distinct total in increasing order, a simulation either
    draws the whole gap with one multinomial (gaps of `jump` or more) or adds
    one pre-drawn feature at a time, refilling the block of `block_size`
    features when it runs out.

    `batch_small` (Cell Ranger 7.2+) adds the features of gaps of at least
    that many in one step; the stream is the same, only the arithmetic of the
    log-likelihood differs.
    """
    totals = np.asarray(totals, dtype=np.int64)
    distinct = np.flatnonzero(np.bincount(totals))
    steps = np.diff(distinct)
    n_profile = len(profile)
    adjusted = scipy_adjusted_profile(profile)
    log_p = np.log(profile)
    simulated = np.empty((len(distinct), num_sims))

    block = state.choice(n_profile, size=block_size, p=profile, replace=True)
    position = 0

    def take(n: int) -> np.ndarray:
        nonlocal block, position
        out = np.empty(n, dtype=np.int64)
        filled = 0
        while filled < n:
            available = min(n - filled, block_size - position)
            out[filled : filled + available] = block[position : position + available]
            filled += available
            position += available
            if position >= block_size:
                block = state.choice(n_profile, size=block_size, p=profile, replace=True)
                position = 0
        return out

    small = steps < jump
    for sim in range(num_sims):
        counts = state.multinomial(distinct[0], adjusted).astype(np.int64)
        loglik = multinomial_logpmf_dense(counts, distinct[0], adjusted)
        simulated[0, sim] = loglik
        i = 1
        while i < len(distinct):
            if not small[i - 1]:
                counts += state.multinomial(steps[i - 1], adjusted)
                loglik = multinomial_logpmf_dense(counts, distinct[i], adjusted)
                simulated[i, sim] = loglik
                i += 1
                continue
            j = i
            while j < len(distinct) and small[j - 1]:
                j += 1
            run_steps = steps[i - 1 : j - 1]
            features = take(int(run_steps.sum()))
            n_values = np.arange(distinct[i - 1] + 1, distinct[j - 1] + 1, dtype=np.float64)
            after = occurrence_rank(features, offset=counts)
            if batch_small is None:
                terms = log_p[features] + np.log(n_values / after)
                path = np.cumsum(np.concatenate(([loglik], terms)))[1:]
                ends = np.cumsum(run_steps) - 1
                simulated[i:j, sim] = path[ends]
                loglik = float(path[-1])
            else:
                loglik = _batched_run(loglik, counts, features, distinct, i, j, run_steps, log_p, batch_small, simulated[:, sim])
            np.add.at(counts, features, 1)
            i = j
    return _ranks_lower(distinct, simulated, totals, observed)


def _batched_run(loglik, counts, features, distinct, i, j, run_steps, log_p, batch_small, column) -> float:
    """Cell Ranger 7.2+: gaps of `batch_small` or more are added in one step."""
    from scipy.special import gammaln

    working = counts.copy()
    start = 0
    for k, step in enumerate(run_steps):
        chunk = features[start : start + step]
        n_before = distinct[i - 1 + k]
        if step >= batch_small:
            added = np.bincount(chunk, minlength=working.size)
            touched = np.flatnonzero(added)
            c = added[touched].astype(np.float64)
            x = working[touched].astype(np.float64)
            loglik += float(
                np.sum(c * log_p[touched])
                + gammaln(n_before + step + 1.0)
                - gammaln(n_before + 1.0)
                - np.sum(gammaln(x + c + 1.0) - gammaln(x + 1.0))
            )
            working[touched] += added[touched]
        else:
            for offset, feature in enumerate(chunk):
                working[feature] += 1
                loglik += log_p[feature] + np.log(float(n_before + offset + 1) / working[feature])
        column[i + k] = loglik
        start += step
    return loglik


def _ranks_lower(distinct: np.ndarray, simulated: np.ndarray, totals: np.ndarray, observed: np.ndarray) -> np.ndarray:
    index = np.searchsorted(distinct, totals)
    n_lower = np.zeros(len(totals), dtype=np.int64)
    for row in np.unique(index):
        members = np.flatnonzero(index == row)
        n_lower[members] = _count_lower(simulated[row], observed[members])
    return n_lower


# ---------------------------------------------------------------------------
# Cell Ranger 9.0+
# ---------------------------------------------------------------------------


def cellranger_v9(
    profile: np.ndarray,
    totals: np.ndarray,
    observed: np.ndarray,
    num_sims: int,
    generator: np.random.Generator,
    *,
    order: Optional[np.ndarray] = None,
    chunk: int = 500,
) -> np.ndarray:
    """Cell Ranger 9.0+: every simulation draws `max(totals)` features from `generator`.

    A feature is the first whose cumulative probability is >= a uniform draw.
    The log-likelihood after t draws accumulates log t - log(occurrences of
    the drawn feature so far) + log p. `order` permutes the profile before
    the cumulative sum (10.0+ sorts features by id for this).
    """
    totals = np.asarray(totals, dtype=np.int64)
    observed = np.asarray(observed, dtype=np.float64)
    p = profile if order is None else profile[order]
    cumulative = np.cumsum(p)
    log_p = np.log(p)
    max_n = int(totals.max())
    log_t = np.log(np.arange(1, max_n + 1, dtype=np.float64))
    n_lower = np.zeros(len(totals), dtype=np.int64)
    columns = totals - 1
    done = 0
    while done < num_sims:
        k = min(chunk, num_sims - done)
        uniforms = generator.random((k, max_n))
        draws = np.searchsorted(cumulative, uniforms, side="left")
        # occurrence of each draw within its own simulation
        keyed = draws + (np.arange(k)[:, None] * len(p))
        rank = occurrence_rank(keyed.ravel()).reshape(k, max_n)
        terms = (log_t[None, :] - np.log(rank)) + log_p[draws]
        paths = np.cumsum(terms, axis=1)
        values = paths[:, columns]  # k × candidates
        n_lower += (values < observed[None, :]).sum(axis=0)
        done += k
    return n_lower


# ---------------------------------------------------------------------------
# STAR
# ---------------------------------------------------------------------------


def starsolo(
    profile_nonzero: np.ndarray,
    totals: np.ndarray,
    observed: np.ndarray,
    num_sims: int,
    *,
    seed_base: int = 19760110,
    chunk: int = 250,
    progress: Optional[Callable[[int], None]] = None,
) -> np.ndarray:
    """STAR's EmptyDrops_CR simulations, then the p-value counts.

    Simulation i walks n = 1 .. max(totals), drawing one feature per step
    from `std::mt19937(seed_base·(i+1))`, and accumulates
    ((ll + log p_g) + log n) - log(count_g): the order STAR adds them in.
    """
    totals = np.asarray(totals, dtype=np.int64)
    observed = np.asarray(observed, dtype=np.float64)
    distribution = LibstdcxxDiscrete(profile_nonzero)
    # STAR passes the normalised profile and its logarithms (ambProfileLogPnon0).
    log_p = np.log(np.asarray(profile_nonzero, dtype=np.float64))
    max_n = int(totals.max())
    log_n = np.log(np.arange(1, max_n + 1, dtype=np.float64))
    n_lower = np.zeros(len(totals), dtype=np.int64)
    columns = totals - 1
    done = 0
    while done < num_sims:
        k = min(chunk, num_sims - done)
        draws = np.empty((k, max_n), dtype=np.int64)
        for row in range(k):
            raw = mt19937_raw(seed_base * (done + row + 1), 2 * max_n)
            draws[row] = distribution.sample_from_raw(raw)
        keyed = draws + (np.arange(k)[:, None] * len(profile_nonzero))
        count = occurrence_rank(keyed.ravel()).reshape(k, max_n)
        paths = np.empty((k, max_n))
        acc = np.zeros(k)
        for n in range(max_n):
            acc = acc + log_p[draws[:, n]]
            acc = acc + log_n[n]
            acc = acc - np.log(count[:, n])
            paths[:, n] = acc
        values = paths[:, columns]
        n_lower += (values < observed[None, :]).sum(axis=0)
        done += k
        if progress is not None:
            progress(done)
    return n_lower
