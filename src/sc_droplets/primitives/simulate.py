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

import math
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
    pre-drawn features from a block of `block_size`. The feature stream is
    continuous; what differs between releases is when the block is refilled,
    and so where the random stream goes:

    - one feature at a time (3.x, and 7.2+ for gaps below `batch_small`):
      the moment the last feature of the block is used, a whole new block is
      drawn;
    - a batched gap (7.2+, `batch_small` <= gap < `jump`): if the gap would
      reach or pass the end of the block, the block becomes its unused tail
      followed by as many fresh features as were used, before the gap is
      taken.

    The log-likelihood is accumulated feature by feature either way; 7.2+
    sums a batched gap with log-gamma terms, which is the same quantity to the
    last bits of a double.
    """
    totals = np.asarray(totals, dtype=np.int64)
    distinct = np.flatnonzero(np.bincount(totals))
    steps = np.diff(distinct)
    n_profile = len(profile)
    adjusted = scipy_adjusted_profile(profile)
    log_p = np.log(profile)
    simulated = np.empty((len(distinct), num_sims))

    def fresh(size: int) -> np.ndarray:
        return state.choice(n_profile, size=size, p=profile, replace=True)

    block = fresh(block_size)
    position = 0

    def take(run_steps: np.ndarray) -> np.ndarray:
        """The features a run of small gaps uses, refilling the block as Cell Ranger does."""
        nonlocal block, position
        total = int(run_steps.sum())
        if position + total < block_size:
            # The common case: no refill can happen within this run.
            out = block[position : position + total]
            position += total
            return out
        pieces = []
        for step in run_steps.tolist():
            if batch_small is not None and step >= batch_small:
                if position + step >= block_size:
                    block = np.concatenate((block[position:], fresh(position)))
                    position = 0
                pieces.append(block[position : position + step])
                position += step
                continue
            remaining = step
            while remaining:
                available = min(remaining, block_size - position)
                pieces.append(block[position : position + available])
                position += available
                remaining -= available
                if position >= block_size:
                    block = fresh(block_size)
                    position = 0
        return np.concatenate(pieces) if pieces else np.zeros(0, dtype=np.int64)

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
            features = take(run_steps)
            n_values = np.arange(distinct[i - 1] + 1, distinct[j - 1] + 1, dtype=np.float64)
            after = occurrence_rank(features, offset=counts)
            terms = log_p[features] + np.log(n_values / after)
            path = np.cumsum(np.concatenate(([loglik], terms)))[1:]
            ends = np.cumsum(run_steps) - 1
            simulated[i:j, sim] = path[ends]
            loglik = float(path[-1])
            np.add.at(counts, features, 1)
            i = j
    return _ranks_lower(distinct, simulated, totals, observed)


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


# ---------------------------------------------------------------------------
# DropletUtils
# ---------------------------------------------------------------------------


def _first_occurrence_rank(draws: np.ndarray, size: int) -> np.ndarray:
    """`occurrence_rank` with numpy's radix sort when the values fit 16 bits."""
    if size <= 1 << 16:
        return occurrence_rank(draws.astype(np.uint16))
    return occurrence_rank(draws)


def dropletutils(
    ambient: np.ndarray,
    totals: np.ndarray,
    observed: np.ndarray,
    niters: int,
    seed: int,
    *,
    alpha: float = math.inf,
    sampler: str = "alias",
    chunk: int = 1000,
    progress: Optional[Callable[[int], None]] = None,
) -> np.ndarray:
    """DropletUtils' `montecarlo_pval`: per tested droplet, how many simulations are <= it.

    `seed` is what R's `set.seed` was given; `generateSeedVectors` draws one
    64-bit seed per iteration from it, and iteration i (1-based) uses
    `pcg32(seed_i, i)`. With finite `alpha` each iteration first draws a
    Dirichlet profile, one `gamma(alpha·ambient_g, 1)` per feature in order.
    A path of max(totals) features is then drawn, from Boost's
    `discrete_distribution` (`sampler="alias"`, 1.6 on) or by inverting the
    cumulative profile (`sampler="inverse"`, 1.4), and its log-probability is
    accumulated draw by draw: + log a_g (or log(alpha·a_g + c_g)), then
    - log(c_g + 1), c_g being how often g was drawn before. The count for a
    droplet with total T is the number of iterations whose path at T is <= its
    observed value (C++ `lower_bound` over the sorted observed values).

    `sampler="r"` is 1.0-1.2: one R stream for every iteration, u = unif_rand()
    times the profile's sum, inverted through the cumulative profile.
    """
    from sc_droplets.primitives import boostrandom as br
    from sc_droplets.primitives.rmath import clog
    from sc_droplets.primitives.rrandom import RMersenneTwister, generate_seed_vectors

    ambient = np.asarray(ambient, dtype=np.float64)
    totals = np.asarray(totals, dtype=np.int64)
    observed = np.asarray(observed, dtype=np.float64)
    n_above = np.zeros(totals.size, dtype=np.int64)
    if totals.size == 0 or ambient.size == 0 or niters <= 0:
        return n_above
    distinct, where = np.unique(totals, return_inverse=True)
    checkpoints = 2 * distinct - 1  # the path value after T draws sits at 2T - 1
    longest = int(distinct[-1])
    size = ambient.size
    use_alpha = math.isfinite(alpha)
    # std::log, as the C library computes it (numpy's own log can differ in the last bit)
    log_ambient = clog(ambient)  # a zero probability's log is -inf, as in C
    log_count = np.array([0.0] + [math.log(k) for k in range(1, longest + 1)])
    scaled = ambient * alpha if use_alpha else None
    rng = RMersenneTwister(seed)

    def path_counts(draws: np.ndarray) -> np.ndarray:
        rank = _first_occurrence_rank(draws, size)
        terms = np.empty(2 * draws.size)
        if use_alpha:
            terms[0::2] = np.log(scaled[draws] + (rank - 1))
        else:
            terms[0::2] = log_ambient[draws]
        terms[1::2] = -log_count[rank]
        path = np.cumsum(terms)
        return (path[checkpoints][where] <= observed).astype(np.int64)

    if sampler == "r":
        if use_alpha:
            raise NotImplementedError("DropletUtils 1.2's Dirichlet-multinomial draws (R::rgamma) are not reproduced")
        cumulative = np.cumsum(ambient)  # std::partial_sum
        total = cumulative[-1]
        for it in range(niters):
            u = rng.unif_rand(longest) * total
            draws = np.searchsorted(cumulative, u, side="left")
            np.minimum(draws, size - 1, out=draws)
            n_above += path_counts(draws)
            if progress is not None and (it + 1) % chunk == 0:
                progress(it + 1)
        return n_above

    seeds = generate_seed_vectors(rng, niters)
    states, incs = br.pcg32_start(seeds, np.arange(1, niters + 1, dtype=np.uint64))
    jumps = br.Pcg32Jumps(2 * longest + 2)
    if not use_alpha:
        if sampler == "alias":
            threshold, alias_of = br.alias_table(ambient)
        else:
            cumulative = np.cumsum(ambient)
    for start in range(0, niters, chunk):
        stop = min(niters, start + chunk)
        if use_alpha:
            engine = br.Lockstep(states[start:stop], incs[start:stop])
            profiles = br.gamma_profile(engine, scaled)
            after_gamma = engine.state
        for row, it in enumerate(range(start, stop)):
            state, inc = (after_gamma[row] if use_alpha else states[it]), incs[it]
            if use_alpha:
                if sampler == "alias":
                    threshold, alias_of = br.alias_table_fast(profiles[row])
                else:
                    cumulative = np.cumsum(profiles[row])
            if sampler == "alias":
                draws, _ = br.sample_alias(jumps, state, inc, longest, threshold, alias_of)
            else:
                draws, _ = br.sample_inverse(jumps, state, inc, longest, cumulative)
            n_above += path_counts(draws)
        if progress is not None:
            progress(stop)
    return n_above
