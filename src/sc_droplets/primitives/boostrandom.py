"""PCG32 and Boost.Random, as DropletUtils' Monte Carlo draws them.

From 1.4 on, `emptyDrops` gives every iteration its own `pcg32` (O'Neill's
`setseq_xsh_rr_64_32`) seeded with a 64-bit seed from R and the iteration
number as its stream, and draws from it with Boost.Random:

- `pcg32(seed, stream)`: increment `inc = 2·stream + 1`, state
  `(seed + inc)·M + inc` (mod 2**64) with M = 6364136223846793005; each draw
  returns the XSH-RR output of the current state, then steps it.
- `uniform_int_distribution(0, n - 1)` with a 32-bit engine: a bucket of
  `(2**32 - 1) // n` values (one more if `(2**32 - 1) % n == n - 1`), the
  raw value divided by it, and a redraw while the quotient is n or more.
  With n = 1 it returns 0 without drawing.
- `uniform_01<double>`: raw · 2**-32.
- `uniform_real_distribution(0, s)` (1.4): raw / 2**32 · s.
- `discrete_distribution<int, double>`: Walker's alias table, built as
  Boost builds it (weights over their mean, then the "below" and "above"
  lists in index order, each "above" weight absorbing deficits until it falls
  below 1 and is processed next); a draw takes `uniform_int` for the column,
  then `uniform_01` for the test against its threshold.
- `gamma_distribution(a, 1)`: for a < 1, Ahrens and Dieter's GS algorithm
  with an exponential draw; for a > 1, the Cauchy-envelope rejection; a = 1
  is one exponential draw. Exponential draws use Boost's 256-layer ziggurat,
  whose integer and float parts come from two engine outputs.

Written from Boost 1.90's `boost/random` headers (Boost Software License
1.0) and `pcg_random.hpp` as shipped by dqrng (Apache-2.0 or MIT, Copyright
2014-2017 Melissa O'Neill). Checked against DropletUtils' own
`montecarlo_pval` (tests/test_primitives_dropletutils.py).
"""

from __future__ import annotations

import math
from typing import Optional, Tuple

import numpy as np

from sc_droplets.primitives._boost_tables import TABLE_X, TABLE_Y

MULTIPLIER = np.uint64(6364136223846793005)
_U32 = np.uint64(0xFFFFFFFF)
TWO_M32 = 1.0 / 4294967296.0
_TX = np.array(TABLE_X, dtype=np.float64)
_TY = np.array(TABLE_Y, dtype=np.float64)
_PI = 3.14159265358979323846
_E = math.exp(1.0)


# ---------------------------------------------------------------------------
# PCG32
# ---------------------------------------------------------------------------


def pcg32_start(seeds: np.ndarray, streams: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """The state and increment of `pcg32(seed, stream)` for each pair."""
    seeds = np.atleast_1d(np.asarray(seeds, dtype=np.uint64))
    streams = np.atleast_1d(np.asarray(streams, dtype=np.uint64))
    with np.errstate(over="ignore"):
        inc = (streams << np.uint64(1)) | np.uint64(1)
        state = (seeds + inc) * MULTIPLIER + inc
    return state, inc


def pcg32_output(state: np.ndarray) -> np.ndarray:
    """XSH-RR output (as uint64 holding 32 bits) of each state."""
    xorshifted = (((state >> np.uint64(18)) ^ state) >> np.uint64(27)) & _U32
    rot = state >> np.uint64(59)
    left = (np.uint64(32) - rot) & np.uint64(31)
    return ((xorshifted >> rot) | (xorshifted << left)) & _U32


class Pcg32Jumps:
    """Jump-ahead coefficients: after k steps, state = A[k]·state + B[k]·inc."""

    def __init__(self, size: int = 1 << 12):
        self.a = np.ones(1, dtype=np.uint64)
        self.b = np.zeros(1, dtype=np.uint64)
        self._grow(size)

    def _grow(self, size: int) -> None:
        if size < self.a.size:
            return
        n = max(size + 1, 2 * self.a.size)
        with np.errstate(over="ignore"):
            powers = np.full(n, MULTIPLIER, dtype=np.uint64)
            powers[0] = 1
            a = np.multiply.accumulate(powers)
            b = np.empty(n, dtype=np.uint64)
            b[0] = 0
            b[1:] = np.add.accumulate(a[:-1])
        self.a, self.b = a, b

    def states(self, state: np.uint64, inc: np.uint64, n: int) -> np.ndarray:
        """The n states starting from `state` (the first is `state` itself)."""
        self._grow(n)
        with np.errstate(over="ignore"):
            return self.a[:n] * state + self.b[:n] * inc

    def raw(self, state: np.uint64, inc: np.uint64, n: int) -> np.ndarray:
        """The next n outputs of one stream."""
        return pcg32_output(self.states(state, inc, n))

    def advance(self, state: np.uint64, inc: np.uint64, n: int) -> np.uint64:
        """The state after n draws."""
        self._grow(n)
        with np.errstate(over="ignore"):
            return self.a[n] * state + self.b[n] * inc


class Lockstep:
    """Many PCG32 streams stepped together; each call draws from the rows asked for."""

    def __init__(self, state: np.ndarray, inc: np.ndarray):
        self.state = np.array(state, dtype=np.uint64)
        self.inc = np.array(inc, dtype=np.uint64)

    def draw(self, rows: np.ndarray) -> np.ndarray:
        current = self.state[rows]
        with np.errstate(over="ignore"):
            self.state[rows] = current * MULTIPLIER + self.inc[rows]
        return pcg32_output(current)

    def uniform01(self, rows: np.ndarray) -> np.ndarray:
        return self.draw(rows).astype(np.float64) * TWO_M32


# ---------------------------------------------------------------------------
# uniform_int and discrete_distribution
# ---------------------------------------------------------------------------


def int_bucket(n: int) -> int:
    """Bucket size of `uniform_int_distribution(0, n - 1)` with a 32-bit engine (n >= 2)."""
    brange = 0xFFFFFFFF
    bucket = brange // n
    if brange % n == n - 1:
        bucket += 1
    return bucket


def alias_table(weights: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Boost's alias table for `discrete_distribution<int, double>`, built step by step.

    Returns (threshold, alias): column r is kept when the test is below
    threshold[r], else alias[r] is returned.
    """
    weights = np.asarray(weights, dtype=np.float64)
    size = weights.size
    average = float(np.add.accumulate(weights)[-1]) / size  # std::accumulate, then / size
    normalised = (weights / average).tolist()
    below = [[v, i] for i, v in enumerate(normalised) if v < 1.0]
    above = [[v, i] for i, v in enumerate(normalised) if not v < 1.0]
    threshold = [0.0] * size
    alias = [0] * size
    b = a = 0
    while b < len(below) and a < len(above):
        threshold[below[b][1]] = below[b][0]
        alias[below[b][1]] = above[a][1]
        above[a][0] -= 1.0 - below[b][0]
        if above[a][0] < 1.0:
            below[b] = above[a]
            a += 1
        else:
            b += 1
    for k in range(b, len(below)):
        threshold[below[k][1]] = 1.0
    for k in range(a, len(above)):
        threshold[above[k][1]] = 1.0
    return np.array(threshold), np.array(alias, dtype=np.int64)


def alias_table_fast(weights: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """`alias_table` without the loop, for a new table per iteration.

    The pairing is the same; it follows from cumulative deficits (D) and
    excesses (E): the "above" weight in use while the k-th "below" weight is
    processed is the first j with E_j >= D_(k-1), and weight j drops below 1
    right after the first k with D_k > E_j. Thresholds of the weights that
    drop are computed from these sums, not by Boost's running subtraction, so
    they can differ from Boost's in the last bits; a draw only changes if a
    test value falls between the two, which is a chance of about 1e-6 per run.
    """
    weights = np.asarray(weights, dtype=np.float64)
    size = weights.size
    average = np.add.accumulate(weights)[-1] / size
    w = weights / average
    is_below = w < 1.0
    below = np.flatnonzero(is_below)
    above = np.flatnonzero(~is_below)
    threshold = np.ones(size)
    alias = np.zeros(size, dtype=np.int64)
    if below.size == 0 or above.size == 0:
        return threshold, alias
    deficit = np.add.accumulate(1.0 - w[below])
    excess = np.add.accumulate(w[above] - 1.0)
    before = np.concatenate(([0.0], deficit[:-1]))
    current = np.searchsorted(excess, before, side="left")
    done = current < above.size
    threshold[below[done]] = w[below[done]]
    alias[below[done]] = above[current[done]]
    # Above weights that fall below 1: those whose capacity the deficits exceed.
    drops = np.flatnonzero(excess < deficit[-1])
    after = np.searchsorted(deficit, excess[drops], side="right")
    previous = np.concatenate(([0.0], excess[:-1]))
    absorbed = deficit[after] - previous[drops]
    remaining = w[above[drops]] - absorbed
    nxt = drops + 1
    paired = nxt < above.size
    threshold[above[drops[paired]]] = remaining[paired]
    alias[above[drops[paired]]] = above[nxt[paired]]
    return threshold, alias


def sample_alias(
    jumps: Pcg32Jumps,
    state: np.uint64,
    inc: np.uint64,
    count: int,
    threshold: np.ndarray,
    alias: np.ndarray,
) -> Tuple[np.ndarray, np.uint64]:
    """`count` draws of `discrete_distribution` from one stream, and the state after them."""
    size = threshold.size
    out = np.empty(count, dtype=np.int64)
    if count == 0:
        return out, state
    if size == 1:
        # uniform_int(0, 0) does not draw; only the test value is drawn.
        out[:] = 0
        return out, jumps.advance(state, inc, count)
    bucket = np.uint64(int_bucket(size))
    filled = 0
    while filled < count:
        need = count - filled
        raw = jumps.raw(state, inc, 2 * need)
        column = raw[0::2] // bucket
        rejected = np.flatnonzero(column >= size)
        take = need if rejected.size == 0 else int(rejected[0])
        col = column[:take].astype(np.int64)
        test = raw[1 : 2 * take : 2].astype(np.float64) * TWO_M32
        out[filled : filled + take] = np.where(test < threshold[col], col, alias[col])
        filled += take
        state = jumps.advance(state, inc, 2 * take + (1 if rejected.size else 0))
    return out, state


def sample_inverse(
    jumps: Pcg32Jumps,
    state: np.uint64,
    inc: np.uint64,
    count: int,
    cumulative: np.ndarray,
) -> Tuple[np.ndarray, np.uint64]:
    """1.4's draws: `uniform_real_distribution(0, total)`, then the first cumulative weight >= it."""
    raw = jumps.raw(state, inc, count)
    u = raw.astype(np.float64) / 4294967296.0 * cumulative[-1]
    chosen = np.searchsorted(cumulative, u, side="left")
    np.minimum(chosen, cumulative.size - 1, out=chosen)
    return chosen, jumps.advance(state, inc, count)


# ---------------------------------------------------------------------------
# Exponential and gamma, one draw per row
# ---------------------------------------------------------------------------


def exponential(engine: Lockstep, rows: np.ndarray) -> np.ndarray:
    """One `exponential_distribution<double>(1)` draw for each row."""
    result = np.empty(rows.size)
    shift = np.zeros(rows.size)
    pending = np.arange(rows.size)
    while pending.size:
        r = rows[pending]
        first = engine.draw(r)
        second = engine.draw(r)
        layer = (first & np.uint64(255)).astype(np.int64)
        # generate_int_float_pair<double, 8>: 24 bits of the first output, 29 of the second.
        unit = (first >> np.uint64(8)).astype(np.float64) * (1.0 / 16777216.0)
        unit = unit + (second & np.uint64((1 << 29) - 1)).astype(np.float64)
        unit = unit * (0.5 / 268435456.0)
        x = unit * _TX[layer]
        inside = x < _TX[layer + 1]
        result[pending[inside]] = shift[pending[inside]] + x[inside]
        base = ~inside & (layer == 0)
        shift[pending[base]] += _TX[1]
        edge = np.flatnonzero(~inside & (layer != 0))
        accepted = np.zeros(pending.size, dtype=bool)
        if edge.size:
            i = layer[edge]
            xe = x[edge]
            y01 = engine.uniform01(r[edge])
            y = _TY[i] + y01 * (_TY[i + 1] - _TY[i])
            above_upper = (_TX[i] - _TX[i + 1]) * y01 - (_TX[i] - xe)
            above_lower = y - (_TY[i + 1] + (_TX[i + 1] - xe) * _TY[i + 1])
            keep = (above_upper < 0) & ((above_lower < 0) | (y < np.exp(-xe)))
            hit = edge[keep]
            result[pending[hit]] = xe[keep] + shift[pending[hit]]
            accepted[hit] = True
        pending = pending[~inside & ~accepted]
    return result


def gamma(engine: Lockstep, shape: float, rows: np.ndarray) -> np.ndarray:
    """One `gamma_distribution<double>(shape, 1)` draw for each row."""
    if shape == 1.0:
        return exponential(engine, rows)
    result = np.empty(rows.size)
    pending = np.arange(rows.size)
    if shape > 1.0:
        root = math.sqrt(2.0 * shape - 1.0)
        while pending.size:
            r = rows[pending]
            y = np.tan(_PI * engine.uniform01(r))
            x = root * y + shape - 1.0
            positive = x > 0.0
            ok = np.zeros(pending.size, dtype=bool)
            idx = np.flatnonzero(positive)
            if idx.size:
                u = engine.uniform01(r[idx])
                yy = y[idx]
                bound = (1.0 + yy * yy) * np.exp((shape - 1.0) * np.log(x[idx] / (shape - 1.0)) - root * yy)
                good = ~(u > bound)
                result[pending[idx[good]]] = x[idx[good]]
                ok[idx[good]] = True
            pending = pending[~ok]
        return result
    p = _E / (shape + _E)
    while pending.size:
        r = rows[pending]
        u = engine.uniform01(r)
        y = exponential(engine, r)
        small = u < p
        x = np.where(small, np.exp(-y / shape), 1.0 + y)
        with np.errstate(over="ignore", divide="ignore"):
            q = np.where(small, p * np.exp(-x), p + (1.0 - p) * np.power(x, shape - 1.0))
        ok = u < q
        result[pending[ok]] = x[ok]
        pending = pending[~ok]
    return result


def gamma_profile(engine: Lockstep, shapes: np.ndarray, rows: Optional[np.ndarray] = None) -> np.ndarray:
    """For each row (stream), one gamma draw per shape, in order: rows × shapes."""
    if rows is None:
        rows = np.arange(engine.state.size)
    out = np.empty((rows.size, shapes.size))
    for g, shape in enumerate(shapes.tolist()):
        out[:, g] = gamma(engine, shape, rows)
    return out
