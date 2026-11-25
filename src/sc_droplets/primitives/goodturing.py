"""Simple Good-Turing (Gale & Sampson 1995) as Cell Ranger and STAR compute it.

Both smooth the ambient profile of empty droplets, and both differ from the
textbook in small ways that move the profile and so the p-values:

| | Cell Ranger (`sgt_cellranger`) | STAR (`SimpleGoodTuringStar`) |
|---|---|---|
| Fewest distinct frequencies | 10, else error | 5, else no estimate |
| Turing vs LGT switch | z = 1.65, index by position | z = 1.96, index by frequency |
| Slope check | error if the log-log slope > -1 | none |
| P0 | n of the smallest frequency / N | n_1 / N (0 without singletons) |

`sgt_cellranger` is a port of Cell Ranger 3.0.2's `lib/python/cellranger/sgt.py`
(MIT licence, Copyright (c) 2018 10X Genomics), itself a port of Gale and
Sampson's S code; the same file is unchanged up to Cell Ranger 10.1.
`SimpleGoodTuringStar` is a port of STAR's `source/SimpleGoodTuring/sgt.h`
(STAR: MIT licence, Copyright (c) 2019 Alexander Dobin; sgt.h: Copyright (c)
David Elworthy 2004, free to use with this notice).

Reference: W. A. Gale & G. Sampson (1995), Good-Turing frequency estimation
without tears, Journal of Quantitative Linguistics 2:3, 217-237,
doi:10.1080/09296179508590051.
"""

from __future__ import annotations

import math
from typing import Dict, Tuple

import numpy as np
from scipy import stats


class GoodTuringError(ValueError):
    """The estimator does not apply to these frequencies."""


# ---------------------------------------------------------------------------
# Cell Ranger
# ---------------------------------------------------------------------------


def _averaging_transform(r: np.ndarray, nr: np.ndarray) -> np.ndarray:
    d = np.concatenate((np.ones(1, dtype=int), np.diff(r)))
    dr = np.concatenate((0.5 * (d[1:] + d[0:-1]), np.array((d[-1],), dtype=float)))
    return nr.astype(float) / dr


def simple_good_turing_cellranger(xr: np.ndarray, xnr: np.ndarray, clamp_slope: bool = False) -> Tuple[np.ndarray, float]:
    """Cell Ranger's `simple_good_turing`: adjusted frequencies r* and P0.

    `xr` are the distinct nonzero frequencies, ascending; `xnr` how many items
    have each. A log-log slope above -1 is an error up to Cell Ranger 8.x;
    from 9.0 (`clamp_slope`) a slope of -1 or more is set to -1.
    """
    xr = xr.astype(float)
    xnr = xnr.astype(float)
    xN = np.sum(xr * xnr)
    xnrz = _averaging_transform(xr, xnr)
    slope, _intercept, _, _, _ = stats.linregress(np.log(xr), np.log(xnrz))
    if clamp_slope:
        if slope >= -1:
            slope = -1.0
    elif slope > -1:
        raise GoodTuringError(
            f"the log-log slope is {slope:.3f} > -1; Simple Good-Turing does not apply to these data"
        )
    xrst = xr * np.power(1 + 1 / xr, 1 + slope)
    xrstrel = xrst / xr
    # Turing estimates exist only where the next frequency is r + 1.
    xrtry = xr == np.concatenate((xr[1:] - 1, np.zeros(1)))
    xrstarel = np.zeros(len(xr))
    xrstarel[xrtry] = (xr[xrtry] + 1) / xr[xrtry] * np.concatenate((xnr[1:], np.zeros(1)))[xrtry] / xnr[xrtry]
    tursd = np.ones(len(xr))
    for i in range(len(xr)):
        if xrtry[i]:
            # i is a position, not a frequency (as in the S code Cell Ranger ported)
            tursd[i] = float(i + 2) / xnr[i] * np.sqrt(xnr[i + 1] * (1 + xnr[i + 1] / xnr[i]))
    xrstcmbrel = np.zeros(len(xr))
    useturing = True
    for r in range(len(xr)):
        if not useturing:
            xrstcmbrel[r] = xrstrel[r]
        elif np.abs(xrstrel[r] - xrstarel[r]) * (1 + r) / tursd[r] > 1.65:
            xrstcmbrel[r] = xrstarel[r]
        else:
            useturing = False
            xrstcmbrel[r] = xrstrel[r]
    sumpraw = np.sum(xrstcmbrel * xr * xnr / xN)
    xrstcmbrel = xrstcmbrel * (1 - xnr[0] / xN) / sumpraw
    p0 = xnr[0] / xN
    return xr * xrstcmbrel, float(p0)


def sgt_cellranger(frequencies: np.ndarray, clamp_slope: bool = False) -> Tuple[np.ndarray, float]:
    """Cell Ranger's `sgt_proportions`: a proportion for each item, and P0.

    `frequencies` are the items' counts, all positive. Raises GoodTuringError
    with fewer than 10 distinct counts (Cell Ranger then skips its
    EmptyDrops step).
    """
    frequencies = np.asarray(frequencies, dtype=np.int64)
    if len(frequencies) == 0:
        raise GoodTuringError("no frequencies")
    if np.count_nonzero(frequencies) != len(frequencies):
        raise ValueError("frequencies must be greater than zero")
    freqfreqs = np.bincount(frequencies)
    use_freqs = np.flatnonzero(freqfreqs)
    if len(use_freqs) < 10:
        raise GoodTuringError(f"too few distinct nonzero frequencies ({len(use_freqs)}); Cell Ranger needs 10")
    rstar, p0 = simple_good_turing_cellranger(use_freqs, freqfreqs[use_freqs], clamp_slope=clamp_slope)
    rstar_by_freq: Dict[int, float] = dict(zip(use_freqs.tolist(), rstar.tolist()))
    rstar_sum = np.sum(freqfreqs[use_freqs] * rstar)
    rstar_i = np.fromiter((rstar_by_freq[f] for f in frequencies.tolist()), dtype=float, count=len(frequencies))
    pstar = (1 - p0) * (rstar_i / rstar_sum)
    return pstar, p0


# ---------------------------------------------------------------------------
# STAR
# ---------------------------------------------------------------------------


class SimpleGoodTuringStar:
    """STAR's `SGT<uint32>`: add (frequency, count) pairs, `analyse`, then `estimate`.

    `analyse` returns False with fewer than 5 distinct frequencies, leaving
    every estimate 0 (STAR does not check, and its P0 is then undefined).
    """

    MIN_INPUT = 5

    def __init__(self) -> None:
        self._freq: Dict[int, int] = {}
        self._estimate: Dict[int, float] = {}
        self.p_zero: float = float("nan")

    def add(self, observation: int, frequency: int) -> None:
        self._freq[int(observation)] = self._freq.get(int(observation), 0) + int(frequency)

    def analyse(self) -> bool:
        if len(self._freq) < self.MIN_INPUT:
            return False
        obs = sorted(self._freq)
        freq = [self._freq[o] for o in obs]
        rows = len(obs)
        big_n = sum(o * f for o, f in zip(obs, freq))
        self.p_zero = (freq[0] / float(big_n)) if obs[0] == 1 else 0.0
        log_obs = [0.0] * rows
        log_z = [0.0] * rows
        mean_x = 0.0
        mean_y = 0.0
        prev = 0
        for r in range(rows):
            k = float(2 * obs[r] - prev) if r == rows - 1 else float(obs[r + 1])
            z = 2 * freq[r] / (k - prev)
            log_obs[r] = math.log(float(obs[r]))
            log_z[r] = math.log(z)
            mean_x += log_obs[r]
            mean_y += log_z[r]
            prev = obs[r]
        mean_x /= rows
        mean_y /= rows
        xys = 0.0
        xsquares = 0.0
        for r in range(rows):
            xys += (log_obs[r] - mean_x) * (log_z[r] - mean_y)
            xsquares += (log_obs[r] - mean_x) * (log_obs[r] - mean_x)
        slope = xys / xsquares
        intercept = mean_y - slope * mean_x

        def smoothed(i: int) -> float:
            return math.exp(intercept + slope * math.log(float(i)))

        index = {o: r for r, o in enumerate(obs)}
        r_star = [0.0] * rows
        indiff_seen = False
        for r in range(rows):
            obs1 = obs[r] + 1
            y = obs1 * smoothed(obs1) / smoothed(obs[r])
            nxt = index.get(obs1)
            if nxt is None:
                indiff_seen = True
            elif not indiff_seen:
                next_n = freq[nxt]
                f = freq[r]
                x = (obs1 * next_n) / float(f)
                bound = 1.96 * math.sqrt(float(obs1) * obs1 * next_n / (float(f) * f) * (1 + next_n / float(f)))
                if abs(x - y) <= bound:
                    indiff_seen = True
                else:
                    r_star[r] = x
            if indiff_seen:
                r_star[r] = y
        big_n_prime = 0.0
        for r in range(rows):
            big_n_prime += freq[r] * r_star[r]
        for r in range(rows):
            self._estimate[obs[r]] = (1 - self.p_zero) * r_star[r] / big_n_prime
        return True

    def estimate(self, observation: int) -> float:
        """The probability of one item seen `observation` times (P0 for 0)."""
        if observation == 0:
            return self.p_zero
        return self._estimate.get(int(observation), 0.0)


# ---------------------------------------------------------------------------
# edgeR
# ---------------------------------------------------------------------------


def good_turing_edger(counts: np.ndarray, conf: float = 1.96) -> Tuple[np.ndarray, np.ndarray, float, int]:
    """edgeR's `goodTuring(x, conf)`: (count, proportion, P0, n0).

    `count` holds the distinct positive counts (after `as.integer`, which
    truncates), `proportion` the probability of one item seen that often, P0
    the probability of the unseen, n0 the number of zeros. Port of edgeR's
    `R/goodTuring.R` and `src/good_turing.c` (GPL-2 or later, Copyright (C)
    Aaron Lun and Gordon Smyth), which differs from Gale and Sampson's code
    in two ways: the last Z_r uses 2(r - r_prev) as its denominator, and once
    the smoothed estimate is chosen it is kept for every larger count.
    """
    x = np.trunc(np.asarray(counts, dtype=np.float64)).astype(np.int64)
    values, freqs = np.unique(x, return_counts=True)
    n0 = int(freqs[0]) if values.size and values[0] == 0 else 0
    if n0:
        values, freqs = values[1:], freqs[1:]
    obs = values.tolist()
    freq = freqs.tolist()
    rows = len(obs)
    out = [0.0] * rows
    if rows == 0:
        return values, np.zeros(0), 0.0, n0
    big_n = 0.0
    xys = mean_x = mean_y = xsquares = 0.0
    log_obs = [math.log(o) for o in obs]
    last = rows - 1
    for i in range(rows):
        big_n += float(obs[i] * freq[i])
        prev = 0 if i == 0 else obs[i - 1]
        xx = 2 * (obs[i] - prev) if i == last else obs[i + 1] - prev
        log_o = log_obs[i]
        log_z = math.log(2 * freq[i]) - math.log(xx)
        mean_x += log_o
        mean_y += log_z
        xys += log_o * log_z
        xsquares += log_o * log_o
    mean_x /= rows
    mean_y /= rows
    xys -= mean_x * mean_y * rows
    xsquares -= mean_x * mean_x * rows
    slope = xys / xsquares if xsquares != 0 else float("nan")
    p_zero = 0.0 if obs[0] != 1 else freq[0] / big_n
    big_n_prime = 0.0
    indiff = False
    for i in range(rows):
        next_obs = obs[i] + 1
        y = next_obs * math.exp(slope * (math.log(next_obs) - log_obs[i]))
        if i == last or obs[i + 1] != next_obs:
            indiff = True
        if not indiff:
            xv = float(next_obs) * freq[i + 1] / freq[i]
            if abs(xv - y) <= conf * xv * math.sqrt(1.0 / freq[i + 1] + 1.0 / freq[i]):
                indiff = True
            else:
                out[i] = xv
        if indiff:
            out[i] = y
        big_n_prime += out[i] * freq[i]
    factor = (1.0 - p_zero) / big_n_prime
    return values, np.array(out) * factor, p_zero, n0


def good_turing_proportions_edger(profile: np.ndarray) -> np.ndarray:
    """edgeR's `goodTuringProportions` for one column: each zero gets P0 / n0."""
    profile = np.asarray(profile, dtype=np.float64)
    values, proportion, p_zero, n0 = good_turing_edger(profile)
    out = np.empty(profile.size)
    zero = profile == 0
    out[zero] = p_zero / n0 if n0 else float("nan")
    as_int = np.trunc(profile[~zero]).astype(np.int64)
    out[~zero] = proportion[np.searchsorted(values, as_int)]
    return out
