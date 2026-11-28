"""R's own arithmetic where DropletUtils relies on it.

Ports of R's `stats` and `base` functions (GPL-2 or later, Copyright (C)
1995-2024 The R Core Team): `quantile(type = 7)`, `p.adjust("BH")`,
`optimize` (Brent's `fmin`, `src/library/stats/src/optimize.c`), `mean`,
`all.equal` for numbers, and `seq.int(from, to, length.out)`.
"""

from __future__ import annotations

import math
from typing import Callable

import numpy as np

DOUBLE_EPS = 2.220446049250313e-16


def quantile7(x: np.ndarray, prob: float) -> float:
    """`quantile(x, prob, type = 7, names = FALSE)` for one probability.

    When the two order statistics around the position are equal, R returns
    that value without interpolating (so without (1 - h)·q + h·q's rounding);
    R 3.5.1 to 4.5.3 all do.
    """
    x = np.sort(np.asarray(x, dtype=np.float64))
    n = x.size
    index = 1 + max(n - 1, 0) * prob
    lo = math.floor(index)
    hi = math.ceil(index)
    q = float(x[lo - 1])
    if index > lo and x[hi - 1] != q:
        h = index - lo
        q = (1 - h) * q + h * float(x[hi - 1])
    return q


def p_adjust_bh(p: np.ndarray) -> np.ndarray:
    """`p.adjust(p, "BH")`: NA stays NA and does not count towards n."""
    p = np.asarray(p, dtype=np.float64)
    out = np.full(p.size, np.nan)
    ok = ~np.isnan(p)
    values = p[ok]
    n = values.size
    if n == 0:
        return out
    i = np.arange(n, 0, -1, dtype=np.float64)
    order = np.argsort(-values, kind="stable")  # order(p, decreasing = TRUE)
    adjusted = np.minimum(1.0, np.minimum.accumulate(n / i * values[order]))
    back = np.empty(n)
    back[order] = adjusted
    out[ok] = back
    return out


def brent_fmin(f: Callable[[float], float], ax: float, bx: float, tol: float) -> float:
    """R's `Brent_fmin`: the minimiser of f on (ax, bx), as `optimize` finds it."""
    c = (3.0 - math.sqrt(5.0)) * 0.5
    eps = math.sqrt(DOUBLE_EPS)
    a, b = ax, bx
    v = a + c * (b - a)
    w = x = v
    d = e = 0.0
    fx = f(x)
    fv = fw = fx
    tol3 = tol / 3.0
    while True:
        xm = (a + b) * 0.5
        tol1 = eps * abs(x) + tol3
        t2 = tol1 * 2.0
        if abs(x - xm) <= t2 - (b - a) * 0.5:
            break
        p = q = r = 0.0
        if abs(e) > tol1:
            r = (x - w) * (fx - fv)
            q = (x - v) * (fx - fw)
            p = (x - v) * q - (x - w) * r
            q = (q - r) * 2.0
            if q > 0.0:
                p = -p
            else:
                q = -q
            r = e
            e = d
        if abs(p) >= abs(q * 0.5 * r) or p <= q * (a - x) or p >= q * (b - x):
            e = (b - x) if x < xm else (a - x)
            d = c * e
        else:
            d = p / q
            u = x + d
            if u - a < t2 or b - u < t2:
                d = tol1
                if x >= xm:
                    d = -d
        if abs(d) >= tol1:
            u = x + d
        elif d > 0.0:
            u = x + tol1
        else:
            u = x - tol1
        fu = f(u)
        if fu <= fx:
            if u < x:
                b = x
            else:
                a = x
            v, w, x = w, x, u
            fv, fw, fx = fw, fx, fu
        else:
            if u < x:
                a = u
            else:
                b = u
            if fu <= fw or w == x:
                v, fv = w, fw
                w, fw = u, fu
            elif fu <= fv or v == x or v == w:
                v, fv = u, fu
    return x


def optimize_max(f: Callable[[float], float], lower: float, upper: float) -> float:
    """`optimize(f, c(lower, upper), maximum = TRUE)$maximum` with R's default tolerance."""
    return brent_fmin(lambda t: -f(t), lower, upper, DOUBLE_EPS**0.25)


def r_sum(values: np.ndarray) -> float:
    """R's `sum` of doubles, which accumulates in long double: approximated by an exact sum."""
    return math.fsum(np.asarray(values, dtype=np.float64).tolist())


def r_mean(values: np.ndarray) -> float:
    """R's `mean`: a long-double mean, refined by the mean of the residuals."""
    values = np.asarray(values, dtype=np.float64)
    n = values.size
    s = math.fsum(values.tolist()) / n
    if math.isfinite(s):
        s += math.fsum((values - s).tolist()) / n
    return s


def r_median(values: np.ndarray) -> float:
    """R's `median` of doubles: the middle value, or R's `mean` of the middle two."""
    values = np.sort(np.asarray(values, dtype=np.float64).ravel())
    n = values.size
    if n == 0:
        return math.nan
    half = (n + 1) // 2
    if n % 2 == 1:
        return float(values[half - 1])
    return r_mean(values[half - 1:half + 1])


def r_mad(values: np.ndarray, center: float, constant: float = 1.4826) -> float:
    """R's `mad(x, center)`: `constant * median(abs(x - center))`."""
    return constant * r_median(np.abs(np.asarray(values, dtype=np.float64) - center))


def all_equal(target: np.ndarray, current: np.ndarray, tolerance: float = 1.5e-8) -> bool:
    """`isTRUE(all.equal(target, current))` for two numeric vectors of one length."""
    target = np.asarray(target, dtype=np.float64)
    current = np.asarray(current, dtype=np.float64)
    if target.size == 0:
        return True
    xy = np.mean(np.abs(target - current))
    xn = np.mean(np.abs(target))
    if math.isfinite(xn) and xn > tolerance:
        xy = xy / xn
    return bool(xy <= tolerance)


def seq_length_out(start: float, end: float, length: int) -> np.ndarray:
    """`seq.int(start, end, length.out = length)` as doubles: start + i·by, the last set to `end`.

    Checked against R 3.5.1 to 4.5.3.
    """
    out = np.empty(length)
    if length > 0:
        out[0] = start
    if length > 1:
        out[-1] = end
    if length > 2:
        by = (end - start) / (length - 1)
        for i in range(1, length - 1):
            out[i] = start + i * by
    return out
