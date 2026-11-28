"""R's `smooth.spline(x, y, df = ...)` and `predict(fit, deriv = ...)`, operation for operation.

DropletUtils' `barcodeRanks` (1.0 to 1.26) finds the knee as the point of
least signed curvature of a smoothing spline fitted with `df = 20`. The knee
is a data point, so only the order of the curvatures matters, but it moves
when the spline does; this port keeps R's arithmetic in R's order:

- the data: x made unique after rounding (x - mean)/tol, tol = 1e-6·IQR(x);
  x scaled to [0, 1]; `.nknots.smspl(nx)` inner knots at `seq.int(1, nx,
  length.out = nknots)` (truncated indices), each end knot repeated 3 times;
- `sbart.c`: the Gram matrix of second derivatives (`sgram.f`), X'WX and X'Wy
  (`stxwx.f`), `ratio = tr(X'WX)/tr(Sigma)` over columns 3 .. nk-3, and a
  golden-section / parabolic search of spar in (-1.5, 1.5) (tol 1e-4, eps
  2e-8, at most 500 steps) for the spar whose fit has `sum(lev) = df`,
  minimising 3 + (df - sum(lev))**2, with lambda = ratio·16**(6·spar - 2);
- each fit (`sslvrg.f`): the banded system X'WX + lambda·Sigma solved by
  LINPACK's `dpbfa`/`dpbsl`, fitted values by `bvalue`, leverages from the
  band of the inverse (`sinerp.f`);
- the fit R returns is the one of the *last* spar the search evaluated, which
  is not always the spar it reports;
- `predict(fit, deriv = d)` at the fitted x: `bvalue` with d derivatives,
  divided by range**d.

Ported from R 4.5's `src/library/stats/R/smspline.R`, `src/sbart.c`,
`sgram.f`, `stxwx.f`, `sslvrg.f`, `sinerp.f`, `bvalue.f`, `bsplvd.f` and
LINPACK's `dpbfa.f`/`dpbsl.f` (GPL-2 or later; Copyright (C) 1995-2024 The R
Core Team; the Fortran from Hastie and Tibshirani's GAMFIT and de Boor's
"A Practical Guide to Splines").
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from decimal import Decimal, localcontext
from typing import List, Optional, Tuple

import numpy as np

from sc_droplets.primitives.rstats import quantile7, r_mean, seq_length_out


def nknots_smspl(n: int) -> int:
    """`.nknots.smspl(n)`: the number of inner knots for n unique x."""
    if n < 50:
        return n
    a1, a2, a3, a4 = math.log2(50), math.log2(100), math.log2(140), math.log2(200)
    if n < 200:
        value = 2 ** (a1 + (a2 - a1) * (n - 50) / 150)
    elif n < 800:
        value = 2 ** (a2 + (a3 - a2) * (n - 200) / 600)
    elif n < 3200:
        value = 2 ** (a3 + (a4 - a3) * (n - 800) / 2400)
    else:
        value = 200 + (n - 3200) ** 0.2
    return int(math.trunc(value))


def r_pow16(y: float) -> float:
    """16**y as R computes it on x86-64 (`powl`, rounded): correctly rounded."""
    with localcontext() as ctx:
        ctx.prec = 60
        return float(Decimal(16) ** Decimal(y))


# ---------------------------------------------------------------------------
# B-splines (1-based arrays: index 0 unused, as in the Fortran)
# ---------------------------------------------------------------------------


def _interv(xt: List[float], n: int, x: float) -> Tuple[int, int]:
    """R's `interv(xt, n, x, 0, 0, ...)`: (ileft, mflag), 1-based."""
    if x < xt[1]:
        return 0, -1
    if x >= xt[n]:
        return n, 1
    lo, hi = 1, n
    while hi - lo > 1:  # xt[lo] <= x < xt[hi]
        mid = (lo + hi) // 2
        if x >= xt[mid]:
            lo = mid
        else:
            hi = mid
    return lo, 0


def _bsplvd(t: List[float], k: int, x: float, left: int, nderiv: int) -> List[List[float]]:
    """de Boor's `bsplvd`: db[i][m] = (m-1)-th derivative of the i-th B-spline at x."""
    db = [[0.0] * (nderiv + 1) for _ in range(k + 1)]
    deltal = [0.0] * 21
    deltar = [0.0] * 21
    state = {"j": 1}

    def bsplvb(jhigh: int, index: int) -> None:
        if index != 2:
            state["j"] = 1
            db[1][1] = 1.0
            if state["j"] >= jhigh:
                return
        while True:
            j = state["j"]
            jp1 = j + 1
            deltar[j] = t[left + j] - x
            deltal[j] = x - t[left + 1 - j]
            saved = 0.0
            for i in range(1, j + 1):
                term = db[i][1] / (deltar[i] + deltal[jp1 - i])
                db[i][1] = saved + deltar[i] * term
                saved = deltal[jp1 - i] * term
            db[jp1][1] = saved
            state["j"] = jp1
            if not jp1 < jhigh:
                return

    mhigh = max(min(nderiv, k), 1)
    kp1 = k + 1
    bsplvb(kp1 - mhigh, 1)
    if mhigh == 1:
        return db
    ideriv = mhigh
    for _m in range(2, mhigh + 1):
        jp1mid = 1
        for j in range(ideriv, k + 1):
            db[j][ideriv] = db[jp1mid][1]
            jp1mid += 1
        ideriv -= 1
        bsplvb(kp1 - ideriv, 2)
    a = [[0.0] * (k + 1) for _ in range(k + 1)]
    jlow = 1
    for i in range(1, k + 1):
        for j in range(jlow, k + 1):
            a[j][i] = 0.0
        jlow = i
        a[i][i] = 1.0
    for m in range(2, mhigh + 1):
        kp1mm = kp1 - m
        fkp1mm = float(kp1mm)
        il = left
        i = k
        for _ in range(kp1mm):
            factor = fkp1mm / (t[il + kp1mm] - t[il])
            for j in range(1, i + 1):
                a[i][j] = (a[i][j] - a[i - 1][j]) * factor
            il -= 1
            i -= 1
        for i in range(1, k + 1):
            total = 0.0
            for j in range(max(i, m), k + 1):
                total = a[j][i] * db[j][m] + total
            db[i][m] = total
    return db


def _bvalue(t: List[float], bcoef: List[float], n: int, k: int, x: float, jderiv: int) -> float:
    """de Boor's `bvalue` as R ships it: the jderiv-th derivative of the spline at x."""
    if jderiv >= k:
        return 0.0
    if x != t[n + 1] or t[n + 1] != t[n + k]:
        i, mflag = _interv(t, n + k, x)
        if mflag != 0:
            return 0.0
    else:
        i = n
    km1 = k - 1
    if km1 <= 0:
        return bcoef[i]
    aj = [0.0] * 21
    dm = [0.0] * 21
    dp = [0.0] * 21
    jcmin = 1
    imk = i - k
    if imk >= 0:
        for j in range(1, km1 + 1):
            dm[j] = x - t[i + 1 - j]
    else:
        jcmin = 1 - imk
        for j in range(1, i + 1):
            dm[j] = x - t[i + 1 - j]
        for j in range(i, km1 + 1):
            aj[k - j] = 0.0
            dm[j] = dm[i]
    jcmax = k
    nmi = n - i
    if nmi >= 0:
        for j in range(1, km1 + 1):
            dp[j] = t[i + j] - x
    else:
        jcmax = k + nmi
        for j in range(1, jcmax + 1):
            dp[j] = t[i + j] - x
        for j in range(jcmax, km1 + 1):
            aj[j + 1] = 0.0
            dp[j] = dp[jcmax]
    for jc in range(jcmin, jcmax + 1):
        aj[jc] = bcoef[imk + jc]
    for j in range(1, jderiv + 1):
        kmj = k - j
        fkmj = float(kmj)
        ilo = kmj
        for jj in range(1, kmj + 1):
            aj[jj] = ((aj[jj + 1] - aj[jj]) / (dm[ilo] + dp[jj])) * fkmj
            ilo -= 1
    if jderiv != km1:
        for j in range(jderiv + 1, km1 + 1):
            kmj = k - j
            ilo = kmj
            for jj in range(1, kmj + 1):
                aj[jj] = (aj[jj + 1] * dm[ilo] + aj[jj] * dp[jj]) / (dm[ilo] + dp[jj])
                ilo -= 1
    return aj[1]


# ---------------------------------------------------------------------------
# The penalised regression
# ---------------------------------------------------------------------------


def _sgram(tb: List[float], nb: int) -> Tuple[List[float], ...]:
    sg = [[0.0] * (nb + 1) for _ in range(4)]
    ileft = 1
    for i in range(1, nb + 1):
        ileft, _ = _interv(tb, nb + 1, tb[i])
        v = _bsplvd(tb, 4, tb[i], ileft, 3)
        yw1 = [0.0] + [v[ii][3] for ii in range(1, 5)]
        v = _bsplvd(tb, 4, tb[i + 1], ileft, 3)
        yw2 = [0.0] + [v[ii][3] - yw1[ii] for ii in range(1, 5)]
        wpt = tb[i + 1] - tb[i]

        def term(ii: int, jj: int) -> float:
            return wpt * (yw1[ii] * yw1[jj] + (yw2[ii] * yw1[jj] + yw2[jj] * yw1[ii]) * 0.5 + yw2[ii] * yw2[jj] * 0.3330)

        if ileft >= 4:
            width, base = 4, ileft - 4
        elif ileft in (1, 2, 3):
            width, base = ileft, 0
        else:
            continue
        for ii in range(1, width + 1):
            for band in range(4):
                jj = ii + band
                if jj <= width:
                    sg[band][base + ii] += term(ii, jj)
    return sg[0], sg[1], sg[2], sg[3]


def _stxwx(x: List[float], z: List[float], w: List[float], k: int, knot: List[float], n: int):
    y = [0.0] * (n + 1)
    hs = [[0.0] * (n + 1) for _ in range(4)]
    eps = 0.1e-9
    for i in range(1, k + 1):
        ileft, mflag = _interv(knot, n + 1, x[i])
        if mflag == 1:
            if x[i] <= knot[ileft] + eps:
                ileft -= 1
            else:
                return y, hs
        v = _bsplvd(knot, 4, x[i], ileft, 1)
        ww = w[i] * w[i]
        b = [0.0] + [v[r][1] for r in range(1, 5)]
        for r in range(1, 5):
            j = ileft - 4 + r
            y[j] = y[j] + ww * z[i] * b[r]
            hs[0][j] = hs[0][j] + ww * (b[r] * b[r])
            for band in range(1, 5 - r):
                hs[band][j] = hs[band][j] + ww * b[r] * b[r + band]
    return y, hs


def _ddot(xs: List[float], ys: List[float]) -> float:
    total = 0.0
    for a, b in zip(xs, ys):
        total = total + a * b
    return total


def _dpbfa(abd: List[List[float]], n: int, m: int) -> int:
    """LINPACK `dpbfa`: Cholesky factor of a banded positive definite matrix, in place."""
    for j in range(1, n + 1):
        s = 0.0
        ik = m + 1
        jk = max(j - m, 1)
        mu = max(m + 2 - j, 1)
        for kk in range(mu, m + 1):
            length = kk - mu
            t = abd[kk][j] - _ddot([abd[ik + r][jk] for r in range(length)], [abd[mu + r][j] for r in range(length)])
            t = t / abd[m + 1][jk]
            abd[kk][j] = t
            s = s + t * t
            ik -= 1
            jk += 1
        s = abd[m + 1][j] - s
        if s <= 0.0:
            return j
        abd[m + 1][j] = math.sqrt(s)
    return 0


def _dpbsl(abd: List[List[float]], n: int, m: int, b: List[float]) -> None:
    """LINPACK `dpbsl`: solve with the factor from `_dpbfa`, in place."""
    for k in range(1, n + 1):
        lm = min(k - 1, m)
        la = m + 1 - lm
        lb = k - lm
        t = _ddot([abd[la + r][k] for r in range(lm)], [b[lb + r] for r in range(lm)])
        b[k] = (b[k] - t) / abd[m + 1][k]
    for kb in range(1, n + 1):
        k = n + 1 - kb
        lm = min(k - 1, m)
        la = m + 1 - lm
        lb = k - lm
        b[k] = b[k] / abd[m + 1][k]
        t = -b[k]
        for r in range(lm):
            b[lb + r] = b[lb + r] + t * abd[la + r][k]


def _sinerp(abd: List[List[float]], nk: int) -> List[List[float]]:
    """`sinerp` (flag 0): the band of the inverse needed for the leverages."""
    p1ip = [[0.0] * (nk + 1) for _ in range(5)]
    wjm3 = [0.0, 0.0, 0.0, 0.0]
    wjm2 = [0.0, 0.0, 0.0]
    wjm1 = [0.0, 0.0]
    c1 = c2 = c3 = 0.0
    for i in range(1, nk + 1):
        j = nk - i + 1
        c0 = 1.0 / abd[4][j]
        if j <= nk - 3:
            c1 = abd[1][j + 3] * c0
            c2 = abd[2][j + 2] * c0
            c3 = abd[3][j + 1] * c0
        elif j == nk - 2:
            c1 = 0.0
            c2 = abd[2][j + 2] * c0
            c3 = abd[3][j + 1] * c0
        elif j == nk - 1:
            c1 = c2 = 0.0
            c3 = abd[3][j + 1] * c0
        else:
            c1 = c2 = c3 = 0.0
        p1ip[1][j] = 0.0 - (c1 * wjm3[1] + c2 * wjm3[2] + c3 * wjm3[3])
        p1ip[2][j] = 0.0 - (c1 * wjm3[2] + c2 * wjm2[1] + c3 * wjm2[2])
        p1ip[3][j] = 0.0 - (c1 * wjm3[3] + c2 * wjm2[2] + c3 * wjm1[1])
        p1ip[4][j] = (
            c0 * c0
            + c1 * c1 * wjm3[1]
            + 2.0 * c1 * c2 * wjm3[2]
            + 2.0 * c1 * c3 * wjm3[3]
            + c2 * c2 * wjm2[1]
            + 2.0 * c2 * c3 * wjm2[2]
            + c3 * c3 * wjm1[1]
        )
        wjm3[1] = wjm2[1]
        wjm3[2] = wjm2[2]
        wjm3[3] = p1ip[2][j]
        wjm2[1] = wjm1[1]
        wjm2[2] = p1ip[3][j]
        wjm1[1] = p1ip[4][j]
    return p1ip


@dataclass
class _Fit:
    coef: List[float]
    fitted: List[float]
    lev: List[float]
    crit: float
    info: int


def _sslvrg(dofoff, x, y, w, ssw, n, knot, nk, icrit, lam, xwy, hs, sg) -> _Fit:
    coef = [0.0] + [xwy[i] for i in range(1, nk + 1)]
    abd = [[0.0] * (nk + 1) for _ in range(5)]
    for i in range(1, nk + 1):
        abd[4][i] = hs[0][i] + lam * sg[0][i]
    for i in range(1, nk):
        abd[3][i + 1] = hs[1][i] + lam * sg[1][i]
    for i in range(1, nk - 1):
        abd[2][i + 2] = hs[2][i] + lam * sg[2][i]
    for i in range(1, nk - 2):
        abd[1][i + 3] = hs[3][i] + lam * sg[3][i]
    info = _dpbfa(abd, nk, 3)
    if info != 0:
        return _Fit(coef, [0.0] * (n + 1), [0.0] * (n + 1), math.nan, info)
    _dpbsl(abd, nk, 3, coef)
    sz = [0.0] * (n + 1)
    for i in range(1, n + 1):
        sz[i] = _bvalue(knot, coef, nk, 4, x[i], 0)
    lev = [0.0] * (n + 1)
    crit = math.nan
    if icrit >= 1:
        p1ip = _sinerp(abd, nk)
        eps = 1e-11
        for i in range(1, n + 1):
            xv = x[i]
            ileft, mflag = _interv(knot, nk + 1, xv)
            if mflag == -1:
                ileft = 4
                xv = knot[4] + eps
            elif mflag == 1:
                ileft = nk
                xv = knot[nk + 1] - eps
            j = ileft - 3
            v = _bsplvd(knot, 4, xv, ileft, 1)
            b0, b1, b2, b3 = v[1][1], v[2][1], v[3][1], v[4][1]
            lev[i] = (
                p1ip[4][j] * (b0 * b0)
                + 2.0 * p1ip[3][j] * b0 * b1
                + 2.0 * p1ip[2][j] * b0 * b2
                + 2.0 * p1ip[1][j] * b0 * b3
                + p1ip[4][j + 1] * (b1 * b1)
                + 2.0 * p1ip[3][j + 1] * b1 * b2
                + 2.0 * p1ip[2][j + 1] * b1 * b3
                + p1ip[4][j + 2] * (b2 * b2)
                + 2.0 * p1ip[3][j + 2] * b2 * b3
                + p1ip[4][j + 3] * (b3 * b3)
            ) * (w[i] * w[i])
        df = 0.0
        if icrit == 1:
            rss = ssw
            sumw = 0.0
            for i in range(1, n + 1):
                resid = (y[i] - sz[i]) * w[i]
                rss = rss + resid * resid
                df = df + lev[i]
                sumw = sumw + w[i] * w[i]
            shrink = 1.0 - (dofoff + 1.0 * df) / sumw
            crit = (rss / sumw) / (shrink * shrink)
        else:
            for i in range(1, n + 1):
                df = df + lev[i]
            gap = dofoff - df
            crit = 3 + gap * gap
    return _Fit(coef, sz, lev, crit, 0)


@dataclass
class SmoothSpline:
    """What `predict` needs from a `smooth.spline` fit, and what the search did."""

    x: np.ndarray  # unique x, original scale
    fitted: np.ndarray
    knot: List[float]  # 1-based
    coef: List[float]  # 1-based
    nk: int
    min: float
    range: float
    spar: float
    df: float
    iterations: int
    criterion: str

    def predict(self, x: Optional[np.ndarray] = None, deriv: int = 0) -> np.ndarray:
        """`predict(fit, x, deriv)$y` for x within the data range."""
        x = self.x if x is None else np.asarray(x, dtype=np.float64)
        xs = (x - self.min) / self.range
        if np.any((xs < 0) | (xs > 1)):
            raise ValueError("predict outside the fitted range is not implemented")
        out = np.array([_bvalue(self.knot, self.coef, self.nk, 4, float(v), deriv) for v in xs])
        if deriv > 0:
            out = out / (self.range if deriv == 1 else self.range * self.range if deriv == 2 else self.range**deriv)
        return out


def smooth_spline(x: np.ndarray, y: np.ndarray, df: Optional[float] = None) -> SmoothSpline:
    """`smooth.spline(x, y, df = df)` with unit weights and R's default controls."""
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    n = x.size
    iqr = quantile7(x, 0.75) - quantile7(x, 0.25)
    tol = 1e-6 * iqr
    if not (math.isfinite(tol) and tol > 0):
        raise ValueError("'tol' must be strictly positive and finite")
    xx = np.rint((x - r_mean(x)) / tol)
    order = np.argsort(x, kind="stable")
    xxs = xx[order]
    nd = np.concatenate(([True], xxs[:-1] < xxs[1:]))
    ux = x[order][nd]
    nx = ux.size
    if nx <= 3:
        raise ValueError("need at least four unique 'x' values")
    if nx == n:
        wbar = np.ones(nx)
        ybar = y[order]
        yssw = 0.0
    else:
        group = np.cumsum(nd) - 1
        ys = y[order]
        wbar = np.bincount(group).astype(np.float64)
        sums = np.array([math.fsum(ys[group == g].tolist()) for g in range(nx)])
        squares = np.array([math.fsum((ys[group == g] ** 2).tolist()) for g in range(nx)])
        ybar = sums / wbar
        yssw = math.fsum((squares - wbar * ybar**2).tolist())
    r_ux = ux[-1] - ux[0]
    xbar = (ux - ux[0]) / r_ux
    nknots = nknots_smspl(nx)
    picks = np.trunc(seq_length_out(1, nx, nknots)).astype(np.int64) - 1
    knot_values = [xbar[0]] * 3 + xbar[picks].tolist() + [xbar[-1]] * 3
    nk = nknots + 2
    icrit, dofoff = 1, 0.0
    if df is not None and 1 < df <= nx:
        icrit, dofoff = 3, float(df)

    knot = [0.0] + [float(v) for v in knot_values]
    xs = [0.0] + xbar.tolist()
    ys = [0.0] + ybar.tolist()
    ws = [0.0] + [math.sqrt(v) if v > 0 else v for v in wbar.tolist()]
    sg = _sgram(knot, nk)
    xwy, hs = _stxwx(xs, ys, ws, nx, knot, nk)
    t1 = t2 = 0.0
    for i in range(3, nk - 2):  # C: i = 2 .. nk - 4, 0-based
        t1 += hs[0][i]
        t2 += sg[0][i]
    ratio = t1 / t2

    last: List[_Fit] = []

    def evaluate(spar: float) -> float:
        lam = ratio * r_pow16(spar * 6.0 - 2.0)
        fit = _sslvrg(dofoff, xs, ys, ws, yssw, nx, knot, nk, icrit, lam, xwy, hs, sg)
        last[:] = [fit]
        return fit.crit

    spar, iterations = _sbart_search(evaluate, lambda: last[0].info)
    fit = last[0]
    return SmoothSpline(
        x=ux,
        fitted=np.array(fit.fitted[1:]),
        knot=knot,
        coef=fit.coef,
        nk=nk,
        min=float(ux[0]),
        range=float(r_ux),
        spar=spar,
        df=float(math.fsum(fit.lev[1:])),
        iterations=iterations,
        criterion="df" if icrit == 3 else "GCV",
    )


def _fsign(x: float, y: float) -> float:
    """R's `fsign`: |x| with the sign of y, where y = -0.0 counts as positive."""
    return abs(x) if y >= 0 else -abs(x)


def _sbart_search(evaluate, error, low: float = -1.5, high: float = 1.5, tol: float = 1e-4, eps: float = 2e-8, maxit: int = 500):
    """The spar search of `sbart.c` (Forsythe, Malcolm and Moler's fmin, R's variant)."""
    gold = 0.381966011250105151795413165634
    big = 1e100
    a, b = low, high
    v = a + gold * (b - a)
    w = x = v
    d = e = 0.0
    fx = evaluate(x)
    fv = fw = fx
    iterations = 0
    while error() == 0:
        xm = (a + b) * 0.5
        tol1 = eps * abs(x) + tol / 3.0
        tol2 = tol1 * 2.0
        iterations += 1
        if abs(x - xm) <= tol2 - (b - a) * 0.5 or iterations > maxit:
            break
        golden = abs(e) <= tol1 or fx >= big or fv >= big or fw >= big
        if not golden:
            r = (x - w) * (fx - fv)
            q = (x - v) * (fx - fw)
            p = (x - v) * q - (x - w) * r
            q = (q - r) * 2.0
            if q > 0.0:
                p = -p
            q = abs(q)
            r = e
            e = d
            if abs(p) >= abs(0.5 * q * r) or q == 0.0 or p <= q * (a - x) or p >= q * (b - x):
                golden = True
            else:
                d = p / q
                u = x + d
                if u - a < tol2 or b - u < tol2:
                    d = _fsign(tol1, xm - x)
        if golden:
            e = (a - x) if x >= xm else (b - x)
            d = gold * e
        u = x + (d if abs(d) >= tol1 else _fsign(tol1, d))
        fu = evaluate(u)
        if not math.isfinite(fu):
            fu = 2.0 * big
        if fu <= fx:
            if u >= x:
                a = x
            else:
                b = x
            v, fv = w, fw
            w, fw = x, fx
            x, fx = u, fu
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
    return x, iterations
