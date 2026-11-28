"""The 1-D kernel estimators of R's ks package that DropletQC uses.

DropletQC (`identify_empty_drops`) calls `ks::kdde(x, deriv.order = 1)` with every default: a
binned estimate of the density's first derivative on 401 grid points, its bandwidth from
`ks::hpi(x, nstage = 2, binned = TRUE, deriv.order = 1)`, Wand and Jones' two-stage plug-in
selector, whose pilot functionals are themselves binned estimates. This ports those pieces,
operation for operation:

- `R/normal.R`: `dnorm.deriv`, `dnorm.deriv.sum`, `psins.1d`;
- `R/kfe.R`: `kfe.1d`;
- `R/selector.R`: `hpi` (derivatives, two stages);
- `R/binning.R`: `binning`, `linbin.ks`, `symconv.1d`, and `src/ks.c`: `massdist1d`;
- `R/kdde.R`: `kdde`, `kdde.binned`, `kdde.binned.1d`.

ks is GPL-2 | GPL-3, Copyright (C) Tarn Duong; read at 1.13.2 – 1.15.3 (CRAN). Two lines changed
in 1.14.3, and both are kept (`legacy=True` is ks before 1.14.3):

- the kernel's half-width in grid points, `ceiling((4+r)*h*(M-1)/(b-a))` before, `ceiling((4+r)*
  h/delta)` from 1.14.3: the same number, rounded on another path;
- where the counts sit in the zero-padded convolution, `(N+1)/2` before, `round(N/2)+1` from
  1.14.3, which with R's halves-to-even rounding is one cell later when (N-1)/2 is odd.

Arithmetic follows R: `dnorm` as R's C code computes it (its second formula beyond five standard
deviations); `x^y` as R on x86-64 (`powl`, so correctly rounded; `x^2` as `x*x`); `sum` in long
double; `seq(a, b, length = 401)`. The convolution is an FFT in both, numpy's here and R's there,
so estimates agree to the last bits, not in them.
"""

from __future__ import annotations

import math
from decimal import Decimal, localcontext
from typing import NamedTuple

import numpy as np

from sc_droplets.primitives.rmath import ld_sum
from sc_droplets.primitives.rstats import r_mean, seq_length_out

#: ks's default grid size for one dimension (`default.gridsize(1)`, `default.bgridsize(1)`).
GRID_SIZE = 401
_M_1_SQRT_2PI = 0.398942280401432677939946059934
#: R's cut-off beyond which `dnorm` is 0: sqrt(-2 ln2 (DBL_MIN_EXP + 1 - DBL_MANT_DIG)).
_DNORM_ZERO = math.sqrt(-2 * math.log(2) * (-1021 + 1 - 53))


def r_pow(x: float, y: float) -> float:
    """R's `x ^ y` for doubles: `x*x` for y = 2, else `powl` rounded to double (x86-64)."""
    if y == 2.0:
        return x * x
    if x == 1.0 or y == 0.0:
        return 1.0
    if x < 0:
        if y != math.floor(y):
            return math.nan
        magnitude = r_pow(-x, y)
        return -magnitude if int(y) % 2 else magnitude
    if x == 0.0:
        return 0.0 if y > 0 else math.inf
    with localcontext() as ctx:
        ctx.prec = 60
        return float(Decimal(x) ** Decimal(y))


def r_dnorm(x: float, sigma: float) -> float:
    """R's `dnorm(x, 0, sigma)` (nmath/dnorm.c), for finite x and sigma > 0."""
    z = abs(x / sigma)
    if z < 5:
        return _M_1_SQRT_2PI * math.exp(-0.5 * z * z) / sigma
    if z > _DNORM_ZERO:
        return 0.0
    # R_forceint is nearbyint: halves to even, as Python's round.
    x1 = math.ldexp(round(math.ldexp(z, 16)), -16)
    x2 = z - x1
    return _M_1_SQRT_2PI / sigma * (math.exp(-0.5 * x1 * x1) * math.exp((-0.5 * x2 - x1) * x2))


def dnorm_deriv(x: np.ndarray, sigma: float, r: int) -> np.ndarray:
    """`dnorm.deriv(x, mu = 0, sigma, deriv.order = r)`: the r-th derivative of the normal density,
    by the Hermite recursion."""
    x = np.atleast_1d(np.asarray(x, dtype=np.float64))
    sign = r_pow(-1.0, r)
    scale = r_pow(sigma, r)
    out = np.empty(x.size)
    for i, value in enumerate(x.tolist()):
        phi = r_dnorm(value, sigma)
        arg = value / sigma
        hmold0, hmold1, hmnew = 1.0, arg, 1.0
        if r == 1:
            hmnew = hmold1
        for k in range(2, r + 1):
            hmnew = arg * hmold1 - (k - 1) * hmold0
            hmold0, hmold1 = hmold1, hmnew
        out[i] = sign * phi * hmnew / scale
    return out


def psins_1d(r: int, sigma: float) -> float:
    """`psins.1d(r, sigma)`: the normal-scale density functional psi_r."""
    if r % 2:
        return 0.0
    numerator = r_pow(-1.0, r // 2) * math.factorial(r)
    denominator = r_pow(2 * sigma, r + 1) * math.factorial(r // 2) * r_pow(math.pi, 0.5)
    return numerator / denominator


def r_sd(x: np.ndarray) -> float:
    """R's `sd(x)`: the square root of `var`, with R's two-pass mean and long-double sums."""
    x = np.asarray(x, dtype=np.float64)
    mean = r_mean(x)
    deviations = x - mean
    return math.sqrt(ld_sum(deviations * deviations) / (x.size - 1))


class Binned(NamedTuple):
    eval_points: np.ndarray
    counts: np.ndarray


def massdist_1d(x: np.ndarray, a: float, b: float, m: int) -> np.ndarray:
    """ks's C `massdist1d` with unit weights: linear binning onto m points from a to b."""
    est = np.zeros(m)
    delta = (b - a) / (m - 1)
    top = m - 2
    for value in np.asarray(x, dtype=np.float64).tolist():
        if not math.isfinite(value):
            continue
        pos = (value - a) / delta
        ix = math.floor(pos)
        fraction = pos - ix
        if 0 <= ix <= top:
            est[ix] += 1 - fraction
            est[ix + 1] += fraction
        elif ix == -1:
            est[0] += fraction
        elif ix == top + 1:
            est[ix] += 1 - fraction
    return est


def binning(x: np.ndarray, h: float, supp: float, size: int = GRID_SIZE) -> Binned:
    """`binning(x, h = h, supp = supp)` in one dimension: a linear grid from min - supp*h to
    max + supp*h, and the linearly binned counts."""
    x = np.asarray(x, dtype=np.float64)
    a = float(x.min()) - supp * h
    b = float(x.max()) + supp * h
    points = seq_length_out(a, b, size)
    return Binned(points, massdist_1d(x, points[0], points[-1], size))


def symconv_1d(keval: np.ndarray, counts: np.ndarray, legacy: bool) -> np.ndarray:
    """`symconv.1d`: the kernel convolved with the counts by zero-padded FFTs."""
    m, n = counts.size, keval.size
    offset = (n + 1) // 2 if legacy else round(n / 2) + 1  # 1-based, as in R
    p = 2 ** math.ceil(math.log2(m + n))
    kz = np.zeros(p)
    cz = np.zeros(p)
    kz[:n] = keval
    cz[offset - 1:offset - 1 + m] = counts
    values = np.fft.ifft(np.fft.fft(kz) * np.fft.fft(cz)).real
    return values[n - 1:n - 1 + m]


def kdde_binned_1d(h: float, r: int, binned: Binned, legacy: bool) -> np.ndarray:
    """`kdde.binned.1d`: the binned estimate of the r-th density derivative on the grid."""
    points, counts = binned
    n = ld_sum(counts)
    a, b, m = float(points.min()), float(points.max()), points.size
    delta = (b - a) / (m - 1)
    reach = (4 + r) * h * (m - 1) / (b - a) if legacy else (4 + r) * h / delta
    half = min(math.ceil(reach), m - 1)
    grid = np.arange(-(half - 1), half, dtype=np.float64)
    keval = dnorm_deriv(delta * grid, h, r) / n
    return symconv_1d(keval, counts, legacy)


def kfe_1d(x: np.ndarray, g: float, r: int, legacy: bool) -> float:
    """`kfe.1d(x, g, r, inc = 1, binned = TRUE)`: the binned estimate of the functional psi_r."""
    n = np.asarray(x).size
    binned = binning(x, h=g, supp=4 + r)
    # kdde.binned is given H = g^2 and takes h back as sqrt(H).
    est = kdde_binned_1d(math.sqrt(g * g), r, binned, legacy)
    return ld_sum(binned.counts * est * n) / (n * n)


def hpi_deriv(x: np.ndarray, r: int, legacy: bool) -> float:
    """`hpi(x, nstage = 2, binned = TRUE, deriv.order = r)` for r >= 1: Wand and Jones' two-stage
    plug-in bandwidth for the r-th density derivative."""
    x = np.asarray(x, dtype=np.float64)
    n = x.size
    k2r4 = float(dnorm_deriv(0.0, 1.0, 2 * r + 4)[0])
    k2r6 = float(dnorm_deriv(0.0, 1.0, 2 * r + 6)[0])
    m2 = 1.0
    mr = psins_1d(2 * r, 1.0)
    psi2r8 = psins_1d(2 * r + 8, r_sd(x))
    gamse2r6 = r_pow(2 * k2r6 / (-m2 * psi2r8 * n), 1 / (2 * r + 9))
    psi2r6 = kfe_1d(x, gamse2r6, 2 * r + 6, legacy)
    gamse2r4 = r_pow(2 * k2r4 / (-m2 * psi2r6 * n), 1 / (2 * r + 7))
    psi2r4 = kfe_1d(x, gamse2r4, 2 * r + 4, legacy)
    return r_pow((2 * r + 1) * mr / (r_pow(m2, 2) * psi2r4 * n), 1 / (2 * r + 5))


class Derivative(NamedTuple):
    eval_points: np.ndarray
    estimate: np.ndarray
    h: float


def kdde_deriv(x: np.ndarray, r: int, legacy: bool) -> Derivative:
    """`kdde(x, deriv.order = r)` with ks's defaults, for r >= 1 and more than one value: the
    binned estimate on 401 points from min - (3.7 + r)h to max + (3.7 + r)h."""
    x = np.asarray(x, dtype=np.float64)
    if x.size < 2:
        raise ValueError("ks's binned estimate needs more than one value")
    h = hpi_deriv(x, r, legacy)
    # kdde.binned makes H = h^2 and binning takes its width back as sqrt(H).
    binned = binning(x, h=math.sqrt(h * h), supp=3.7 + r)
    return Derivative(binned.eval_points, kdde_binned_1d(h, r, binned, legacy), h)
