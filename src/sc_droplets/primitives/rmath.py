"""R's `lgamma` and R's sums, bit for bit.

DropletUtils' estimate of the Dirichlet-multinomial alpha maximises a
likelihood so flat near its optimum that the last bits of each evaluation
steer R's `optimize`. Reproducing it needs both of these exactly:

- `lgammafn`: R's own log-gamma (nmath `lgamma.c`, `gamma.c`,
  `lgammacor.c`, `chebyshev.c`), not the platform's: for 0 < x <= 10,
  log|gamma(x)| from a Chebyshev series on [1, 2] and the recurrence; above
  10, Stirling's formula with `lgammacor`'s 5-term correction.
- `ld_sum` / `ld_cumsum`: R's `sum` and `cumsum` of doubles accumulate in
  `long double`, which on x86-64 is the x87 80-bit format (64-bit
  mantissa), rounded to double at the end. Where numpy's `longdouble` is that
  format it is used directly; elsewhere (Apple silicon, aarch64 Linux) the
  accumulation is emulated exactly in software. R built for those platforms
  has no extended precision, so this reproduces R on x86-64, the platform of
  Bioconductor's builds and containers.

Ports of R's `src/nmath` (GPL-2 or later; Copyright (C) 1998-2024 Ross Ihaka
and the R Core Team; the Chebyshev series from W. Fullerton's SLATEC code).
"""

from __future__ import annotations

import math
from typing import Tuple

import numpy as np

M_LN_SQRT_2PI = 0.918938533204672741780329736406

_GAMCS = np.array([
    +.8571195590989331421920062399942e-2,
    +.4415381324841006757191315771652e-2,
    +.5685043681599363378632664588789e-1,
    -.4219835396418560501012500186624e-2,
    +.1326808181212460220584006796352e-2,
    -.1893024529798880432523947023886e-3,
    +.3606925327441245256578082217225e-4,
    -.6056761904460864218485548290365e-5,
    +.1055829546302283344731823509093e-5,
    -.1811967365542384048291855891166e-6,
    +.3117724964715322277790254593169e-7,
    -.5354219639019687140874081024347e-8,
    +.9193275519859588946887786825940e-9,
    -.1577941280288339761767423273953e-9,
    +.2707980622934954543266540433089e-10,
    -.4646818653825730144081661058933e-11,
    +.7973350192007419656460767175359e-12,
    -.1368078209830916025799499172309e-12,
    +.2347319486563800657233471771688e-13,
    -.4027432614949066932766570534699e-14,
    +.6910051747372100912138336975257e-15,
    -.1185584500221992907052387126192e-15,
])  # the first ngam = 22 of R's 42 coefficients

_ALGMCS = np.array([
    +.1666389480451863247205729650822e+0,
    -.1384948176067563840732986059135e-4,
    +.9810825646924729426157171547487e-8,
    -.1809129475572494194263306266719e-10,
    +.6221098041892605227126015543416e-13,
])  # nalgm = 5

_XBIG = 94906265.62425156


def _c(function, value: float) -> float:
    """C semantics where Python raises: log(0) is -inf, log of a negative number is NaN."""
    if value > 0:
        return function(value)
    if value == 0:
        return -math.inf
    return math.nan if value < 0 or math.isnan(value) else function(value)


def clog(x: np.ndarray) -> np.ndarray:
    """The C library's `log`, element by element, as R calls it.

    numpy's vectorised log may use its own implementation (AVX-512 on x86-64)
    whose last bit can differ from the C library's.
    """
    x = np.asarray(x, dtype=np.float64)
    values = [_c(math.log, v) for v in x.ravel().tolist()]
    return np.array(values, dtype=np.float64).reshape(x.shape)


def clog10(x: np.ndarray) -> np.ndarray:
    """The C library's `log10`, element by element."""
    x = np.asarray(x, dtype=np.float64)
    values = [_c(math.log10, v) for v in x.ravel().tolist()]
    return np.array(values, dtype=np.float64).reshape(x.shape)


def _chebyshev(x: np.ndarray, coefficients: np.ndarray) -> np.ndarray:
    twox = x * 2
    b0 = np.zeros_like(x)
    b1 = np.zeros_like(x)
    b2 = np.zeros_like(x)
    n = coefficients.size
    for i in range(1, n + 1):
        b2 = b1
        b1 = b0
        b0 = twox * b1 - b2 + coefficients[n - i]
    return (b0 - b2) * 0.5


def _lgammacor(x: np.ndarray) -> np.ndarray:
    out = np.empty_like(x)
    small = x < _XBIG
    tmp = 10 / x[small]
    out[small] = _chebyshev(tmp * tmp * 2 - 1, _ALGMCS) / x[small]
    out[~small] = 1 / (x[~small] * 12)
    return out


def _gamma_small(x: np.ndarray) -> np.ndarray:
    """gammafn(x) for 0 < x <= 10."""
    n = np.trunc(x).astype(np.int64)
    y = x - n
    n = n - 1
    value = _chebyshev(y * 2 - 1, _GAMCS) + 0.9375
    below = n < 0  # 0 < x < 1: one division by x
    value[below] = value[below] / x[below]
    for i in range(1, 10):
        up = n >= i
        value[up] = value[up] * (y[up] + i)
    return value


def lgammafn(x) -> np.ndarray:
    """R's `lgamma(x)` for x > 0."""
    x = np.asarray(x, dtype=np.float64)
    scalar = x.ndim == 0
    x = np.atleast_1d(x)
    if np.any(~(x > 0)):
        raise ValueError("lgammafn here handles x > 0 only")
    out = np.empty_like(x)
    tiny = x < 1e-306
    out[tiny] = -clog(x[tiny])
    small = ~tiny & (x <= 10)
    out[small] = clog(np.abs(_gamma_small(x[small])))
    large = x > 10
    xl = x[large]
    log_xl = clog(xl)
    value = M_LN_SQRT_2PI + (xl - 0.5) * log_xl - xl
    mid = xl <= 4934720.0
    value[mid] = value[mid] + _lgammacor(xl[mid])
    huge = xl > 1e17
    value[huge] = xl[huge] * (log_xl[huge] - 1.0)
    out[large] = value
    return out[0] if scalar else out


def lfactorial(x) -> np.ndarray:
    """R's `lfactorial(x)` = lgamma(x + 1)."""
    return lgammafn(np.asarray(x, dtype=np.float64) + 1)


# ---------------------------------------------------------------------------
# Long-double accumulation
# ---------------------------------------------------------------------------

X87 = np.finfo(np.longdouble).nmant == 63
_TWO63 = float(2**63)
_TWO64 = float(2**64)


def _ld_add(num: int, exp: int, value: float) -> Tuple[int, int]:
    """(num·2**exp) + value, rounded to a 64-bit mantissa (ties to even), normalised."""
    p, q = float(value).as_integer_ratio()
    e = -(q.bit_length() - 1)
    common = min(exp, e)
    num = (num << (exp - common)) + (p << (e - common))
    exp = common
    if num == 0:
        return 0, 0
    length = abs(num).bit_length()
    if length > 64:
        excess = length - 64
        sign = -1 if num < 0 else 1
        quotient, remainder = divmod(abs(num), 1 << excess)
        half = 1 << (excess - 1)
        if remainder > half or (remainder == half and quotient & 1):
            quotient += 1
        num, exp = sign * quotient, exp + excess
        if abs(num).bit_length() > 64:  # rounding carried into a new bit
            num, exp = num // 2, exp + 1
    elif length < 64:
        num, exp = num << (64 - length), exp - (64 - length)
    return num, exp


def _ld_to_float(num: int, exp: int) -> float:
    if num == 0:
        return 0.0
    # int / power of two: Python rounds the exact quotient to nearest (ties to even)
    return num / (1 << -exp) if exp < 0 else float(num << exp)


def ld_cumsum(values) -> np.ndarray:
    """R's `cumsum` of doubles: a long-double running sum, each value rounded to double."""
    values = np.asarray(values, dtype=np.float64)
    if X87:
        return np.cumsum(values, dtype=np.longdouble).astype(np.float64)
    if not np.all(np.isfinite(values)):
        # -inf, inf or NaN take over the sum, whatever the precision
        return np.cumsum(values)
    out = np.empty(values.size)
    num, exp = 0, 0
    for k, value in enumerate(values.tolist()):
        num, exp = _ld_add(num, exp, value)
        out[k] = _ld_to_float(num, exp)
    return out


def ld_sum(values) -> float:
    """R's `sum` of doubles: a long-double running sum, rounded to double once."""
    values = np.asarray(values, dtype=np.float64).ravel()
    if values.size == 0:
        return 0.0
    if X87:
        return float(np.cumsum(values, dtype=np.longdouble)[-1])
    if not np.all(np.isfinite(values)):
        # -inf (a zero probability's log), inf or NaN take over the sum, whatever the precision
        return float(np.cumsum(values)[-1])
    num, exp = 0, 0
    i, n = 0, values.size
    length = 256
    while i < n:
        if num == 0:
            num, exp = _ld_add(num, exp, float(values[i]))
            i += 1
            continue
        # With |s| in [2**63, 2**64)·2**exp, s + v rounds to a multiple of 2**exp, as long as
        # the result stays in that binade: s/2**exp + round(v/2**exp), ties to even.
        chunk = values[i : i + length]
        q = np.ldexp(chunk, -exp)
        if not np.all(np.isfinite(q)) or np.any(np.abs(q) >= 2.0**62):
            stop = int(np.argmax(~np.isfinite(q) | (np.abs(q) >= 2.0**62)))
        else:
            stop = chunk.size
        floor = np.floor(q[:stop])
        tie = (q[:stop] - floor) == 0.5
        steps = np.where(tie, floor, np.rint(q[:stop]))
        approx = float(num) + np.cumsum(steps)
        inside = (np.abs(approx) > _TWO63 * (1 + 1e-9)) & (np.abs(approx) < _TWO64 * (1 - 1e-9))
        if not inside.all():
            stop = int(np.argmin(inside))
        if stop > 0:
            exact = np.cumsum(steps[:stop].astype(np.int64))
            total = num + int(exact[-1])
            correction = 0
            for t in np.flatnonzero(tie[:stop]).tolist():
                before = num + (int(exact[t - 1]) if t else 0) + correction
                if (before + int(floor[t])) & 1:
                    correction += 1
            num = total + correction
            i += stop
        if stop < chunk.size:
            num, exp = _ld_add(num, exp, float(values[i]))
            i += 1
            length = 256
        else:
            length = min(length * 2, 1 << 20)
    return _ld_to_float(num, exp)
