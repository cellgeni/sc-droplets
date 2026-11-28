"""DropletUtils' cell calling: barcodeRanks, emptyDrops, emptyDropsCellRanger, defaultDrops.

Ports of Bioconductor DropletUtils (GPL-3; Copyright (C) Aaron Lun, Jonathan
Griffiths and the DropletUtils authors): `R/barcodeRanks.R`, `R/emptyDrops.R`,
`R/ambientProfileEmpty.R`, `R/emptyDropsCellRanger.R`, `R/defaultDrops.R`,
`R/utils.R` and `src/montecarlo_pval.cpp`, at the releases each variant
names; with edgeR's Good-Turing (`primitives.goodturing`), R's
`smooth.spline` (`primitives.smoothspline`), R's arithmetic
(`primitives.rstats`) and DropletUtils' random numbers
(`primitives.rrandom`, `primitives.boostrandom`). Pages:
docs/methods/dropletutils/.

DropletUtils sees every column of the matrix as a barcode, including those
with no counts, in the matrix's order.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from decimal import Decimal, localcontext
from typing import Any, Dict, Optional, Tuple

import numpy as np
import pandas as pd
import scipy.sparse as sp

from sc_droplets.methods._common import result
from sc_droplets.primitives.goodturing import good_turing_proportions_edger
from sc_droplets.primitives.rmath import clog, clog10, ld_cumsum, ld_sum, lfactorial, lgammafn
from sc_droplets.primitives.rstats import all_equal, optimize_max, p_adjust_bh, quantile7
from sc_droplets.primitives.simulate import dropletutils as simulate_dropletutils
from sc_droplets.primitives.smoothspline import smooth_spline
from sc_droplets.registry import Method, MethodInfo, Reference, Variant, param, register

TOOL = "DropletUtils"
LATEST = "1.34.0"
SEED = 100  # the DropletUtils vignette's set.seed(100); R itself has no default
BROWSE = "https://code.bioconductor.org/browse/DropletUtils/blob/RELEASE_3_22/R/"
REFERENCES = (
    Reference(
        "Lun et al. 2019, EmptyDrops: distinguishing cells from empty droplets in droplet-based "
        "single-cell RNA sequencing data, Genome Biology 20:63",
        doi="10.1186/s13059-019-1662-y",
    ),
)


class DropletUtilsError(ValueError):
    """DropletUtils would stop with an error here."""


# ---------------------------------------------------------------------------
# R arithmetic used here
# ---------------------------------------------------------------------------


def r_pow10(y: float) -> float:
    """`10^y` as R computes it on x86-64 (`powl`, then rounded): correctly rounded."""
    with localcontext() as ctx:
        ctx.prec = 60
        return float(Decimal(10) ** Decimal(y))


def _rounded_counts(counts: sp.csr_matrix) -> sp.csr_matrix:
    """`round(m)`: half to even."""
    if counts.dtype.kind != "f":
        return counts
    out = counts.copy()
    out.data = np.rint(out.data)
    out.eliminate_zeros()
    return out


def _int_totals(counts: sp.csr_matrix) -> np.ndarray:
    """`as.integer(round(colSums(m)))`."""
    return np.rint(np.asarray(counts.sum(axis=1), dtype=np.float64).ravel()).astype(np.int64)


_ld_cumsum = ld_cumsum  # R's cumsum (long double)


# ---------------------------------------------------------------------------
# barcodeRanks
# ---------------------------------------------------------------------------


@dataclass
class Ranks:
    """The barcode-rank curve and the points barcodeRanks finds on it."""

    knee: float
    inflection: float
    x: np.ndarray  # log10 mid-rank of the runs on the curve
    y: np.ndarray  # log10 total
    diagnostics: Dict[str, Any]


def _runs(totals: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Run totals (decreasing) and their mid-ranks, as `rle` over the sorted totals."""
    s = np.sort(totals)[::-1]
    if s.size == 0:
        return s, s.astype(np.float64)
    starts = np.flatnonzero(np.concatenate(([True], s[1:] != s[:-1])))
    lengths = np.diff(np.concatenate((starts, [s.size])))
    rank = np.cumsum(lengths).astype(np.float64) - (lengths - 1) / 2.0
    return s[starts].astype(np.float64), rank


def _curve_bounds(x: np.ndarray, y: np.ndarray, exclude_from: Optional[float]) -> Tuple[int, int]:
    """0-based (left, right) edges: the largest slope up to the least slope."""
    d1n = np.diff(y) / np.diff(x)
    skip = 0
    if exclude_from is not None:
        skip = min(d1n.size - 1, int(np.sum(x <= math.log10(exclude_from))))
        d1n = d1n[skip:]
    right = int(np.argmin(d1n))
    left = int(np.argmax(d1n[: right + 1]))
    return left + skip, right + skip


def _interpolate(target: np.ndarray, cumulative: np.ndarray, stepwise: np.ndarray, x: np.ndarray, y: np.ndarray):
    """`.interpolate_on_curve`: findInterval(left.open = TRUE), then linear interpolation."""
    hi = np.searchsorted(cumulative, target, side="left")
    lo = hi - 1
    prop = (target - cumulative[lo]) / stepwise[hi]
    return x[lo] + prop * (x[hi] - x[lo]), y[lo] + prop * (y[hi] - y[lo])


def barcode_ranks(totals: np.ndarray, settings: Dict[str, Any], *, lower: float, exclude_from=None, df=None,
                  fit_bounds=None, window=None, gradient_threshold=None) -> Ranks:
    """`barcodeRanks()`'s knee and inflection for the given column sums, in one release's form."""
    totals = np.asarray(totals, dtype=np.float64)
    if settings.get("round_totals"):
        totals = np.rint(totals)
    run_totals, run_rank = _runs(totals)
    keep = run_totals > lower
    rule = settings["knee"]
    if rule == "trace":
        keep &= ~(run_rank <= exclude_from)
        if keep.sum() < 2:
            raise DropletUtilsError("insufficient unique points for computing knee/inflection points")
    elif keep.sum() < 3:
        raise DropletUtilsError("insufficient unique points for computing knee/inflection points")
    y = clog10(run_totals[keep])
    x = clog10(run_rank[keep])
    diagnostics: Dict[str, Any] = {"rule": rule, "n_points": int(x.size)}
    if rule == "trace":
        knee, inflection = _trace(x, y, window, gradient_threshold, diagnostics)
        return Ranks(knee, inflection, x, y, diagnostics)

    left, right = _curve_bounds(x, y, exclude_from if settings.get("exclude") else None)
    inflection = r_pow10(float(y[right]))
    if fit_bounds is None:
        chosen = np.arange(left, right + 1)
    else:
        chosen = np.flatnonzero((y > math.log10(fit_bounds[0])) & (y < math.log10(fit_bounds[1])))
    diagnostics.update(left=int(left), right=int(right))
    if rule == "chord":
        if chosen.size >= 4:
            cx, cy = x[chosen], y[chosen]
            gradient = (cy[-1] - cy[0]) / (cx[-1] - cx[0])
            intercept = cy[0] - cx[0] * gradient
            above = np.flatnonzero(cy >= cx * gradient + intercept)
            dist = np.abs(gradient * cx[above] - cy[above] + intercept) / math.sqrt(gradient * gradient + 1)
            knee = r_pow10(float(cy[above[int(np.argmax(dist))]]))
            diagnostics["chord"] = [float(cx[0]), float(cy[0]), float(cx[-1]), float(cy[-1])]
        else:
            knee = r_pow10(float(y[chosen[0]]))
        return Ranks(knee, inflection, x, y, diagnostics)

    # spline curvature
    if chosen.size >= 4 or not settings.get("fallback"):
        fit = smooth_spline(x[chosen], y[chosen], df=df)
        d1 = fit.predict(deriv=1)
        d2 = fit.predict(deriv=2)
        curvature = d2 / (1.0 + d1 * d1) ** 1.5
        at = int(np.nanargmin(curvature))
        if settings.get("index_bug"):
            knee = r_pow10(float(y[at]))  # 1.0 - 1.8 index the full curve with a position in the fit
        else:
            knee = r_pow10(float(y[chosen][at]))
        diagnostics.update(
            fitted_x=fit.x.tolist(), fitted_y=fit.fitted.tolist(), spar=fit.spar, df=fit.df,
            spline_iterations=fit.iterations,
        )
    else:
        knee = r_pow10(float(y[chosen[0]]))
    return Ranks(knee, inflection, x, y, diagnostics)


def _trace(x: np.ndarray, y: np.ndarray, window: float, threshold: float, diagnostics: Dict[str, Any]) -> Tuple[float, float]:
    """1.30's knee and inflection by a window traced along the curve."""
    dx, dy = np.diff(x), np.diff(y)
    stepwise = np.concatenate(([0.0], np.sqrt(dx * dx + dy * dy)))
    cumulative = _ld_cumsum(stepwise)
    rhs = cumulative + window
    scan = rhs <= cumulative[-1]
    if not scan.any():
        value = r_pow10(float(y[-1]))
        return value, value
    left_x, left_y = x[scan], y[scan]
    right_x, right_y = _interpolate(rhs[scan], cumulative, stepwise, x, y)
    gx, gy = left_x - right_x, left_y - right_y
    gap = np.sqrt(gx * gx + gy * gy)
    gradient = (right_y - left_y) / (right_x - left_x)
    intercept = right_y - gradient * right_x
    mid_x, mid_y = _interpolate(cumulative[scan] + window / 2, cumulative, stepwise, x, y)
    mid_above = mid_y > gradient * mid_x + intercept
    elbows = np.flatnonzero(~mid_above & (gradient < threshold))
    if elbows.size == 0:
        inflection_window = int(np.argmin(gradient))
        maybe_knee = np.flatnonzero(mid_above)
    else:
        first = int(elbows[0])
        maybe_knee = np.flatnonzero(mid_above[: first + 1])
        before = int(np.searchsorted(left_x, mid_x[first], side="right"))
        inflection_window = int(np.argmin(gradient[:before]))
    if maybe_knee.size:
        knee_window = int(maybe_knee[int(np.argmin(gap[maybe_knee]))])
    else:
        knee_window = inflection_window
    diagnostics.update(
        knee_window=[float(left_x[knee_window]), float(right_x[knee_window])],
        inflection_window=[float(left_x[inflection_window]), float(right_x[inflection_window])],
        n_windows=int(scan.sum()),
    )
    return r_pow10(float(mid_y[knee_window])), r_pow10(float(mid_y[inflection_window]))


_BR_VARIANTS = (
    Variant("v1.0", since="1.0", until="1.10", summary="spline curvature between the numerical edges; "
            "the knee indexed in the full curve (bug)",
            settings={"knee": "spline", "index_bug": True}, defaults={"df": 20.0}),
    Variant("v1.10", since="1.10", until="1.12", summary="the first 50 barcodes skipped; knee within the fit; "
            "first fitted point below 4 points",
            settings={"knee": "spline", "exclude": True, "fallback": True}, defaults={"df": 20.0, "exclude_from": 50}),
    Variant("v1.12", since="1.12", until="1.28", summary="totals rounded to integers",
            settings={"knee": "spline", "exclude": True, "fallback": True, "round_totals": True},
            defaults={"df": 20.0, "exclude_from": 50}),
    Variant("v1.28", since="1.28", until="1.30", summary="knee: farthest point above the chord between the edges",
            settings={"knee": "chord", "exclude": True, "round_totals": True}, defaults={"exclude_from": 50}),
    Variant("v1.30", since="1.30", summary="knee and inflection by a window traced along the curve",
            settings={"knee": "trace", "round_totals": True},
            defaults={"exclude_from": 50, "window": 1.0, "gradient_threshold": -1.0}),
)


def _br_variant(name: str) -> Variant:
    return next(v for v in _BR_VARIANTS if v.name == name)


def _check_given(method: str, variant: Variant, params: Any, allowed: Dict[str, bool]) -> None:
    for name, ok in allowed.items():
        if not ok and getattr(params, name) is not None:
            flag = "--" + name.replace("_", "-")
            raise ValueError(f"{method}: DropletUtils {variant.versions} has no {flag}")


@register
class BarcodeRanks(Method):
    info = MethodInfo(
        id="dropletutils-barcoderanks",
        title="DropletUtils barcodeRanks",
        tool=TOOL,
        category="knee",
        summary="The knee (or inflection) of the log-log barcode-rank curve: spline curvature (1.0, 1.10), "
        "chord distance (1.28) or a traced window (1.30); cells are the barcodes at or above it.",
        docs="methods/dropletutils/barcoderanks.md",
        variants=_BR_VARIANTS,
        latest_known=LATEST,
        references=REFERENCES,
        source=BROWSE + "barcodeRanks.R",
        original_license="GPL-3",
    )

    @dataclass
    class Params:
        lower: float = param(100.0, "Runs with total at or below this are not on the curve", original="lower", minimum=0)
        exclude_from: Optional[float] = param(None, "1.10+: barcodes ranked this high or higher are skipped (default 50)", original="exclude.from", minimum=0)
        df: Optional[float] = param(None, "Up to 1.26: degrees of freedom of the smoothing spline (default 20)", original="df", minimum=0)
        fit_lower: Optional[float] = param(None, "Up to 1.28: fit the runs with total above this (with --fit-upper)", original="fit.bounds[1]", minimum=0)
        fit_upper: Optional[float] = param(None, "Up to 1.28: ... and below this", original="fit.bounds[2]", minimum=0)
        window: Optional[float] = param(None, "1.30+: window length along the curve (default 1)", original="window", minimum=0)
        gradient_threshold: Optional[float] = param(None, "1.30+: slope below which a window under its line is an elbow (default -1)", original="gradient.threshold")
        use: str = param("knee", "The cut: knee or inflection", choices=("knee", "inflection"))

    def run(self, droplets, params, variant, random_state):
        rule = variant.settings["knee"]
        _check_given(self.info.id, variant, params, {
            "exclude_from": bool(variant.settings.get("exclude")) or rule == "trace",
            "df": rule == "spline",
            "fit_lower": rule != "trace",
            "fit_upper": rule != "trace",
            "window": rule == "trace",
            "gradient_threshold": rule == "trace",
        })
        bounds = _fit_bounds(self.info.id, params)
        ranks = barcode_ranks(
            droplets.totals, variant.settings, lower=params.lower, exclude_from=params.exclude_from, df=params.df,
            fit_bounds=bounds, window=params.window, gradient_threshold=params.gradient_threshold,
        )
        cut = ranks.knee if params.use == "knee" else ranks.inflection
        totals = _br_totals(droplets.totals, variant.settings)
        return result(
            type(self), variant, params, droplets, totals >= cut,
            thresholds={"knee": ranks.knee, "inflection": ranks.inflection, "umis": cut},
            diagnostics=ranks.diagnostics,
        )


def _br_totals(totals: np.ndarray, settings: Dict[str, Any]) -> np.ndarray:
    totals = np.asarray(totals, dtype=np.float64)
    return np.rint(totals) if settings.get("round_totals") else totals


def _fit_bounds(method: str, params: Any) -> Optional[Tuple[float, float]]:
    if (params.fit_lower is None) != (params.fit_upper is None):
        raise ValueError(f"{method}: give both --fit-lower and --fit-upper, or neither")
    return None if params.fit_lower is None else (params.fit_lower, params.fit_upper)


# ---------------------------------------------------------------------------
# emptyDrops
# ---------------------------------------------------------------------------


_ED_CLASSES = ["retained", "rescued", "tested", "ignored", "ambient"]
_CR_CLASSES = ["initial", "rescued", "candidate", "background"]


def _ambient_profile(counts: sp.csr_matrix, ambient: np.ndarray, safe: bool, stop_empty: bool):
    """Features kept (any count anywhere) and their Good-Turing ambient proportions."""
    features = np.flatnonzero(np.asarray(counts.sum(axis=0)).ravel() > 0)
    profile = np.asarray(counts[ambient][:, features].sum(axis=0), dtype=np.float64).ravel()
    if profile.sum() == 0:
        if stop_empty:
            raise DropletUtilsError("no counts available to estimate the ambient profile")
        raise DropletUtilsError("the ambient droplets hold no counts; edgeR's goodTuring has nothing to fit")
    prop = good_turing_proportions_edger(profile)
    if safe and np.isnan(prop).any():
        # R: `if (any(still.zero))` on NA stops with "missing value where TRUE/FALSE needed".
        raise DropletUtilsError("edgeR's Good-Turing gave no proportions (too few distinct ambient counts)")
    if safe:
        still_zero = prop <= 0
        if still_zero.any():
            pseudo = 1.0 / profile.sum()
            prop[still_zero] = pseudo / still_zero.sum()
            prop[~still_zero] = prop[~still_zero] * (1.0 - pseudo)
    return features, profile, prop


def _estimate_alpha(ambient_counts: sp.csr_matrix, prop: np.ndarray, totals: np.ndarray) -> float:
    """`.estimate_alpha`: the Dirichlet-multinomial MLE of alpha on (0.01, 10000).

    The likelihood is flat enough near its optimum that R's rounding steers
    `optimize`, so it is evaluated as R evaluates it: R's lgamma, and sums in
    long double in R's order (barcodes in order, features ascending).
    """
    block = ambient_counts.tocsr()
    block.sort_indices()
    per_prop = prop[block.indices]
    x = block.data.astype(np.float64)
    distinct, where = np.unique(totals.astype(np.float64), return_inverse=True)
    n = totals.size

    def loglik(alpha: float) -> float:
        prop_alpha = per_prop * alpha
        by_total = lgammafn(distinct + alpha)[where]
        return (
            float(lgammafn(alpha)) * n
            - ld_sum(by_total)
            + ld_sum(lgammafn(x + prop_alpha))
            - ld_sum(lgammafn(prop_alpha))
        )

    return optimize_max(loglik, 0.01, 10000.0)


def _observed(counts: sp.csr_matrix, rows: np.ndarray, prop: np.ndarray, alpha: float) -> np.ndarray:
    """The data-dependent log-probability of each tested droplet (`.compute_multinom_prob_data`)."""
    block = counts[rows].tocsr()
    block.sort_indices()
    x = block.data.astype(np.float64)
    g = block.indices
    if math.isinf(alpha):
        terms = x * clog(prop)[g] - lfactorial(x)
    else:
        alpha_prop = alpha * prop[g]
        terms = lgammafn(alpha_prop + x) - lfactorial(x) - lgammafn(alpha_prop)
    indptr = block.indptr
    # tapply(p.n0, j, sum): R's sum, droplet by droplet
    return np.array([ld_sum(terms[indptr[i] : indptr[i + 1]]) for i in range(rows.size)])


def _rest(totals: np.ndarray, alpha: float) -> np.ndarray:
    t = totals.astype(np.float64)
    if math.isinf(alpha):
        return lfactorial(t)
    return lfactorial(t) + float(lgammafn(alpha)) - lgammafn(t + alpha)


def _alpha_value(method: str, text: Optional[str]) -> Optional[float]:
    """None for "estimate", else the number (inf for the multinomial)."""
    word = str(text).strip().lower()
    if word == "estimate":
        return None
    if word in ("inf", "infinity", "multinomial"):
        return math.inf
    try:
        value = float(word)
    except ValueError:
        raise ValueError(f"{method}: --alpha must be estimate, inf or a positive number, got {text!r}") from None
    if not value > 0:
        raise ValueError(f"{method}: --alpha must be positive, got {text!r}")
    return value


def _round_for(settings: Dict[str, Any], counts: sp.csr_matrix, round_counts: Optional[bool]) -> sp.csr_matrix:
    mode = settings.get("round", "never")
    if mode == "never" or counts.dtype.kind != "f" or round_counts is False:
        return counts
    if mode == "if-needed":
        cs = np.asarray(counts.sum(axis=1)).ravel()
        rs = np.asarray(counts.sum(axis=0)).ravel()
        if all_equal(cs, np.rint(cs)) and all_equal(rs, np.rint(rs)):
            return counts
    return _rounded_counts(counts)


@dataclass
class _Test:
    pvalue: np.ndarray
    logprob: np.ndarray
    limited: np.ndarray
    tested: np.ndarray
    alpha: float
    diagnostics: Dict[str, Any]


def _test_empty_drops(counts, totals, ambient, keep, *, niters, alpha, seed, sampler, safe, stop_empty) -> _Test:
    """`.test_empty_drops`: profile, alpha, observed values and Monte Carlo p-values."""
    features, profile, prop = _ambient_profile(counts, ambient, safe, stop_empty)
    kept = counts[:, features].tocsr()
    if alpha is None:
        alpha = _estimate_alpha(kept[ambient], prop, totals[ambient])
    rows = np.flatnonzero(keep)
    observed = _observed(kept, rows, prop, alpha)
    n_above = simulate_dropletutils(prop, totals[rows], observed, niters, seed, alpha=alpha, sampler=sampler)
    n = totals.size
    pvalue = np.full(n, np.nan)
    logprob = np.full(n, np.nan)
    limited = np.zeros(n, dtype=bool)
    pvalue[rows] = (n_above + 1) / (niters + 1)
    logprob[rows] = observed + _rest(totals[rows], alpha)
    limited[rows] = n_above == 0
    top = np.argsort(-prop, kind="stable")[:20]
    diagnostics = {
        "alpha": alpha if math.isfinite(alpha) else "Inf",
        "niters": niters,
        "seed": seed,
        "n_tested": int(rows.size),
        "n_ambient": int(ambient.sum()),
        "n_features": int(features.size),
        "ambient_total": float(profile.sum()),
        "profile_top": [int(features[i]) for i in top],
    }
    return _Test(pvalue, logprob, limited, rows, alpha, diagnostics)


_ED_VARIANTS = (
    Variant("v1.0", since="1.0", until="1.2", summary="multinomial; R's random numbers",
            settings={"sampler": "r", "alpha": False, "knee": "v1.0"}),
    Variant("v1.2", since="1.2", until="1.4", summary="Dirichlet-multinomial with estimated alpha",
            settings={"sampler": "r", "alpha": True, "knee": "v1.0"}, defaults={"alpha": "estimate"}),
    Variant("v1.4", since="1.4", until="1.6", summary="PCG32 per iteration, inverse-CDF draws",
            settings={"sampler": "inverse", "alpha": True, "stop_empty": True, "knee": "v1.0"},
            defaults={"alpha": "estimate"}),
    Variant("v1.6", since="1.6", until="1.8", summary="Boost alias-table draws; non-integer matrices rounded",
            settings={"sampler": "alias", "alpha": True, "stop_empty": True, "round": "if-needed", "knee": "v1.0"},
            defaults={"alpha": "estimate"}),
    Variant("v1.8", since="1.8", until="1.10", summary="round=; pseudo-count for zero probabilities",
            settings={"sampler": "alias", "alpha": True, "stop_empty": True, "round": "if-needed",
                      "round_param": True, "safe": True, "knee": "v1.0"},
            defaults={"alpha": "estimate", "round_counts": True}),
    Variant("v1.10", since="1.10", until="1.12", summary="by.rank (retain still from the given lower); knee v1.10",
            settings={"sampler": "alias", "alpha": True, "stop_empty": True, "round": "if-needed",
                      "round_param": True, "safe": True, "by_rank": "lower", "knee": "v1.10"},
            defaults={"alpha": "estimate", "round_counts": True}),
    Variant("v1.12", since="1.12", until="1.20", summary="always rounds; test.ambient FDR rule; knee v1.12",
            settings={"sampler": "alias", "alpha": True, "stop_empty": True, "round": "always", "round_param": True,
                      "safe": True, "by_rank": "lower", "fdr_na": True, "retain_from_lower": True, "knee": "v1.12"},
            defaults={"alpha": "estimate", "round_counts": True}),
    Variant("v1.20", since="1.20", until="1.28", summary="by.rank by exact rank",
            settings={"sampler": "alias", "alpha": True, "stop_empty": True, "round": "always", "round_param": True,
                      "safe": True, "by_rank": "first", "fdr_na": True, "retain_from_lower": True, "knee": "v1.12"},
            defaults={"alpha": "estimate", "round_counts": True}),
    Variant("v1.28", since="1.28", until="1.30", summary="multinomial again (alpha = Inf); knee v1.28",
            settings={"sampler": "alias", "alpha": True, "stop_empty": True, "round": "always", "round_param": True,
                      "safe": True, "by_rank": "first", "fdr_na": True, "retain_from_lower": True, "knee": "v1.28"},
            defaults={"alpha": "inf", "round_counts": True}),
    Variant("v1.30", since="1.30", summary="knee v1.30",
            settings={"sampler": "alias", "alpha": True, "stop_empty": True, "round": "always", "round_param": True,
                      "safe": True, "by_rank": "first", "fdr_na": True, "retain_from_lower": True, "knee": "v1.30"},
            defaults={"alpha": "inf", "round_counts": True}),
)


@register
class EmptyDrops(Method):
    info = MethodInfo(
        id="dropletutils-emptydrops",
        title="DropletUtils emptyDrops",
        tool=TOOL,
        category="statistical-test",
        summary="Every droplet above `lower` tested against the Good-Turing ambient profile by Monte Carlo "
        "(multinomial, or Dirichlet-multinomial 1.2-1.26); droplets at or above the knee always kept; BH FDR.",
        docs="methods/dropletutils/emptydrops.md",
        variants=_ED_VARIANTS,
        latest_known=LATEST,
        references=REFERENCES,
        source=BROWSE + "emptyDrops.R",
        original_license="GPL-3",
        feature_types=("Gene Expression",),
        default_seed=SEED,
    )

    @dataclass
    class Params:
        lower: float = param(100.0, "Droplets with this total or less are ambient", original="lower", minimum=0)
        niters: int = param(10000, "Monte Carlo iterations", original="niters", minimum=1)
        test_ambient: str = param("false", "Test the ambient droplets too (true); na (1.12+) also keeps them in the FDR", original="test.ambient", choices=("false", "true", "na"))
        ignore: Optional[float] = param(None, "Droplets with this total or less are not tested", original="ignore", minimum=0)
        alpha: Optional[str] = param(None, "1.2+: estimate, inf (multinomial) or a number (default: the release's)", original="alpha")
        retain: Optional[float] = param(None, "Droplets with at least this total are kept untested (default: the knee)", original="retain", minimum=0)
        by_rank: Optional[int] = param(None, "1.10+: the droplets ranked below this many are ambient", original="by.rank", minimum=0)
        round_counts: Optional[bool] = param(None, "1.8+: round non-integer counts (default on)", original="round")
        fdr: float = param(0.001, "Largest FDR called a cell", minimum=0, maximum=1)

    def run(self, droplets, params, variant, random_state):
        s = variant.settings
        _check_given(self.info.id, variant, params, {
            "alpha": bool(s.get("alpha")),
            "by_rank": bool(s.get("by_rank")),
            "round_counts": bool(s.get("round_param")),
        })
        if params.test_ambient == "na" and not s.get("fdr_na"):
            raise ValueError(f"{self.info.id}: --test-ambient na needs DropletUtils 1.12 or later")
        alpha = math.inf if params.alpha is None else _alpha_value(self.info.id, params.alpha)
        if s["sampler"] == "r" and (alpha is None or math.isfinite(alpha)):
            raise NotImplementedError(
                f"{self.info.id}: DropletUtils {variant.versions} draws its Dirichlet profiles with R's rgamma "
                "in one R stream, which is not reproduced; pass --alpha inf for the multinomial"
            )
        seed = SEED if random_state is None else int(random_state)
        counts = _round_for(s, droplets.counts, params.round_counts)
        totals = _int_totals(counts)
        n = totals.size
        lower = params.lower
        if params.by_rank is not None:
            if params.by_rank >= n:
                raise DropletUtilsError("not have enough columns for supplied 'by.rank'")
            order = np.argsort(-totals, kind="stable")
            if s["by_rank"] == "lower":
                lower = float(totals[order[params.by_rank]])
                ambient = totals <= lower
            else:
                ambient = np.ones(n, dtype=bool)
                ambient[order[: params.by_rank]] = False
                lower = float(totals[ambient].max())
        else:
            ambient = totals <= lower
        test_ambient = params.test_ambient != "false"
        keep = (totals > 0) if test_ambient else ~ambient
        if params.ignore is not None:
            keep &= totals > params.ignore
        test = _test_empty_drops(
            counts, totals, ambient, keep, niters=params.niters, alpha=alpha, seed=seed, sampler=s["sampler"],
            safe=bool(s.get("safe")), stop_empty=bool(s.get("stop_empty")),
        )
        retain = params.retain
        knee_info: Dict[str, Any] = {}
        if retain is None:
            knee_lower = lower if s.get("retain_from_lower") else params.lower
            br = _br_variant(s["knee"])
            ranks = barcode_ranks(
                np.asarray(counts.sum(axis=1), dtype=np.float64).ravel(), br.settings, lower=knee_lower,
                exclude_from=br.defaults.get("exclude_from"), df=br.defaults.get("df"),
                window=br.defaults.get("window"), gradient_threshold=br.defaults.get("gradient_threshold"),
            )
            retain = ranks.knee
            knee_info = {"knee_variant": br.name, "inflection": ranks.inflection}
        tmp = test.pvalue.copy()
        always = totals >= retain
        tmp[always] = 0.0
        if test_ambient and s.get("fdr_na") and params.test_ambient == "true":
            tmp[totals <= lower] = np.nan
        fdr = p_adjust_bh(tmp)
        cells = np.nan_to_num(fdr, nan=np.inf) <= params.fdr
        classes = np.full(n, "ambient", dtype=object)
        if params.ignore is not None:
            classes[~ambient & (totals <= params.ignore)] = "ignored"
        classes[test.tested] = "tested"
        classes[cells & ~always] = "rescued"
        classes[always] = "retained"
        diagnostics = dict(test.diagnostics, lower=lower, retain=retain, **knee_info)
        return result(
            type(self), variant, params, droplets, cells,
            classes=pd.Categorical(classes, categories=_ED_CLASSES),
            columns={"pvalue": test.pvalue, "fdr": fdr, "logprob": test.logprob, "limited": test.limited},
            thresholds={"retain": float(retain), "lower": float(lower), "fdr": params.fdr},
            diagnostics=diagnostics,
        )


# ---------------------------------------------------------------------------
# emptyDropsCellRanger
# ---------------------------------------------------------------------------


@register
class EmptyDropsCellRanger(Method):
    info = MethodInfo(
        id="dropletutils-emptydropscellranger",
        title="DropletUtils emptyDropsCellRanger",
        tool=TOOL,
        category="statistical-test",
        summary="STARsolo's EmptyDrops_CR rewritten in R: robust-maximum cut, ambient ranks 45,000-90,000, "
        "candidates just below the cut tested with emptyDrops' multinomial Monte Carlo.",
        docs="methods/dropletutils/emptydropscellranger.md",
        variants=(
            Variant("v1.14", since="1.14", until="1.18.1", summary="candidates' minimum from the unsorted totals (bug)",
                    settings={"median_bug": True}),
            Variant("v1.18.1", since="1.18.1", summary="candidates' minimum from the sorted totals"),
        ),
        latest_known=LATEST,
        references=REFERENCES,
        source=BROWSE + "emptyDropsCellRanger.R",
        original_license="GPL-3",
        default_seed=SEED,
    )

    @dataclass
    class Params:
        expected_cells: float = param(3000.0, "Expected number of cells", original="n.expected.cells", minimum=0)
        max_percentile: float = param(0.99, "Robust maximum percentile", original="max.percentile", minimum=0, maximum=1)
        max_min_ratio: float = param(10.0, "Robust maximum ÷ this is the initial cut", original="max.min.ratio", minimum=0)
        umi_min: float = param(500.0, "Fewest UMIs a candidate needs", original="umi.min", minimum=0)
        umi_min_frac_median: float = param(0.01, "Fraction of the middle initial cell's total a candidate needs", original="umi.min.frac.median", minimum=0)
        cand_max_n: int = param(20000, "Most candidates", original="cand.max.n", minimum=0)
        ind_min: int = param(45000, "First ambient position (1-based)", original="ind.min", minimum=1)
        ind_max: int = param(90000, "Last ambient position", original="ind.max", minimum=1)
        round_counts: bool = param(True, "Round non-integer counts", original="round")
        niters: int = param(10000, "Monte Carlo iterations", original="niters", minimum=1)
        fdr: float = param(0.01, "Largest FDR called a cell", minimum=0, maximum=1)

    def run(self, droplets, params, variant, random_state):
        seed = SEED if random_state is None else int(random_state)
        counts = droplets.counts
        if params.round_counts:
            counts = _rounded_counts(counts)
        totals = _int_totals(counts)
        n = totals.size
        order = np.argsort(-totals, kind="stable")
        max_ind = int(np.rint(params.expected_cells * (1 - params.max_percentile)))
        position = min(n, max_ind)
        # R: an index of 0 selects nothing, and max(numeric(0), 1) is 1.
        n_umi_max = float(totals[order[position - 1]]) if position >= 1 else math.nan
        retain = 1.0 if position < 1 else max(float(np.rint(n_umi_max / params.max_min_ratio)), 1.0)
        n_simple = int(np.sum(totals >= retain))
        middle = int(n_simple / 2)  # R truncates a fractional index
        if middle >= 1:
            at = totals[middle - 1] if variant.settings.get("median_bug") else totals[order[middle - 1]]
            min_umi = max(params.umi_min, float(np.rint(params.umi_min_frac_median * at)))
        else:
            min_umi = params.umi_min
        cand_last = min(n_simple + params.cand_max_n, int(np.sum(totals > min_umi)))
        ambient = np.zeros(n, dtype=bool)
        ambient[order[min(params.ind_min, n) - 1 : min(params.ind_max, n)]] = True
        keep = np.zeros(n, dtype=bool)
        first = n_simple + 1
        positions = np.arange(first, cand_last + 1) if cand_last >= first else np.arange(first, cand_last - 1, -1)
        positions = positions[(positions >= 1) & (positions <= n)]
        keep[order[positions - 1]] = True
        test = _test_empty_drops(
            counts, totals, ambient, keep, niters=params.niters, alpha=math.inf, seed=seed, sampler="alias",
            safe=True, stop_empty=True,
        )
        tmp = test.pvalue.copy()
        always = totals >= retain
        tmp[always] = 0.0
        fdr = p_adjust_bh(tmp)
        cells = np.nan_to_num(fdr, nan=np.inf) <= params.fdr
        classes = np.full(n, "background", dtype=object)
        classes[test.tested] = "candidate"
        classes[cells & ~always] = "rescued"
        classes[always] = "initial"
        diagnostics = dict(test.diagnostics, robust_max=n_umi_max, retain=retain, n_initial=n_simple,
                           min_umi=min_umi, n_candidates=int(keep.sum()))
        return result(
            type(self), variant, params, droplets, cells,
            classes=pd.Categorical(classes, categories=_CR_CLASSES),
            columns={"pvalue": test.pvalue, "fdr": fdr, "logprob": test.logprob, "limited": test.limited},
            thresholds={"robust_max": n_umi_max, "umis": retain, "candidate_umis": min_umi, "fdr": params.fdr},
            diagnostics=diagnostics,
        )


# ---------------------------------------------------------------------------
# defaultDrops
# ---------------------------------------------------------------------------


@register
class DefaultDrops(Method):
    info = MethodInfo(
        id="dropletutils-defaultdrops",
        title="DropletUtils defaultDrops",
        tool=TOOL,
        category="rank-threshold",
        summary="Cell Ranger 2's rule in R: 10 % of the 99th percentile (type 7) of the top `expected` totals; "
        "cells strictly above it.",
        docs="methods/dropletutils/defaultdrops.md",
        variants=(Variant("v1.0", since="1.0", summary="unchanged since 1.0"),),
        latest_known=LATEST,
        references=(
            Reference(
                "Zheng et al. 2017, Massively parallel digital transcriptional profiling of single cells, "
                "Nature Communications 8:14049",
                doi="10.1038/ncomms14049",
            ),
        ),
        source=BROWSE + "defaultDrops.R",
        original_license="GPL-3",
    )

    @dataclass
    class Params:
        expected_cells: int = param(3000, "How many top barcodes the quantile is taken over", original="expected", minimum=0)
        upper_quant: float = param(0.99, "The quantile", original="upper.quant", minimum=0, maximum=1)
        lower_prop: float = param(0.1, "Fraction of the quantile that is the threshold", original="lower.prop", minimum=0, maximum=1)

    def run(self, droplets, params, variant, random_state):
        libs = np.asarray(droplets.totals, dtype=np.float64)
        order = np.argsort(-libs, kind="stable")
        top = libs[order[: params.expected_cells]]
        if top.size == 0:
            raise DropletUtilsError("no barcodes to take the quantile of")
        quantile = quantile7(top, params.upper_quant)
        threshold = quantile * params.lower_prop
        return result(
            type(self), variant, params, droplets, libs > threshold,
            thresholds={"quantile": quantile, "umis": threshold},
        )
