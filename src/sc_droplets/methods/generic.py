"""Generic rules that belong to no single tool: top N, a fixed threshold, MAD outliers, Otsu and
multi-Otsu thresholds, and a kernel-density valley.

Pages: docs/methods/generic/{top-n,threshold,mad,otsu,multiotsu,kde-valley}.md. MAD follows
scuttle's `isOutlier` (R's `log2`, `median` and `mad`); Otsu and multi-Otsu are scikit-image's
(`primitives/otsu.py`).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional

import numpy as np
import pandas as pd

from sc_droplets.methods._common import result, top_by_rank
from sc_droplets.primitives.otsu import histogram, threshold_multiotsu, threshold_otsu
from sc_droplets.primitives.rmath import clog2, clog10
from sc_droplets.primitives.rstats import r_mad, r_median
from sc_droplets.registry import Method, MethodInfo, Reference, param, register

_GENERIC = "Generic"
_UNITS = ("umis", "reads", "genes")
_UNIT_LABEL = {"umis": "UMIs", "reads": "reads", "genes": "genes"}
_TRANSFORMS = ("log10", "arcsinh", "none")
#: Points of the kernel density's grid, as UMI-tools evaluates it.
GRID_POINTS = 10000


class RuleError(ValueError):
    """Where a generic rule cannot place its cut."""


def candidates(values: np.ndarray, min_counts: float) -> np.ndarray:
    """Indices of the barcodes with at least `min_counts`; the others are empty."""
    index = np.flatnonzero(values >= min_counts)
    if index.size == 0:
        raise RuleError(f"no barcode has at least --min-counts {min_counts:g}")
    return index


def transform(values: np.ndarray, how: str) -> np.ndarray:
    """log10 or arcsinh with the C library's functions, element by element, or the values."""
    if how == "log10":
        return clog10(values)
    if how == "arcsinh":
        return np.array([math.asinh(v) for v in np.asarray(values, dtype=np.float64).tolist()], dtype=np.float64)
    return np.asarray(values, dtype=np.float64)


def untransform(value: float, how: str) -> float:
    """A cut on the transformed scale, back in counts."""
    if how == "log10":
        return 10.0 ** value
    if how == "arcsinh":
        return math.sinh(value)
    return float(value)


def _scale_label(params) -> str:
    unit = _UNIT_LABEL[params.unit]
    return unit if params.transform == "none" else f"{params.transform} {unit}"


def _histogram_diagnostics(x: np.ndarray, n_bins: int, label: str, cuts) -> dict:
    counts, centers = histogram(x, n_bins)
    half = (centers[1] - centers[0]) / 2 if centers.size > 1 else 0.5
    return {
        "histogram_counts": counts.astype(np.int64),
        "histogram_edges": np.concatenate([centers - half, centers[-1:] + half]),
        "histogram_label": label,
        "cuts": np.asarray([c for c, _ in cuts], dtype=np.float64),
        "cut_labels": [str(t) for _, t in cuts],
    }


@register
class TopN(Method):
    info = MethodInfo(
        id="top-n",
        title="Top N barcodes",
        tool=_GENERIC,
        category="rank-threshold",
        summary="The N barcodes with the most counts are cells.",
        docs="methods/generic/top-n.md",
        counts="either",
        references=(
            Reference(
                "The rule behind Cell Ranger --force-cells, STARsolo TopCells, UMI-tools "
                "--set-cell-number and alevin --forceCells; each breaks ties its own way (see their pages)."
            ),
        ),
    )

    @dataclass
    class Params:
        n_cells: int = param(help="How many barcodes to call", original="--force-cells", minimum=0)
        ties: str = param(
            "first",
            "first: exactly N, ties broken by input order; include: every barcode with at least the N-th total",
            choices=("first", "include"),
        )
        unit: str = param("umis", "What to rank by", choices=("umis", "reads", "genes"))

    def run(self, droplets, params, variant, random_state):
        values = droplets.values(params.unit)
        n = min(params.n_cells, values.shape[0])
        if params.ties == "first":
            is_cell = top_by_rank(values, n)
        else:
            nth = np.sort(values)[::-1][n - 1] if n > 0 else np.inf
            is_cell = values >= nth if n > 0 else np.zeros(values.shape[0], dtype=bool)
        nth_value = float(np.sort(values)[::-1][n - 1]) if n > 0 else float("nan")
        return result(
            type(self),
            variant,
            params,
            droplets,
            is_cell,
            thresholds={params.unit: nth_value},
            diagnostics={"rank": n},
        )


@register
class Threshold(Method):
    info = MethodInfo(
        id="threshold",
        title="Fixed threshold",
        tool=_GENERIC,
        category="rank-threshold",
        summary="Barcodes with at least a fixed number of UMIs, reads or genes are cells.",
        docs="methods/generic/threshold.md",
        counts="either",
        references=(
            Reference(
                "Used as Drop-seq SelectCellsByNumTranscripts MIN_TRANSCRIPTS_PER_CELL, "
                "alevin-fry --min-reads, and the 'more than 200 genes' rule of the scanpy and Seurat tutorials."
            ),
        ),
    )

    @dataclass
    class Params:
        min_counts: float = param(help="The threshold", original="MIN_TRANSCRIPTS_PER_CELL", minimum=0)
        unit: str = param("umis", "What the threshold applies to", choices=("umis", "reads", "genes"))
        strict: bool = param(False, "Require more than the threshold, not at least it")

    def run(self, droplets, params, variant, random_state):
        values = droplets.values(params.unit)
        is_cell = values > params.min_counts if params.strict else values >= params.min_counts
        return result(
            type(self),
            variant,
            params,
            droplets,
            is_cell,
            thresholds={params.unit: float(params.min_counts)},
        )


@register
class Mad(Method):
    info = MethodInfo(
        id="mad",
        title="MAD outliers",
        tool=_GENERIC,
        category="density",
        summary="Among the barcodes above a floor, those whose log total is more than k MADs below the median "
                "are empty (scuttle's isOutlier, lower tail); the rest are cells.",
        docs="methods/generic/mad.md",
        counts="either",
        references=(
            Reference(
                "McCarthy, Campbell, Lun & Wills 2017, Scater: pre-processing, quality control, normalization "
                "and visualization of single-cell RNA-seq data in R, Bioinformatics 33:1179",
                doi="10.1093/bioinformatics/btw777",
            ),
        ),
        source="https://github.com/LTLA/scuttle/blob/master/R/isOutlier.R",
        original_license="GPL-3.0",
    )

    @dataclass
    class Params:
        nmads: float = param(3.0, "How many MADs below the median a barcode may be", original="nmads", minimum=0)
        min_counts: float = param(100.0, "Candidates have at least this value; the others are empty", minimum=0)
        log: bool = param(True, "Work on log2 values", original="log")
        min_diff: Optional[float] = param(None, "The least distance from the median to the limit", original="min.diff",
                                          minimum=0)
        unit: str = param("umis", "What the rule applies to", choices=_UNITS)

    def run(self, droplets, params, variant, random_state):
        values = droplets.values(params.unit)
        index = candidates(values, params.min_counts)
        metric = clog2(values[index]) if params.log else values[index].astype(np.float64)
        median = r_median(metric)
        mad = r_mad(metric, median)
        spread = params.nmads * mad
        if params.min_diff is not None and not math.isnan(spread):
            spread = max(params.min_diff, spread)
        elif params.min_diff is not None:
            spread = params.min_diff
        lower = median - spread
        is_cell = np.zeros(values.size, dtype=bool)
        is_cell[index[~(metric < lower)]] = True
        back = (lambda v: 2.0 ** v) if params.log else float
        label = ("log2 " if params.log else "") + _UNIT_LABEL[params.unit]
        return result(
            type(self), variant, params, droplets, is_cell,
            thresholds={params.unit: back(lower), "median": back(median), "mad": float(mad)},
            diagnostics=_histogram_diagnostics(metric, 100, label, [(median, "median"), (lower, f"limit: {params.nmads:g} MADs")]),
        )


@register
class Otsu(Method):
    info = MethodInfo(
        id="otsu",
        title="Otsu threshold",
        tool=_GENERIC,
        category="density",
        summary="Otsu's threshold (scikit-image's threshold_otsu) on the log totals of the barcodes above a "
                "floor; the upper class is the cells.",
        docs="methods/generic/otsu.md",
        counts="either",
        references=(
            Reference(
                "Otsu 1979, A threshold selection method from gray-level histograms, IEEE Transactions on "
                "Systems, Man, and Cybernetics 9:62",
                doi="10.1109/TSMC.1979.4310076",
            ),
        ),
        source="https://github.com/scikit-image/scikit-image/blob/v0.26.0/src/skimage/filters/thresholding.py",
        original_license="BSD-3-Clause",
    )

    @dataclass
    class Params:
        min_counts: float = param(100.0, "Candidates have at least this value; the others are empty", minimum=0)
        transform: str = param("log10", "The scale the histogram is on", choices=_TRANSFORMS)
        bins: int = param(256, "Histogram bins", original="nbins", minimum=2)
        unit: str = param("umis", "What the rule applies to", choices=_UNITS)

    def run(self, droplets, params, variant, random_state):
        values = droplets.values(params.unit)
        index = candidates(values, params.min_counts)
        x = transform(values[index], params.transform)
        cut = threshold_otsu(x, params.bins)
        is_cell = np.zeros(values.size, dtype=bool)
        is_cell[index[x > cut]] = True
        return result(
            type(self), variant, params, droplets, is_cell,
            thresholds={params.unit: untransform(cut, params.transform), "transformed": float(cut)},
            diagnostics=_histogram_diagnostics(x, params.bins, _scale_label(params), [(cut, "Otsu's cut")]),
        )


@register
class MultiOtsu(Method):
    info = MethodInfo(
        id="multiotsu",
        title="Multi-Otsu thresholds",
        tool=_GENERIC,
        category="density",
        summary="Otsu's method with several classes (scikit-image's threshold_multiotsu) on the log totals of "
                "every barcode with a count; the top class is the cells.",
        docs="methods/generic/multiotsu.md",
        counts="either",
        references=(
            Reference(
                "Liao, Chen & Chung 2001, A fast algorithm for multilevel thresholding, Journal of Information "
                "Science and Engineering 17:713",
                doi="10.6688/JISE.2001.17.5.1",
            ),
            Reference(
                "Heiser, Wang, Chen, Hughes & Lau 2021, Automated quality control and cell identification of "
                "droplet-based single-cell data using dropkick, Genome Research 31:1742",
                doi="10.1101/gr.271908.120",
            ),
        ),
        source="https://github.com/scikit-image/scikit-image/blob/v0.26.0/src/skimage/filters/_multiotsu.pyx",
        original_license="BSD-3-Clause",
    )

    @dataclass
    class Params:
        classes: int = param(3, "How many classes", original="classes", minimum=2, maximum=4)
        min_counts: float = param(1.0, "Candidates have at least this value; the others are empty", minimum=0)
        transform: str = param("log10", "The scale the histogram is on", choices=_TRANSFORMS)
        bins: int = param(256, "Histogram bins", original="nbins", minimum=2)
        unit: str = param("umis", "What the rule applies to", choices=_UNITS)

    def run(self, droplets, params, variant, random_state):
        values = droplets.values(params.unit)
        index = candidates(values, params.min_counts)
        x = transform(values[index], params.transform)
        cuts = threshold_multiotsu(x, params.classes, params.bins)
        level = np.digitize(x, bins=cuts)
        names = ["below floor"] + [f"class {i + 1}" for i in range(params.classes)]
        codes = np.zeros(values.size, dtype=np.int8)
        codes[index] = level + 1
        is_cell = codes == params.classes
        thresholds = {params.unit: untransform(float(cuts[-1]), params.transform)}
        thresholds.update({f"threshold_{i + 1}": float(c) for i, c in enumerate(cuts)})
        return result(
            type(self), variant, params, droplets, is_cell,
            classes=pd.Categorical.from_codes(codes, categories=names),
            thresholds=thresholds,
            diagnostics=_histogram_diagnostics(x, params.bins, _scale_label(params),
                                               [(c, f"threshold {i + 1}") for i, c in enumerate(cuts)]),
        )


def bandwidth_method(text: str):
    """`--bandwidth`: scott, silverman, or a positive factor."""
    if text in ("scott", "silverman"):
        return text
    try:
        factor = float(text)
    except ValueError:
        raise RuleError(f"--bandwidth must be scott, silverman or a number, not {text!r}") from None
    if not factor > 0:
        raise RuleError(f"--bandwidth must be positive, not {text}")
    return factor


def density_valley(x: np.ndarray, bandwidth) -> dict:
    """The kernel density of x on a 10,000-point grid, its peaks, and the lowest point between the two
    highest (the first such point)."""
    from scipy.signal import argrelextrema
    from scipy.stats import gaussian_kde

    if x.size < 2 or np.all(x == x[0]):
        raise RuleError("a kernel density needs at least two different values")
    kde = gaussian_kde(x, bw_method=bandwidth)
    grid = np.linspace(x.min(), x.max(), GRID_POINTS)
    density = kde(grid)
    peaks = argrelextrema(density, np.greater)[0]
    if peaks.size < 2:
        raise RuleError(f"the density has {peaks.size} peak: no valley between two (try another --bandwidth or --min-counts)")
    highest = peaks[np.lexsort((peaks, -density[peaks]))[:2]]
    low, high = int(highest.min()), int(highest.max())
    valley = low + int(np.argmin(density[low:high + 1]))
    return {"grid": grid, "density": density, "peaks": peaks, "bounds": (low, high), "valley": valley}


@register
class KdeValley(Method):
    info = MethodInfo(
        id="kde-valley",
        title="Kernel-density valley",
        tool=_GENERIC,
        category="density",
        summary="The lowest point of a kernel density of log totals between its two highest peaks; barcodes "
                "above it are cells.",
        docs="methods/generic/kde-valley.md",
        counts="either",
        references=(
            Reference(
                "Smith, Heger & Sudbery 2017, UMI-tools: modeling sequencing errors in Unique Molecular "
                "Identifiers to improve quantification accuracy, Genome Research 27:491",
                doi="10.1101/gr.209601.116",
            ),
        ),
    )

    @dataclass
    class Params:
        min_counts: float = param(100.0, "Candidates have at least this value; the others are empty", minimum=0)
        bandwidth: str = param("scott", "Bandwidth: scott, silverman, or a factor such as 0.1", original="bw_method")
        unit: str = param("umis", "What the rule applies to", choices=_UNITS)

    def run(self, droplets, params, variant, random_state):
        bandwidth = bandwidth_method(params.bandwidth)
        values = droplets.values(params.unit)
        index = candidates(values, params.min_counts)
        x = clog10(values[index])
        found = density_valley(x, bandwidth)
        cut_log = float(found["grid"][found["valley"]])
        cut = 10.0 ** cut_log
        is_cell = values > cut
        thin = slice(None, None, GRID_POINTS // 1000)
        n_cells = int(is_cell.sum())
        return result(
            type(self), variant, params, droplets, is_cell,
            thresholds={params.unit: cut, "log10": cut_log},
            diagnostics={
                "grid": found["grid"][thin], "density": found["density"][thin],
                "peaks": found["grid"][found["peaks"]],
                "valleys": np.array([cut_log]), "valley_passing": np.array([n_cells]), "chosen": cut_log,
                "density_title": f"Density of log10 {_UNIT_LABEL[params.unit]} (bandwidth {params.bandwidth})",
            },
        )
