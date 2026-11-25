"""STARsolo's cell filters (`--soloCellFilter`) and its `soloBasicCellFilter.awk`.

Ports of STAR (github.com/alexdobin/STAR, MIT licence, Copyright (c) 2019
Alexander Dobin): `source/SoloFeature_cellFiltering.cpp`,
`source/SoloFeature_emptyDrops_CR.cpp`, `source/SimpleGoodTuring/sgt.h`
(Copyright (c) David Elworthy 2004) and `extras/scripts/soloBasicCellFilter.awk`,
at the releases each variant names.

STAR counts only barcodes with at least one count ("detected", nCB), in the
order of its barcode list, and rounds each matrix entry to an integer (half
away from zero) when it loads a matrix. Pages: docs/methods/starsolo/.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from typing import Optional

import numpy as np
import pandas as pd
import scipy.sparse as sp

from sc_droplets.methods._common import result
from sc_droplets.primitives.goodturing import SimpleGoodTuringStar
from sc_droplets.primitives.pvalues import bh_starsolo, monte_carlo_pvalues
from sc_droplets.primitives.rounding import round_half_away
from sc_droplets.primitives.simulate import starsolo as simulate_starsolo
from sc_droplets.registry import Method, MethodInfo, Reference, Variant, param, register

TOOL = "STARsolo"
SOURCE = "https://github.com/alexdobin/STAR/blob/2.7.11b/source/SoloFeature_cellFiltering.cpp"
LATEST = "2.7.11b"
REFERENCES = (
    Reference(
        "Kaminow, Yunusov & Dobin 2021, STARsolo: accurate, fast and versatile mapping/quantification "
        "of single-cell and single-nucleus RNA-seq data, bioRxiv",
        doi="10.1101/2021.05.05.442755",
    ),
)
CLASSES = ["initial", "rescued", "candidate", "background"]


class _Detected:
    """STAR's view of a matrix: detected barcodes, their integer totals, sorted totals."""

    def __init__(self, droplets):
        counts = droplets.counts
        if counts.dtype.kind == "f":
            counts = counts.copy()
            counts.data = np.floor(counts.data + 0.5)  # std::round on non-negative values
        totals = np.asarray(counts.sum(axis=1)).ravel()
        self.index = np.flatnonzero(totals > 0)
        self.counts = sp.csr_matrix(counts[self.index])
        self.totals = totals[self.index].astype(np.int64)
        self.sorted = np.sort(self.totals)[::-1]

    @property
    def n(self) -> int:
        return int(self.index.size)

    def to_all(self, mask: np.ndarray, n_all: int) -> np.ndarray:
        out = np.zeros(n_all, dtype=bool)
        out[self.index[mask]] = True
        return out


def _robust_max_threshold(detected: _Detected, expected: float, percentile: float, ratio: float, rounding: str, floor: bool):
    """The CellRanger2.2 cut: the robust maximum (31st largest by default) ÷ ratio."""
    if rounding == "add-half":
        index = int(expected * (1.0 - percentile) + 0.5)  # (uint64)(x + 0.5)
        robust_max = int(detected.sorted[min(detected.n - 1, index)])
        minimum = int(robust_max / ratio + 0.5)
    else:
        index = int(round_half_away(expected * (1.0 - percentile)))
        robust_max = int(detected.sorted[min(detected.n - 1, index)])
        minimum = int(round_half_away(robust_max / ratio))
    if floor:
        minimum = max(minimum, 1)
    return robust_max, minimum


_CR22_VARIANTS = (
    Variant(
        "v1",
        since="2.7.3a",
        until="2.7.5a",
        summary="robust max ÷ ratio, rounded by adding 0.5; no floor, so a cut of 0 keeps every detected barcode",
        settings={"rounding": "add-half", "floor": False},
    ),
    Variant(
        "v2",
        since="2.7.5a",
        until="2.7.8a",
        summary="as v1, with the cut floored at 1",
        settings={"rounding": "add-half", "floor": True},
    ),
    Variant(
        "v3",
        since="2.7.8a",
        summary="as v2, rounding with C round(); the three numbers become optional",
        settings={"rounding": "round", "floor": True},
    ),
)


@register
class CellRanger22(Method):
    info = MethodInfo(
        id="starsolo-cellranger2.2",
        title="STARsolo CellRanger2.2",
        tool=TOOL,
        category="rank-threshold",
        summary="STARsolo's copy of Cell Ranger 2's OrdMag without the bootstrap: barcodes with at least "
        "the 31st-largest total ÷ 10 are cells. STARsolo's default filter.",
        docs="methods/starsolo/cellranger2-2.md",
        variants=_CR22_VARIANTS,
        latest_known=LATEST,
        references=REFERENCES,
        source=SOURCE,
        original_license="MIT",
        feature_types=None,
    )

    @dataclass
    class Params:
        expected_cells: float = param(3000.0, "Expected number of cells", original="nExpectedCells", minimum=0)
        max_percentile: float = param(0.99, "Percentile of the top expected cells taken as the robust maximum", original="maxPercentile", minimum=0, maximum=1)
        max_min_ratio: float = param(10.0, "Robust maximum ÷ this is the cut", original="maxMinRatio", minimum=0)

    def run(self, droplets, params, variant, random_state):
        detected = _Detected(droplets)
        if detected.n == 0:
            return result(type(self), variant, params, droplets, np.zeros(droplets.n_barcodes, bool), notes=["no barcode has counts"])
        robust_max, minimum = _robust_max_threshold(
            detected, params.expected_cells, params.max_percentile, params.max_min_ratio, **variant.settings
        )
        cells = detected.to_all(detected.totals >= minimum, droplets.n_barcodes)
        return result(
            type(self),
            variant,
            params,
            droplets,
            cells,
            thresholds={"robust_max": float(robust_max), "umis": float(minimum)},
            diagnostics={"robust_max": robust_max, "cut": minimum, "n_detected": detected.n},
        )


@register
class TopCells(Method):
    info = MethodInfo(
        id="starsolo-topcells",
        title="STARsolo TopCells",
        tool=TOOL,
        category="rank-threshold",
        summary="Barcodes with at least the total of the (N+1)-th detected barcode are cells: "
        "N + 1 cells, plus ties, as STAR indexes from 0.",
        docs="methods/starsolo/topcells.md",
        variants=(
            Variant(
                "v1",
                since="2.7.3a",
                until="2.7.5a",
                summary="indexes max(nCB - 1, N): keeps every detected barcode while N < nCB - 1; no floor",
                settings={"buggy": True, "floor": False},
            ),
            Variant(
                "v2",
                since="2.7.5a",
                until="2.7.6a",
                summary="as v1, with the cut floored at 1",
                settings={"buggy": True, "floor": True},
            ),
            Variant(
                "v3",
                since="2.7.6a",
                summary="indexes min(nCB - 1, N): the (N+1)-th largest total and every barcode with as many",
                settings={"buggy": False, "floor": True},
            ),
        ),
        latest_known=LATEST,
        references=REFERENCES,
        source=SOURCE,
        original_license="MIT",
        feature_types=None,
    )

    @dataclass
    class Params:
        n_cells: int = param(help="Number of cells", original="TopCells <N>", minimum=0)

    def run(self, droplets, params, variant, random_state):
        detected = _Detected(droplets)
        if detected.n == 0:
            return result(type(self), variant, params, droplets, np.zeros(droplets.n_barcodes, bool), notes=["no barcode has counts"])
        if variant.settings["buggy"]:
            index = max(detected.n - 1, params.n_cells)
            if index >= detected.n:
                raise ValueError(
                    f"STAR {variant.versions} reads past the end of its sorted totals here "
                    f"(index {index} of {detected.n}); its result is undefined"
                )
        else:
            index = min(detected.n - 1, params.n_cells)
        minimum = int(detected.sorted[index])
        if variant.settings["floor"]:
            minimum = max(minimum, 1)
        cells = detected.to_all(detected.totals >= minimum, droplets.n_barcodes)
        return result(type(self), variant, params, droplets, cells, thresholds={"umis": float(minimum)}, diagnostics={"index": index})


@register
class EmptyDropsCR(Method):
    info = MethodInfo(
        id="starsolo-emptydrops-cr",
        title="STARsolo EmptyDrops_CR",
        tool=TOOL,
        category="statistical-test",
        summary="CellRanger2.2's cut, then STAR's own EmptyDrops: barcodes below the cut that differ from the "
        "ambient profile of ranks 45,000 to 90,000 are rescued.",
        docs="methods/starsolo/emptydrops-cr.md",
        deterministic=False,
        variants=(
            Variant(
                "v1",
                since="2.7.8a",
                summary="CellRanger2.2 (v3), then a multinomial test with STAR's Good-Turing profile and one "
                "std::mt19937 per simulation",
            ),
        ),
        latest_known=LATEST,
        references=REFERENCES,
        source="https://github.com/alexdobin/STAR/blob/2.7.11b/source/SoloFeature_emptyDrops_CR.cpp",
        original_license="MIT",
        feature_types=None,
        default_seed=19760110,
    )

    @dataclass
    class Params:
        expected_cells: float = param(3000.0, "Expected number of cells", original="nExpectedCells", minimum=0)
        max_percentile: float = param(0.99, "Percentile of the top expected cells taken as the robust maximum", original="maxPercentile", minimum=0, maximum=1)
        max_min_ratio: float = param(10.0, "Robust maximum ÷ this is the initial cut", original="maxMinRatio", minimum=0)
        ind_min: int = param(45000, "First rank (0-based) of the ambient droplets", original="indMin", minimum=0)
        ind_max: int = param(90000, "Rank after the last ambient droplet", original="indMax", minimum=0)
        umi_min: int = param(500, "Fewest UMIs a candidate may have", original="umiMin", minimum=0)
        umi_min_frac_median: float = param(0.01, "Fraction of the initial cells' (upper) median a candidate needs", original="umiMinFracMedian", minimum=0)
        cand_max_n: int = param(20000, "Most candidates tested", original="candMaxN", minimum=0)
        fdr: float = param(0.01, "Largest adjusted p-value called a cell", original="FDR", minimum=0)
        sim_n: int = param(10000, "Number of simulations", original="simN", minimum=1)

    def run(self, droplets, params, variant, random_state):
        detected = _Detected(droplets)
        n_all = droplets.n_barcodes
        classes = np.full(n_all, "background", dtype=object)
        notes = []
        if detected.n == 0:
            return result(type(self), variant, params, droplets, np.zeros(n_all, bool), notes=["no barcode has counts"])
        robust_max, minimum = _robust_max_threshold(
            detected, params.expected_cells, params.max_percentile, params.max_min_ratio, rounding="round", floor=True
        )
        initial = detected.totals >= minimum
        n_simple = int(initial.sum())
        called = initial.copy()
        pvalue = np.full(n_all, np.nan)
        padj = np.full(n_all, np.nan)
        loglik = np.full(n_all, np.nan)
        thresholds = {"robust_max": float(robust_max), "umis": float(minimum), "fdr": params.fdr}
        diagnostics = {"robust_max": robust_max, "cut": minimum, "n_initial": n_simple, "n_detected": detected.n}

        tested = self._test(detected, n_simple, params, random_state, notes, diagnostics)
        if tested is not None:
            cand, observed, p, adjusted = tested
            called[cand[adjusted <= params.fdr]] = True
            where = detected.index[cand]
            pvalue[where], padj[where], loglik[where] = p, adjusted, observed
            thresholds["candidate_umis"] = float(diagnostics["min_umi"])
            classes[where] = "candidate"
        classes[detected.index[initial]] = "initial"
        classes[detected.index[called & ~initial]] = "rescued"
        cells = detected.to_all(called, n_all)
        return result(
            type(self),
            variant,
            params,
            droplets,
            cells,
            classes=pd.Categorical(classes, categories=CLASSES),
            columns={"pvalue": pvalue, "fdr": padj, "loglik": loglik},
            thresholds=thresholds,
            diagnostics=diagnostics,
            notes=notes,
        )

    @staticmethod
    def _test(detected, n_simple, params, random_state, notes, diagnostics):
        """STAR's emptyDrops_CR(): None when it stops early, else the candidates' results."""
        n_cb = detected.n
        if n_cb <= params.ind_min:
            notes.append(
                f"only {n_cb} detected barcodes, not more than ind_min={params.ind_min}: STAR skips EmptyDrops_CR"
            )
            return None
        counts = detected.counts
        n_features = counts.shape[1]
        detected_features = np.zeros(n_features, dtype=bool)
        detected_features[np.unique(counts.indices[counts.data > 0])] = True
        # descending by count, ascending by index
        order = np.lexsort((np.arange(n_cb), -detected.totals))
        ambient_rows = order[params.ind_min : min(n_cb, params.ind_max)]
        ambient = np.asarray(counts[ambient_rows].sum(axis=0)).ravel().astype(np.int64)
        frequency = Counter(ambient.tolist())
        if len(frequency) <= 1:
            notes.append("the ambient droplets hold no counts: STAR skips EmptyDrops_CR")
            return None
        frequency[0] -= n_features - int(detected_features.sum())
        max_freq = max(frequency)
        sgt = SimpleGoodTuringStar()
        for value, count in sorted(frequency.items()):
            if value != 0:
                sgt.add(value, count)
        if not sgt.analyse():
            notes.append("fewer than 5 distinct ambient counts: STAR's Good-Turing gives no estimate and its result is undefined")
            return None
        by_frequency = np.array([sgt.estimate(f) for f in range(max_freq + 1)], dtype=np.float64)
        with np.errstate(divide="ignore", invalid="ignore"):
            by_frequency[0] = by_frequency[0] / frequency[0]
        profile = np.zeros(n_features)
        profile[detected_features] = by_frequency[ambient[detected_features]]
        norm = np.cumsum(profile)[-1]
        nonzero = profile > 0
        profile[nonzero] = profile[nonzero] / norm
        log_profile = np.zeros(n_features)
        log_profile[nonzero] = np.log(profile[nonzero])

        sorted_totals = detected.sorted
        min_umi = int(params.umi_min_frac_median * sorted_totals[n_simple // 2])
        min_umi = max(params.umi_min, min_umi)
        first = n_simple
        last = first
        stop = min(first + params.cand_max_n, n_cb)
        while last < stop and detected.totals[order[last]] >= min_umi:
            last += 1
        diagnostics["min_umi"] = min_umi
        diagnostics["ambient_ranks"] = [params.ind_min, min(n_cb, params.ind_max)]
        if last == first:
            notes.append("no candidate barcodes")
            return None
        cand = order[first:last]
        cand_totals = detected.totals[cand]
        max_count = int(cand_totals[0])
        log_factorial = np.concatenate(([0.0, 0.0], np.cumsum(np.log(np.arange(2, max_count + 1, dtype=np.float64)))))
        rows = counts[cand]
        entry_rows = np.repeat(np.arange(rows.shape[0]), np.diff(rows.indptr))
        data = rows.data.astype(np.int64)
        observed = (
            log_factorial[cand_totals]
            - np.bincount(entry_rows, weights=log_factorial[data], minlength=rows.shape[0])
            + np.bincount(entry_rows, weights=log_profile[rows.indices] * data, minlength=rows.shape[0])
        )
        seed = 19760110 if random_state is None else int(random_state)
        n_lower = simulate_starsolo(profile[nonzero], cand_totals, observed, params.sim_n, seed_base=seed)
        p = monte_carlo_pvalues(n_lower, params.sim_n)
        adjusted = bh_starsolo(p)
        diagnostics["n_candidates"] = int(cand.size)
        diagnostics["max_candidate_umis"] = max_count
        diagnostics["profile_top"] = _top_features(profile, n=20)
        return cand, observed, p, adjusted


def _top_features(profile: np.ndarray, n: int) -> list:
    top = np.argsort(-profile, kind="stable")[:n]
    return [int(i) for i in top]


@register
class BasicAwk(Method):
    info = MethodInfo(
        id="starsolo-basic-awk",
        title="STARsolo soloBasicCellFilter.awk",
        tool=TOOL,
        category="rank-threshold",
        summary="The awk script shipped with STAR: the (maxCells·(1-maxPercentile)+1)-th largest total ÷ "
        "ratio, unrounded, or the top N with exactCells.",
        docs="methods/starsolo/basic-awk.md",
        variants=(Variant("v1", since="2.7.2a", summary="extras/scripts/soloBasicCellFilter.awk, unchanged since 2.7.2a"),),
        latest_known=LATEST,
        references=REFERENCES,
        source="https://github.com/alexdobin/STAR/blob/2.7.11b/extras/scripts/soloBasicCellFilter.awk",
        original_license="MIT",
        feature_types=None,
    )

    @dataclass
    class Params:
        max_cells: float = param(3000.0, "Expected number of cells", original="maxCells", minimum=0)
        max_percentile: float = param(0.99, "Percentile taken as the robust maximum", original="maxPercentile", minimum=0, maximum=1)
        max_min_ratio: float = param(10.0, "Robust maximum ÷ this is the cut", original="maxMinRatio", minimum=0)
        exact_cells: int = param(0, "If above 0, keep this many cells instead (plus one, plus ties)", original="exactCells", minimum=0)

    def run(self, droplets, params, variant, random_state):
        detected = _Detected(droplets)
        ascending = detected.sorted[::-1]
        length = detected.n

        def at(position: int) -> float:  # awk arrays are 1-based; a missing element is 0
            return float(ascending[position - 1]) if 1 <= position <= length else 0.0

        robust_max = at(-int((1 - params.max_percentile) * params.max_cells) + length)
        minimum = robust_max / params.max_min_ratio
        if params.exact_cells > 0:
            minimum = at(1) if length < params.exact_cells else at(length - params.exact_cells)
        cells = detected.to_all(detected.totals >= minimum, droplets.n_barcodes)
        if minimum <= 0:
            # awk compares a barcode with no counts as 0, which passes a cut of 0
            cells[:] = True
        return result(
            type(self),
            variant,
            params,
            droplets,
            cells,
            thresholds={"robust_max": robust_max, "umis": float(minimum)},
        )
