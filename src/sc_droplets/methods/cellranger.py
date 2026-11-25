"""Cell Ranger's cell calling, release by release.

Two methods:

- `cellranger`: what `cellranger count` / `multi` called by default. 2.x: OrdMag.
  3.0 onwards: OrdMag, then an EmptyDrops-style test of the barcodes below it.
  `--force-cells` switches to the top N, as Cell Ranger does.
- `cellranger-ordmag`: OrdMag alone (Cell Ranger's `ordmag` mode, and what it
  runs for antibody-only libraries).

Provenance, by release:

- 2.x and 3.x: a port of Cell Ranger 3.0.2 (`lib/python/cellranger/cell_calling.py`,
  `stats.py`, `sgt.py`, `mro/stages/counter/filter_barcodes`), which 10x
  released under the MIT licence (Copyright (c) 2015-2018 10X Genomics). Its
  OrdMag is the 2.x algorithm.
- 4.0 to 7.1: no public source; refused (documented gap).
- 7.2 onwards: written from the specification on the method pages, not from
  10x's code, which is under a licence that does not allow copying.

Pages: docs/methods/cellranger/.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import scipy.sparse as sp

from sc_droplets.methods._common import result
from sc_droplets.primitives.goodturing import GoodTuringError, sgt_cellranger
from sc_droplets.primitives.legacysort import legacy_argsort
from sc_droplets.primitives.multinomial import multinomial_loglik_rows, scipy_adjusted_profile
from sc_droplets.primitives.pvalues import bh_cellranger, monte_carlo_pvalues
from sc_droplets.primitives.rounding import round_half_away
from sc_droplets.primitives.simulate import cellranger_v1
from sc_droplets.registry import Gap, Method, MethodInfo, Reference, Variant, param, register

TOOL = "Cell Ranger"
LATEST = "10.1.0"
REFERENCES = (
    Reference(
        "Zheng et al. 2017, Massively parallel digital transcriptional profiling of single cells, Nat Commun 8:14049",
        doi="10.1038/ncomms14049",
    ),
    Reference(
        "Lun et al. 2019, EmptyDrops: distinguishing cells from empty droplets in droplet-based "
        "single-cell RNA sequencing data, Genome Biol 20:63",
        doi="10.1186/s13059-019-1662-y",
    ),
    Reference(
        "10x Genomics, Cell Ranger: calling cell barcodes",
        url="https://www.10xgenomics.com/support/software/cell-ranger/latest/algorithms-overview/cr-gex-algorithm",
    ),
)
GAPS = (
    Gap(
        "4.0",
        "7.2",
        "Cell Ranger 4.0 to 7.1 published no source; their release notes describe changes "
        "(targeted calling, per-species EmptyDrops, changed parameters) that cannot be reproduced exactly.",
    ),
)
CLASSES = ["initial", "rescued", "candidate", "background"]

# Constants of Cell Ranger 2.x / 3.x (cellranger/constants.py, cell_calling.py)
ORDMAG_QUANTILE = 0.99
ORDMAG_BOOTSTRAPS = 100
DEFAULT_RECOVERED_CELLS = 3000
N_CANDIDATE_BARCODES = 20000
N_PARTITIONS = 90000
MIN_UMIS = 500
MIN_UMI_FRAC_OF_MEDIAN = 0.01
MAX_ADJ_PVALUE = 0.01
NUM_SIMS = 10000


# ---------------------------------------------------------------------------
# Matrix structure: GEM groups and genomes
# ---------------------------------------------------------------------------


def gem_groups(barcodes: np.ndarray) -> np.ndarray:
    """The GEM group of each barcode: the integer after its last '-', or 1."""
    groups = np.ones(len(barcodes), dtype=np.int64)
    for i, barcode in enumerate(barcodes):
        head, sep, tail = str(barcode).rpartition("-")
        if sep and tail.isdigit():
            groups[i] = int(tail)
    return groups


def genomes(features: pd.DataFrame) -> List[Tuple[str, np.ndarray]]:
    """Each genome and its features, in order of first appearance (Cell Ranger's order)."""
    if "genome" not in features.columns:
        return [("", np.ones(len(features), dtype=bool))]
    values = features["genome"].astype(str).to_numpy()
    names = list(dict.fromkeys(values.tolist()))
    return [(name, values == name) for name in names]


def _row_sums(matrix: sp.csr_matrix) -> np.ndarray:
    return np.asarray(matrix.sum(axis=1), dtype=np.int64).ravel()


# ---------------------------------------------------------------------------
# OrdMag, Cell Ranger 2.x / 3.x
# ---------------------------------------------------------------------------


def _find_within_ordmag_v2(x: np.ndarray, baseline_index: int) -> int:
    ascending = np.sort(x)
    baseline = ascending[-baseline_index]  # the b-th largest; b == 0 wraps to the smallest
    cutoff = max(1, round_half_away(0.1 * baseline))
    return int(len(x) - np.searchsorted(ascending, cutoff))


def ordmag_v2(counts: np.ndarray, recovered_cells: int, state: np.random.RandomState) -> Tuple[np.ndarray, Dict]:
    """Cell Ranger 2.x/3.x `filter_cellular_barcodes_ordmag`: indices of the cells, and metrics."""
    nonzero = counts[counts > 0]
    if nonzero.size == 0:
        return np.zeros(0, dtype=np.int64), {"n": 0}
    baseline_index = int(round_half_away(float(recovered_cells) * (1 - ORDMAG_QUANTILE)))
    baseline_index = min(baseline_index, len(nonzero) - 1)
    top = np.array(
        [_find_within_ordmag_v2(state.choice(nonzero, len(nonzero)), baseline_index) for _ in range(ORDMAG_BOOTSTRAPS)]
    )
    n = int(round_half_away(float(np.mean(top))))
    chosen = np.sort(legacy_argsort(counts)[::-1][0:n])
    ascending = np.sort(nonzero)
    baseline = int(ascending[-baseline_index])
    return chosen, {
        "n": n,
        "bootstrap": top,
        "recovered_cells": recovered_cells,
        "baseline_rank": baseline_index,
        "baseline": baseline,
        "cutoff": max(1, int(round_half_away(0.1 * baseline))),
    }


def top_n_v2(counts: np.ndarray, force: int) -> np.ndarray:
    """Cell Ranger 2.x/3.x `filter_cellular_barcodes_fixed_cutoff`."""
    top = min(force, int((counts > 0).sum()))
    return np.sort(legacy_argsort(counts)[::-1][0:top])


# ---------------------------------------------------------------------------
# EmptyDrops step, Cell Ranger 3.x
# ---------------------------------------------------------------------------


def _profile_v3(counts: sp.csr_matrix, ambient_rows: np.ndarray, features: np.ndarray) -> np.ndarray:
    """`estimate_profile_sgt`: Good-Turing proportions, P0 shared by the zero features."""
    profile = np.asarray(counts[ambient_rows][:, features].sum(axis=0), dtype=np.int64).ravel()
    zero = np.flatnonzero(profile == 0)
    nonzero = np.flatnonzero(profile)
    smoothed, p0 = sgt_cellranger(profile[nonzero])
    if zero.size == 0:
        # Cell Ranger 3.0.2 divides P0 by zero here and then fails its sum-to-one check.
        raise GoodTuringError("every feature has ambient counts, which Cell Ranger 3.0.2 fails on")
    out = np.repeat(p0 / len(zero), len(features))
    out[nonzero] = smoothed
    return out


def nonambient_v3(
    counts: sp.csr_matrix, initial: np.ndarray, state: np.random.RandomState, num_sims: int = NUM_SIMS
) -> Optional[Dict]:
    """Cell Ranger 3.x `find_nonambient_barcodes` on one GEM group (rows) and all GEX features.

    Returns None where Cell Ranger stops early (no ambient barcodes, Good-Turing
    not applicable, no candidates).
    """
    umis = _row_sums(counts)
    order = legacy_argsort(umis)
    empty = np.sort(order[::-1][N_PARTITIONS // 2 : N_PARTITIONS])
    use = np.intersect1d(empty, np.flatnonzero(umis), assume_unique=True)
    if use.size == 0:
        return None
    features = np.flatnonzero(np.asarray(counts.sum(axis=0)).ravel())
    try:
        profile = _profile_v3(counts, use, features)
    except GoodTuringError as error:
        return {"stopped": str(error)}
    if not initial.any():
        return None
    median = float(np.median(umis[initial]))
    min_umis = int(max(MIN_UMIS, round_half_away(float(np.ceil(median * MIN_UMI_FRAC_OF_MEDIAN)))))
    masked = initial | (umis < min_umis)
    n_unmasked = int((~masked).sum())
    filled = np.where(masked, np.iinfo(np.int64).max, umis)  # numpy.ma sorts masked entries last
    candidates = legacy_argsort(filled)[0:n_unmasked][-N_CANDIDATE_BARCODES:]
    if candidates.size == 0:
        return None
    adjusted = scipy_adjusted_profile(profile)
    rows = counts[candidates][:, features]
    with np.errstate(divide="ignore"):
        observed = multinomial_loglik_rows(rows, np.log(adjusted))
    totals = umis[candidates]
    n_lower = cellranger_v1(profile, totals, observed, num_sims, state)
    pvalues = monte_carlo_pvalues(n_lower, num_sims)
    adjusted_p = bh_cellranger(pvalues)
    return {
        "candidates": candidates,
        "loglik": observed,
        "pvalue": pvalues,
        "fdr": adjusted_p,
        "nonambient": adjusted_p <= MAX_ADJ_PVALUE,
        "min_umis": min_umis,
        "median_initial": median,
        "ambient_rows": use,
        "profile": profile,
        "features": features,
    }


# ---------------------------------------------------------------------------
# Methods
# ---------------------------------------------------------------------------

_V2 = Variant(
    "v2",
    since="2.0",
    until="4.0",
    summary="OrdMag: 100 bootstraps of the count of barcodes within 10x of the 30th-largest (for 3,000 "
    "expected cells); legacy numpy RNG seeded 0; Python 2 rounding",
)


@register
class OrdMag(Method):
    info = MethodInfo(
        id="cellranger-ordmag",
        title="Cell Ranger OrdMag",
        tool=TOOL,
        category="rank-threshold",
        summary="Barcodes within an order of magnitude of a robust maximum: the 99th percentile of the top "
        "expected cells ÷ 10, averaged over 100 bootstraps.",
        docs="methods/cellranger/ordmag.md",
        deterministic=True,
        variants=(_V2,),
        gaps=GAPS,
        latest_known=LATEST,
        references=REFERENCES,
        source="https://github.com/lfperales/cellranger/blob/3.0.2/lib/python/cellranger/stats.py",
        original_license="MIT (3.0.2); 10x licence (7.2+)",
        default_seed=0,
    )

    @dataclass
    class Params:
        expect_cells: Optional[int] = param(None, "Expected number of cells (default 3,000 per GEM group)", original="--expect-cells", minimum=1)

    def run(self, droplets, params, variant, random_state):
        state = np.random.RandomState(0 if random_state is None else random_state)
        groups = gem_groups(droplets.barcodes)
        unique = sorted(set(groups.tolist()))
        recovered = int(float(params.expect_cells) / len(unique)) if params.expect_cells else DEFAULT_RECOVERED_CELLS
        cells = np.zeros(droplets.n_barcodes, dtype=bool)
        diagnostics: Dict = {"groups": {}}
        for genome, feature_mask in genomes(droplets.features):
            sub = droplets.counts[:, np.flatnonzero(feature_mask)]
            for group in unique:
                rows = np.flatnonzero(groups == group)
                chosen, info = ordmag_v2(_row_sums(sub[rows]), recovered, state)
                cells[rows[chosen]] = True
                diagnostics["groups"][f"{genome or 'all'}-{group}"] = {k: v for k, v in info.items()}
        first = next(iter(diagnostics["groups"].values()))
        return result(
            type(self),
            variant,
            params,
            droplets,
            cells,
            thresholds={"umis": float(first.get("cutoff", 0))},
            diagnostics=diagnostics,
        )


@register
class CellRanger(Method):
    info = MethodInfo(
        id="cellranger",
        title="Cell Ranger cell calling",
        tool=TOOL,
        category="statistical-test",
        summary="What cellranger count called: OrdMag (2.x), then from 3.0 an EmptyDrops-style test that "
        "rescues barcodes below the OrdMag cut whose counts differ from the ambient profile.",
        docs="methods/cellranger/cellranger.md",
        deterministic=False,
        variants=(
            Variant(
                "v2",
                since="2.0",
                until="3.0",
                summary="OrdMag only, per genome and GEM group",
                settings={"nonambient": False},
            ),
            Variant(
                "v3",
                since="3.0",
                until="4.0",
                summary="OrdMag, then EmptyDrops on ranks 45,000–90,000 as ambient, the top 20,000 candidates "
                "with at least max(500, 1% of the median cell), 10,000 simulations, BH FDR ≤ 0.01",
                settings={"nonambient": True},
            ),
        ),
        gaps=GAPS,
        latest_known=LATEST,
        references=REFERENCES,
        source="https://github.com/lfperales/cellranger/blob/3.0.2/lib/python/cellranger/cell_calling.py",
        original_license="MIT (3.0.2); 10x licence (7.2+)",
        default_seed=0,
    )

    @dataclass
    class Params:
        expect_cells: Optional[int] = param(None, "Expected number of cells (default 3,000 per GEM group)", original="--expect-cells", minimum=1)
        force_cells: Optional[int] = param(None, "Call exactly this many cells (the top N) instead", original="--force-cells", minimum=1)

    def run(self, droplets, params, variant, random_state):
        state = np.random.RandomState(0 if random_state is None else random_state)
        notes: List[str] = []
        n_all = droplets.n_barcodes
        if "Antibody Capture" in set(droplets.obs.columns) or (
            droplets.feature_types is not None and "Antibody Capture" in droplets.feature_types
        ):
            notes.append("Cell Ranger 3.x first removed antibody aggregates using read counts a matrix does not hold")
        groups = gem_groups(droplets.barcodes)
        unique = sorted(set(groups.tolist()))
        recovered = int(float(params.expect_cells) / len(unique)) if params.expect_cells else DEFAULT_RECOVERED_CELLS
        initial = np.zeros(n_all, dtype=bool)
        diagnostics: Dict = {"ordmag": {}}
        for genome, feature_mask in genomes(droplets.features):
            sub = droplets.counts[:, np.flatnonzero(feature_mask)]
            for group in unique:
                rows = np.flatnonzero(groups == group)
                bc_counts = _row_sums(sub[rows])
                if params.force_cells:
                    chosen = top_n_v2(bc_counts, int(float(params.force_cells) / len(unique)))
                    info = {"n": int(chosen.size)}
                else:
                    chosen, info = ordmag_v2(bc_counts, recovered, state)
                initial[rows[chosen]] = True
                diagnostics["ordmag"][f"{genome or 'all'}-{group}"] = info

        rescued = np.zeros(n_all, dtype=bool)
        classes = np.full(n_all, "background", dtype=object)
        pvalue = np.full(n_all, np.nan)
        fdr = np.full(n_all, np.nan)
        loglik = np.full(n_all, np.nan)
        thresholds: Dict[str, float] = {}
        if variant.settings["nonambient"] and not params.force_cells:
            for group in unique:
                rows = np.flatnonzero(groups == group)
                found = nonambient_v3(droplets.counts[rows], initial[rows], state)
                if found is None:
                    notes.append(f"GEM group {group}: no barcodes to test against the ambient profile")
                    continue
                if "stopped" in found:
                    notes.append(f"GEM group {group}: EmptyDrops step skipped: {found['stopped']}")
                    continue
                where = rows[found["candidates"]]
                classes[where] = "candidate"
                pvalue[where], fdr[where], loglik[where] = found["pvalue"], found["fdr"], found["loglik"]
                rescued[where[found["nonambient"]]] = True
                thresholds["candidate_umis"] = float(found["min_umis"])
                diagnostics.setdefault("nonambient", {})[str(group)] = {
                    "min_umis": found["min_umis"],
                    "median_initial": found["median_initial"],
                    "n_candidates": int(found["candidates"].size),
                    "n_ambient": int(found["ambient_rows"].size),
                    "ambient_ranks": [N_PARTITIONS // 2, N_PARTITIONS],
                }
            thresholds["fdr"] = MAX_ADJ_PVALUE
        classes[initial] = "initial"
        classes[rescued & ~initial] = "rescued"
        cells = initial | rescued
        if initial.any():
            thresholds["umis"] = float(droplets.totals[initial].min())
        return result(
            type(self),
            variant,
            params,
            droplets,
            cells,
            classes=pd.Categorical(classes, categories=CLASSES),
            columns={"pvalue": pvalue, "fdr": fdr, "loglik": loglik} if variant.settings["nonambient"] else {},
            thresholds=thresholds,
            diagnostics=diagnostics,
            notes=notes,
        )
