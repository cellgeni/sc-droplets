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
- 4.0 to 5.x: no public source. 10x's own outputs of these releases are reproduced
  exactly by the 3.x port with two changes seen in them (and in 7.2's specification): the
  simulation seeds a stream of its own, and from 5.0 antibody aggregates are removed first
  (the method page lists the datasets).
- 6.0 to 7.1 (and 1.x): no public source, and 10x's outputs match no rule exactly; refused
  (documented gaps).
- 7.2 to 10.1: written from the behavioural specification in
  docs/_research/cellranger-7.2-plus-spec.md and the method pages, not from
  10x's code, which is under a licence that does not allow copying.

Pages: docs/methods/cellranger/.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import scipy.sparse as sp

from sc_droplets.methods._common import result
from sc_droplets.primitives.goodturing import GoodTuringError, sgt_cellranger
from sc_droplets.primitives.legacysort import legacy_argsort
from sc_droplets.primitives.multinomial import multinomial_loglik_rows, scipy_adjusted_profile
from sc_droplets.primitives.pvalues import bh_cellranger, monte_carlo_pvalues
from sc_droplets.primitives.rounding import round_half_away
from sc_droplets.primitives.simulate import cellranger_v1, cellranger_v9
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
        "1.0",
        "2.0",
        "Cell Ranger 1.x published no source. Its paper describes OrdMag (10% of the top 1% of the "
        "expected cells), but not whether 1.x bootstrapped it as 2.x does.",
    ),
    Gap(
        "6.0",
        "7.2",
        "Cell Ranger 6.0 to 7.1 published no source, and no rule reproduces 10x's own outputs of these "
        "releases exactly: the closest differ by a few barcodes at the FDR cut.",
    ),
)
CLASSES = ["initial", "rescued", "candidate", "removed", "background"]
GEX = "Gene Expression"
ANTIBODY = "Antibody Capture"
ANTIGEN = "Antigen Capture"

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

# 7.2+ (from the specification)
Q = 1 - 0.99  # 0.010000000000000009, as Cell Ranger computes it
MIN_RECOVERED = 50
MAX_RECOVERED = 1 << 18
GRID_POINTS = 2000

# The 13 protein-coding mitochondrial genes of human and of mouse (8.0+ mito filter).
MITO_IDS = frozenset(
    [
        "ENSG00000198888", "ENSG00000198763", "ENSG00000198804", "ENSG00000198712", "ENSG00000228253",
        "ENSG00000198899", "ENSG00000198938", "ENSG00000198840", "ENSG00000212907", "ENSG00000198886",
        "ENSG00000198786", "ENSG00000198695", "ENSG00000198727",
        "ENSMUSG00000064341", "ENSMUSG00000064345", "ENSMUSG00000064351", "ENSMUSG00000064354",
        "ENSMUSG00000064356", "ENSMUSG00000064357", "ENSMUSG00000064358", "ENSMUSG00000064360",
        "ENSMUSG00000065947", "ENSMUSG00000064363", "ENSMUSG00000064367", "ENSMUSG00000064368",
        "ENSMUSG00000064370",
    ]
)


# ---------------------------------------------------------------------------
# Chemistry: ambient range and FDR (7.2+)
# ---------------------------------------------------------------------------

#: Chemistry name -> kit family. Cell Ranger matches descriptions, which differ
#: between releases for the same kit (SC3Pv3 became SC3Pv3-polyA in 9.0); the
#: family is what the per-release tables below are written in.
CHEMISTRY_FAMILY: Dict[str, str] = {
    **{n: "3p-v1v2" for n in ("SC3Pv1", "SC3Pv2")},
    **{n: "3p-v3" for n in ("SC3Pv3", "SC3Pv3-polyA", "SC3Pv3-CS1")},
    **{n: "3p-v3-oh" for n in ("SC3Pv3-OH", "SC3Pv3-polyA-OCM", "SC3Pv3-CS1-OCM")},
    **{n: "3p-v3-ht" for n in ("SC3Pv3HT", "SC3Pv3HT-polyA", "SC3Pv3HT-CS1")},
    "SC3Pv3LT": "3p-v3-lt",
    **{n: "3p-v4" for n in ("SC3Pv4", "SC3Pv4-polyA", "SC3Pv4-CS1")},
    **{n: "3p-v4-oh" for n in ("SC3Pv4-OH", "SC3Pv4-polyA-OCM", "SC3Pv4-CS1-OCM")},
    "SC3Pv4HT": "3p-v4-ht",
    **{n: "5p-v2" for n in ("SC5P-PE", "SC5P-R2", "SC5P-R1")},
    **{n: "5p-v2-oh" for n in ("SC5P-R2-OH", "SC5P-R2-OCM")},
    "SC5PHT": "5p-ht",
    "SC5PHT-v3": "5p-v3-ht",
    **{n: "5p-v3" for n in ("SC5P-PE-v3", "SC5P-R2-v3", "SC5P-R1-v3")},
    **{n: "5p-v3-oh" for n in ("SC5P-R2-OH-v3", "SC5P-R2-OCM-v3")},
    **{n: "5p-v3-ocm" for n in ("SC5P-PE-OCM-v3", "SC5P-R1-OCM-v3")},
    **{
        n: "flex-v1"
        for n in (
            "SFRP", "MFRP", "MFRP-R1", "MFRP-RNA", "MFRP-Ab", "MFRP-RNA-R1", "MFRP-Ab-R1", "MFRP-47",
            "MFRP-R1-48-uncollapsed", "MFRP-uncollapsed", "MFRP-Ab-R2pos50", "MFRP-CRISPR", "MFRP-96-R1",
            "MFRP-96-RNA-R2", "SFRP-no-trim-R2", "MFRP-R1-no-trim-R2",
        )
    },
}

#: Chemistry descriptions (as written in the HDF5 `chemistry_description`) -> name, every release.
DESCRIPTIONS: Dict[str, str] = {
    "Single Cell 3' v1": "SC3Pv1", "Single Cell 3' v2": "SC3Pv2", "Single Cell 3' v3": "SC3Pv3",
    "Single Cell 3' v3 (polyA)": "SC3Pv3-polyA", "Single Cell 3' v3 (CS1)": "SC3Pv3-CS1",
    "Single Cell 3' v3 OH": "SC3Pv3-OH", "Single Cell 3' v3 (polyA) OCM": "SC3Pv3-polyA-OCM",
    "Single Cell 3' v3 (CS1) OCM": "SC3Pv3-CS1-OCM", "Single Cell 3' v3 HT": "SC3Pv3HT",
    "Single Cell 3' v3 HT (polyA)": "SC3Pv3HT-polyA", "Single Cell 3' v3 HT (CS1)": "SC3Pv3HT-CS1",
    "Single Cell 3' v3 LT": "SC3Pv3LT", "Single Cell 3' v4": "SC3Pv4", "Single Cell 3' v4 (polyA)": "SC3Pv4-polyA",
    "Single Cell 3' v4 (CS1)": "SC3Pv4-CS1", "Single Cell 3' v4 OH": "SC3Pv4-OH",
    "Single Cell 3' v4 (polyA) OCM": "SC3Pv4-polyA-OCM", "Single Cell 3' v4 (CS1) OCM": "SC3Pv4-CS1-OCM",
    "Single Cell 3' v4 HT": "SC3Pv4HT", "Single Cell 5' PE": "SC5P-PE", "Single Cell 5' R2-only": "SC5P-R2",
    "Single Cell 5' R1-only": "SC5P-R1", "Single Cell 5' R2-only OH": "SC5P-R2-OH",
    "Single Cell 5' R2-only OCM": "SC5P-R2-OCM", "Single Cell 5' HT": "SC5PHT", "Single Cell 5' HT v3": "SC5PHT-v3",
    "Single Cell 5' PE v3": "SC5P-PE-v3", "Single Cell 5' R2-only v3": "SC5P-R2-v3",
    "Single Cell 5' R1-only v3": "SC5P-R1-v3", "Single Cell 5' R2-only OH v3": "SC5P-R2-OH-v3",
    "Single Cell 5' R2-only OCM v3": "SC5P-R2-OCM-v3", "Single Cell 5' PE OCM v3": "SC5P-PE-OCM-v3",
    "Single Cell 5' R1-only OCM v3": "SC5P-R1-OCM-v3", "Fixed RNA Profiling": "MFRP",
    "Fixed RNA Profiling (Singleplex)": "SFRP", "Fixed RNA Profiling (probe barcode on R1)": "MFRP-R1",
    "Flex Gene Expression": "MFRP-RNA", "Flex Gene Expression (probe barcode on R1)": "MFRP-RNA-R1",
    "GEM-X Flex v2 Gene Expression": "Flex-v2-RNA-R2", "GEM-X Flex v2": "Flex-v2-R1",
    "GEM-X Flex v2 singleplex": "Flex-v2-singleplex",
}

#: (range rule, FDR) per kit family and release; families missing from a
#: release take "default" (the release's code matches no group for them).
_D7, _D9 = "D7", "D9"
RANGE_RULES: Dict[str, Dict[str, Tuple[str, float]]] = {
    "7.2": {"3p-v3-lt": ("LT", 0.01), "default": (_D7, 0.01)},
    "8.0.0": {
        "3p-v3-lt": ("LT", 0.01),
        **{f: ("V8", 0.001) for f in ("3p-v4", "3p-v4-oh", "5p-v3", "5p-v3-oh")},
        "default": (_D7, 0.01),
    },
    "8.0.1": {
        "3p-v3-lt": ("LT", 0.01),
        "3p-v3-ht": ("HT8", 0.01),
        "5p-ht": ("HT8", 0.01),
        **{f: ("V8", 0.001) for f in ("3p-v4", "3p-v4-oh", "5p-v3", "5p-v3-oh")},
        "default": (_D7, 0.01),
    },
    "9.0": {
        "3p-v3-lt": ("LT", 0.01),
        **{f: ("R9", 0.001) for f in ("3p-v3-ht", "3p-v4", "3p-v4-oh", "5p-ht", "5p-v3", "5p-v3-oh")},
        **{f: (_D9, 0.001) for f in ("3p-v3", "5p-v2")},
        "default": (_D9, 0.01),
    },
}
RANGE_RULES["10.0"] = {**RANGE_RULES["9.0"], "flex-v2": ("F2", 0.01)}
RANGE_RULES["10.1"] = {k: v for k, v in RANGE_RULES["10.0"].items() if k != "3p-v3-lt"}


def chemistry_name(chemistry: Optional[str]) -> Optional[str]:
    """A chemistry given as a name or as a description, as its name (None if unknown)."""
    if chemistry is None:
        return None
    if chemistry in CHEMISTRY_FAMILY or chemistry.startswith("Flex-v2"):
        return chemistry
    return DESCRIPTIONS.get(chemistry)


def chemistry_family(chemistry: Optional[str]) -> str:
    name = chemistry_name(chemistry)
    if name is None:
        return "other"
    if name.startswith("Flex-v2"):
        return "flex-v2"
    return CHEMISTRY_FAMILY.get(name, "other")


def ambient_range(release: str, family: str, n_tags: Optional[int]) -> Tuple[Tuple[int, int], float, str]:
    """((low, high), FDR, rule) for a release, kit family and number of probe barcodes / OCM tags."""
    rules = RANGE_RULES[release]
    rule, fdr = rules.get(family, rules["default"])
    m = n_tags
    if rule == "LT":
        n = 9000
    elif rule == "D7":
        n = 45000 * m if m is not None and m > 1 else 90000
    elif rule == "V8":
        n = 80000 * m if m is not None and m > 1 else 160000
    elif rule == "HT8":
        n = 160000
    elif rule == "R9":
        n = 160000 if m is None else 40000 * m
    elif rule == "D9":
        n = 90000 if m is None else 90000 * m
    elif rule == "F2":
        n = 160000 * (m if m else 1)
    else:  # pragma: no cover - the tables hold no other rule
        raise ValueError(rule)
    return (n // 2, n), fdr, rule


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


def genomes(features: pd.DataFrame, sort: bool = False) -> List[Tuple[str, np.ndarray]]:
    """Each genome and its features: in order of first appearance, or sorted (10.0+)."""
    if "genome" not in features.columns:
        return [("", np.ones(len(features), dtype=bool))]
    values = features["genome"].astype(str).to_numpy()
    names = list(dict.fromkeys(values.tolist()))
    if sort:
        names = sorted(names)
    return [(name, values == name) for name in names]


def _row_sums(matrix: sp.csr_matrix) -> np.ndarray:
    return np.asarray(matrix.sum(axis=1), dtype=np.int64).ravel()


def _feature_types(droplets) -> np.ndarray:
    from sc_droplets.counts import feature_type_column

    column = feature_type_column(droplets.features)
    if column is None:
        return np.full(droplets.n_features, GEX, dtype=object)
    return droplets.features[column].astype(str).to_numpy()


# ---------------------------------------------------------------------------
# OrdMag
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
    size = len(counts)
    chosen = np.sort(legacy_argsort(counts, window=(size - n, size))[::-1][0:n])
    baseline = int(np.sort(nonzero)[-baseline_index])
    return chosen, {
        "n": n,
        "bootstrap": top,
        "recovered_cells": recovered_cells,
        "baseline_rank": baseline_index,
        "baseline": baseline,
        "cutoff": max(1, int(round_half_away(0.1 * baseline))),
    }


def _within_ordmag_v7(ascending: np.ndarray, index):
    """Number of values >= max(1, round(0.1 × the (index+1)-th largest)), for one index or an array."""
    baseline = ascending[-(np.asarray(index) + 1)]
    cutoff = np.maximum(1, np.round(0.1 * baseline)).astype(np.int64)
    return len(ascending) - np.searchsorted(ascending, cutoff, side="left")


def ordmag_v7(
    counts: np.ndarray, expected: Optional[int], cap: int, tie_rule: str
) -> Tuple[np.ndarray, Dict]:
    """Cell Ranger 7.2+ OrdMag (see the page): estimated expected cells, tie extension."""
    nonzero = counts[counts > 0]
    if nonzero.size == 0:
        return np.zeros(0, dtype=np.int64), {"n": 0}
    state = np.random.RandomState(0)  # private, fresh per genome and GEM group
    info: Dict = {}
    if expected is None:
        grid = np.unique(np.round(2 ** np.linspace(1, np.log2(cap), GRID_POINTS)).astype(np.int64))
        index = np.minimum(np.round(grid * Q).astype(np.int64), len(nonzero) - 1)
        picks = np.empty(ORDMAG_BOOTSTRAPS, dtype=np.int64)
        for b in range(ORDMAG_BOOTSTRAPS):
            sample = np.sort(state.choice(nonzero, len(nonzero)))
            found = _within_ordmag_v7(sample, index)
            loss = (found - grid) ** 2 / grid
            picks[b] = grid[int(np.argmin(loss))]
        expected = max(int(np.round(np.mean(picks))), MIN_RECOVERED)
        info["estimated_cells"] = expected
        info["estimates"] = picks
    else:
        expected = max(expected, MIN_RECOVERED)
    index = min(int(np.round(expected * Q)), len(nonzero) - 1)
    top = np.array(
        [int(_within_ordmag_v7(np.sort(state.choice(nonzero, len(nonzero))), index)) for _ in range(ORDMAG_BOOTSTRAPS)]
    )
    n0 = int(np.round(np.mean(top)))
    n = _extend_ties(np.sort(nonzero)[::-1], n0, tie_rule) if n0 > 0 else 0
    order = np.argsort(counts, kind="stable")[::-1]  # ties: higher index first
    chosen = np.sort(order[:n])
    ascending = np.sort(nonzero)
    baseline = int(ascending[-(index + 1)])
    info.update(
        {
            "n": n,
            "n_before_ties": n0,
            "bootstrap": top,
            "recovered_cells": expected,
            "baseline_rank": index + 1,
            "baseline": baseline,
            "cutoff": int(max(1, np.round(0.1 * baseline))),
        }
    )
    return chosen, info


def _extend_ties(descending: np.ndarray, n0: int, rule: str) -> int:
    """The tie extension as written: it always takes one barcode past the tie run."""
    tie = descending[n0 - 1]
    j = n0 - 1
    accepted = n0
    while j + 1 < len(descending) and descending[j] == tie:
        j += 1
        if (j + 1 - n0) > 0.2 * n0:
            return n0 if rule == "revert" else accepted
        accepted = j + 1
    return j + 1 if rule == "revert" else accepted


def top_n_v2(counts: np.ndarray, force: int) -> np.ndarray:
    """Cell Ranger 2.x/3.x `filter_cellular_barcodes_fixed_cutoff`."""
    top = min(force, int((counts > 0).sum()))
    size = len(counts)
    return np.sort(legacy_argsort(counts, window=(size - top, size))[::-1][0:top])


def top_n_v7(counts: np.ndarray, force: int) -> np.ndarray:
    """7.2+: the same count, in stable descending order (ties: higher index first)."""
    top = min(force, int((counts > 0).sum()))
    return np.sort(np.argsort(counts, kind="stable")[::-1][0:top])


# ---------------------------------------------------------------------------
# EmptyDrops step
# ---------------------------------------------------------------------------


def _profile(counts: sp.csr_matrix, ambient_rows: np.ndarray, features: np.ndarray, release: str) -> np.ndarray:
    """Good-Turing ambient profile over `features`; P0 shared by the zero features."""
    profile = np.asarray(counts[ambient_rows][:, features].sum(axis=0), dtype=np.int64).ravel()
    zero = np.flatnonzero(profile == 0)
    nonzero = np.flatnonzero(profile)
    smoothed, p0 = sgt_cellranger(profile[nonzero], clamp_slope=release not in ("3.x", "7.2", "8.0.0", "8.0.1"))
    if zero.size == 0:
        if release == "3.x":
            # Cell Ranger 3.0.2 divides P0 by zero here and then fails its sum-to-one check.
            raise GoodTuringError("every feature has ambient counts, which Cell Ranger 3.0.2 fails on")
        out = np.empty(len(features))
        out[nonzero] = smoothed / smoothed.sum()
        return out
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
    size = len(umis)
    order = legacy_argsort(umis, window=(size - N_PARTITIONS, size - N_PARTITIONS // 2))
    empty = np.sort(order[::-1][N_PARTITIONS // 2 : N_PARTITIONS])
    use = np.intersect1d(empty, np.flatnonzero(umis), assume_unique=True)
    if use.size == 0:
        return None
    features = np.flatnonzero(np.asarray(counts.sum(axis=0)).ravel())
    try:
        profile = _profile(counts, use, features, "3.x")
    except GoodTuringError as error:
        return {"stopped": str(error)}
    if not initial.any():
        return None
    median = float(np.median(umis[initial]))
    min_umis = int(max(MIN_UMIS, round_half_away(float(np.ceil(median * MIN_UMI_FRAC_OF_MEDIAN)))))
    masked = initial | (umis < min_umis)
    n_unmasked = int((~masked).sum())
    filled = np.where(masked, np.iinfo(np.int64).max, umis)  # numpy.ma sorts masked entries last
    candidates = legacy_argsort(filled, window=(n_unmasked - N_CANDIDATE_BARCODES, n_unmasked))[0:n_unmasked][-N_CANDIDATE_BARCODES:]
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
        "ambient_ranks": [N_PARTITIONS // 2, N_PARTITIONS],
        "n_ambient": int(use.size),
    }


def nonambient_v7(
    counts: sp.csr_matrix,
    feature_ids: np.ndarray,
    originals: np.ndarray,
    settings: Dict,
    ambient: Tuple[int, int],
    fdr: float,
    ed_minimum_umis: int,
    generator: Optional[np.random.Generator],
) -> Optional[Dict]:
    """Cell Ranger 7.2+ EmptyDrops on one genome's GEX features and one GEM group (rows)."""
    release = settings["release"]
    umis = _row_sums(counts)
    low, high = ambient
    size = len(umis)
    descending = legacy_argsort(umis, window=(size - high, size - low))[::-1]
    empty = descending[low:high] if len(umis) >= low else np.zeros(0, dtype=np.int64)
    use = np.sort(empty[umis[empty] > 0])
    if use.size == 0:
        return None
    if np.flatnonzero(originals).sum() == 0:
        return None  # Cell Ranger tests the sum of the original cells' indices
    if settings["candidates"] == "median":
        median = float(np.median(umis[originals]))
        min_umis = int(max(ed_minimum_umis, np.ceil(median * 0.01)))
    else:
        min_umis = max(ed_minimum_umis, 1 + (int(umis[empty].max()) if empty.size else 0))
    candidates = np.flatnonzero(~originals & (umis >= min_umis))
    if candidates.size == 0:
        return None
    features = np.flatnonzero(np.asarray(counts.sum(axis=0)).ravel())
    try:
        profile = _profile(counts, use, features, release)
    except GoodTuringError as error:
        return {"stopped": str(error)}
    rows = counts[candidates][:, features]
    totals = umis[candidates]
    if settings["simulator"] == "v1":
        with np.errstate(divide="ignore"):
            observed = multinomial_loglik_rows(rows, np.log(scipy_adjusted_profile(profile)))
        state = np.random.RandomState(0)  # Cell Ranger seeds the global stream with 0 per call
        n_lower = cellranger_v1(profile, totals, observed, settings["num_sims"], state, batch_small=20)
    else:
        observed = multinomial_loglik_rows(rows, np.log(profile))
        order = None
        if settings["feature_order"] == "id":
            ids = np.array([str(i).encode("utf-8") for i in feature_ids[features]])
            order = np.argsort(ids)
        n_lower = cellranger_v9(profile, totals, observed, settings["num_sims"], generator, order=order)
    pvalues = monte_carlo_pvalues(n_lower, settings["num_sims"])
    adjusted_p = bh_cellranger(pvalues)
    return {
        "candidates": candidates,
        "loglik": observed,
        "pvalue": pvalues,
        "fdr": adjusted_p,
        "nonambient": adjusted_p <= fdr,
        "min_umis": min_umis,
        "ambient_ranks": [low, high],
        "n_ambient": int(use.size),
    }


# ---------------------------------------------------------------------------
# Filters around calling (7.2+)
# ---------------------------------------------------------------------------


def antibody_aggregates(counts: sp.csr_matrix, n_tags: Optional[int]) -> Tuple[np.ndarray, Dict]:
    """The count-only antibody aggregate rule (5.0+, as in 7.2–10.1): a mask of aggregates."""
    dense = counts.toarray()
    per_antibody = dense.sum(axis=0)
    signal = np.flatnonzero(per_antibody >= 1000)
    n = int(signal.size)
    aggregates = np.zeros(counts.shape[0], dtype=bool)
    if n < 5:
        return aggregates, {"n_signal_antibodies": n}
    k = (n_tags or 1) * 25
    totals = dense[:, signal].sum(axis=1)
    size = counts.shape[0]
    candidates = legacy_argsort(totals, window=(size - k, size))[-k:]
    hits = np.zeros(size, dtype=np.int64)
    for antibody in signal:
        hits[legacy_argsort(dense[:, antibody], window=(size - k, size))[-k:]] += 1
    f = 0.6 if n > 26 else (-0.02 * n) + 1.1
    needed = int(np.round(n * f))
    aggregates[candidates[hits[candidates] >= needed]] = True
    return aggregates, {"n_signal_antibodies": n, "hits_needed": needed}


def high_occupancy_gems(barcodes: Sequence[str], tool_version: str = "10.1", gem_length: int = 16) -> Tuple[List[str], Dict]:
    """Cell Ranger's high-occupancy GEM filter for multiplexed Flex (7.2+).

    `barcodes` are the cells called in *every* sample of the run (Cell Ranger
    applies this after merging the samples): GEM barcode (`gem_length`
    nucleotides) + probe barcode + "-" + GEM group. Returns the barcodes to
    remove (those in GEMs holding more cells than a simulation of random GEM
    occupancy allows) and what the filter computed.
    """
    from sc_droplets.registry import ToolVersion

    newer = ToolVersion(tool_version) >= ToolVersion("10.0")
    called = sorted(set(map(str, barcodes)))
    gems = [b[:gem_length] for b in called]
    labels = [b[gem_length:] for b in called]
    per_gem = pd.Series(gems).value_counts()
    partitions = 115000 * (1 / 1.65)
    lam = len(called) / max(int(partitions), int(per_gem.size)) if called else 0.0
    info: Dict = {"lambda": lam, "n_gems": int(per_gem.size)}
    if lam == 0:
        threshold = 0
    else:
        order = list(dict.fromkeys(labels))  # first appearance
        frequency = pd.Series(labels).value_counts()
        probabilities = np.array([frequency[label] for label in order], dtype=np.float64) / len(labels)
        state = np.random.RandomState(0)
        if newer:
            occupancy = state.negative_binomial(lam, 0.5, 1_000_000)
        else:
            occupancy = state.poisson(lam, 1_000_000)
        occupancy = occupancy[occupancy > 0]
        width = int(occupancy.max())
        draws = state.choice(len(order), size=(occupancy.size, width), p=probabilities)
        valid = np.arange(width)[None, :] < occupancy[:, None]
        ordered = np.sort(np.where(valid, draws, len(order)), axis=1)
        new_value = np.ones_like(ordered, dtype=bool)
        new_value[:, 1:] = ordered[:, 1:] != ordered[:, :-1]
        distinct = (new_value & (ordered < len(order))).sum(axis=1)
        if newer:
            weights = np.bincount(distinct) * np.arange(distinct.max() + 1)
            threshold = int(np.flatnonzero(np.cumsum(weights) / weights.sum() >= 0.999)[0])
        else:
            threshold = int(np.ceil(np.quantile(distinct, 0.999)))
    crowded = set(per_gem.index[per_gem.to_numpy() > threshold])
    removed = [b for b, g in zip(called, gems) if g in crowded]
    info.update({"threshold": threshold, "n_removed": len(removed)})
    return removed, info


def _read_barcode_list(path: str) -> set:
    """Barcodes from Cell Ranger's aggregate_barcodes.csv (a 'barcode' column) or a plain list."""
    table = pd.read_csv(path, sep=None, engine="python")
    column = "barcode" if "barcode" in table.columns else table.columns[0]
    values = table[column].astype(str).tolist()
    if column == table.columns[0] and column != "barcode":
        values = [column] + values  # no header: the first line was a barcode
    return set(values)


# ---------------------------------------------------------------------------
# Methods
# ---------------------------------------------------------------------------

_V2 = Variant(
    "v2",
    since="2.0",
    until="6.0",
    summary="OrdMag: 100 bootstraps of the count of barcodes within 10x of the 30th-largest (for 3,000 "
    "expected cells); legacy numpy RNG seeded 0; Python 2 rounding",
    settings={"release": "3.x"},
)


def _v7_settings(release: str) -> Dict:
    newer = release not in ("7.2", "8.0.0", "8.0.1")
    return {
        "release": release,
        "tie_rule": "revert" if release == "7.2" else "truncate",
        "candidates": "median" if release == "7.2" else "above_ambient",
        "simulator": "v9" if newer else "v1",
        "num_sims": 100000 if newer else 10000,
        "feature_order": "id" if release in ("10.0", "10.1") else "index",
        "loop": "group_first" if newer else "genome_first",
        "post_filters": release != "7.2",
        "sorted_genomes": release in ("10.0", "10.1"),
    }


_V7_VARIANTS = (
    Variant("v7.2", since="7.2", until="8.0", summary="OrdMag with estimated expected cells and tie extension; "
            "EmptyDrops on chemistry-dependent ranks, candidates ≥ max(500, 1% of the median), no cap",
            settings=_v7_settings("7.2")),
    Variant("v8.0", since="8.0", until="8.0.1", summary="candidates must lie above the ambient range; GEM-X "
            "ranges and FDR 0.001; global-minimum-UMI and mito filters", settings=_v7_settings("8.0.0")),
    Variant("v8.0.1", since="8.0.1", until="9.0", summary="as 8.0, with the HT range for HT chemistries set by hand",
            settings=_v7_settings("8.0.1")),
    Variant("v9.0", since="9.0", until="10.0", summary="100,000 simulations from a PCG64 stream seeded 42; "
            "FDR 0.001 for 3′ v3/v4, 5′ and HT; Good-Turing slope clamped", settings=_v7_settings("9.0")),
    Variant("v10.0", since="10.0", until="10.1", summary="Flex v2 ranges; features in id order in the simulation",
            settings=_v7_settings("10.0")),
    Variant("v10.1", since="10.1", summary="the LT range removed", settings=_v7_settings("10.1")),
)


_ORDMAG_V7_VARIANTS = (
    Variant("v7.2", since="7.2", until="8.0", summary="baseline index + 1, numpy rounding, a private "
            "RandomState(0); expected cells estimated when not given; the tie extension reverted when over 20%",
            settings=_v7_settings("7.2")),
    Variant("v8.0", since="8.0", until="8.0.1", summary="the tie extension truncated at 20% instead",
            settings=_v7_settings("8.0.0")),
    Variant("v8.0.1", since="8.0.1", until="9.0", summary="as 8.0; the estimate's cap follows 8.0.1's ranges",
            settings=_v7_settings("8.0.1")),
    Variant("v9.0", since="9.0", until="10.0", summary="as 8.0; the estimate's cap follows 9.0's ranges",
            settings=_v7_settings("9.0")),
    Variant("v10.0", since="10.0", until="10.1", summary="as 8.0; the estimate's cap follows 10.0's ranges",
            settings=_v7_settings("10.0")),
    Variant("v10.1", since="10.1", summary="as 8.0; the estimate's cap follows 10.1's ranges",
            settings=_v7_settings("10.1")),
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
        variants=(_V2,) + _ORDMAG_V7_VARIANTS,
        gaps=GAPS,
        latest_known=LATEST,
        references=REFERENCES,
        source="https://github.com/lfperales/cellranger/blob/3.0.2/lib/python/cellranger/stats.py",
        original_license="MIT (3.0.2); 10x licence (7.2+)",
        default_seed=0,
        feature_types=None,
    )

    @dataclass
    class Params:
        expect_cells: Optional[int] = param(None, "Expected number of cells (2.x–5.x: 3,000 per GEM group if not given; 7.2+: estimated)", original="--expect-cells", minimum=1)
        chemistry: Optional[str] = param(None, "7.2+: chemistry name or description, for the estimate's cap (default: the input file's)", original="--chemistry")
        n_probe_barcodes: Optional[int] = param(None, "7.2+: probe barcodes / OCM tags of the sample, for the estimate's cap", original="probe_barcode_ids", minimum=1)

    def run(self, droplets, params, variant, random_state):
        legacy = variant.name == "v2"
        if legacy and (params.chemistry or params.n_probe_barcodes):
            raise ValueError("cellranger-ordmag before 7.2 takes no chemistry or probe barcodes")
        gex, types = _gex_or_antibody(droplets)
        groups = gem_groups(droplets.barcodes)
        unique = sorted(set(groups.tolist()))
        cells = np.zeros(droplets.n_barcodes, dtype=bool)
        diagnostics: Dict = {"groups": {}}
        notes: List[str] = []
        state = np.random.RandomState(0 if random_state is None else random_state)
        if legacy:
            recovered = int(float(params.expect_cells) / len(unique)) if params.expect_cells else DEFAULT_RECOVERED_CELLS
        else:
            chemistry = params.chemistry or droplets.source.get("attrs", {}).get("chemistry_description")
            (low, _), _, _ = ambient_range(variant.settings["release"], chemistry_family(chemistry), params.n_probe_barcodes)
            cap = min(low, MAX_RECOVERED)
            expected = _expected_v7(params.expect_cells, len(unique))
        for genome, feature_mask in genomes(droplets.features):
            columns = np.flatnonzero(feature_mask & gex)
            sub = droplets.counts[:, columns]
            for group in unique:
                rows = np.flatnonzero(groups == group)
                if legacy:
                    chosen, info = ordmag_v2(_row_sums(sub[rows]), recovered, state)
                else:
                    chosen, info = ordmag_v7(_row_sums(sub[rows]), expected, cap, variant.settings["tie_rule"])
                cells[rows[chosen]] = True
                diagnostics["groups"][f"{genome or 'all'}-{group}"] = info
        first = next(iter(diagnostics["groups"].values()))
        return result(
            type(self), variant, params, droplets, cells,
            thresholds={"umis": float(first.get("cutoff", 0))}, diagnostics=diagnostics, notes=notes,
        )


def _expected_v7(expect: Optional[int], n_groups: int) -> Optional[int]:
    """`int(float(E) / G)`, then `// G` again if nonzero (as Cell Ranger 7.2+ does)."""
    if not expect:
        return None
    staged = int(float(expect) / n_groups)
    return staged // n_groups if staged else None


def _remove_aggregates(droplets, params, types: np.ndarray, notes: List[str], diagnostics: Dict) -> np.ndarray:
    """The barcodes kept for calling (5.0+): without Cell Ranger's list of antibody aggregates, or
    without those of its count rule."""
    kept = np.ones(droplets.n_barcodes, dtype=bool)
    if params.aggregate_barcodes:
        listed = _read_barcode_list(params.aggregate_barcodes)
        kept &= ~np.isin(droplets.barcodes.astype(str), list(listed))
        diagnostics["aggregates_removed"] = int((~kept).sum())
    elif params.filter_aggregates and (types == ANTIBODY).any():
        # Antibody-only runs too: Cell Ranger filters aggregates whenever there is an antibody library.
        aggregates, info = antibody_aggregates(droplets.counts[:, np.flatnonzero(types == ANTIBODY)], params.n_probe_barcodes)
        kept &= ~aggregates
        diagnostics["aggregates_removed"] = int(aggregates.sum())
        diagnostics["aggregates"] = info
        notes.append(
            "antibody aggregates removed by Cell Ranger's count rule; its read-based rule needs read counts "
            "a matrix does not hold (give aggregate_barcodes= for Cell Ranger's own list)"
        )
    return kept


def _gex_or_antibody(droplets) -> Tuple[np.ndarray, np.ndarray]:
    """The features Cell Ranger calls on: Gene Expression, or Antibody Capture in an antibody-only run."""
    types = _feature_types(droplets)
    gex = types == GEX
    if not gex.any():
        gex = types == ANTIBODY
        if not gex.any():
            raise ValueError("Cell Ranger calls cells on Gene Expression (or Antibody Capture) features; there are none")
    return gex, types


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
        deterministic=True,
        variants=(
            Variant("v2", since="2.0", until="3.0", summary="OrdMag only, per genome and GEM group",
                    settings={"release": "2.x"}),
            Variant(
                "v3",
                since="3.0",
                until="4.0",
                summary="OrdMag, then EmptyDrops on ranks 45,000–90,000 as ambient, the top 20,000 candidates "
                "with at least max(500, 1% of the median cell), 10,000 simulations, BH FDR ≤ 0.01",
                settings={"release": "3.x"},
            ),
            Variant("v4", since="4.0", until="5.0", summary="as 3.x, but the simulation draws from a stream "
                    "of its own, seeded 0", settings={"release": "4.x"}),
            Variant("v5", since="5.0", until="6.0", summary="as 4.x, after removing antibody aggregates "
                    "(the count rule of 7.2+)", settings={"release": "5.x"}),
        )
        + _V7_VARIANTS,
        gaps=GAPS,
        latest_known=LATEST,
        references=REFERENCES,
        source="https://github.com/lfperales/cellranger/blob/3.0.2/lib/python/cellranger/cell_calling.py",
        original_license="MIT (3.0.2); 10x licence (7.2+)",
        default_seed=0,
        feature_types=None,
    )

    @dataclass
    class Params:
        expect_cells: Optional[int] = param(None, "Expected number of cells (2.x–5.x: 3,000 per GEM group if not given; 7.2+: estimated)", original="--expect-cells", minimum=1)
        force_cells: Optional[int] = param(None, "Call exactly this many cells (the top N) instead", original="--force-cells", minimum=1)
        chemistry: Optional[str] = param(None, "7.2+: chemistry name or description (default: the input file's)", original="--chemistry")
        n_probe_barcodes: Optional[int] = param(None, "7.2+: number of probe barcodes / OCM tags of this sample (Flex, OCM)", original="probe_barcode_ids", minimum=1)
        emptydrops_minimum_umis: int = param(500, "7.2+: fewest UMIs a candidate may have", original="emptydrops-minimum-umis", minimum=0)
        global_minimum_umis: int = param(0, "8.0+: remove called cells with fewer UMIs", original="global-minimum-umis", minimum=0)
        max_mito_percent: float = param(100.0, "8.0+: remove called cells with a higher mitochondrial %", original="max-mito-percent", minimum=0)
        filter_aggregates: bool = param(True, "5.0+: remove antibody aggregates before calling", original="filter-aggregates")
        aggregate_barcodes: Optional[str] = param(None, "5.0+: Cell Ranger's aggregate_barcodes.csv, removed before calling instead of the count rule", original="aggregate_barcodes.csv")

    def run(self, droplets, params, variant, random_state):
        release = variant.settings["release"]
        if release in ("2.x", "3.x", "4.x", "5.x"):
            unused = [n for n in ("chemistry", "n_probe_barcodes") if getattr(params, n)]
            unused += [n for n, d in (("emptydrops_minimum_umis", 500), ("global_minimum_umis", 0), ("max_mito_percent", 100.0)) if getattr(params, n) != d]
            if unused:
                raise ValueError(f"cellranger before 7.2 takes no {', '.join(unused)}")
            if params.aggregate_barcodes and release != "5.x":
                raise ValueError("cellranger before 5.0 takes no aggregate_barcodes")
            return self._run_v3(droplets, params, variant, random_state)
        return self._run_v7(droplets, params, variant, random_state)

    # -- 2.x to 5.x ------------------------------------------------------------

    def _run_v3(self, droplets, params, variant, random_state):
        release = variant.settings["release"]
        seed = 0 if random_state is None else random_state
        state = np.random.RandomState(seed)
        notes: List[str] = []
        n_all = droplets.n_barcodes
        gex, types = _gex_or_antibody(droplets)
        diagnostics: Dict = {"ordmag": {}}
        if release == "5.x":
            kept = _remove_aggregates(droplets, params, types, notes, diagnostics)
        else:
            kept = np.ones(n_all, dtype=bool)
            if (types == ANTIBODY).any():
                notes.append(f"Cell Ranger {release} first removed antibody aggregates using read counts a matrix does not hold")
        index = np.flatnonzero(kept)
        n_kept = index.size
        counts = droplets.counts[index][:, np.flatnonzero(gex)]
        features = droplets.features.iloc[np.flatnonzero(gex)]
        groups = gem_groups(droplets.barcodes[index])
        unique = sorted(set(groups.tolist()))
        recovered = int(float(params.expect_cells) / len(unique)) if params.expect_cells else DEFAULT_RECOVERED_CELLS
        initial = np.zeros(n_kept, dtype=bool)
        for genome, feature_mask in genomes(features):
            sub = counts[:, np.flatnonzero(feature_mask)]
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
        calls = _Calls(n_kept)
        calls.initial = initial
        with_test = release != "2.x"
        if with_test and not params.force_cells:
            for group in unique:
                rows = np.flatnonzero(groups == group)
                # 3.x simulates on the stage's stream, after OrdMag's bootstraps; from 4.0 the
                # simulation seeds a stream of its own (as 7.2+ does), so OrdMag no longer moves it.
                stream = state if release == "3.x" else np.random.RandomState(seed)
                found = nonambient_v3(counts[rows], initial[rows], stream)
                calls.add_test(found, rows, f"GEM group {group}", notes, diagnostics)
            calls.thresholds["fdr"] = MAX_ADJ_PVALUE
        out = _Calls(n_all)
        out.absorb(calls, index)
        return out.result(self, variant, params, droplets, diagnostics, notes, with_test=with_test)

    # -- 7.2+ --------------------------------------------------------------------

    def _run_v7(self, droplets, params, variant, random_state):
        settings = variant.settings
        release = settings["release"]
        n_all = droplets.n_barcodes
        notes: List[str] = []
        gex, types = _gex_or_antibody(droplets)
        antibody_only = not (types == GEX).any()
        chemistry = params.chemistry or droplets.source.get("attrs", {}).get("chemistry_description")
        if chemistry is None:
            notes.append("no chemistry given or recorded in the input: using the default ambient range and FDR")
        family = chemistry_family(chemistry)
        (low, high), fdr, rule = ambient_range(release, family, params.n_probe_barcodes)
        diagnostics: Dict = {"chemistry": chemistry or "", "kit": family, "range_rule": rule, "ordmag": {}}

        kept = _remove_aggregates(droplets, params, types, notes, diagnostics)
        if (types == ANTIGEN).any():
            notes.append("Cell Ranger's antigen (BEAM) aggregate rule is not reproduced")

        index = np.flatnonzero(kept)
        counts = droplets.counts[index][:, np.flatnonzero(gex)]
        features = droplets.features.iloc[np.flatnonzero(gex)]
        ids = features["gene_ids"].astype(str).to_numpy() if "gene_ids" in features.columns else features.index.astype(str).to_numpy()
        groups = gem_groups(droplets.barcodes[index])
        unique = sorted(set(groups.tolist()))
        genome_list = genomes(features, sort=settings["sorted_genomes"])
        n_kept = index.size

        # Initial calls per (genome, GEM group)
        per_call: Dict[Tuple[str, int], np.ndarray] = {}
        cap = min(low, MAX_RECOVERED)
        expected = _expected_v7(params.expect_cells, len(unique))
        for genome, mask in genome_list:
            sub = counts[:, np.flatnonzero(mask)]
            for group in unique:
                rows = np.flatnonzero(groups == group)
                bc = _row_sums(sub[rows])
                called = np.zeros(n_kept, dtype=bool)
                if params.force_cells:
                    chosen = top_n_v7(bc, params.force_cells // len(unique))
                    info = {"n": int(chosen.size)}
                else:
                    chosen, info = ordmag_v7(bc, expected, cap, settings["tie_rule"])
                called[rows[chosen]] = True
                per_call[(genome, group)] = called
                diagnostics["ordmag"][f"{genome or 'all'}-{group}"] = info
        initial = np.zeros(n_kept, dtype=bool)
        for called in per_call.values():
            initial |= called

        tested = _Calls(n_kept)
        tested.initial = initial.copy()
        if not params.force_cells and not antibody_only:
            generator = np.random.default_rng(42) if settings["simulator"] == "v9" else None

            def run_one(genome, mask, group, originals):
                rows = np.flatnonzero(groups == group)
                sub = counts[rows][:, np.flatnonzero(mask)]
                found = nonambient_v7(sub, ids[np.flatnonzero(mask)], originals[rows], settings, (low, high), fdr,
                                      params.emptydrops_minimum_umis, generator)
                label = f"GEM group {group}" + (f", genome {genome}" if genome else "")
                tested.add_test(found, rows, label, notes, diagnostics)
                if found is not None and "candidates" in found:
                    per_call[(genome, group)][rows[found["candidates"][found["nonambient"]]]] = True

            if settings["loop"] == "genome_first":
                for genome, mask in genome_list:
                    for group in unique:
                        originals = np.zeros(n_kept, dtype=bool)
                        for (g, grp), called in per_call.items():
                            if grp == group:
                                originals |= called
                        run_one(genome, mask, group, originals)
            else:
                for group in unique:
                    originals = initial.copy()
                    for genome, mask in genome_list:
                        run_one(genome, mask, group, originals)
            tested.thresholds["fdr"] = fdr

        # 8.0+ filters on the called cells, per genome and GEM group. They count the genome's
        # Gene Expression features only, whatever the cells were called on: in an antibody-only
        # run every total is 0, so a positive global minimum removes every cell (as Cell Ranger does).
        if settings["post_filters"] and not params.force_cells:
            removed_min = removed_mito = 0
            kept_counts = droplets.counts[index]
            all_ids = (
                droplets.features["gene_ids"].astype(str).to_numpy()
                if "gene_ids" in droplets.features.columns
                else droplets.features.index.astype(str).to_numpy()
            )
            all_genomes = (
                droplets.features["genome"].astype(str).to_numpy() if "genome" in droplets.features.columns else None
            )
            for genome, _ in genome_list:
                in_genome = types == GEX
                if all_genomes is not None:
                    in_genome &= all_genomes == genome
                columns = np.flatnonzero(in_genome)
                totals = _row_sums(kept_counts[:, columns]) if columns.size else np.zeros(n_kept, dtype=np.int64)
                mito_columns = columns[[str(i).rpartition("_")[2] in MITO_IDS for i in all_ids[columns]]]
                mito = _row_sums(kept_counts[:, mito_columns]) if mito_columns.size else np.zeros(n_kept, dtype=np.int64)
                for (g, group), called in per_call.items():
                    if g != genome:
                        continue
                    low_umis = called & (totals < params.global_minimum_umis)
                    removed_min += int(low_umis.sum())
                    called &= ~low_umis
                    with np.errstate(divide="ignore", invalid="ignore"):
                        percent = 100.0 * mito / totals
                    high_mito = called & (percent > params.max_mito_percent)
                    removed_mito += int(high_mito.sum())
                    called &= ~high_mito
            diagnostics["removed_global_minimum"] = removed_min
            diagnostics["removed_mito"] = removed_mito

        final = np.zeros(n_kept, dtype=bool)
        for called in per_call.values():
            final |= called
        tested.final = final
        diagnostics["ambient_range"] = [low, high]
        out = _Calls(n_all)
        out.absorb(tested, index)
        return out.result(self, variant, params, droplets, diagnostics, notes, with_test=not antibody_only)


class _Calls:
    """Initial calls, test results and final cells, assembled into a CallResult."""

    def __init__(self, n: int):
        self.n = n
        self.initial = np.zeros(n, dtype=bool)
        self.final: Optional[np.ndarray] = None
        self.rescued = np.zeros(n, dtype=bool)
        self.candidate = np.zeros(n, dtype=bool)
        self.removed = np.zeros(n, dtype=bool)
        self.pvalue = np.full(n, np.nan)
        self.fdr = np.full(n, np.nan)
        self.loglik = np.full(n, np.nan)
        self.thresholds: Dict[str, float] = {}

    def add_test(self, found, rows, label, notes, diagnostics) -> None:
        if found is None:
            notes.append(f"{label}: no barcodes to test against the ambient profile")
            return
        if "stopped" in found:
            notes.append(f"{label}: EmptyDrops step skipped: {found['stopped']}")
            return
        where = rows[found["candidates"]]
        self.candidate[where] = True
        self.pvalue[where], self.fdr[where], self.loglik[where] = found["pvalue"], found["fdr"], found["loglik"]
        self.rescued[where[found["nonambient"]]] = True
        self.thresholds["candidate_umis"] = float(found["min_umis"])
        diagnostics.setdefault("nonambient", {})[label] = {
            "min_umis": found["min_umis"],
            "n_candidates": int(found["candidates"].size),
            "n_ambient": found["n_ambient"],
            "ambient_ranks": found["ambient_ranks"],
        }

    def absorb(self, other: "_Calls", index: np.ndarray) -> None:
        """Place calls made on a subset of the barcodes (`index`) into this full set; the rest were removed."""
        for name in ("initial", "rescued", "candidate", "pvalue", "fdr", "loglik"):
            getattr(self, name)[index] = getattr(other, name)
        self.final = np.zeros(self.n, dtype=bool)
        self.final[index] = other.final if other.final is not None else (other.initial | other.rescued)
        self.removed[:] = True
        self.removed[index] = False
        self.thresholds = dict(other.thresholds)

    def result(self, method, variant, params, droplets, diagnostics, notes, with_test):
        cells = self.final if self.final is not None else (self.initial | self.rescued)
        classes = np.full(self.n, "background", dtype=object)
        classes[self.candidate] = "candidate"
        classes[self.initial & cells] = "initial"
        classes[self.rescued & ~self.initial & cells] = "rescued"
        # aggregates, and cells the 8.0+ filters took out
        classes[self.removed | ((self.initial | self.rescued) & ~cells)] = "removed"
        thresholds = dict(self.thresholds)
        if (self.initial & cells).any():
            thresholds["umis"] = float(droplets.totals[self.initial & cells].min())
        return result(
            type(method),
            variant,
            params,
            droplets,
            cells,
            classes=pd.Categorical(classes, categories=CLASSES),
            columns={"pvalue": self.pvalue, "fdr": self.fdr, "loglik": self.loglik} if with_test else {},
            thresholds=thresholds,
            diagnostics=diagnostics,
            notes=notes,
        )
