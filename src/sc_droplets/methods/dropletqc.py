"""DropletQC's empty-droplet call: a valley in the nuclear fraction.

Page: docs/methods/dropletqc/empty-drops.md. A port of `identify_empty_drops` (DropletQC v1.0.0,
`R/identify_empty_drops.R`, MIT licence, Copyright (c) 2021 Walter Muskovic), with the ks
estimators it calls in `primitives/ks.py`.

DropletQC estimates the density derivative of the barcodes' nuclear fractions, cuts at the grid
point where its sign changes the second time, calls the barcodes below the cut empty droplets, and
rescues those with both a nuclear fraction and a UMI count above two floors.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from packaging.version import Version

from sc_droplets.methods._common import result
from sc_droplets.primitives.ks import kdde_deriv
from sc_droplets.registry import Method, MethodInfo, Reference, Variant, param, register

#: ks's first release whose binned estimator differs (see primitives/ks.py).
KS_CHANGE = Version("1.14.3")
#: The oldest ks checked against.
KS_OLDEST = Version("1.13.2")
#: The fraction columns `fraction = "auto"` looks for, in order.
AUTO_COLUMNS = ("nuclear_fraction", "unspliced_fraction")
CLASSES = ["cell", "rescued", "empty_droplet", "not given"]
#: Points of the scatter kept in the diagnostics, at most.
MAX_POINTS = 5000


class DropletQCError(ValueError):
    """No nuclear fraction, or nothing to call."""


def ks_behaviour(version: str) -> bool:
    """Whether a ks release estimates as ks before 1.14.3 did (True) or as 1.14.3 and later."""
    parsed = Version(version)
    if parsed < KS_OLDEST:
        raise DropletQCError(f"ks {version} is older than any release checked ({KS_OLDEST})")
    return parsed < KS_CHANGE


def nuclear_fraction(droplets, fraction: str):
    """(the fraction per barcode, where it came from)."""
    obs = droplets.obs
    if fraction == "auto":
        for column in AUTO_COLUMNS:
            if column in obs.columns:
                return obs[column].to_numpy(dtype=np.float64), column
        raise DropletQCError(
            "needs each barcode's nuclear fraction: an obs column 'nuclear_fraction' (DropletQC's, from a BAM), "
            "spliced and unspliced counts for the unspliced fraction (--splicing, or scd.add_metrics), or "
            "--fraction naming a column or a table of barcode and fraction"
        )
    if fraction in obs.columns:
        return obs[fraction].to_numpy(dtype=np.float64), fraction
    path = Path(fraction)
    if not path.is_file():
        raise DropletQCError(f"--fraction {fraction}: neither an obs column nor a file")
    from sc_droplets.metrics import _matched

    sep = "\t" if any(s in path.name for s in (".tsv", ".txt")) else ","
    table = pd.read_csv(path, sep=sep, header=None, comment="#")
    if table.shape[1] < 2:
        raise DropletQCError(f"{path}: expected two columns, barcode and nuclear fraction")
    if not _is_number(table.iloc[0, 1]):
        table = table.iloc[1:]
    values = pd.Series(table.iloc[:, 1].astype(float).to_numpy(), index=table.iloc[:, 0].astype(str))
    return _matched(pd.Index(droplets.barcodes.astype(str)), values), path.name


def _is_number(value) -> bool:
    try:
        float(value)
    except (TypeError, ValueError):
        return False
    return True


def identify_empty_drops(nf: np.ndarray, umi: np.ndarray, nf_rescue: float, umi_rescue: float, legacy: bool):
    """`identify_empty_drops` on the barcodes given: (status per barcode, the cut, the derivative,
    the number of sign runs)."""
    derivative = kdde_deriv(nf, 1, legacy)
    positive = derivative.estimate > 0
    starts = np.flatnonzero(np.diff(positive.astype(np.int8))) + 1
    runs = starts.size + 1
    # R: eval.points[sum(rle(estimate > 0)$lengths[1:2])], NA with a single run.
    cut = float(derivative.eval_points[starts[1] - 1]) if runs >= 3 else (
        float(derivative.eval_points[-1]) if runs == 2 else float("nan"))
    status = np.full(nf.size, "cell", dtype=object)
    status[nf < cut] = "empty_droplet"  # a NaN cut marks nothing, as R's NA index
    status[(nf > nf_rescue) & (umi > umi_rescue) & (status == "empty_droplet")] = "rescued"
    return status, cut, derivative, runs


@register
class DropletQCEmpty(Method):
    info = MethodInfo(
        id="dropletqc",
        title="DropletQC empty droplets",
        tool="DropletQC",
        category="density",
        summary="DropletQC's identify_empty_drops: the barcodes below the valley after the first peak of the "
                "nuclear fraction's density are empty droplets, unless both their nuclear fraction and UMI "
                "count clear two floors.",
        docs="methods/dropletqc/empty-drops.md",
        counts="umis",
        variants=(Variant("v1.0", since="1.0.0", summary="identify_empty_drops of v1.0.0, the only release"),),
        references=(
            Reference(
                "Muskovic & Powell 2021, DropletQC: improved identification of empty droplets and damaged cells "
                "in single-cell RNA-seq data, Genome Biology 22:329",
                doi="10.1186/s13059-021-02547-0",
            ),
        ),
        source="https://github.com/powellgenomicslab/DropletQC/blob/v1.0.0/R/identify_empty_drops.R",
        original_license="MIT",
    )

    @dataclass
    class Params:
        fraction: str = param(
            "auto",
            "Each barcode's nuclear fraction: an obs column, a table of barcode and fraction, or 'auto' "
            "(nuclear_fraction, else unspliced_fraction)",
        )
        among: Optional[str] = param(
            None,
            "The barcodes DropletQC is given: a stored call or boolean obs column (DropletQC runs on a "
            "caller's cells); without it, every barcode with at least --min-counts",
        )
        min_counts: float = param(500.0, "Without --among, the barcodes given have at least this many UMIs",
                                  minimum=0)
        nf_rescue: float = param(0.05, "Rescued if the nuclear fraction is above this ...", original="nf_rescue",
                                 minimum=0)
        umi_rescue: float = param(1000.0, "... and the UMI count above this", original="umi_rescue", minimum=0)
        ks_version: str = param("1.15.3", "The installed ks, whose estimator the call follows (1.14.3 changed it)")

    def run(self, droplets, params, variant, random_state):
        legacy = ks_behaviour(params.ks_version)
        nf_all, source = nuclear_fraction(droplets, params.fraction)
        umi = droplets.totals
        if params.among is not None:
            if params.among not in droplets.obs.columns:
                raise DropletQCError(f"--among {params.among}: no such stored call or obs column")
            given = droplets.obs[params.among].to_numpy(dtype=bool)
        else:
            given = umi >= params.min_counts
        notes = []
        if source != "nuclear_fraction" and source in AUTO_COLUMNS:
            notes.append(f"the nuclear fraction is the {source.replace('_', ' ')}, of UMIs: an approximation of "
                         "DropletQC's, which counts reads in a BAM")
        nf = nf_all[given]
        missing = np.isnan(nf)
        if missing.any():
            # DropletQC sets 0/0 to 0; a barcode with no spliced or unspliced count is the same.
            notes.append(f"{int(missing.sum()):,} barcodes given have no fraction: set to 0, as DropletQC's 0/0")
            nf = np.where(missing, 0.0, nf)
        if nf.size < 2:
            raise DropletQCError(f"DropletQC needs more than one barcode; {nf.size} given")
        status, cut, derivative, runs = identify_empty_drops(nf, umi[given], params.nf_rescue, params.umi_rescue,
                                                             legacy)
        if runs < 4:
            notes.append("one peak only in the nuclear fraction: there may be no empty droplets (DropletQC's warning)")
        full = np.full(umi.size, "not given", dtype=object)
        full[given] = status
        nf_column = np.full(umi.size, np.nan)
        nf_column[given] = nf
        is_cell = np.isin(full, ("cell", "rescued"))
        index = np.flatnonzero(given)
        if index.size > MAX_POINTS:
            index = np.sort(np.random.default_rng(0).choice(index, MAX_POINTS, replace=False))
        return result(
            type(self), variant, params, droplets, is_cell,
            classes=pd.Categorical(full, categories=CLASSES),
            columns={"nuclear_fraction": nf_column},
            thresholds={"nuclear_fraction": cut, "bandwidth": derivative.h, "nf_rescue": params.nf_rescue,
                        "umi_rescue": params.umi_rescue, "given": float(given.sum())},
            diagnostics={
                "rule": "dropletqc",
                "grid": derivative.eval_points.tolist(),
                "derivative": derivative.estimate.tolist(),
                "points_nf": nf_column[index].tolist(),
                "points_umi": umi[index].tolist(),
                "points_class": full[index].tolist(),
                "fraction_source": source,
            },
            notes=notes,
        )
