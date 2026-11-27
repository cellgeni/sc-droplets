"""What the plots read: totals, ranks and the stored calls to draw."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, List, Optional, Sequence, Union

import numpy as np

from sc_droplets.counts import Droplets
from sc_droplets.registry import parse_spec
from sc_droplets.result import UNS_KEY


@dataclass
class Curve:
    """The barcode-rank curve: totals sorted decreasing, zeros dropped."""

    totals: np.ndarray  # per barcode, input order
    order: np.ndarray  # barcode indices by decreasing total (stable)
    sorted: np.ndarray  # totals[order], > 0 only
    unit: str

    @property
    def n(self) -> int:
        return int(self.sorted.shape[0])

    def rank_of(self) -> np.ndarray:
        """1-based rank of every barcode on the curve (0 for zero totals)."""
        ranks = np.zeros(self.totals.shape[0], dtype=np.int64)
        ranks[self.order[: self.n]] = np.arange(1, self.n + 1)
        return ranks

    def steps(self):
        """(rank, total) points tracing the curve: each run of equal totals as a flat step."""
        values = self.sorted
        if values.size == 0:
            return np.array([]), np.array([])
        change = np.flatnonzero(np.diff(values)) + 1
        starts = np.concatenate([[0], change])
        ends = np.concatenate([change, [values.size]])
        x = np.empty(2 * starts.size)
        y = np.empty(2 * starts.size)
        x[0::2] = starts + 1
        x[1::2] = ends
        y[0::2] = values[starts]
        y[1::2] = values[starts]
        return x, y


def agreement_groups(cells: np.ndarray) -> np.ndarray:
    """Each barcode's agreement across calls (one column per call): `every call`, `some calls` or
    `no call`; with one call, `cell` or `not a cell`."""
    n = cells.shape[1]
    hits = cells.sum(axis=1)
    if n == 1:
        return np.where(hits == 1, "cell", "not a cell").astype(object)
    return np.where(hits == n, "every call", np.where(hits > 0, "some calls", "no call")).astype(object)


def curve_from_totals(totals: np.ndarray, unit: str = "umis") -> Curve:
    totals = np.asarray(totals, dtype=np.float64)
    order = np.argsort(-totals, kind="stable")
    ordered = totals[order]
    return Curve(totals=totals, order=order, sorted=ordered[ordered > 0], unit=unit)


def droplets_for(
    adata: Any,
    keys: Sequence[str],
    *,
    layer: Optional[str] = None,
    feature_types: Union[None, str, Sequence[str]] = None,
) -> Droplets:
    """What the calls counted: their layer, feature types and genome, unless given. A stored list
    of feature types comes back from h5ad as an array; it is taken as a list, never truth-tested."""
    record = adata.uns.get(UNS_KEY, {}).get(keys[0], {}) if keys else {}
    counting = record.get("counting", {}) if hasattr(record, "get") else {}
    if layer is None:
        layer = counting.get("layer")
    if feature_types is None:
        stored = counting.get("feature_types")
        if stored is None:
            feature_types = "Gene Expression" if _has_type(adata, "Gene Expression") else "all"
        else:
            feature_types = stored if isinstance(stored, str) else [str(s) for s in np.atleast_1d(stored)]
    return Droplets.from_anndata(adata, layer=layer, feature_types=feature_types, genome=counting.get("genome"))


def curve_for(
    adata: Any,
    keys: Sequence[str],
    *,
    layer: Optional[str] = None,
    feature_types: Union[None, str, Sequence[str]] = None,
) -> Curve:
    """The curve the calls were made on: their layer, feature types and genome, unless given."""
    return curve_from_totals(droplets_for(adata, keys, layer=layer, feature_types=feature_types).totals)


def _has_type(adata: Any, kind: str) -> bool:
    from sc_droplets.counts import feature_type_column

    column = feature_type_column(adata.var)
    return column is None or kind in set(adata.var[column].astype(str))


def resolve_calls(adata: Any, calls: Union[None, str, Sequence[str]], run_missing: bool = True) -> List[str]:
    """Keys of the calls to draw: stored keys as given, specs run first."""
    from sc_droplets.api import call

    stored = adata.uns.get(UNS_KEY, {})
    if calls is None:
        return list(stored)
    if isinstance(calls, str):
        calls = [calls]
    keys = []
    for entry in calls:
        if entry in stored:
            keys.append(entry)
            continue
        if not run_missing:
            raise KeyError(f"no stored call {entry!r}; stored: {', '.join(stored) or 'none'}")
        parse_spec(entry)
        keys.append(call(adata, entry).key)
    return keys
