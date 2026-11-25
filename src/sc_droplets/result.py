"""The outcome of one call, and how it is stored on an AnnData.

On the AnnData a call keyed `k` writes:

- `obs[k]`: bool, the barcode is a cell;
- `obs[k + "_class"]`: categorical, when the method has more than two states
  (CellBender's priors: cell / transition / empty);
- `obs[k + "_" + column]` for per-barcode values (p-values, FDR, scores);
- `uns["sc_droplets"][k]`: the method, the tool version asked for, the variant
  it resolved to, every parameter, the thresholds, the number of cells, the
  package version, what was counted, and small diagnostics for the plots.

Everything in `uns` is plain: strings, numbers, bools, arrays and dicts of
them, so it survives `write_h5ad` and a plot can be redrawn from the file.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from typing import Any, Dict, List, Mapping, Optional

import numpy as np
import pandas as pd

UNS_KEY = "sc_droplets"


@dataclass
class CallResult:
    method: str
    title: str
    variant: str
    tool_version: Optional[str]
    params: Dict[str, Any]
    barcodes: np.ndarray
    totals: np.ndarray
    is_cell: np.ndarray
    classes: Optional[pd.Categorical] = None
    columns: Dict[str, np.ndarray] = field(default_factory=dict)
    thresholds: Dict[str, float] = field(default_factory=dict)
    diagnostics: Dict[str, Any] = field(default_factory=dict)
    notes: List[str] = field(default_factory=list)
    counting: Dict[str, Any] = field(default_factory=dict)
    key: Optional[str] = None

    def __post_init__(self) -> None:
        self.is_cell = np.asarray(self.is_cell, dtype=bool)
        n = len(self.barcodes)
        if self.is_cell.shape != (n,):
            raise ValueError(f"is_cell has shape {self.is_cell.shape}, expected ({n},)")
        for name, values in self.columns.items():
            if len(values) != n:
                raise ValueError(f"column {name!r} has {len(values)} values, expected {n}")
        if self.classes is not None and len(self.classes) != n:
            raise ValueError(f"classes has {len(self.classes)} values, expected {n}")

    @property
    def n_cells(self) -> int:
        return int(self.is_cell.sum())

    @property
    def min_cell_total(self) -> Optional[float]:
        """The smallest total among the called cells (None with no cells)."""
        if not self.is_cell.any():
            return None
        return float(self.totals[self.is_cell].min())

    @property
    def table(self) -> pd.DataFrame:
        """One row per barcode: total, call, class and every extra column."""
        frame = pd.DataFrame({"total": self.totals, "is_cell": self.is_cell}, index=pd.Index(self.barcodes, name="barcode"))
        if self.classes is not None:
            frame["class"] = self.classes
        for name, values in self.columns.items():
            frame[name] = values
        return frame

    def summary(self) -> Dict[str, Any]:
        """A JSON-ready description of the call (no per-barcode data)."""
        from sc_droplets import __version__

        return {
            "method": self.method,
            "title": self.title,
            "tool_version": self.tool_version,
            "variant": self.variant,
            "params": dict(self.params),
            "n_barcodes": int(len(self.barcodes)),
            "n_cells": self.n_cells,
            "min_cell_total": self.min_cell_total,
            "thresholds": dict(self.thresholds),
            "counting": dict(self.counting),
            "notes": list(self.notes),
            "package_version": __version__,
        }

    def to_json(self) -> str:
        return json.dumps(_jsonable(self.summary()), indent=2)

    def store(self, adata: Any, key: str) -> None:
        """Write the call into `adata.obs` and `adata.uns["sc_droplets"]`."""
        if list(adata.obs_names) != list(self.barcodes):
            raise ValueError("this call was made on different barcodes from this AnnData")
        # A call stored again under the same key replaces the old one whole:
        # its columns go, even those the new call does not write.
        for name in _owned(adata, key):
            if f"{key}_{name}" in adata.obs.columns:
                del adata.obs[f"{key}_{name}"]
        adata.obs[key] = self.is_cell
        if self.classes is not None:
            adata.obs[key + "_class"] = self.classes
        for name, values in self.columns.items():
            adata.obs[f"{key}_{name}"] = values
        record = self.summary()
        record["columns"] = (["class"] if self.classes is not None else []) + list(self.columns)
        record["diagnostics"] = self.diagnostics
        if UNS_KEY not in adata.uns:
            adata.uns[UNS_KEY] = {}
        adata.uns[UNS_KEY][key] = h5ad_safe(record)
        self.key = key


def _owned(adata: Any, key: str) -> List[str]:
    """The obs columns an earlier call under `key` wrote."""
    previous = adata.uns.get(UNS_KEY, {}).get(key, {})
    return [str(c) for c in np.atleast_1d(previous.get("columns", []))]


def stored_calls(adata: Any) -> Dict[str, Mapping[str, Any]]:
    """Every call stored on an AnnData, by key."""
    return dict(adata.uns.get(UNS_KEY, {}))


def h5ad_safe(value: Any) -> Any:
    """`value` with None dropped and sequences as arrays, so h5ad can store it."""
    if isinstance(value, Mapping):
        return {str(k): h5ad_safe(v) for k, v in value.items() if v is not None}
    if isinstance(value, (list, tuple)):
        items = [h5ad_safe(v) for v in value if v is not None]
        if not items:
            return np.array([], dtype=object).astype(str)
        if all(isinstance(v, str) for v in items):
            return np.array(items, dtype=object)
        if all(isinstance(v, (bool, np.bool_)) for v in items):
            return np.array(items, dtype=bool)
        if all(isinstance(v, (int, float, np.integer, np.floating)) and not isinstance(v, bool) for v in items):
            return np.asarray(items)
        return {str(i): v for i, v in enumerate(items)}
    if isinstance(value, np.generic):
        return value.item()
    return value


def _jsonable(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(v) for v in value]
    if isinstance(value, np.ndarray):
        return [_jsonable(v) for v in value.tolist()]
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value
