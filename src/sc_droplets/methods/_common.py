"""Helpers every method module uses."""

from __future__ import annotations

import dataclasses
from typing import Any, Dict, Optional, Type

import numpy as np
import pandas as pd

from sc_droplets.counts import Droplets
from sc_droplets.registry import Method, Variant
from sc_droplets.result import CallResult


def result(
    cls: Type[Method],
    variant: Variant,
    params: Any,
    droplets: Droplets,
    is_cell: np.ndarray,
    *,
    classes: Optional[pd.Categorical] = None,
    columns: Optional[Dict[str, np.ndarray]] = None,
    thresholds: Optional[Dict[str, float]] = None,
    diagnostics: Optional[Dict[str, Any]] = None,
    notes: Optional[list] = None,
) -> CallResult:
    """A CallResult for `droplets`, filled with the method's identity."""
    return CallResult(
        method=cls.info.id,
        title=cls.info.title,
        variant=variant.name,
        tool_version=None,
        params=dataclasses.asdict(params),
        barcodes=droplets.barcodes,
        totals=droplets.totals,
        is_cell=is_cell,
        classes=classes,
        columns=columns or {},
        thresholds=thresholds or {},
        diagnostics=diagnostics or {},
        notes=notes or [],
    )


def top_by_rank(values: np.ndarray, n: int) -> np.ndarray:
    """A mask of the `n` largest values; ties broken by input order (stable)."""
    mask = np.zeros(values.shape[0], dtype=bool)
    if n > 0:
        mask[np.argsort(-values, kind="stable")[:n]] = True
    return mask
