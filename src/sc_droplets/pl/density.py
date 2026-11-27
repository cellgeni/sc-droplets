"""The distribution of log totals, with each call's cut: `scd.pl.density`."""

from __future__ import annotations

from typing import Any, Optional, Sequence, Union

import numpy as np

from sc_droplets.pl._data import curve_for, resolve_calls
from sc_droplets.pl.knee import CURVE_COLOUR, PALETTE

_UNIT_LABEL = {"umis": "UMIs", "reads": "reads", "genes": "genes"}


def density(
    adata: Any,
    calls: Union[None, str, Sequence[str]] = None,
    *,
    layer: Optional[str] = None,
    feature_types: Union[None, str, Sequence[str]] = None,
    min_counts: float = 1,
    bins: int = 100,
    weighted: bool = False,
    ax: Any = None,
    title: Optional[str] = None,
    legend: bool = True,
):
    """Draw the histogram of log10 totals and each call's cut; return the Axes.

    `calls`: stored keys or method specs (a spec not yet stored is run first);
    None draws every stored call. Each call gets a dashed line at its smallest
    cell total, labelled with its number of cells.

    Barcodes below `min_counts` are left out (log10 of zero is undefined). The
    y axis is logarithmic: empty droplets outnumber cells by orders of
    magnitude. `weighted` weights each barcode by its total, so each bar is the
    share of all counts its barcodes hold: the cells' peak then stands out.
    """
    import matplotlib.pyplot as plt

    keys = resolve_calls(adata, calls)
    curve = curve_for(adata, keys, layer=layer, feature_types=feature_types)
    if ax is None:
        _, ax = plt.subplots(figsize=(6.4, 4.8))
    totals = curve.totals
    shown = totals[totals >= max(min_counts, np.nextafter(0, 1))]
    unit = _UNIT_LABEL.get(curve.unit, curve.unit)
    if shown.size:
        x = np.log10(shown)
        weights = shown / shown.sum() if weighted else None
        counts, edges = np.histogram(x, bins=bins, weights=weights)
        ax.stairs(counts, edges, fill=True, color=CURVE_COLOUR, alpha=0.35, label=f"{shown.size:,} barcodes")
        ax.stairs(counts, edges, color=CURVE_COLOUR, lw=0.8)
    ax.set_yscale("log")
    ax.set_xlabel(f"log10 {unit} per barcode")
    ax.set_ylabel(f"share of {unit}" if weighted else "barcodes")
    for i, key in enumerate(keys):
        colour = PALETTE[i % len(PALETTE)]
        cells = adata.obs[key].to_numpy(dtype=bool) & (totals > 0)
        n = int(cells.sum())
        label = f"{key}: {n:,} cells"
        if n == 0:
            ax.plot([], [], color=colour, label=label + " (none)")
            continue
        ax.axvline(np.log10(totals[cells].min()), color=colour, ls="--", lw=1.2, label=label)
    if title is not None:
        ax.set_title(title)
    if legend and (keys or shown.size):
        ax.legend(frameon=False, fontsize="small")
    ax.grid(True, alpha=0.3)
    return ax
