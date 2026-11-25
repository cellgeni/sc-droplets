"""The barcode-rank (knee) plot, with each call's cut-off on it."""

from __future__ import annotations

from typing import Any, Optional, Sequence, Union

import numpy as np

from sc_droplets.pl._data import curve_for, resolve_calls

#: Okabe–Ito, readable for colour-blind viewers and in print.
PALETTE = ("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#000000")
CURVE_COLOUR = "#7f7f7f"
#: Above this many called barcodes, the single-call overlay draws a sample.
MAX_POINTS = 20000

_UNIT_LABEL = {"umis": "UMIs", "reads": "reads", "genes": "genes"}


def knee(
    adata: Any,
    calls: Union[None, str, Sequence[str]] = None,
    *,
    layer: Optional[str] = None,
    feature_types: Union[None, str, Sequence[str]] = None,
    rotated: bool = False,
    cumulative: bool = False,
    ax: Any = None,
    title: Optional[str] = None,
    legend: bool = True,
):
    """Draw the barcode-rank curve and the calls on it; return the Axes.

    `calls`: stored keys or method specs (a spec not yet stored is run first);
    None draws every stored call. With one call, its cells are marked on the
    curve; with several, each gets a dashed line at its smallest cell total
    and a dot where its cells end.

    `rotated` puts totals on x and ranks on y (Pachter's "rotated knee").
    `cumulative` plots the cumulative share of all counts against rank, the
    view of Drop-seq and UMI-tools' distance method.
    """
    import matplotlib.pyplot as plt

    keys = resolve_calls(adata, calls)
    curve = curve_for(adata, keys, layer=layer, feature_types=feature_types)
    if ax is None:
        _, ax = plt.subplots(figsize=(6.4, 4.8))
    ranks = curve.rank_of()
    unit = _UNIT_LABEL.get(curve.unit, curve.unit)

    if cumulative:
        share = np.cumsum(curve.sorted) / curve.sorted.sum() if curve.n else np.array([])
        x = np.arange(1, curve.n + 1)
        ax.plot(x, share, color=CURVE_COLOUR, lw=1.5, label="all barcodes")
        ax.set_xscale("log")
        ax.set_xlabel("Barcode rank")
        ax.set_ylabel(f"Cumulative share of {unit}")
        for i, key in enumerate(keys):
            colour = PALETTE[i % len(PALETTE)]
            cells = adata.obs[key].to_numpy(dtype=bool)
            if cells.any():
                last = int(ranks[cells].max())
                ax.axvline(last, color=colour, ls="--", lw=1, label=f"{key}: {int(cells.sum()):,} cells")
    else:
        x, y = curve.steps()
        if rotated:
            x, y = y, x
        ax.plot(x, y, color=CURVE_COLOUR, lw=1.5, label="all barcodes", zorder=1)
        ax.set_xscale("log")
        ax.set_yscale("log")
        rank_label, total_label = "Barcode rank", f"{unit} per barcode"
        ax.set_xlabel(total_label if rotated else rank_label)
        ax.set_ylabel(rank_label if rotated else total_label)
        for i, key in enumerate(keys):
            colour = PALETTE[i % len(PALETTE)]
            cells = adata.obs[key].to_numpy(dtype=bool) & (ranks > 0)
            n = int(cells.sum())
            label = f"{key}: {n:,} cells"
            if n == 0:
                ax.plot([], [], color=colour, label=label + " (none)")
                continue
            cut = float(curve.totals[cells].min())
            last = int(ranks[cells].max())
            if len(keys) == 1:
                index = np.flatnonzero(cells)
                if index.size > MAX_POINTS:
                    index = np.random.default_rng(0).choice(index, MAX_POINTS, replace=False)
                px, py = ranks[index], curve.totals[index]
                if rotated:
                    px, py = py, px
                ax.scatter(px, py, s=4, color=colour, label=label, zorder=2, linewidths=0)
            line = ax.axvline if rotated else ax.axhline
            line(cut, color=colour, ls="--", lw=1, zorder=0, label=None if len(keys) == 1 else label)
            point = (cut, last) if rotated else (last, cut)
            ax.plot(*point, marker="o", ms=5, color=colour, zorder=3)
    if title is not None:
        ax.set_title(title)
    if legend and (keys or curve.n):
        ax.legend(frameon=False, fontsize="small")
    ax.grid(True, which="major", alpha=0.3)
    return ax
