"""The unspliced fraction against the calls: `scd.pl.splicing`."""

from __future__ import annotations

from typing import Any, Optional, Sequence, Union

import numpy as np

from sc_droplets.pl._data import curve_for, resolve_calls
from sc_droplets.pl.knee import PALETTE

#: Points drawn in the scatter panel, at most.
MAX_POINTS = 20000
GROUP_COLOURS = {"every call": PALETTE[0], "some calls": PALETTE[1], "no call": "#999999",
                 "cell": PALETTE[0], "not a cell": "#999999"}


def _groups(cells: np.ndarray) -> np.ndarray:
    """Each barcode's agreement: called by every call, by some, by none (one call: cell or not)."""
    n = cells.shape[1]
    hits = cells.sum(axis=1)
    if n == 1:
        return np.where(hits == 1, "cell", "not a cell").astype(object)
    return np.where(hits == n, "every call", np.where(hits > 0, "some calls", "no call")).astype(object)


def splicing(
    adata: Any,
    calls: Union[None, str, Sequence[str]] = None,
    *,
    metric: str = "unspliced_fraction",
    min_counts: float = 10,
    layer: Optional[str] = None,
    feature_types: Union[None, str, Sequence[str]] = None,
    fig: Any = None,
):
    """Draw a per-barcode measure (the unspliced fraction by default) against the calls; return the
    Figure.

    Panels: the measure against log10 total, each barcode coloured by how many calls took it; for
    each call, the measure in its cells and in the barcodes it left out (the AUROC between them in
    the label); and, with several calls, the measure's distribution in barcodes every call took,
    some took and none took. Barcodes below `min_counts` are left out: too few molecules for a
    fraction. `obs[metric]` must exist (`scd.add_metrics`).
    """
    import matplotlib.pyplot as plt

    from sc_droplets.metrics import auroc

    if metric not in adata.obs.columns:
        raise KeyError(f"no obs column {metric!r}: run scd.add_metrics (with splicing= for the unspliced fraction)")
    keys = resolve_calls(adata, calls)
    if not keys:
        raise ValueError("no calls to draw: give method specs or stored keys")
    curve = curve_for(adata, keys, layer=layer, feature_types=feature_types)
    totals = curve.totals
    values = adata.obs[metric].to_numpy(dtype=np.float64)
    shown = (totals >= max(min_counts, 1)) & ~np.isnan(values)
    cells = np.column_stack([adata.obs[k].to_numpy(dtype=bool) for k in keys])
    groups = _groups(cells)

    several = len(keys) > 1
    if fig is None:
        fig = plt.figure(figsize=(16, 4.8) if several else (11.5, 4.8))
    axes = fig.subplots(1, 3 if several else 2)
    label = metric.replace("_", " ")

    ax = axes[0]
    index = np.flatnonzero(shown)
    if index.size > MAX_POINTS:
        index = np.sort(np.random.default_rng(0).choice(index, MAX_POINTS, replace=False))
    order = ["no call", "not a cell", "some calls", "every call", "cell"]
    for name in order:
        pick = index[groups[index] == name]
        if pick.size:
            ax.scatter(np.log10(totals[pick]), values[pick], s=3, linewidths=0, color=GROUP_COLOURS[name],
                       label=f"{name}: {int(((groups == name) & shown).sum()):,}", zorder=2 if name != "no call" else 1)
    ax.set_xlabel("log10 count per barcode")
    ax.set_ylabel(label)
    ax.set_title(f"{label} against count", fontsize="medium")
    ax.legend(frameon=False, fontsize="small", markerscale=3)
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    data, ticks, positions = [], [], []
    for i, key in enumerate(keys):
        inside = values[cells[:, i] & shown]
        outside = values[~cells[:, i] & shown]
        data += [inside, outside]
        positions += [3 * i + 1, 3 * i + 2]
        ticks.append((3 * i + 1.5, f"{key}\nAUROC {auroc(inside, outside):.2f}"))
    parts = ax.boxplot([d if d.size else np.array([np.nan]) for d in data], positions=positions, widths=0.8,
                       showfliers=False, patch_artist=True)
    for j, box in enumerate(parts["boxes"]):
        box.set_facecolor(PALETTE[0] if j % 2 == 0 else "#cccccc")
    ax.set_xticks([t for t, _ in ticks], [name for _, name in ticks], fontsize="small")
    ax.set_ylabel(label)
    ax.set_title(f"cells (blue) and the barcodes left out, ≥ {min_counts:g} counts", fontsize="medium")
    ax.grid(True, axis="y", alpha=0.3)

    if several:
        ax = axes[2]
        bins = np.linspace(np.nanmin(values[shown]), np.nanmax(values[shown]), 41) if shown.any() else 10
        for name in ("every call", "some calls", "no call"):
            pick = values[(groups == name) & shown]
            if pick.size:
                ax.hist(pick, bins=bins, density=True, histtype="step", lw=1.5, color=GROUP_COLOURS[name],
                        label=f"{name}: median {np.median(pick):.2f}")
        ax.set_xlabel(label)
        ax.set_ylabel("density")
        ax.set_title("where the calls disagree", fontsize="medium")
        ax.legend(frameon=False, fontsize="small")
        ax.grid(True, alpha=0.3)
    fig.tight_layout()
    return fig
