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
    color: Optional[str] = None,
    track: Union[None, str, Sequence[str]] = None,
    bins: int = 100,
):
    """Draw the barcode-rank curve and the calls on it; return the Axes.

    `calls`: stored keys or method specs (a spec not yet stored is run first);
    None draws every stored call. With one call, its cells are marked on the
    curve; with several, each gets a dashed line at its smallest cell total
    and a dot where its cells end.

    `rotated` puts totals on x and ranks on y (Pachter's "rotated knee").
    `cumulative` plots the cumulative share of all counts against rank, the
    view of Drop-seq and UMI-tools' distance method.

    `color` colours every barcode on the curve by an obs column (the unspliced
    fraction, say). `track` adds a panel under the curve for each obs column
    named, with its median over `bins` log-spaced bins of rank, so a change in
    composition at the knee shows (CellBender's view, Fleming et al. 2023,
    Fig. 1d). With `track` the function makes its own figure: leave `ax` out.
    """
    import matplotlib.pyplot as plt

    tracks = [track] if isinstance(track, str) else list(track or [])
    if (color or tracks) and cumulative:
        raise ValueError("color= and track= draw on the rank curve, not the cumulative one")
    if tracks and rotated:
        raise ValueError("track= shares the rank axis: it cannot be rotated")
    for column in ([color] if color else []) + tracks:
        if column not in adata.obs.columns:
            raise KeyError(f"no obs column {column!r} to colour or track by (scd.add_metrics adds the unspliced and mitochondrial fractions)")
    keys = resolve_calls(adata, calls)
    curve = curve_for(adata, keys, layer=layer, feature_types=feature_types)
    track_axes = []
    if tracks:
        if ax is not None:
            raise ValueError("track= draws its own figure: leave ax= out")
        figure = plt.figure(figsize=(6.4, 4.8 + 1.6 * len(tracks)))
        grid = figure.add_gridspec(1 + len(tracks), 1, height_ratios=[3] + [1] * len(tracks), hspace=0.08)
        ax = figure.add_subplot(grid[0])
        track_axes = [figure.add_subplot(grid[i + 1], sharex=ax) for i in range(len(tracks))]
    elif ax is None:
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
        if color:
            points = _colour_by(ax, adata, curve, color, rotated)
            # One bar beside every panel, so the tracks keep the curve's rank axis.
            ax.figure.colorbar(points, ax=[ax] + track_axes, label=color.replace("_", " "), pad=0.01)
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
            if len(keys) == 1 and not color:  # with color=, the cells keep the measure's colours
                index = np.flatnonzero(cells)
                if index.size > MAX_POINTS:
                    index = np.random.default_rng(0).choice(index, MAX_POINTS, replace=False)
                px, py = ranks[index], curve.totals[index]
                if rotated:
                    px, py = py, px
                ax.scatter(px, py, s=4, color=colour, label=label, zorder=2, linewidths=0)
            line = ax.axvline if rotated else ax.axhline
            line(cut, color=colour, ls="--", lw=1, zorder=0, label=None if len(keys) == 1 and not color else label)
            point = (cut, last) if rotated else (last, cut)
            ax.plot(*point, marker="o", ms=5, color=colour, zorder=3)
    if title is not None:
        ax.set_title(title)
    if legend and (keys or curve.n):
        ax.legend(frameon=False, fontsize="small")
    ax.grid(True, which="major", alpha=0.3)
    for column, panel in zip(tracks, track_axes):
        _track(panel, adata, curve, column, bins)
    if track_axes:
        ax.tick_params(labelbottom=False)
        ax.set_xlabel("")
        track_axes[-1].set_xlabel("Barcode rank")
    return ax


def _colour_by(ax, adata, curve, column: str, rotated: bool):
    """Every barcode on the curve, coloured by an obs column: the top 5,000 ranks and a sample of
    the rest, so the knee keeps its detail."""
    values = adata.obs[column].to_numpy(dtype=np.float64)
    on_curve = curve.order[: curve.n]
    head, tail = on_curve[:5000], on_curve[5000:]
    if tail.size > MAX_POINTS:
        tail = np.random.default_rng(0).choice(tail, MAX_POINTS, replace=False)
    index = np.concatenate([head, tail])
    index = index[~np.isnan(values[index])]
    ranks = curve.rank_of()[index]
    x, y = (curve.totals[index], ranks) if rotated else (ranks, curve.totals[index])
    return ax.scatter(x, y, c=values[index], s=3, linewidths=0, cmap="viridis", zorder=2)


def track_medians(values: np.ndarray, curve, bins: int = 100):
    """Medians of `values` along the rank curve in log-spaced bins of rank: (bin edges in rank, the
    median of each bin; NaN where a bin has no value)."""
    n = curve.n
    edges = np.unique(np.round(np.logspace(0, np.log10(max(n, 1)), bins + 1)).astype(np.int64))
    edges[-1] = n
    ordered = values[curve.order[:n]]
    medians = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        chunk = ordered[lo - 1:hi]
        chunk = chunk[~np.isnan(chunk)]
        medians.append(float(np.median(chunk)) if chunk.size else np.nan)
    return edges, np.array(medians)


def _track(ax, adata, curve, column: str, bins: int) -> None:
    edges, medians = track_medians(adata.obs[column].to_numpy(dtype=np.float64), curve, bins)
    ax.stairs(medians, edges, color=PALETTE[2], lw=1.4, baseline=None)
    ax.set_ylabel(column.replace("_", " "), fontsize="small")
    ax.grid(True, which="major", alpha=0.3)
