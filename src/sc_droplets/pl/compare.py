"""Several calls side by side: `scd.pl.compare`."""

from __future__ import annotations

from typing import Any, Optional, Sequence, Union

import numpy as np

from sc_droplets.pl._data import curve_for, resolve_calls
from sc_droplets.pl.knee import CURVE_COLOUR, PALETTE

#: Log-spaced rank bins of the call strips.
STRIP_BINS = 200


def strip_shares(cells: np.ndarray, curve, bins: int = STRIP_BINS):
    """For one call, the share of each log-spaced bin of rank that it called: (bin edges in rank, the
    shares)."""
    n = curve.n
    edges = np.unique(np.round(np.logspace(0, np.log10(max(n, 1)), bins + 1)).astype(np.int64))
    edges[-1] = n
    ordered = cells[curve.order[:n]].astype(np.float64)
    cumulative = np.concatenate([[0.0], np.cumsum(ordered)])
    called = cumulative[edges[1:]] - cumulative[edges[:-1] - 1]
    sizes = edges[1:] - edges[:-1] + 1
    return edges, called / sizes


def compare(
    adata: Any,
    calls: Union[None, str, Sequence[str]] = None,
    *,
    layer: Optional[str] = None,
    feature_types: Union[None, str, Sequence[str]] = None,
    fig: Any = None,
):
    """Draw several calls side by side; return the Figure.

    Left: the barcode-rank curve with each call's cut, and under it one strip per call, shaded by
    the share of each log-spaced bin of rank the call took, so where the calls part is plain.
    Right: the number of cells per call, and the pairwise Jaccard index of their cell sets.
    """
    import matplotlib.pyplot as plt

    from sc_droplets.api import agreement

    keys = resolve_calls(adata, calls)
    if not keys:
        raise ValueError("no calls to compare: give method specs or stored keys")
    curve = curve_for(adata, keys, layer=layer, feature_types=feature_types)
    ranks = curve.rank_of()
    cells = {k: adata.obs[k].to_numpy(dtype=bool) for k in keys}
    if fig is None:
        fig = plt.figure(figsize=(14, 5.6 + 0.3 * len(keys)), layout="constrained")
    grid = fig.add_gridspec(2, 2, width_ratios=[1.5, 1], height_ratios=[2.2, 0.6 + 0.25 * len(keys)])
    knee_ax = fig.add_subplot(grid[0, 0])
    strip_ax = fig.add_subplot(grid[1, 0], sharex=knee_ax)
    bar_ax = fig.add_subplot(grid[0, 1])
    jaccard_ax = fig.add_subplot(grid[1, 1])

    x, y = curve.steps()
    knee_ax.plot(x, y, color=CURVE_COLOUR, lw=1.5, label="all barcodes")
    knee_ax.set_xscale("log")
    knee_ax.set_yscale("log")
    knee_ax.set_ylabel("count per barcode")
    for i, key in enumerate(keys):
        on_curve = cells[key] & (ranks > 0)
        if on_curve.any():
            knee_ax.axhline(float(curve.totals[on_curve].min()), color=PALETTE[i % len(PALETTE)], ls="--", lw=1,
                            label=f"{key}: {int(cells[key].sum()):,}")
    knee_ax.legend(frameon=False, fontsize="small")
    knee_ax.grid(True, which="major", alpha=0.3)
    knee_ax.tick_params(labelbottom=False)

    for i, key in enumerate(keys):
        edges, shares = strip_shares(cells[key], curve)
        strip_ax.pcolormesh(edges, [i, i + 0.85], shares[np.newaxis, :], cmap="Blues", vmin=0, vmax=1, shading="flat")
    strip_ax.set_yticks(np.arange(len(keys)) + 0.425, keys, fontsize="small")
    strip_ax.set_ylim(len(keys), 0)
    strip_ax.set_xlabel("Barcode rank")
    strip_ax.set_title("share of each rank bin called", fontsize="small", loc="left")

    counts = [int(cells[k].sum()) for k in keys]
    bars = bar_ax.barh(np.arange(len(keys)), counts, color=[PALETTE[i % len(PALETTE)] for i in range(len(keys))])
    bar_ax.bar_label(bars, labels=[f"{c:,}" for c in counts], fontsize="small", padding=3)
    bar_ax.set_yticks(np.arange(len(keys)), keys, fontsize="small")
    bar_ax.invert_yaxis()
    bar_ax.set_xlabel("cells")
    bar_ax.set_xlim(0, max(counts + [1]) * 1.2)
    bar_ax.grid(True, axis="x", alpha=0.3)

    matrix = agreement(adata, keys).to_numpy()
    jaccard_ax.imshow(matrix, cmap="viridis", vmin=0, vmax=1)
    for (a, b), value in np.ndenumerate(matrix):
        jaccard_ax.text(b, a, f"{value:.2f}", ha="center", va="center", fontsize="x-small",
                        color="white" if value < 0.6 else "black")
    jaccard_ax.set_xticks(np.arange(len(keys)), keys, fontsize="x-small", rotation=45, ha="right")
    jaccard_ax.set_yticks(np.arange(len(keys)), keys, fontsize="x-small")
    jaccard_ax.set_title("Jaccard index of the cell sets", fontsize="small")
    return fig
