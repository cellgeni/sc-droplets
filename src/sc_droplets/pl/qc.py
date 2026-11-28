"""Quality measures against the calls: `scd.pl.qc`."""

from __future__ import annotations

from typing import Any, Optional, Sequence, Union

import numpy as np

from sc_droplets.pl._data import agreement_groups, droplets_for, resolve_calls
from sc_droplets.pl.knee import GROUP_COLOURS

#: Points drawn in each scatter panel, at most.
MAX_POINTS = 20000
_DRAW_ORDER = ("no call", "not a cell", "some calls", "every call", "cell")


def disagreement_band(adata: Any, keys: Sequence[str], totals: np.ndarray):
    """The range of totals where the calls disagree: from the lowest of their smallest cell totals
    to the highest. None when they agree on every barcode."""
    smallest = [float(totals[adata.obs[k].to_numpy(dtype=bool)].min()) for k in keys if adata.obs[k].any()]
    if len(smallest) < 2 or min(smallest) == max(smallest):
        return None
    return min(smallest), max(smallest)


def _scatter(ax, x, y, groups, shown) -> None:
    index = np.flatnonzero(shown)
    if index.size > MAX_POINTS:
        index = np.sort(np.random.default_rng(0).choice(index, MAX_POINTS, replace=False))
    for name in _DRAW_ORDER:
        pick = index[groups[index] == name]
        if pick.size:
            ax.scatter(x[pick], y[pick], s=3, linewidths=0, color=GROUP_COLOURS[name], zorder=2 if name != "no call" else 1,
                       label=f"{name}: {int(((groups == name) & shown).sum()):,}")
    ax.legend(frameon=False, fontsize="small", markerscale=3)
    ax.grid(True, alpha=0.3)


def qc(
    adata: Any,
    calls: Union[None, str, Sequence[str]] = None,
    *,
    min_counts: float = 10,
    layer: Optional[str] = None,
    feature_types: Union[None, str, Sequence[str]] = None,
    fig: Any = None,
):
    """Draw the barcodes' quality measures against the calls; return the Figure.

    Panels: log10 genes detected against log10 count, each barcode coloured by how many calls took
    it; the mitochondrial fraction against count, when `obs` has it (`scd.add_metrics`); and, with
    several calls that disagree, the genes per count (log10 genes − log10 count) of the barcodes in
    the barcodes some calls took, beside those every call took and those none took in the band of
    totals where the calls disagree. Barcodes below `min_counts` are left out.
    """
    import matplotlib.pyplot as plt

    keys = resolve_calls(adata, calls)
    if not keys:
        raise ValueError("no calls to draw: give method specs or stored keys")
    droplets = droplets_for(adata, keys, layer=layer, feature_types=feature_types)
    genes = droplets.values("genes")
    totals = droplets.totals
    shown = totals >= max(min_counts, 1)
    cells = np.column_stack([adata.obs[k].to_numpy(dtype=bool) for k in keys])
    groups = agreement_groups(cells)
    log_total = np.log10(np.where(totals > 0, totals, 1))
    log_genes = np.log10(np.where(genes > 0, genes, 1))
    band = disagreement_band(adata, keys, totals)
    has_mito = "mito_fraction" in adata.obs.columns
    panels = 1 + has_mito + (band is not None)
    if fig is None:
        fig = plt.figure(figsize=(5.4 * panels, 4.8), layout="constrained")
    axes = np.atleast_1d(fig.subplots(1, panels))

    ax = axes[0]
    _scatter(ax, log_total, log_genes, groups, shown)
    ax.set_xlabel("log10 count per barcode")
    ax.set_ylabel("log10 genes detected")
    ax.set_title("genes against count", fontsize="medium")
    if band is not None:
        for edge in band:
            ax.axvline(np.log10(edge), color="#555555", ls=":", lw=1)

    i = 1
    if has_mito:
        ax = axes[i]
        mito = adata.obs["mito_fraction"].to_numpy(dtype=np.float64)
        _scatter(ax, log_total, mito, groups, shown & ~np.isnan(mito))
        ax.set_xlabel("log10 count per barcode")
        ax.set_ylabel("mitochondrial fraction")
        ax.set_title("mitochondrial fraction against count", fontsize="medium")
        i += 1

    if band is not None:
        ax = axes[i]
        inside = shown & (totals >= band[0]) & (totals <= band[1])
        complexity = log_genes - log_total
        # The disputed barcodes beside the agreed cells and the barcodes no call took in the band.
        picks = {"every call": shown & (groups == "every call"), "some calls": shown & (groups == "some calls"),
                 "no call": inside & (groups == "no call")}
        used = np.logical_or.reduce(list(picks.values()))
        bins = np.linspace(complexity[used].min(), complexity[used].max(), 31) if used.any() else 10
        for name, mask in picks.items():
            pick = complexity[mask]
            if pick.size:
                where = " in the band" if name == "no call" else ""
                ax.hist(pick, bins=bins, density=True, histtype="step", lw=1.5, color=GROUP_COLOURS[name],
                        label=f"{name}{where}: {pick.size:,}, median {np.median(pick):.2f}")
        ax.set_xlabel("log10 genes − log10 count")
        ax.set_ylabel("density")
        ax.set_title(f"the barcodes the calls disagree on ({band[0]:g} – {band[1]:g} counts)", fontsize="medium")
        ax.legend(frameon=False, fontsize="small")
        ax.grid(True, alpha=0.3)
    return fig
