"""Each method's own view of a call, drawn from what the call stored.

Threshold methods (OrdMag, CellRanger2.2, TopCells, the awk script, top N,
fixed threshold, DropletUtils barcodeRanks and defaultDrops): the rank curve,
the cells, and the cut, with the robust maximum it came from and, for OrdMag,
the spread of its bootstraps; for barcodeRanks, the knee, the inflection and
what found them (the fitted spline, the chord, or the knee's window).

Test methods (Cell Ranger 3+, STARsolo EmptyDrops_CR, DropletUtils emptyDrops
and emptyDropsCellRanger): three panels. The rank curve with each barcode's
class (initial or retained, rescued, tested and kept empty) and the ambient
ranks shaded; the candidates' log-likelihood under the ambient profile against
their total, called ones marked; and the candidates' p-values with the FDR cut.
"""

from __future__ import annotations

from typing import Any, Dict, Optional

import numpy as np

from sc_droplets.pl._data import curve_for, resolve_calls
from sc_droplets.pl.knee import CURVE_COLOUR, PALETTE
from sc_droplets.result import UNS_KEY

CLASS_COLOURS = {"initial": PALETTE[0], "retained": PALETTE[0], "rescued": PALETTE[1], "candidate": "#999999", "tested": "#999999"}
KEPT_EMPTY = {"candidate": "tested, kept empty", "tested": "tested, kept empty"}


def diagnostics(adata: Any, call: str, *, fig: Any = None):
    """Draw one call's diagnostic view and return the matplotlib Figure.

    `call` is a stored key or a method spec (run first if not stored).
    """
    import matplotlib.pyplot as plt

    key = resolve_calls(adata, [call])[0]
    record = adata.uns[UNS_KEY][key]
    has_test = f"{key}_pvalue" in adata.obs.columns
    if fig is None:
        fig = plt.figure(figsize=(13, 4.2) if has_test else (6.4, 4.8))
    curve = curve_for(adata, [key])
    if has_test:
        axes = fig.subplots(1, 3)
        _classes_on_curve(axes[0], adata, key, record, curve)
        _loglik_panel(axes[1], adata, key)
        _pvalue_panel(axes[2], adata, key, record)
    else:
        ax = fig.subplots(1, 1)
        _threshold_on_curve(ax, adata, key, record, curve)
    fig.suptitle(f"{record.get('title', key)}: {int(record.get('n_cells', 0)):,} cells\n{key}", fontsize="medium")
    fig.tight_layout()
    return fig


def _curve(ax, curve) -> None:
    x, y = curve.steps()
    ax.plot(x, y, color=CURVE_COLOUR, lw=1.5, zorder=1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Barcode rank")
    ax.set_ylabel("UMIs per barcode")
    ax.grid(True, which="major", alpha=0.3)


def _threshold_on_curve(ax, adata, key, record, curve) -> None:
    _curve(ax, curve)
    cells = adata.obs[key].to_numpy(dtype=bool)
    ranks = curve.rank_of()
    on_curve = cells & (ranks > 0)
    if on_curve.any():
        index = np.flatnonzero(on_curve)
        if index.size > 20000:
            index = np.random.default_rng(0).choice(index, 20000, replace=False)
        ax.scatter(ranks[index], curve.totals[index], s=4, color=PALETTE[0], label="cells", zorder=2, linewidths=0)
    thresholds: Dict[str, float] = dict(record.get("thresholds", {}))
    diag: Dict[str, Any] = dict(record.get("diagnostics", {}))
    if "umis" in thresholds:
        ax.axhline(thresholds["umis"], color=PALETTE[1], ls="--", lw=1, label=f"cut: {thresholds['umis']:g}")
    robust = thresholds.get("robust_max")
    group = _first_group(diag)
    if robust is None and group is not None and "baseline" in group:
        robust = float(group["baseline"])
    if robust is not None:
        ax.axhline(robust, color=PALETTE[2], ls=":", lw=1, label=f"robust maximum: {robust:g}")
    if group is not None and "bootstrap" in group:
        boot = np.asarray(group["bootstrap"], dtype=float)
        ax.axvspan(boot.min(), boot.max(), color=PALETTE[0], alpha=0.12, label=f"OrdMag bootstraps: {int(boot.min())}–{int(boot.max())} cells")
    _knee_geometry(ax, thresholds, diag)
    ax.legend(frameon=False, fontsize="small")


def _knee_geometry(ax, thresholds: Dict[str, float], diag: Dict[str, Any]) -> None:
    """barcodeRanks: knee and inflection, and the spline, chord or window that found the knee."""
    for name, colour in (("knee", PALETTE[3]), ("inflection", PALETTE[4])):
        if name in thresholds and thresholds.get("umis") != thresholds[name]:
            ax.axhline(thresholds[name], color=colour, ls=":", lw=1, label=f"{name}: {thresholds[name]:.6g}")
    if "fitted_x" in diag:
        ax.plot(10 ** np.asarray(diag["fitted_x"], dtype=float), 10 ** np.asarray(diag["fitted_y"], dtype=float),
                color=PALETTE[3], lw=1.2, label="smoothing spline (df 20)", zorder=3)
    if "chord" in diag:
        x0, y0, x1, y1 = (float(v) for v in diag["chord"])
        ax.plot([10**x0, 10**x1], [10**y0, 10**y1], color=PALETTE[3], lw=1, ls="-.", label="chord", zorder=3)
    if "knee_window" in diag:
        left, right = (10 ** float(v) for v in diag["knee_window"])
        ax.axvspan(left, right, color=PALETTE[3], alpha=0.15, label="knee's window")


def _first_group(diag: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    for name in ("groups", "ordmag"):
        groups = diag.get(name)
        if groups:
            return dict(next(iter(groups.values())))
    return None


def _classes_on_curve(ax, adata, key, record, curve) -> None:
    _curve(ax, curve)
    ranks = curve.rank_of()
    classes = adata.obs[f"{key}_class"].astype(str).to_numpy()
    for name, colour in CLASS_COLOURS.items():
        members = np.flatnonzero((classes == name) & (ranks > 0))
        if members.size > 20000:
            members = np.random.default_rng(0).choice(members, 20000, replace=False)
        if members.size == 0 and name in ("retained", "tested", "initial", "candidate"):
            continue
        label = KEPT_EMPTY.get(name, name)
        ax.scatter(ranks[members], curve.totals[members], s=5, color=colour, label=f"{label} ({members.size:,})", zorder=2, linewidths=0)
    ambient = _ambient_ranks(record)
    if ambient is not None:
        ax.axvspan(ambient[0] + 1, ambient[1], color="#cccccc", alpha=0.4, label="ambient ranks")
    thresholds = dict(record.get("thresholds", {}))
    if "retain" in thresholds:
        ax.axhline(thresholds["retain"], color=PALETTE[0], ls="--", lw=1, label=f"retain (knee): {thresholds['retain']:.6g}")
    if "lower" in thresholds:
        ax.axhline(thresholds["lower"], color="#777777", ls=":", lw=1, label=f"lower: {thresholds['lower']:g}")
    ax.legend(frameon=False, fontsize="small")
    ax.set_title("Classes on the rank curve", fontsize="medium")


def _ambient_ranks(record) -> Optional[np.ndarray]:
    diag = dict(record.get("diagnostics", {}))
    if "ambient_ranks" in diag:
        return np.asarray(diag["ambient_ranks"], dtype=float)
    groups = diag.get("nonambient")
    if groups:
        first = dict(next(iter(groups.values())))
        if "ambient_ranks" in first:
            return np.asarray(first["ambient_ranks"], dtype=float)
    return None


def _loglik_panel(ax, adata, key) -> None:
    column = f"{key}_loglik" if f"{key}_loglik" in adata.obs.columns else f"{key}_logprob"
    loglik = adata.obs[column].to_numpy(dtype=float)
    tested = np.isfinite(loglik)
    totals = curve_for(adata, [key]).totals
    called = adata.obs[key].to_numpy(dtype=bool)
    for mask, colour, label in ((tested & ~called, CLASS_COLOURS["candidate"], "kept empty"), (tested & called, CLASS_COLOURS["rescued"], "called")):
        ax.scatter(totals[mask], loglik[mask], s=6, color=colour, label=f"{label} ({int(mask.sum()):,})", linewidths=0)
    ax.set_xscale("log")
    ax.set_xlabel("UMIs per candidate")
    ax.set_ylabel("log-likelihood under the ambient profile")
    ax.set_title("Candidates", fontsize="medium")
    ax.legend(frameon=False, fontsize="small")
    ax.grid(True, alpha=0.3)


def _pvalue_panel(ax, adata, key, record) -> None:
    pvalues = adata.obs[f"{key}_pvalue"].to_numpy(dtype=float)
    pvalues = pvalues[np.isfinite(pvalues)]
    if pvalues.size:
        ax.hist(pvalues, bins=50, range=(0, 1), color=PALETTE[5], edgecolor="white")
    fdr = dict(record.get("thresholds", {})).get("fdr")
    title = "p-values of the candidates"
    if fdr is not None:
        title += f" (cells: adjusted p ≤ {fdr:g})"
    ax.set_title(title, fontsize="medium")
    ax.set_xlabel("Monte Carlo p-value")
    ax.set_ylabel("candidates")
