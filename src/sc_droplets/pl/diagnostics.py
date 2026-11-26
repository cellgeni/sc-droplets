"""Each method's own view of a call, drawn from what the call stored.

Threshold methods (OrdMag, CellRanger2.2, TopCells, the awk script, top N,
fixed threshold, DropletUtils barcodeRanks and defaultDrops, UMI-tools): the
rank curve, the cells, and the cut, with the robust maximum it came from and,
for OrdMag, the spread of its bootstraps; for barcodeRanks, the knee, the
inflection and what found them (the fitted spline, the chord, or the knee's
window). UMI-tools adds a second panel: for distance, the cumulative curve of
the last pass with its chord and each pass's knee; for density, the kernel
density of log10 counts with every valley and the chosen one. Kneedle adds its
difference curve with the maxima, their thresholds and the knee; cellsweep its
rule's own curve (the second derivative, the Kneedle difference curve, or the
two Gaussian components and their entropy).

CellBender's priors: the rank curve with the ranks of its classes (probable cells, the other
analysed droplets, the surely-empty droplets) and the priors, and the distribution the priors
came from (the rounded log counts whose mode is the empty prior; or, from 0.3, the log counts,
Otsu's cut and the last kernel density).

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
    diag = dict(record.get("diagnostics", {}))
    if diag.get("prior_rule") in ("mode-median", "otsu-kde"):
        if fig is None:
            fig = plt.figure(figsize=(12.8, 4.8))
        axes = fig.subplots(1, 2)
        _cellbender_curve(axes[0], record)
        _cellbender_prior_panel(axes[1], record)
        fig.suptitle(f"{record.get('title', key)}: {int(record.get('n_cells', 0)):,} probable cells\n{key}", fontsize="medium")
        fig.tight_layout()
        return fig
    second = (_cellsweep_panel if diag.get("rule") in ("threshold", "kneed", "mx_filter")
              else _cumulative_panel if "passes" in diag else _density_panel if "density" in diag
              else _difference_panel if "difference_x" in diag else None)
    if fig is None:
        fig = plt.figure(figsize=(13, 4.2) if has_test else (12.8, 4.8) if second else (6.4, 4.8))
    curve = curve_for(adata, [key])
    if has_test:
        axes = fig.subplots(1, 3)
        _classes_on_curve(axes[0], adata, key, record, curve)
        _loglik_panel(axes[1], adata, key)
        _pvalue_panel(axes[2], adata, key, record)
    elif second is not None:
        axes = fig.subplots(1, 2)
        _threshold_on_curve(axes[0], adata, key, record, curve)
        second(axes[1], record, curve)
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


def _cumulative_panel(ax, record, curve) -> None:
    """UMI-tools distance: the last pass's cumulative curve, its chord, and each pass's knee."""
    diag = dict(record.get("diagnostics", {}))
    passes = np.asarray(diag["passes"], dtype=np.int64).reshape(-1, 2)
    ax.set_title(diag.get("curve_title", "Cumulative counts (UMI-tools distance)"), fontsize="medium")
    ax.set_xlabel("Barcode rank (0-based)")
    ax.set_ylabel(f"Cumulative {diag.get('unit', 'umis')}")
    # Drawn from the stored curve, not the AnnData's barcodes: those may be only the cells.
    x = np.asarray(diag["cumulative_x"], dtype=np.int64)
    y = np.asarray(diag["cumulative_y"], dtype=float)
    at = dict(zip(x.tolist(), y.tolist()))
    points = int(passes[-1, 0])
    ax.plot(x, y, color=CURVE_COLOUR, lw=1.5, label=f"last pass: first {points:,} barcodes")
    ax.plot([0, points - 1], [at[0], at[points - 1]], color=PALETTE[3], ls="-.", lw=1, label="chord")
    for i, (_, knee) in enumerate(passes):
        if knee in at:
            last = i == passes.shape[0] - 1
            ax.plot([knee], [at[int(knee)]], marker="o", ms=7 if last else 4, color=PALETTE[1] if last else "#999999",
                    ls="none", label=f"knee: {knee + 1:,} cells" if last else None, zorder=3)
    ax.legend(frameon=False, fontsize="small")
    ax.grid(True, alpha=0.3)


def _density_panel(ax, record, curve) -> None:
    """UMI-tools density: the kernel density of log10 counts, its valleys, and the chosen one."""
    diag = dict(record.get("diagnostics", {}))
    grid = np.asarray(diag["grid"], dtype=float)
    ax.plot(grid, np.asarray(diag["density"], dtype=float), color=CURVE_COLOUR, lw=1.5, label="kernel density")
    valleys = np.asarray(diag.get("valleys", []), dtype=float)
    passing = np.asarray(diag.get("valley_passing", []), dtype=np.int64)
    chosen = diag.get("chosen")
    other_labelled = False
    for t, n_pass in zip(valleys, passing):
        picked = chosen is not None and t == float(chosen)
        if picked:
            label = f"chosen valley: {n_pass:,} cells"
        else:
            label = None if other_labelled else "other valleys"
            other_labelled = True
        ax.axvline(t, color=PALETTE[1] if picked else "#999999", ls="--", lw=1.2 if picked else 0.8, label=label)
    ax.set_title("Density of log10 counts (UMI-tools density)", fontsize="medium")
    ax.set_xlabel("log10 count per barcode")
    ax.set_ylabel("density")
    ax.legend(frameon=False, fontsize="small")
    ax.grid(True, alpha=0.3)


def _difference_panel(ax, record, curve) -> None:
    """Kneedle: the difference curve, its maxima with their thresholds, and the knee's maximum."""
    diag = dict(record.get("diagnostics", {}))
    thresholds = dict(record.get("thresholds", {}))
    ax.plot(np.asarray(diag["difference_x"], dtype=float), np.asarray(diag["difference_y"], dtype=float),
            color=CURVE_COLOUR, lw=1.2, label="difference curve")
    mx = np.asarray(diag.get("maxima_x", []), dtype=float)
    if mx.size:
        ax.scatter(mx, np.asarray(diag["maxima_y"], dtype=float), s=10, color="#999999", zorder=2, label="local maxima")
        ax.scatter(mx, np.asarray(diag["maxima_threshold"], dtype=float), s=6, marker="_", color=PALETTE[3], zorder=2,
                   label=f"thresholds (S = {thresholds.get('sensitivity', float('nan')):g})")
    ax.plot([float(diag["knee_x"])], [float(diag["knee_y"])], marker="o", ms=7, color=PALETTE[1], ls="none", zorder=3,
            label=f"knee: {int(thresholds.get('knee_rank', 0)):,} cells")
    ax.set_title("Kneedle difference curve", fontsize="medium")
    ax.set_xlabel("x, scaled to [0, 1]")
    ax.set_ylabel("y − x, scaled")
    ax.legend(frameon=False, fontsize="small")
    ax.grid(True, alpha=0.3)


def _cellsweep_panel(ax, record, curve) -> None:
    """cellsweep: the curve its automatic rule cut on, and the cut."""
    diag = dict(record.get("diagnostics", {}))
    cutoff = float(dict(record.get("thresholds", {})).get("umis", float("nan")))
    rule = str(diag["rule"])
    if rule == "threshold":
        ranks = np.asarray(diag["curve_index"], dtype=float) + 1
        d2 = np.asarray(diag["second_derivative"], dtype=float)
        ax.plot(ranks, d2, color=CURVE_COLOUR, lw=1.2, label="second derivative of smoothed log10 counts")
        totals = np.asarray(diag["rank_curve"], dtype=float)
        at = ranks[np.argmin(np.abs(totals - cutoff))]
        ax.axvline(at, color=PALETTE[1], ls="--", lw=1, label=f"cut: {cutoff:g}")
        ax.set_xscale("log")
        ax.set_xlabel("Barcode rank")
        ax.set_ylabel("second derivative")
    elif rule == "kneed":
        ax.plot(np.asarray(diag["difference_x"], dtype=float), np.asarray(diag["difference_y"], dtype=float),
                color=CURVE_COLOUR, lw=1.2, label="Kneedle difference curve")
        ax.plot([float(diag["knee_x"])], [float(diag["knee_y"])], marker="o", ms=7, color=PALETTE[1], ls="none",
                label=f"knee: cut {cutoff:g}", zorder=3)
        ax.set_xlabel("rank, scaled to [0, 1]")
        ax.set_ylabel("y − x, scaled")
    else:
        from scipy.stats import norm

        means = np.asarray(diag["mixture_means"], dtype=float)
        sds = np.asarray(diag["mixture_sds"], dtype=float)
        weights = np.asarray(diag["mixture_weights"], dtype=float)
        v = np.log1p(np.asarray(diag["mixture_totals"], dtype=float))
        grid = np.linspace(v.min(), v.max(), 400)
        for k in range(means.size):
            ax.plot(grid, weights[k] * norm.pdf(grid, means[k], sds[k]), lw=1.2, color=PALETTE[2 + k], label=f"component {k + 1}")
        twin = ax.twinx()
        twin.plot(v, np.asarray(diag["mixture_entropy"], dtype=float), color="#999999", lw=0.8)
        twin.set_ylabel("entropy of the posterior", color="#777777")
        ax.axvline(np.log1p(cutoff), color=PALETTE[1], ls="--", lw=1, label=f"cut (highest entropy): {cutoff:g}")
        ax.set_xlabel("log(1 + count)")
        ax.set_ylabel("weighted density")
    ax.set_title(f"cellsweep: {rule}", fontsize="medium")
    ax.legend(frameon=False, fontsize="small")
    ax.grid(True, alpha=0.3)


CELLBENDER_BANDS = (("cell", PALETTE[0]), ("transition", PALETTE[3]), ("empty", "#999999"))


def _cellbender_curve(ax, record) -> None:
    """CellBender: the rank curve it ranks on, with the ranks of its classes and the priors."""
    diag = dict(record.get("diagnostics", {}))
    thresholds = dict(record.get("thresholds", {}))
    # The curve and the class ranks are stored, not read from the AnnData: it may hold only the
    # cells, and from 0.3 CellBender ranks on the totals over the features it keeps.
    ax.plot(np.asarray(diag["rank_curve_x"], dtype=float), np.asarray(diag["rank_curve_y"], dtype=float),
            color=CURVE_COLOUR, lw=1.5, drawstyle="steps-post", zorder=1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Barcode rank")
    trimmed = "n_features" in thresholds
    ax.set_ylabel("UMIs per barcode, over the features kept" if trimmed else "UMIs per barcode")
    ax.grid(True, which="major", alpha=0.3)
    ends = [0] + [int(v) for v in np.asarray(diag["class_ranks"], dtype=np.int64)]
    for (name, colour), start, end in zip(CELLBENDER_BANDS, ends[:-1], ends[1:]):
        if end > start:
            ax.axvspan(start + 1, end, color=colour, alpha=0.15, label=f"{name}: {end - start:,}")
    for name, label, colour, style in (("cell_counts", "cell prior", PALETTE[0], "--"),
                                       ("empty_counts", "empty prior", PALETTE[1], "--"),
                                       ("low_count_cutoff", "low cut", "#777777", ":")):
        if name in thresholds:
            ax.axhline(thresholds[name], color=colour, ls=style, lw=1, label=f"{label}: {thresholds[name]:.6g}")
    ax.legend(frameon=False, fontsize="small")


def _cellbender_prior_panel(ax, record) -> None:
    """CellBender: the distribution the priors came from (0.1, 0.2: the rounded log counts whose
    mode is the empty prior; 0.3: log counts, the last kernel density, Otsu's cut)."""
    diag = dict(record.get("diagnostics", {}))
    thresholds = dict(record.get("thresholds", {}))
    if diag["prior_rule"] == "mode-median":
        ax.bar(np.asarray(diag["mode_values"], dtype=float), np.asarray(diag["mode_counts"], dtype=float),
               width=0.08, color=CURVE_COLOUR, label="barcodes above the low-count threshold")
        lines = (("empty_counts", "empty prior (mode)", PALETTE[1]), ("cell_counts", "cell prior", PALETTE[0]))
        ax.set_xlabel("log(1 + count), rounded to 0.1")
        for name, label, colour in lines:
            ax.axvline(np.log1p(thresholds[name]), color=colour, ls="--", lw=1, label=f"{label}: {thresholds[name]:g}")
        ax.set_title("CellBender 0.1–0.2: mode and median", fontsize="medium")
    else:
        x = np.asarray(diag["log_histogram_x"], dtype=float)
        ax.step(x, np.asarray(diag["log_histogram_y"], dtype=float), where="mid", color=CURVE_COLOUR,
                label="barcodes above the low-count threshold")
        twin = ax.twinx()
        twin.plot(np.asarray(diag["kde_x"], dtype=float), np.asarray(diag["kde_density"], dtype=float),
                  color=PALETTE[2], lw=1.2)
        twin.set_ylabel("kernel density of the resampled empty droplets", color=PALETTE[2])
        for name, label, colour, style in (("empty_counts", "empty prior (density peak)", PALETTE[1], "--"),
                                           ("empty_count_upper_limit", "its upper limit", PALETTE[1], ":"),
                                           ("cell_count_low_limit", "Otsu's cut", PALETTE[3], "-."),
                                           ("cell_counts", "cell prior", PALETTE[0], "--")):
            if name in thresholds:
                ax.axvline(np.log(thresholds[name]), color=colour, ls=style, lw=1, label=f"{label}: {thresholds[name]:.6g}")
        ax.set_xlabel("log count")
        ax.set_title("CellBender 0.3+: Otsu and kernel density", fontsize="medium")
    ax.set_ylabel("barcodes")
    ax.legend(frameon=False, fontsize="small")
    ax.grid(True, alpha=0.3)


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
