"""Draw the figures of the Explained pages and write the data their widgets read.

    uv run python docs/_scripts/explained.py               # every page
    uv run python docs/_scripts/explained.py --page knee   # one page (see PAGES)

Runs on 10x's raw matrices of PBMC 1k (whole cells) and mouse spinal cord nuclei, which
tests/datasets.py downloads on first use. Writes docs/assets/explained/*.png and
docs/assets/explained/curves/*.json, and prints every number the pages quote, so the prose can be
checked against it.

The pages derive identities: on the cumulative curve, the farthest point from a chord is the last
barcode above a mean. Each is asserted here against the methods' own functions on both matrices.
A change to a method that would make a page wrong stops this script, instead of leaving the page
quietly out of date.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tests"))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import sc_droplets as scd  # noqa: E402
from sc_droplets.counts import Droplets  # noqa: E402
from sc_droplets.methods import alevinfry, dropseq, umitools, zumis  # noqa: E402
from sc_droplets.pl.knee import PALETTE  # noqa: E402
from datasets import filtered_barcodes, raw_path  # noqa: E402

OUT = ROOT / "docs" / "assets" / "explained"
DATA = OUT / "curves"
DPI = 100

BLUE, VERMILLION, GREEN, PURPLE, ORANGE, SKY = PALETTE[:6]
GREY = "#8c8c8c"
LIGHT = "#d0d0d0"
TRUTH_COLOURS = {"cell": BLUE, "empty": ORANGE, "low": GREY}
TRUTH_LABELS = {"cell": "cells", "empty": "empty droplets (ambient RNA)", "low": "error barcodes"}

plt.rcParams.update({
    "font.size": 10,
    "axes.titlesize": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
    "figure.dpi": DPI,
})

#: Numbers the pages quote, printed at the end.
NUMBERS: Dict[str, object] = {}


def note(key: str, value) -> None:
    NUMBERS[key] = value


# ---------------------------------------------------------------------------
# Samples
# ---------------------------------------------------------------------------


@dataclass
class Sample:
    key: str
    label: str
    adata: object
    totals: np.ndarray  # per barcode, in matrix order (UMIs)
    c: np.ndarray  # totals above 0, decreasing, ties in matrix order (int64)
    tenx_ranks: np.ndarray  # 1-based ranks of 10x's cells in `c`
    calls: Dict[str, object] = field(default_factory=dict)

    @property
    def n(self) -> int:
        return int(self.c.size)

    @property
    def cum(self) -> np.ndarray:
        return np.cumsum(self.c)

    @property
    def ranks(self) -> np.ndarray:
        return np.arange(1, self.n + 1, dtype=np.float64)

    def call(self, spec: str):
        if spec not in self.calls:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                self.calls[spec] = scd.call(self.adata, spec, inplace=False)
        return self.calls[spec]


SAMPLES = {"pbmc": ("pbmc_1k_v3_3_0", "PBMC 1k"), "nuclei": ("spinal_cord_10_0", "spinal cord nuclei")}


def load(key: str) -> Sample:
    name, label = SAMPLES[key]
    adata = scd.read(raw_path(name))
    totals = Droplets.from_anndata(adata).values("umis")
    index, counts = umitools._observed(totals)  # the barcodes with a count, as the knee rules see them
    order = np.argsort(-counts, kind="stable")
    rank_of = np.full(totals.size, 0, dtype=np.int64)
    rank_of[index[order]] = np.arange(1, order.size + 1)
    tenx = adata.obs_names.isin(filtered_barcodes(name))
    sample = Sample(key, label, adata, totals, counts[order].astype(np.int64), np.sort(rank_of[tenx]))
    note(f"{key}.barcodes_with_a_count", sample.n)
    note(f"{key}.barcodes", int(totals.size))
    note(f"{key}.tenx_cells", int(tenx.sum()))
    note(f"{key}.largest", int(sample.c[0]))
    return sample


def thin(n: int, most: int = 4000, log: bool = True) -> np.ndarray:
    """About `most` 0-based indices into a curve of n points, spaced evenly in log rank (or rank)."""
    if log:
        picks = np.round(np.logspace(0, math.log10(n), most)).astype(np.int64) - 1
    else:
        picks = np.linspace(0, n - 1, most).round().astype(np.int64)
    return np.unique(np.concatenate([picks, [0, n - 1]]))


def farthest(x: np.ndarray, y: np.ndarray, a: int = 0, b: Optional[int] = None) -> int:
    """The first index in [a, b] farthest from the chord through points a and b (b inclusive)."""
    b = x.size - 1 if b is None else b
    area = np.abs((x[b] - x[a]) * (y[a:b + 1] - y[a]) - (y[b] - y[a]) * (x[a:b + 1] - x[a]))
    return a + int(np.argmax(area))


def save(figure, name: str) -> None:
    figure.savefig(OUT / f"{name}.png", dpi=DPI, bbox_inches="tight")
    plt.close(figure)
    print(f"  {name}.png")


def rank_axes(ax, xlabel: str = "rank (log)", ylabel: str = "UMIs (log)") -> None:
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)


def draw_curve(ax, s: Sample, colour: str = GREY, lw: float = 1.4, label: Optional[str] = None, most: int = 4000):
    picks = thin(s.n, most)
    ax.plot(s.ranks[picks], s.c[picks], color=colour, lw=lw, label=label, zorder=1)


def draw_tenx(ax, s: Sample, label: bool = True):
    r = s.tenx_ranks
    ax.scatter(r, s.c[r - 1], s=3, color=BLUE, linewidths=0, zorder=2,
               label=f"10x's cells ({r.size:,})" if label else None)


# ---------------------------------------------------------------------------
# The identities the pages derive, checked against the methods' own code
# ---------------------------------------------------------------------------


def umitools_pass(c: np.ndarray, cum: np.ndarray, m: int) -> int:
    """UMI-tools' knee on the first m points of the cumulative curve, by the identity: the number
    of barcodes after the first whose count exceeds the mean of c_1 .. c_{m-1}."""
    mean = (cum[m - 1] - cum[0]) / (m - 1)
    return int(np.searchsorted(-c[1:m], -mean, side="left"))


def alevinfry_pass(c: np.ndarray, cum: np.ndarray, m: int) -> int:
    """alevin-fry's knee on the first m points: its line ends at x = 1, one step past the last
    point, so the mean divides by m; ties go to the last barcode (>=)."""
    mean = (cum[m - 1] - cum[0]) / m
    return int(np.searchsorted(-c[1:m], -mean, side="right"))


def iterate(step: Callable[[int], int], start: int, grow: Callable[[int], int]) -> List[Tuple[int, int]]:
    """Passes (window, knee) of a refit loop: knee = step(window), window = grow(knee)."""
    knee = step(start)
    passes = [(start, knee)]
    for _ in range(200):
        window = grow(knee)
        new = step(window)
        passes.append((window, new))
        if new == knee:
            break
        knee = new
    return passes


def check_identities(s: Sample) -> None:
    c, cum = s.c, s.cum
    # dropbead: the cells are the barcodes above the mean of the top 15,000.
    top = c[: dropseq.MAX_CELLS].astype(np.float64)
    mean = top.sum() / top.size
    by_rule = dropseq.estimate_cell_number(c.astype(np.float64))
    assert by_rule == int(np.sum(top > mean)), (s.key, "dropbead", by_rule)
    note(f"{s.key}.dropbead.mean", round(mean, 1))
    # the chord widget joins the window's first and last barcodes, as UMI-tools does
    note(f"{s.key}.chord_over_top_15000", umitools_pass(c, cum, dropseq.MAX_CELLS) + 1)
    note(f"{s.key}.dropbead.cells", by_rule)
    # UMI-tools: every pass is the mean rule; the loop is k -> f(min(3k, n)).
    knee, passes = umitools.distance_knee(c)
    for m, k in passes:
        assert umitools_pass(c, cum, m) == k, (s.key, "umitools", m, k)
    assert iterate(lambda m: umitools_pass(c, cum, m), s.n, lambda k: min(3 * k, s.n))[-1][1] == knee
    note(f"{s.key}.umitools.passes", [(int(m), int(k) + 1) for m, k in passes])
    note(f"{s.key}.umitools.first_mean", round((cum[-1] - cum[0]) / (s.n - 1), 2))
    # alevin-fry: the same with m in the denominator, >=, and windows of min(n - 1, 5k).
    knee, passes = alevinfry.get_knee(c)
    for m, k in passes:
        assert alevinfry_pass(c, cum, m) == k, (s.key, "alevin-fry", m, k)
    note(f"{s.key}.alevinfry.passes", [(int(m), int(k) + 1) for m, k in passes])
    # zUMIs: uik on the cumulative curve of the barcodes with at least 100 is the mean rule too.
    floor = int(np.sum(c >= zumis.N_READS_PER_CELL))
    window = c[:floor].astype(np.float64)
    mean = window[1:].mean()
    cut = s.call("zumis-cell-barcodes").thresholds["umis"]
    assert np.sum(c > mean) == np.sum(c >= cut), (s.key, "zumis")
    note(f"{s.key}.zumis.window", floor)
    note(f"{s.key}.zumis.mean", round(mean, 1))
    note(f"{s.key}.zumis.cells", int(np.sum(c >= cut)))


# ---------------------------------------------------------------------------
# The refit loop as a map: k -> f(c k)
# ---------------------------------------------------------------------------


def refit_map(s: Sample, multiplier: float, ks: np.ndarray, form: str = "umitools") -> np.ndarray:
    c, cum = s.c, s.cum
    if form == "umitools":
        return np.array([umitools_pass(c, cum, max(2, min(int(multiplier * k), s.n))) for k in ks])
    return np.array([alevinfry_pass(c, cum, max(2, min(int(multiplier * k), s.n - 1))) for k in ks])


def fixed_points(ks: np.ndarray, g: np.ndarray, smallest: int = 500) -> List[Tuple[int, str]]:
    """Where g crosses the diagonal (on the grid): stable where it crosses from above."""
    side = np.sign(g - ks)
    found = []
    for i in np.flatnonzero(np.diff(side) != 0):
        if ks[i] >= smallest:
            found.append((int(ks[i + 1] if side[i] > 0 else ks[i]), "stable" if side[i] > 0 else "unstable"))
    return found


def note_refits(s: Sample) -> None:
    """Where the refit loop settles for each whole multiplier, and how wide the map's range is."""
    c, cum = s.c, s.cum
    ends = {}
    for multiplier in range(2, 9):
        passes = iterate(lambda m: umitools_pass(c, cum, m), s.n, lambda k: max(2, min(multiplier * k, s.n)))
        ends[multiplier] = passes[-1][1] + 1
    note(f"{s.key}.refit.cells_by_multiplier", ends)
    windows = np.unique(np.round(np.logspace(math.log10(3000), math.log10(s.n), 300)).astype(np.int64))
    knees = [umitools_pass(c, cum, int(m)) + 1 for m in windows]
    note(f"{s.key}.refit.knee_range_windows_3000_to_all", [min(knees), max(knees)])
    knee, passes = umitools.distance_knee(c)
    note(f"{s.key}.umitools.pass_means", [round(float(cum[m - 1] - cum[0]) / (m - 1), 1) for m, _ in passes])
    last = passes[-1][0]
    note(f"{s.key}.umitools.last_window", int(last))
    note(f"{s.key}.umitools.last_mean", round((cum[last - 1] - cum[0]) / (last - 1), 1))
    note(f"{s.key}.umitools.cut", int(c[knee]))


def rho(a: np.ndarray) -> np.ndarray:
    """On a power-law curve c_i = i^-a, the share of a window that lies above its mean."""
    return (1 - a) ** (1 / a)


# ---------------------------------------------------------------------------
# The barcode-rank curve
# ---------------------------------------------------------------------------


def fig_rank_toy() -> None:
    from sc_droplets.simulate import simulate_droplets

    adata = simulate_droplets(seed=0)
    totals = np.asarray(adata.X.sum(axis=1)).ravel()
    truth = adata.obs["truth"].to_numpy()
    order = np.argsort(-totals, kind="stable")
    ranks = np.arange(1, totals.size + 1)
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.8))
    ax = axes[0]
    for kind in ("low", "empty", "cell"):
        at = truth[order] == kind
        ax.scatter(ranks[at], totals[order][at], s=2, color=TRUTH_COLOURS[kind], linewidths=0,
                   label=TRUTH_LABELS[kind], rasterized=True)
    rank_axes(ax)
    ax.set_title("Simulated: the rank curve, coloured by what each barcode holds")
    ax.legend(markerscale=5, loc="lower left")
    ax = axes[1]
    bins = np.linspace(0, math.log10(totals.max()) + 0.05, 70)
    stacks = [np.log10(totals[truth == k]) for k in ("low", "empty", "cell")]
    ax.hist(stacks, bins=bins, stacked=True, color=[TRUTH_COLOURS[k] for k in ("low", "empty", "cell")])
    ax.set_yscale("log")
    ax.set_xlabel("UMIs (log10)")
    ax.set_ylabel("barcodes")
    ax.set_title("The same totals as a histogram")
    fig.tight_layout()
    save(fig, "rank-toy")
    note("toy.n", {k: int(np.sum(truth == k)) for k in ("cell", "empty", "low")})


def fig_rank_anatomy(p: Sample) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    draw_curve(ax, p)
    draw_tenx(ax, p)
    rank_axes(ax)
    last = int(p.tenx_ranks.max())
    ax.annotate("cells: a plateau", (40, p.c[39]), (3, 3e2), arrowprops={"arrowstyle": "->", "color": "k"})
    knee = int(np.sum(p.c >= 1000))
    ax.annotate("the cliff", (knee, 1000), (5e3, 3e3), arrowprops={"arrowstyle": "->", "color": "k"})
    ax.annotate("empty droplets: a shoulder", (2e4, p.c[19999]), (3e4, 1.2e2),
                arrowprops={"arrowstyle": "->", "color": "k"})
    ax.annotate("error barcodes: the tail", (1.5e5, 1.2), (5e2, 1.5), arrowprops={"arrowstyle": "->", "color": "k"})
    ax.set_title(f"{p.label}: {p.n:,} barcodes with a count")
    ax.legend(loc="upper right", markerscale=4)
    fig.tight_layout()
    save(fig, "rank-anatomy")
    note("pbmc.tenx_last_rank", last)
    note("pbmc.count_at_20000", int(p.c[19999]))


def fig_rank_survival(p: Sample) -> None:
    """The rank curve turned on its side is a survival function; its slope is the density."""
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 5.6), sharex=True, gridspec_kw={"height_ratios": [1.2, 1]})
    picks = thin(p.n)
    logc = np.log10(p.c.astype(np.float64))
    ax = axes[0]
    ax.plot(logc[picks], p.ranks[picks], color=GREY, lw=1.4)
    ax.scatter(logc[p.tenx_ranks - 1], p.tenx_ranks, s=3, color=BLUE, linewidths=0, label="10x's cells")
    ax.set_yscale("log")
    ax.set_ylabel("barcodes with at least\nthis many UMIs (log)")
    ax.set_title(f"{p.label}: the rank curve on its side is a survival function")
    ax.legend(loc="upper right", markerscale=4)
    ax = axes[1]
    bins = np.linspace(0, logc.max() + 0.05, 90)
    ax.hist(logc, bins=bins, color=GREY)
    ax.hist(logc[p.tenx_ranks - 1], bins=bins, color=BLUE)
    ax.set_yscale("log")
    ax.set_xlabel("UMIs (log10)")
    ax.set_ylabel("barcodes per bin (log)")
    ax.set_title("the histogram of the same totals: each steep fall above is a hump here")
    fig.tight_layout()
    save(fig, "rank-survival")


def fig_rank_views(p: Sample) -> None:
    fig, axes = plt.subplots(1, 4, figsize=(15, 3.5))
    r, c = p.ranks, p.c.astype(np.float64)
    share = p.cum / p.cum[-1]
    tenx = p.tenx_ranks - 1
    shown = 2000
    ax = axes[0]
    ax.plot(r[:shown], c[:shown], color=GREY, lw=1.4)
    ax.scatter(r[tenx], c[tenx], s=3, color=BLUE, linewidths=0)
    ax.set_xlim(0, shown)
    ax.set_title(f"linear: the top {shown:,} ranks")
    ax.set_xlabel("rank")
    ax.set_ylabel("UMIs")
    ax = axes[1]
    picks = thin(p.n)
    ax.plot(r[picks], c[picks], color=GREY, lw=1.4)
    ax.scatter(r[tenx], c[tenx], s=3, color=BLUE, linewidths=0)
    rank_axes(ax, ylabel="")
    ax.set_title("log-log: every barcode")
    ax = axes[2]
    ax.plot(r[picks], share[picks], color=GREY, lw=1.4)
    ax.scatter(r[tenx], share[tenx], s=3, color=BLUE, linewidths=0)
    ax.set_xscale("log")
    ax.set_xlabel("rank (log)")
    ax.set_ylabel("share of all UMIs")
    ax.set_title("cumulative: UMIs held by the top r")
    ax = axes[3]
    logc = np.log10(c)
    bins = np.linspace(0, logc.max() + 0.05, 70)
    ax.hist(logc, bins=bins, color=GREY)
    ax.hist(logc[tenx], bins=bins, color=BLUE)
    ax.set_yscale("log")
    ax.set_xlabel("UMIs (log10)")
    ax.set_ylabel("barcodes")
    ax.set_title("density: a histogram of log totals")
    fig.tight_layout()
    save(fig, "rank-views")
    note("pbmc.share_held_by_tenx", round(float(share[p.tenx_ranks.size - 1]), 3))


def fig_rank_nuclei(p: Sample, q: Sample) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.9), sharey=True)
    for ax, s in zip(axes, (p, q)):
        draw_curve(ax, s)
        draw_tenx(ax, s)
        rank_axes(ax)
        ax.set_title(f"{s.label}: {s.tenx_ranks.size:,} cells called by 10x")
        ax.legend(loc="lower left", markerscale=4)
    axes[1].annotate("a gradual bend", (6e3, q.c[5999]), (2.5e2, 4e2), arrowprops={"arrowstyle": "->", "color": "k"})
    axes[1].annotate("a second cliff,\nwhere empty droplets end", (1.6e5, 60), (5e3, 8), arrowprops={"arrowstyle": "->", "color": "k"})
    fig.tight_layout()
    save(fig, "rank-nuclei")
    for s in (p, q):
        note(f"{s.key}.median_cell", int(np.median(s.c[s.tenx_ranks - 1])))
        note(f"{s.key}.ambient_median_rank_2x_to_10x", int(np.median(s.c[2 * s.tenx_ranks.size: 10 * s.tenx_ranks.size])))


# ---------------------------------------------------------------------------
# Knee geometry
# ---------------------------------------------------------------------------


def views(s: Sample):
    r = s.ranks
    c = s.c.astype(np.float64)
    cum = s.cum.astype(np.float64)
    return [
        ("UMIs against rank", r, c, False, False),
        ("both on log scales", np.log10(r), np.log10(c), True, True),
        ("cumulative UMIs against rank", r, cum, False, False),
        ("cumulative UMIs against log rank", np.log10(r), cum, True, False),
    ]


def fig_knee_views(p: Sample, q: Sample) -> None:
    fig, axes = plt.subplots(1, 4, figsize=(15, 3.6))
    for ax, (title, x, y, logx, logy), (_, xq, yq, _, _) in zip(axes, views(p), views(q)):
        at = farthest(x, y)
        at_q = farthest(xq, yq)
        picks = thin(p.n) if logx else np.unique(np.concatenate([thin(p.n, 6000, log=False), np.arange(5000)]))
        ax.plot(x[picks], y[picks], color=GREY, lw=1.4)
        ax.plot([x[0], x[-1]], [y[0], y[-1]], color="k", ls="--", lw=1)
        ax.plot(x[at], y[at], "o", color=VERMILLION, ms=7)
        ax.set_title(f"{title}\nknee at rank {at + 1:,}")
        ax.set_xlabel("log10 rank" if logx else "rank")
        if logy:
            ax.set_ylabel("log10 UMIs")
        if not logx:  # everything happens at the left edge: zoom in on it
            inset = ax.inset_axes([0.4, 0.25, 0.55, 0.55])
            near = np.arange(3000)
            inset.plot(x[near], y[near], color=GREY, lw=1.4)
            inset.plot(x[at], y[at], "o", color=VERMILLION, ms=5)
            slope = (y[-1] - y[0]) / (x[-1] - x[0])
            inset.plot(x[near], y[0] + slope * (x[near] - x[0]), color="k", ls="--", lw=1)
            inset.set_ylim(0, y.max() * 1.05)
            inset.set_title("ranks 1 – 3,000", fontsize=8)
            inset.tick_params(labelsize=7)
        note(f"pbmc.view.{title}", at + 1)
        note(f"nuclei.view.{title}", at_q + 1)
    fig.suptitle(f"{p.label}: the farthest point from the chord, in four drawings of the same totals "
                 f"(10x called {p.tenx_ranks.size:,})", y=1.04)
    fig.tight_layout()
    save(fig, "knee-views")


def fig_knee_parallel() -> None:
    x = np.linspace(0, 1, 400)
    y = 1 - np.exp(-4 * x)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    for ax, (sx, sy, title) in zip(axes, ((1.0, 1.0, "the farthest point is where the curve runs parallel to the chord"),
                                          (1.8, 0.55, "stretch either axis: the distances change, the point does not"))):
        X, Y = x * sx, y * sy
        at = farthest(X, Y)
        A, B, P = np.array([X[0], Y[0]]), np.array([X[-1], Y[-1]]), np.array([X[at], Y[at]])
        ax.fill([A[0], B[0], P[0]], [A[1], B[1], P[1]], color=SKY, alpha=0.25, lw=0)
        ax.plot(X, Y, color=GREY, lw=2)
        ax.plot([A[0], B[0]], [A[1], B[1]], color="k", ls="--", lw=1)
        u = (B - A) / np.linalg.norm(B - A)
        ax.plot([P[0] - 0.35 * u[0], P[0] + 0.35 * u[0]], [P[1] - 0.35 * u[1], P[1] + 0.35 * u[1]], color=VERMILLION, lw=1.5)
        foot = A + np.dot(P - A, u) * u
        ax.plot([P[0], foot[0]], [P[1], foot[1]], color=VERMILLION, lw=1, ls=":")
        ax.plot(*P, "o", color=VERMILLION, ms=7)
        ax.text(*(A + [0.02, 0.06]), "A")
        ax.text(*(B + [-0.06, -0.09]), "B")
        ax.text(*(P + [-0.07, 0.04]), "P")
        ax.set_aspect("equal")
        ax.set_title(title, fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])
    fig.tight_layout()
    save(fig, "knee-parallel")


def fig_knee_mean(p: Sample) -> None:
    """dropbead's window on PBMC: the chord, the distance, and the mean it amounts to."""
    N = dropseq.MAX_CELLS
    top = p.c[:N].astype(np.float64)
    share = np.cumsum(top) / top.sum()
    x = np.arange(1, N + 1)
    distance = share - x / N
    mean = top.sum() / N
    cells = int(np.sum(top > mean))
    fig, axes = plt.subplots(1, 3, figsize=(15, 3.7))
    ax = axes[0]
    ax.plot(x, share, color=GREY, lw=1.6)
    ax.plot([0, N], [0, 1], color="k", ls="--", lw=1)
    ax.plot(cells, share[cells - 1], "o", color=VERMILLION, ms=7)
    ax.set_xlabel("rank")
    ax.set_ylabel("share of the top 15,000's UMIs")
    ax.set_title("dropbead's curve: the top 15,000 barcodes")
    shown = 3000
    ax = axes[1]
    ax.plot(x[:shown], distance[:shown], color=GREY, lw=1.6)
    ax.axvline(cells, color=VERMILLION, lw=1)
    ax.set_xlabel("rank")
    ax.set_ylabel("height above the diagonal")
    ax.set_title("it climbs while each step adds more than 1/N")
    ax = axes[2]
    ax.plot(x[:shown], top[:shown], color=GREY, lw=1.6)
    ax.axhline(mean, color="k", ls="--", lw=1)
    ax.fill_between(x[:cells], mean, top[:cells], color=BLUE, alpha=0.25, lw=0)
    ax.axvline(cells, color=VERMILLION, lw=1)
    ax.set_yscale("log")
    ax.set_xlabel("rank")
    ax.set_ylabel("UMIs (log)")
    ax.text(cells * 1.05, mean * 1.3, f"mean of the top 15,000 = {mean:,.0f}")
    ax.set_title(f"so it stops where counts fall below the mean: {cells:,}")
    fig.tight_layout()
    save(fig, "knee-mean")


def cobweb(ax, s: Sample, multiplier: float, form: str, title: str, low: float = 10,
           also_from: Optional[int] = None) -> Tuple[List[Tuple[int, int]], List[Tuple[int, str]]]:
    ks = np.unique(np.round(np.logspace(math.log10(low), math.log10(s.n), 1500)).astype(np.int64))
    g = refit_map(s, multiplier, ks, form)
    ax.plot(ks, np.maximum(g, 1), color=GREY, lw=1.6, label=f"k → knee of the first {multiplier:g}k")
    ax.plot(ks, ks, color="k", lw=0.8, ls="--", label="k = k")
    c, cum = s.c, s.cum
    if form == "umitools":
        step, grow = (lambda m: umitools_pass(c, cum, m)), (lambda k: min(int(multiplier * k), s.n))
    else:
        step, grow = (lambda m: alevinfry_pass(c, cum, m)), (lambda k: min(int(multiplier * k), s.n - 1))
    passes = iterate(step, s.n, grow)

    def web(knees, style, label):
        px, py = [knees[0]], [knees[0]]
        for k in knees[1:]:
            px += [px[-1], k]
            py += [k, k]
        ax.plot(px, py, color=VERMILLION, lw=1.2, ls=style, label=label)
        ax.plot(knees[0], knees[0], "o", color=VERMILLION, ms=5, mfc="white")

    web([k for _, k in passes], "-", "the passes, from every barcode")
    if also_from is not None:
        other = iterate(step, grow(also_from), grow)
        web([also_from] + [k for _, k in other], ":", f"had it started at {also_from:,}")
        note(f"{s.key}.cobweb.x{multiplier:g}.from_{also_from}", other[-1][1] + 1)
    points = fixed_points(ks, g)
    for k, kind in points:
        ax.plot(k, k, "o", color="k", ms=7, mfc="k" if kind == "stable" else "white", zorder=4)
    ax.axvline(s.tenx_ranks.size, color=BLUE, lw=1, ls=":", label=f"10x's cells ({s.tenx_ranks.size:,})")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("knee k (log)")
    ax.set_ylabel("next knee (log)")
    ax.set_title(title)
    return passes, points


def fig_knee_cobweb(q: Sample) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), sharey=True)
    p3, f3 = cobweb(axes[0], q, 3, "umitools", "", low=1000)
    p5, f5 = cobweb(axes[1], q, 5, "alevinfry", "", low=1000, also_from=50000)
    axes[0].set_title(f"×3 (UMI-tools): the passes end at {p3[-1][1] + 1:,}")
    axes[1].set_title(f"×5 (alevin-fry): the passes end at {p5[-1][1] + 1:,}")
    axes[0].legend(loc="upper left", fontsize=8)
    axes[1].legend(loc="upper left", fontsize=8)
    fig.suptitle(f"{q.label}: the refit loop as a map (filled: stable fixed points; open: unstable)", y=1.02)
    fig.tight_layout()
    save(fig, "knee-cobweb")
    note("nuclei.cobweb.x3.fixed", f3)
    note("nuclei.cobweb.x5.fixed", f5)
    note("nuclei.cobweb.x3.passes", len(p3))
    note("nuclei.cobweb.x5.end", p5[-1][1] + 1)
    for m in (2, 4, 6, 8):
        ks = np.unique(np.round(np.logspace(1, math.log10(q.n), 800)).astype(np.int64))
        note(f"nuclei.cobweb.x{m}.fixed", fixed_points(ks, refit_map(q, m, ks)))


def fig_knee_powerlaw() -> None:
    a = np.linspace(0.005, 0.995, 400)
    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    ax.plot(a, rho(a), color=GREY, lw=2)
    ax.axhline(1 / math.e, color="k", lw=0.8, ls=":")
    ax.text(0.62, 1 / math.e + 0.01, "1/e: a flat curve", fontsize=9)
    for multiplier, colour in ((3, GREEN), (5, PURPLE)):
        star = critical_slope(multiplier)
        ax.axhline(1 / multiplier, color=colour, lw=1, ls="--")
        ax.plot(star, 1 / multiplier, "o", color=colour)
        ax.text(star + 0.02, 1 / multiplier + 0.01, f"×{multiplier}: a* = {star:.2f}", color=colour, fontsize=9)
        note(f"powerlaw.critical.x{multiplier}", round(star, 3))
    ax.set_xlabel("local slope a of the log-log curve (counts ∝ rank^-a)")
    ax.set_ylabel("share of the window above its mean")
    ax.set_ylim(0, 0.42)
    fig.tight_layout()
    save(fig, "knee-powerlaw")


def critical_slope(multiplier: float) -> float:
    lo, hi = 1e-9, 1 - 1e-9
    for _ in range(200):
        mid = (lo + hi) / 2
        if rho(np.array([mid]))[0] > 1 / multiplier:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def fig_knee_linear_loglog(p: Sample) -> None:
    from sc_droplets.methods import dropletutils
    from sc_droplets.methods import kneedle as kneedle_method
    from sc_droplets.primitives import kneed

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    _, x, y = kneedle_method.rank_curve(p.totals, 1000.0, False, False)
    found = kneed.prepare(x, y, curve="convex", direction="decreasing", behaviour="0.9")
    knee_rank = int(p.call("kneedle").thresholds["knee_rank"])
    top = int(farthest(x.astype(float), y.astype(float))) + 1
    ax = axes[0]
    ax.plot(x, y, color=GREY, lw=1.6)
    ax.plot([x[0], x[-1]], [y[0], y[-1]], color="k", ls="--", lw=1)
    ax.plot(top, y[top - 1], "o", color=VERMILLION, ms=7, label=f"farthest from the chord: rank {top}")
    ax.plot(knee_rank, y[knee_rank - 1], "s", color=PURPLE, ms=7, label=f"Kneedle's knee: rank {knee_rank}")
    ax.set_xlabel("rank")
    ax.set_ylabel("UMIs")
    ax.set_title(f"linear: barcodes with at least 1,000 UMIs ({x.size:,})")
    ax.legend()
    note("pbmc.kneedle.points", int(x.size))
    note("pbmc.kneedle.farthest", top)
    note("pbmc.kneedle.knee", knee_rank)
    ax = axes[1]
    ranks = dropletutils.barcode_ranks(p.totals, {"knee": "chord", "exclude": True, "round_totals": True},
                                       lower=100.0, exclude_from=50)
    chord = ranks.diagnostics["chord"]
    ax.plot(10 ** ranks.x, 10 ** ranks.y, color=GREY, lw=1.6)
    ax.plot([10 ** chord[0], 10 ** chord[2]], [10 ** chord[1], 10 ** chord[3]], color="k", ls="--", lw=1)
    cells = int(p.call("dropletutils-barcoderanks@1.28").n_cells)
    ax.plot(cells, ranks.knee, "o", color=VERMILLION, ms=7, label=f"barcodeRanks 1.28: {cells:,} cells")
    rank_axes(ax)
    ax.set_title("log-log: barcodeRanks' chord between the plateau and the cliff")
    ax.legend(loc="lower left")
    fig.tight_layout()
    save(fig, "knee-linear-loglog")
    note("pbmc.barcoderanks128.cells", cells)


def fig_knee_curvature(p: Sample) -> None:
    from sc_droplets.methods import dropletutils
    from sc_droplets.primitives.smoothspline import smooth_spline

    ranks = dropletutils.barcode_ranks(p.totals, {"knee": "spline", "exclude": True, "fallback": True, "round_totals": True},
                                       lower=100.0, exclude_from=50, df=20.0)
    left, right = ranks.diagnostics["left"], ranks.diagnostics["right"]
    x, y = ranks.x[left:right + 1], ranks.y[left:right + 1]
    fit = smooth_spline(x, y, df=20.0)
    d1, d2 = fit.predict(deriv=1), fit.predict(deriv=2)
    kappa = d2 / (1 + d1 * d1) ** 1.5
    # the same from the points themselves: second differences on an uneven grid
    fd1 = np.gradient(y, x)
    fd2 = np.gradient(fd1, x)
    raw = fd2 / (1 + fd1 * fd1) ** 1.5
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    ax = axes[0]
    ax.plot(x, y, ".", color=GREY, ms=3, label="runs of equal totals")
    ax.plot(fit.x, fit.fitted, color=BLUE, lw=1.5, label="smooth.spline, df = 20")
    knee_x = x[int(np.nanargmin(kappa))]
    ax.axvline(knee_x, color=VERMILLION, lw=1)
    ax.set_xlabel("log10 rank")
    ax.set_ylabel("log10 UMIs")
    ax.set_title("barcodeRanks fits a spline between the plateau and the cliff")
    ax.legend(loc="lower left")
    ax = axes[1]
    ax.plot(x, raw, color=LIGHT, lw=1, label="from the points: second differences")
    ax.plot(x, kappa, color=BLUE, lw=1.5, label="from the spline")
    ax.axvline(knee_x, color=VERMILLION, lw=1)
    lim = 4 * np.nanmax(np.abs(kappa))
    ax.set_ylim(-lim, lim)
    ax.set_xlabel("log10 rank")
    ax.set_ylabel("curvature κ")
    ax.set_title("raw differences drown the bend; the spline's least κ is the knee")
    ax.legend(loc="upper left", fontsize=8)
    fig.tight_layout()
    save(fig, "knee-curvature")
    note("pbmc.curvature.raw_range", [round(float(np.nanmin(raw)), 1), round(float(np.nanmax(raw)), 1)])
    note("pbmc.curvature.spline_range", [round(float(np.nanmin(kappa)), 2), round(float(np.nanmax(kappa)), 2)])
    note("pbmc.curvature.points", int(x.size))


#: The members of the family, in the pages' order: (label, spec).
MEMBERS = [
    ("dropbead", "dropseq-dropbead"),
    ("UMI-tools distance", "umitools-distance"),
    ("alevin-fry knee", "alevinfry-permit-list"),
    ("zUMIs", "zumis-cell-barcodes"),
    ("Kneedle (yascp's search)", "kneedle"),
    ("cellsweep 1.0.2 (kneed)", "cellsweep-empty"),
    ("cellsweep 0.1 (second derivative)", "cellsweep-empty@0.1.0"),
    ("barcodeRanks 1.26 (curvature)", "dropletutils-barcoderanks@1.26"),
    ("barcodeRanks 1.28 (chord)", "dropletutils-barcoderanks@1.28"),
    ("barcodeRanks 1.30 (tracing)", "dropletutils-barcoderanks@1.30"),
    ("barcodeRanks 1.30 inflection", "dropletutils-barcoderanks@1.30:use=inflection"),
    ("dropEst (expected)", "dropest-bounds"),
    ("alevin (tangent)", "alevin-initial-whitelist"),
]


def member_cells(s: Sample, spec: str) -> int:
    result = s.call(spec)
    if spec == "alevin-initial-whitelist":
        return int(result.thresholds["knee_rank"])  # the tangent alone, before the density correction
    return int(result.n_cells)


def fig_knee_landing(p: Sample, q: Sample) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), sharey=True)
    labels = [label for label, _ in MEMBERS]
    for ax, s in zip(axes, (p, q)):
        cells = [member_cells(s, spec) for _, spec in MEMBERS]
        ys = np.arange(len(MEMBERS))[::-1]
        ax.scatter(cells, ys, color=VERMILLION, s=28, zorder=3)
        for x, y in zip(cells, ys):
            ax.text(x * 1.15, y, f"{x:,}", va="center", fontsize=8)
        ax.axvline(s.tenx_ranks.size, color=BLUE, lw=1, ls=":")
        ax.text(s.tenx_ranks.size * 1.08, len(MEMBERS) - 0.4, "10x", color=BLUE, fontsize=8)
        ax.set_xscale("log")
        ax.set_xlim(1, s.n * 3)
        ax.set_xlabel("cells called (log)")
        ax.set_title(s.label)
        ax.grid(axis="y", color="#eeeeee")
        for (label, spec), n in zip(MEMBERS, cells):
            note(f"{s.key}.cells.{spec}", n)
    axes[0].set_yticks(np.arange(len(MEMBERS))[::-1], labels)
    fig.suptitle("Where each rule stops", y=1.0)
    fig.tight_layout()
    save(fig, "knee-landing")


# ---------------------------------------------------------------------------
# UMI-tools
# ---------------------------------------------------------------------------


def fig_umitools_passes(p: Sample, q: Sample) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    ax = axes[0]
    for s, colour in ((p, BLUE), (q, VERMILLION)):
        _, passes = umitools.distance_knee(s.c)
        knees = [k + 1 for _, k in passes]
        ax.plot(np.arange(1, len(knees) + 1), knees, "o-", color=colour, ms=4, label=s.label)
        ax.axhline(s.tenx_ranks.size, color=colour, lw=1, ls=":")
    ax.set_yscale("log")
    ax.set_xlabel("pass")
    ax.set_ylabel("cells (k + 1, log)")
    ax.set_title("the knee after each pass (dotted: 10x)")
    ax.legend()
    cobweb(axes[1], p, 3, "umitools", f"{p.label}: ×3 settles at {umitools.distance_knee(p.c)[0] + 1:,}")
    cobweb(axes[2], q, 3, "umitools", f"{q.label}: ×3 slides to {umitools.distance_knee(q.c)[0] + 1:,}")
    axes[1].legend(loc="upper left", fontsize=7)
    fig.tight_layout()
    save(fig, "umitools-passes")


def fig_umitools_window(p: Sample, q: Sample) -> None:
    """Each pass's threshold: the mean count of its window, against the curve."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
    for ax, s in zip(axes, (p, q)):
        draw_curve(ax, s)
        _, passes = umitools.distance_knee(s.c)
        cum = s.cum
        shown = passes if len(passes) <= 8 else passes[:5] + passes[-2:]
        for i, (m, k) in enumerate(shown):
            mean = (cum[m - 1] - cum[0]) / (m - 1)
            ax.plot([1, m], [mean, mean], color=VERMILLION, lw=1, alpha=0.4 + 0.6 * i / max(1, len(shown) - 1))
            ax.plot(k + 1, s.c[k], "o", color=VERMILLION, ms=4)
        ax.axvline(s.tenx_ranks.size, color=BLUE, lw=1, ls=":")
        rank_axes(ax)
        ax.set_title(f"{s.label}: each pass's mean, drawn over its window")
    fig.tight_layout()
    save(fig, "umitools-window")


def fig_umitools_diagnostics(p: Sample) -> None:
    adata = p.adata
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        scd.call(adata, "umitools-distance", key_added="umitools")
    save(scd.pl.diagnostics(adata, "umitools"), "umitools-diagnostics")


# ---------------------------------------------------------------------------
# dropbead
# ---------------------------------------------------------------------------


def dropbead_cells(c: np.ndarray, window: int) -> int:
    """dropbead's number of cells by the identity: the top `window` counts above their mean."""
    top = c[:window].astype(np.float64)
    return int(np.sum(top > top.sum() / window))


def fig_dropbead_window(p: Sample, q: Sample) -> None:
    fig, ax = plt.subplots(figsize=(7.4, 4.2))
    for s, colour in ((p, BLUE), (q, VERMILLION)):
        for window in (2000, 10000, 15000, 50000):
            assert dropseq.estimate_cell_number(s.c.astype(np.float64), window) == dropbead_cells(s.c, window), (s.key, window)
        windows = np.unique(np.round(np.logspace(math.log10(300), math.log10(min(s.n, 10 ** 6)), 400)).astype(np.int64))
        cells = [dropbead_cells(s.c, int(w)) for w in windows]
        ax.plot(windows, cells, color=colour, lw=1.8, label=s.label)
        ax.axhline(s.tenx_ranks.size, color=colour, lw=1, ls=":")
        note(f"{s.key}.dropbead.by_window", {w: dropbead_cells(s.c, w) for w in (1000, 5000, 10000, 15000, 30000, 100000)})
    ax.axvline(15000, color="k", lw=1, ls="--")
    ax.text(15000 * 1.08, 0.04, "max.cells = 15,000", fontsize=9, transform=ax.get_xaxis_transform())
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("barcodes on the curve, N (log)")
    ax.set_ylabel("cells (log)")
    ax.set_title("dropbead's answer as its window grows (dotted: 10x's call)")
    ax.legend(loc="upper left")
    fig.tight_layout()
    save(fig, "dropbead-window")


# ---------------------------------------------------------------------------
# alevin-fry
# ---------------------------------------------------------------------------


def fig_alevinfry_bifurcation(p: Sample, q: Sample) -> None:
    multipliers = np.round(np.arange(2.0, 8.001, 0.05), 2)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4))
    for ax, s in zip(axes, (p, q)):
        ks = np.unique(np.round(np.logspace(2, math.log10(s.n), 500)).astype(np.int64))
        ends = []
        for multiplier in multipliers:
            g = refit_map(s, multiplier, ks)
            for k, kind in fixed_points(ks, g, smallest=100):
                ax.plot(multiplier, k, "o", ms=2.6, color="k" if kind == "stable" else "#bbbbbb", zorder=2)
            passes = iterate(lambda m: umitools_pass(s.c, s.cum, m), s.n, lambda k: max(2, min(int(multiplier * k), s.n)))
            ends.append(passes[-1][1] + 1)
        ax.plot(multipliers, ends, color=VERMILLION, lw=1.6, label="where the loop ends, from every barcode", zorder=3)
        ax.axhline(s.tenx_ranks.size, color=BLUE, lw=1, ls=":", label="10x's call")
        for multiplier in (3, 5):
            ax.axvline(multiplier, color="#cccccc", lw=0.8, zorder=0)
        ax.plot([], [], "o", ms=3, color="k", label="stable fixed points")
        ax.plot([], [], "o", ms=3, color="#bbbbbb", label="unstable fixed points")
        ax.set_yscale("log")
        ax.set_ylim(80, s.n * 1.5)
        ax.set_xlabel("multiplier c (window = c × knee)")
        ax.set_ylabel("knee (log)")
        ax.set_title(s.label)
        note(f"{s.key}.bifurcation.ends", {float(m): e for m, e in zip(multipliers, ends) if m in (2.5, 2.7, 2.8, 3.0, 4.0, 4.5, 4.6, 4.7, 4.8, 5.0, 6.0)})
        note(f"{s.key}.bifurcation.holds_from", float(multipliers[np.flatnonzero(np.array(ends) > 100)[0]]))
        jumps = np.flatnonzero(np.array(ends) > 50000)
        if jumps.size:
            note(f"{s.key}.bifurcation.second_cliff_from", float(multipliers[jumps[0]]))
    axes[0].legend(loc="upper left", fontsize=8)
    fig.suptitle("Fixed points of the refit loop against the multiplier (grey lines: UMI-tools' 3 and alevin-fry's 5)", y=1.02)
    fig.tight_layout()
    save(fig, "alevinfry-bifurcation")


def fig_alevinfry_diagnostics(q: Sample) -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        scd.call(q.adata, "alevinfry-permit-list", key_added="alevinfry")
    save(scd.pl.diagnostics(q.adata, "alevinfry"), "alevinfry-diagnostics")


# ---------------------------------------------------------------------------
# zUMIs
# ---------------------------------------------------------------------------


def zumis_cells(c: np.ndarray, floor: float) -> int:
    window = c[: int(np.sum(c >= floor))].astype(np.float64)
    return int(np.sum(c > window[1:].mean()))


def fig_zumis_floor(p: Sample, q: Sample) -> None:
    fig, ax = plt.subplots(figsize=(7.4, 4.2))
    for s, colour in ((p, BLUE), (q, VERMILLION)):
        for floor in (30, 100, 1000):
            by_rule = int(s.call(f"zumis-cell-barcodes:n_reads_per_cell={floor}").n_cells)
            assert by_rule == zumis_cells(s.c, floor), (s.key, floor, by_rule)
        floors = np.unique(np.round(np.logspace(0, math.log10(3000), 300)).astype(np.int64))
        ax.plot(floors, [zumis_cells(s.c, f) for f in floors], color=colour, lw=1.8, label=s.label)
        ax.axhline(s.tenx_ranks.size, color=colour, lw=1, ls=":")
        note(f"{s.key}.zumis.by_floor", {f: zumis_cells(s.c, f) for f in (1, 10, 20, 50, 100, 200, 300, 500, 1000)})
    ax.axvline(100, color="k", lw=1, ls="--")
    ax.text(108, 0.9, "nReadsperCell = 100", fontsize=9, transform=ax.get_xaxis_transform())
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("floor: barcodes with fewer are left off the curve (log)")
    ax.set_ylabel("cells (log)")
    ax.set_title("zUMIs' answer as its floor moves (dotted: 10x's call)")
    ax.legend(loc="lower left")
    fig.tight_layout()
    save(fig, "zumis-floor")


# ---------------------------------------------------------------------------
# Kneedle
# ---------------------------------------------------------------------------


def kneedle_curve(s: Sample):
    from sc_droplets.methods import kneedle as kneedle_method
    from sc_droplets.primitives import kneed

    _, x, y = kneedle_method.rank_curve(s.totals, 1000.0, False, False)
    return x, y, kneed.prepare(x, y, curve="convex", direction="decreasing", behaviour="0.9")


def kneedle_descents(found) -> List[Tuple[int, float, int]]:
    """Each maximum of the difference curve, its fall to the next minimum, and where that is."""
    d, minima = found.y_difference, found.minima
    out = []
    for top in found.maxima:
        after = minima[minima > top]
        end = int(after[0]) if after.size else d.size - 1
        out.append((int(top), float(d[top] - d[end]), end))
    return out


def fig_kneedle(p: Sample, q: Sample) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12, 7.4))
    for row, s in zip(axes, (p, q)):
        x, y, found = kneedle_curve(s)
        knee = int(s.call("kneedle").thresholds["knee_rank"])
        sensitivity = float(s.call("kneedle").thresholds["sensitivity"])
        descents = kneedle_descents(found)
        longest = max(descents, key=lambda t: t[1])
        # yascp's search lowers S by 100 until a maximum falls by more than S / (n - 1): the first
        # such maximum is the knee. Here one fall dwarfs the others, so it is also the longest.
        first = next(t for t in descents if t[1] > sensitivity / (x.size - 1))
        assert first[0] + 1 == knee, (s.key, first, knee)
        assert first == longest, (s.key, first, longest)
        slope = float(y[0] - y[-1]) / (x.size - 1)
        gaps = -np.diff(y.astype(np.float64))
        assert np.all(gaps[longest[0]:longest[2]] <= slope), s.key  # no neighbour gap above the chord's slope
        ranks = np.arange(1, x.size + 1)
        ax = row[0]
        ax.plot(ranks, found.y_difference, color=GREY, lw=1.2)
        ax.plot(found.maxima + 1, found.y_difference[found.maxima], ".", color="k", ms=3, label="local maxima")
        ax.axvspan(longest[0] + 1, longest[2] + 1, color=SKY, alpha=0.2, lw=0, label="the longest fall")
        top = int(np.argmax(found.y_difference)) + 1
        ax.plot(top, found.y_difference[top - 1], "o", color=VERMILLION, ms=6, label=f"highest point: rank {top}")
        ax.plot(knee, found.y_difference[knee - 1], "s", color=PURPLE, ms=6, label=f"the search's knee: rank {knee:,}")
        ax.set_xscale("log")
        ax.set_xlabel("rank (log)")
        ax.set_ylabel("difference curve")
        ax.set_title(f"{s.label}: {x.size:,} barcodes with at least 1,000 UMIs; S = {sensitivity:,.0f}")
        ax.legend(loc="lower left", fontsize=8)
        ax = row[1]
        ax.plot(ranks[:-1], np.maximum(gaps, 0.5), ".", color=GREY, ms=2)
        ax.axhline(slope, color="k", ls="--", lw=1, label=f"the chord's slope: {slope:.1f} UMIs per rank")
        ax.axvspan(longest[0] + 1, longest[2] + 1, color=SKY, alpha=0.2, lw=0)
        ax.axvline(knee, color=PURPLE, lw=1)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("rank (log)")
        ax.set_ylabel("gap to the next barcode (UMIs, log)")
        ax.set_title("no gap in the shaded stretch exceeds the chord's slope")
        ax.legend(loc="upper left", fontsize=8)
        note(f"{s.key}.kneedle.points", int(x.size))
        note(f"{s.key}.kneedle.chord_slope", round(slope, 1))
        note(f"{s.key}.kneedle.highest", top)
        note(f"{s.key}.kneedle.longest_fall", [longest[0] + 1, longest[2] + 1, round(longest[1], 3)])
        note(f"{s.key}.kneedle.second_longest", round(sorted(t[1] for t in descents)[-2], 3))
        note(f"{s.key}.kneedle.maxima", len(descents))
        note(f"{s.key}.kneedle.sensitivity", sensitivity)
        note(f"{s.key}.kneedle.by_release", {r: int(s.call(f"kneedle@{r}").n_cells) for r in ("0.5.1", "0.7.0", "0.8.6", "0.9.0")})
    fig.tight_layout()
    save(fig, "kneedle-gaps")


# ---------------------------------------------------------------------------
# cellsweep
# ---------------------------------------------------------------------------


def fig_cellsweep(p: Sample, q: Sample) -> None:
    from sc_droplets.methods import cellsweep

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    for ax, s in zip(axes[:2], (p, q)):
        cut, counts, smooth, d2 = cellsweep.threshold_cutoff(s.totals)
        ranks = np.arange(1, d2.size + 1)
        cliff = slice(s.tenx_ranks.size // 2, 3 * s.tenx_ranks.size)
        at = cliff.start + int(np.argmin(d2[cliff]))
        ax.plot(ranks, d2, color=GREY, lw=0.8)
        ax.plot(int(np.argmin(d2)) + 1, d2.min(), "o", color=VERMILLION, ms=6, label=f"most negative: rank {int(np.argmin(d2)) + 1}")
        if s.key == "pbmc":  # the nuclei have no cliff for it to dip at
            ax.plot(at + 1, d2[at], "s", color=BLUE, ms=6, label=f"the cliff's dip: rank {at + 1:,}")
        ax.set_xscale("log")
        ax.set_xlabel("rank (log)")
        ax.set_ylabel("second derivative of smoothed log10 UMIs")
        ax.set_title(f"{s.label}: cellsweep to 1.0.1")
        ax.legend(loc="lower right", fontsize=8)
        from scipy.ndimage import gaussian_filter1d  # the same rule against log rank, for the page's aside

        log_rank = np.log10(np.arange(1, counts.size + 1))
        d2_log = np.gradient(np.gradient(gaussian_filter1d(np.log10(counts), 2), log_rank), log_rank)
        note(f"{s.key}.cellsweep.d2_against_log_rank", {"at": int(np.argmin(d2_log)) + 1, "count": float(counts[np.argmin(d2_log)])})
        note(f"{s.key}.cellsweep.v0.1", {"cells": int(s.call("cellsweep-empty@0.1.0").n_cells),
                                         "cut": float(s.call("cellsweep-empty@0.1.0").thresholds["umis"])})
        note(f"{s.key}.cellsweep.d2", {"min": float(f"{d2.min():.3g}"), "at": int(np.argmin(d2)) + 1,
                                       "cliff_min": float(f"{d2[at]:.3g}"), "cliff_at": at + 1})
    ax = axes[2]
    for s, colour in ((p, BLUE), (q, VERMILLION)):
        t = s.totals[s.totals > cellsweep.MIN_COUNTS]
        umis, n_barcodes = np.unique(t, return_counts=True)
        ranks = np.cumsum(n_barcodes[::-1])
        ax.plot(ranks, umis[::-1], ".-", color=colour, ms=2, lw=0.8, label=s.label)
        result = s.call("cellsweep-empty")
        cut, found, index = cellsweep.kneed_cutoff(s.totals, "0.9.0")
        assert index == int(found.maxima[-1]), s.key  # online, S = 1: the walk ends at the last maximum
        note(f"{s.key}.cellsweep.kneed_highest_is_last", int(np.argmax(found.y_difference)) == index)
        ax.plot(int(result.n_cells), result.thresholds["umis"], "o", color="k", ms=6)
        ax.text(int(result.n_cells) * 1.3, result.thresholds["umis"] * 1.3, f"{int(result.n_cells):,} cells", fontsize=8)
        note(f"{s.key}.cellsweep.kneed", {"cells": int(result.n_cells), "cut": float(result.thresholds["umis"]),
                                          "distinct_totals": int(umis.size),
                                          "by_kneed": {v: int(s.call(f"cellsweep-empty:kneed_version={v}").n_cells) for v in ("0.8.5", "0.8.6", "0.9.0")}})
    rank_axes(ax)
    ax.set_title("1.0.2: Kneedle on one point per distinct total")
    ax.legend(loc="lower left", fontsize=8)
    fig.tight_layout()
    save(fig, "cellsweep-rules")


# ---------------------------------------------------------------------------
# barcodeRanks
# ---------------------------------------------------------------------------

BR_KNEES = [
    ("1.0 – 1.8: curvature, mis-indexed", "dropletutils-barcoderanks@1.0", "v", PURPLE),
    ("1.10 – 1.26: spline curvature", "dropletutils-barcoderanks@1.26", "o", VERMILLION),
    ("1.28: chord", "dropletutils-barcoderanks@1.28", "s", GREEN),
    ("1.30: traced window", "dropletutils-barcoderanks@1.30", "D", ORANGE),
    ("inflection, to 1.28", "dropletutils-barcoderanks@1.26:use=inflection", "^", SKY),
    ("inflection, 1.30", "dropletutils-barcoderanks@1.30:use=inflection", "P", BLUE),
]


def fig_barcoderanks(p: Sample, q: Sample) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), sharey=True)
    for ax, s in zip(axes, (p, q)):
        draw_curve(ax, s, colour=LIGHT)
        above = int(np.sum(s.c > 100))
        picks = thin(above)
        ax.plot(s.ranks[picks], s.c[picks], color=GREY, lw=1.6)
        found = {}
        for label, spec, marker, colour in BR_KNEES:
            result = s.call(spec)
            cells, cut = int(result.n_cells), float(result.thresholds["umis"])
            ax.plot(cells, cut, marker, color=colour, ms=7, label=f"{label}: {cells:,}", zorder=3)
            found[spec] = (cells, round(cut, 1))
        ax.axvline(s.tenx_ranks.size, color=BLUE, lw=1, ls=":")
        rank_axes(ax)
        ax.set_title(f"{s.label} (dotted: 10x's {s.tenx_ranks.size:,})")
        ax.legend(loc="lower left", fontsize=8)
        note(f"{s.key}.barcoderanks", found)
    fig.tight_layout()
    save(fig, "barcoderanks-knees")


def fig_barcoderanks_trace(p: Sample) -> None:
    from sc_droplets.methods import dropletutils

    settings = {"knee": "trace", "round_totals": True}
    ranks = dropletutils.barcode_ranks(p.totals, settings, lower=100.0, exclude_from=50, window=1.0, gradient_threshold=-1.0)
    x, y = ranks.x, ranks.y
    fig, ax = plt.subplots(figsize=(5.2, 5.6))
    ax.plot(x, y, color=GREY, lw=1.6)
    for label, key, colour in (("the knee's window", "knee_window", ORANGE), ("the inflection's window", "inflection_window", BLUE)):
        left, right = ranks.diagnostics[key]
        ends = np.interp([left, right], x, y)
        ax.plot([left, right], ends, color=colour, lw=2, label=label)
        inside = (x >= left) & (x <= right)
        ax.plot(x[inside], y[inside], color=colour, lw=3, alpha=0.35)
    ax.plot(math.log10(int(p.call("dropletutils-barcoderanks@1.30").n_cells)), math.log10(ranks.knee), "o", color=ORANGE, ms=7)
    ax.set_aspect("equal")
    ax.set_xlabel("log10 rank")
    ax.set_ylabel("log10 UMIs")
    ax.set_title("1.30: windows one unit long, traced along the curve")
    ax.legend(loc="lower left")
    fig.tight_layout()
    save(fig, "barcoderanks-trace")
    note("pbmc.barcoderanks.trace", {"windows": ranks.diagnostics["n_windows"], "knee": round(ranks.knee, 1), "inflection": round(ranks.inflection, 1)})


# ---------------------------------------------------------------------------
# dropEst
# ---------------------------------------------------------------------------


def dropest_slopes(c: np.ndarray, lag: float = 0.05):
    """dropestr's lagged slope of log count against log rank, placed at rank i + lag."""
    from sc_droplets.primitives.rounding import round_half_even

    n = c.size
    lag_n = int(round_half_even(n * lag))
    log_u = np.log(c.astype(np.float64))
    log_rank = np.log(np.arange(1, n + 1, dtype=np.float64))
    slope = (log_u[lag_n:] - log_u[:-lag_n]) / (log_rank[lag_n:] - log_rank[:-lag_n])
    return np.arange(lag_n + 1, n + 1), slope, lag_n


def fig_dropest(p: Sample, q: Sample) -> None:
    from sc_droplets.methods import dropest

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
    ax = axes[0]
    for floor, colour in ((1, VERMILLION), (20, BLUE)):
        sub = p.c[p.c >= floor]
        at, slope, lag = dropest_slopes(sub)
        bounds = dropest.estimate_cells_number(sub)
        note(f"pbmc.dropest.floor_{floor}", {"barcodes": int(sub.size), "lag": lag, "expected": bounds["expected"]})
        ax.plot(at, slope, color=colour, lw=1.2, label=f"{'every barcode' if floor == 1 else f'at least {floor} UMIs'}: lag {lag:,}, expected {bounds['expected']:,}")
        ax.axvline(bounds["expected"], color=colour, lw=1, ls="--")
    ax.axvline(p.tenx_ranks.size, color="k", lw=1, ls=":")
    ax.set_xscale("log")
    ax.set_xlabel("rank (log)")
    ax.set_ylabel("slope of log UMIs over the lag")
    ax.set_title(f"{p.label}: the lagged slope (dotted: 10x)")
    ax.legend(loc="lower left", fontsize=8)
    ax = axes[1]
    for s, colour in ((p, BLUE), (q, VERMILLION)):
        assert dropest.estimate_cells_number(s.c)["expected"] == int(s.call("dropest-bounds").n_cells)
        floors = np.unique(np.round(np.logspace(0, math.log10(1000), 120)).astype(np.int64))
        expected = []
        for floor in floors:
            try:
                expected.append(dropest.estimate_cells_number(s.c[s.c >= floor])["expected"])
            except Exception:  # noqa: BLE001 - too few barcodes for the lag: no estimate
                expected.append(np.nan)
        ax.plot(floors, expected, color=colour, lw=1.6, label=s.label)
        ax.axhline(s.tenx_ranks.size, color=colour, lw=1, ls=":")
        note(f"{s.key}.dropest.by_floor", {int(f): dropest.estimate_cells_number(s.c[s.c >= f])["expected"] for f in (1, 2, 5, 10, 20, 50, 100, 200, 500)})
        note(f"{s.key}.dropest.lag", dropest.estimate_cells_number(s.c)["lag"])
        note(f"{s.key}.dropest.barcodes_with_2_or_more", int(np.sum(s.c >= 2)))
        note(f"{s.key}.dropest.upper_bound_from_100", dropest.estimate_cells_number(s.c[s.c >= 100])["max"])
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("barcodes given: those with at least this many UMIs (log)")
    ax.set_ylabel("expected cells (log)")
    ax.set_title("dropEst's estimate against what it is given (dotted: 10x)")
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    save(fig, "dropest-lag")


# ---------------------------------------------------------------------------
# alevin
# ---------------------------------------------------------------------------


def alevin_curve(s: Sample, top: int = 100000):
    """alevin's log cumulative curve of the top M, summed from the smallest up, and the tangent."""
    m = min(s.n, top)
    values = s.c[:m][::-1].astype(np.float64)
    total = np.cumsum(values)
    curve = np.log(total)
    i = np.arange(m, dtype=np.float64)
    ratio = np.full(m, np.inf)
    ratio[1:] = curve[1:] / i[1:]
    return m, total, curve, int(np.argmin(ratio))


def fig_alevin(p: Sample, q: Sample) -> None:
    from sc_droplets.methods import alevin

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    m, total, curve, j = alevin_curve(p)
    ax = axes[0]
    i = np.arange(m)
    ax.plot(i, curve, color=GREY, lw=1.6)
    ax.plot([0, m], [0, curve[j] / j * m], color=VERMILLION, lw=1)
    ax.plot(j, curve[j], "o", color=VERMILLION, ms=6)
    ax.set_xlabel("barcodes summed, from the smallest of the top 100,000")
    ax.set_ylabel("ln(cumulative UMIs)")
    ax.set_title(f"{p.label}: the tangent from the origin touches at {m - j:,} from the top")
    inset = ax.inset_axes([0.42, 0.12, 0.5, 0.5])
    near = slice(j - 3000, m)
    inset.plot(i[near], curve[near], color=GREY, lw=1.4)
    inset.plot(i[near], curve[j] / j * i[near], color=VERMILLION, lw=1)
    inset.plot(j, curve[j], "o", color=VERMILLION, ms=4)
    inset.tick_params(labelsize=7)
    for ax, s in zip(axes[1:], (p, q)):
        m, total, curve, j = alevin_curve(s)
        knee = alevin.left_boundary(s.c, m)
        r = np.arange(1, m)  # rank r: the barcodes below it among the top m are ranks r + 1 .. m
        below = total[::-1][1:]  # sum of the counts below rank r
        threshold = below / (m - r) * np.log(below)
        picks = thin(m - 1)
        ax.plot(r[picks], s.c[: m - 1][picks], color=GREY, lw=1.6, label="UMIs")
        ax.plot(r[picks], threshold[picks], color=VERMILLION, lw=1.2, ls="--", label="ln(sum below) × mean below")
        ax.axvline(knee, color=VERMILLION, lw=1)
        ax.axvline(s.tenx_ranks.size, color=BLUE, lw=1, ls=":")
        rank_axes(ax)
        ax.set_title(f"{s.label}: the tangent's {knee:,} cells")
        ax.legend(loc="lower left", fontsize=8)
        mean_below = float(total[j] / (j + 1))
        note(f"{s.key}.alevin", {"cells": knee, "count_at_knee": int(s.c[knee - 1]), "mean_below": round(mean_below, 1),
                                 "ln_sum_below": round(float(curve[j]), 2), "product": round(float(mean_below * curve[j]), 1),
                                 "v0.11": int(s.call("alevin-initial-whitelist@0.11.0").thresholds["knee_rank"]),
                                 "valley": int(s.call("alevin-initial-whitelist").thresholds["gauss_rank"]),
                                 "after_correction": int(s.call("alevin-initial-whitelist").n_cells)})
    fig.tight_layout()
    save(fig, "alevin-tangent")


# ---------------------------------------------------------------------------
# Widget data
# ---------------------------------------------------------------------------


def write_data(s: Sample, specs: List[Tuple[str, str]]) -> None:
    values, starts = np.unique(-s.c, return_index=True)
    runs = np.diff(np.concatenate([np.sort(starts), [s.n]]))
    data = {
        "label": s.label,
        "barcodes": s.n,
        "tenx": int(s.tenx_ranks.size),
        "values": (-values).astype(int).tolist(),
        "runs": runs.astype(int).tolist(),
        "calls": {label: member_cells(s, spec) for label, spec in specs},
    }
    assert sum(data["runs"]) == s.n
    path = DATA / f"{s.key}.json"
    path.write_text(json.dumps(data, separators=(",", ":")) + "\n")
    print(f"  curves/{path.name} ({path.stat().st_size / 1000:.0f} kB)")


# ---------------------------------------------------------------------------


PAGES = ("rank-curve", "knee", "dropbead", "umitools", "alevinfry", "zumis", "kneedle", "cellsweep", "barcoderanks",
         "dropest", "alevin", "data")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--page", choices=PAGES, action="append", help="draw only these pages (repeatable)")
    pages = set(parser.parse_args().page or PAGES)
    OUT.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)
    p, q = load("pbmc"), load("nuclei")
    for s in (p, q):
        check_identities(s)
        note_refits(s)
    print("identities hold on both matrices")
    if "rank-curve" in pages:
        fig_rank_toy()
        fig_rank_anatomy(p)
        fig_rank_survival(p)
        fig_rank_views(p)
        fig_rank_nuclei(p, q)
    if "knee" in pages:
        fig_knee_views(p, q)
        fig_knee_parallel()
        fig_knee_mean(p)
        fig_knee_cobweb(q)
        fig_knee_powerlaw()
        fig_knee_linear_loglog(p)
        fig_knee_curvature(p)
        fig_knee_landing(p, q)
    if "dropbead" in pages:
        fig_dropbead_window(p, q)
    if "umitools" in pages:
        fig_umitools_passes(p, q)
        fig_umitools_window(p, q)
        fig_umitools_diagnostics(p)
    if "alevinfry" in pages:
        fig_alevinfry_bifurcation(p, q)
        fig_alevinfry_diagnostics(q)
    if "zumis" in pages:
        fig_zumis_floor(p, q)
    if "kneedle" in pages:
        fig_kneedle(p, q)
    if "cellsweep" in pages:
        fig_cellsweep(p, q)
    if "barcoderanks" in pages:
        fig_barcoderanks(p, q)
        fig_barcoderanks_trace(p)
    if "dropest" in pages:
        fig_dropest(p, q)
    if "alevin" in pages:
        fig_alevin(p, q)
    if "data" in pages:
        specs = [("dropbead", "dropseq-dropbead"), ("UMI-tools", "umitools-distance"), ("alevin-fry", "alevinfry-permit-list"),
                 ("zUMIs", "zumis-cell-barcodes"), ("Kneedle", "kneedle"), ("barcodeRanks 1.28", "dropletutils-barcoderanks@1.28")]
        for s in (p, q):
            write_data(s, specs)
    print()
    for key in sorted(NUMBERS):
        print(f"{key}: {NUMBERS[key]}")


if __name__ == "__main__":
    main()
