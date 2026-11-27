"""Which rule did a Cell Ranger 4.0-7.1 run follow? Compare 10x's own cells with combinations of
the parts of the 3.x and 7.2 rules (see cellranger-4.0-7.1-outputs.md).

    uv run python docs/_research/cellranger-4.0-7.1-survey.py DATASET [--expect N] [--lt]
        [--aggregates] [--per-genome] [--cap 262144]

DATASET is a name in tests/datasets.py. Single GEM group. Every combination of

    ordmag   v2 (3.x: one stream seeded 0 for the stage) | v7 (7.2: a private RandomState(0),
             the (b+1)-th largest, the tie extension; cells estimated when --expect is not given)
    stream   shared (3.x: the simulation continues the stage's stream after OrdMag) |
             own (a fresh RandomState(0) for each test)
    refill   3.x (the block is refilled when its last feature is used) |
             7.2 (gaps of 20-999 refill it early)
    cand     cap (3.x: the 20,000 with the most UMIs) | all (7.2)

is printed with its disagreement with 10x's filtered barcodes and the adjusted p-values of the
barcodes that differ. --lt takes ranks 4,500-9,000 as ambient (90,000 otherwise), --per-genome
tests each genome on its own features after the one before (6.1+), --aggregates removes antibody
aggregates by their counts first (5.0+).
"""

from __future__ import annotations

import argparse
import itertools
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tests"))
from datasets import filtered_barcodes, raw_path  # noqa: E402

import sc_droplets as scd  # noqa: E402
from sc_droplets.methods import cellranger as cr  # noqa: E402
from sc_droplets.primitives.legacysort import legacy_argsort  # noqa: E402
from sc_droplets.primitives.multinomial import multinomial_loglik_rows, scipy_adjusted_profile  # noqa: E402
from sc_droplets.primitives.pvalues import bh_cellranger, monte_carlo_pvalues  # noqa: E402
from sc_droplets.primitives.simulate import cellranger_v1  # noqa: E402


def emptydrops(counts, called, state, cap, batch_small, high):
    """One EmptyDrops test (Cell Ranger 3.x's, with the options above): candidates and adjusted p."""
    umis = cr._row_sums(counts)
    size = len(umis)
    order = legacy_argsort(umis, window=(size - high, size - high // 2))
    empty = np.sort(order[::-1][high // 2 : high])
    use = empty[umis[empty] > 0]
    features = np.flatnonzero(np.asarray(counts.sum(axis=0)).ravel())
    profile = cr._profile(counts, use, features, "7.2")
    min_umis = int(max(cr.MIN_UMIS, np.ceil(float(np.median(umis[called])) * 0.01)))
    if cap:
        masked = called | (umis < min_umis)
        n_unmasked = int((~masked).sum())
        filled = np.where(masked, np.iinfo(np.int64).max, umis)
        window = (n_unmasked - cr.N_CANDIDATE_BARCODES, n_unmasked)
        candidates = legacy_argsort(filled, window=window)[0:n_unmasked][-cr.N_CANDIDATE_BARCODES :]
    else:
        candidates = np.flatnonzero(~called & (umis >= min_umis))
    rows = counts[candidates][:, features]
    with np.errstate(divide="ignore"):
        observed = multinomial_loglik_rows(rows, np.log(scipy_adjusted_profile(profile)))
    n_lower = cellranger_v1(profile, umis[candidates], observed, cr.NUM_SIMS, state, batch_small=batch_small)
    return candidates, bh_cellranger(monte_carlo_pvalues(n_lower, cr.NUM_SIMS))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("dataset")
    parser.add_argument("--expect", type=int, help="--expect-cells of the run (3.x default 3,000; 7.2 estimates)")
    parser.add_argument("--lt", action="store_true", help="ambient ranks 4,500-9,000")
    parser.add_argument("--aggregates", action="store_true", help="remove antibody aggregates first")
    parser.add_argument("--per-genome", action="store_true", help="EmptyDrops per genome, one after another")
    parser.add_argument("--cap", type=int, default=45000, help="7.2's cap on the estimated cells (7.0: 262144)")
    args = parser.parse_args()

    adata = scd.read(raw_path(args.dataset))
    barcodes = np.asarray(adata.obs_names).astype(str)
    types = adata.var["feature_types"].astype(str).to_numpy()
    counts = adata.X.tocsr()
    keep = np.ones(adata.n_obs, dtype=bool)
    if args.aggregates:
        aggregates, _ = cr.antibody_aggregates(counts[:, np.flatnonzero(types == cr.ANTIBODY)], None)
        keep &= ~aggregates
    gex = counts[np.flatnonzero(keep)][:, np.flatnonzero(types == cr.GEX)]
    barcodes = barcodes[keep]
    genome_list = cr.genomes(adata.var.iloc[np.flatnonzero(types == cr.GEX)])
    expected = set(filtered_barcodes(args.dataset))
    high = 9000 if args.lt else cr.N_PARTITIONS

    stage = np.random.RandomState(0)
    initial = {"v2": np.zeros(gex.shape[0], bool), "v7": np.zeros(gex.shape[0], bool)}
    for _, mask in genome_list:
        totals = cr._row_sums(gex[:, np.flatnonzero(mask)])
        initial["v2"][cr.ordmag_v2(totals, args.expect or cr.DEFAULT_RECOVERED_CELLS, stage)[0]] = True
        initial["v7"][cr.ordmag_v7(totals, args.expect, args.cap, "revert")[0]] = True
    after_ordmag = stage.get_state()
    print(f"{args.dataset}: 10x {len(expected)} cells; OrdMag v2 {int(initial['v2'].sum())}, v7 {int(initial['v7'].sum())}")

    for ordmag, stream, refill, cap in itertools.product(["v2", "v7"], ["shared", "own"], ["3.x", "7.2"], [True, False]):
        if ordmag == "v7" and stream == "shared":
            continue  # 7.2's OrdMag draws from a private stream
        called = initial[ordmag].copy()
        padj = np.full(gex.shape[0], np.nan)
        parts = [(None, None)] if not args.per_genome else genome_list
        for _, mask in parts:
            sub = gex if mask is None else gex[:, np.flatnonzero(mask)]
            if stream == "shared":
                state = np.random.RandomState()
                state.set_state(after_ordmag)
            else:
                state = np.random.RandomState(0)
            candidates, p = emptydrops(sub, called.copy(), state, cap, 20 if refill == "7.2" else None, high)
            padj[candidates] = p
            called[candidates[p <= 0.01]] = True
        ours = set(barcodes[called])
        different = np.isin(barcodes, list((expected - ours) | (ours - expected)))
        p = np.sort(padj[different][~np.isnan(padj[different])])
        print(f"  ordmag={ordmag} stream={stream} refill={refill} cand={'cap' if cap else 'all'}: "
              f"missing {len(expected - ours)}, extra {len(ours - expected)}"
              + (f", adjusted p {p[0]:.4f}-{p[-1]:.4f}" if p.size else ""), flush=True)


if __name__ == "__main__":
    main()
