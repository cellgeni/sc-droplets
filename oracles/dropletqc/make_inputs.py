"""Write the (nuclear fraction, UMIs) tables DropletQC's oracle runs `identify_empty_drops` on.

    uv run python oracles/dropletqc/make_inputs.py OUT_DIR

One gzipped table per case, `barcode<TAB>nf<TAB>umi`, in the matrix's order:

- pbmc_1k_v3_filtered, pbmc_1k_v3_500: 10x's pbmc_1k_v3 (Cell Ranger 3.0), the barcodes of 10x's
  filtered matrix (as DropletQC's paper applies it, after a caller) and every barcode with at least
  500 UMIs; the fraction is kb's nascent / (mature + ambiguous + nascent), from the committed
  `tests/oracles/splicing` counts, 0 where a barcode has none; UMIs are the Gene Expression totals.
- sim_*: simulated fractions and UMIs (numpy, seeded): two peaks, one peak, a small sample, and
  fractions of small counts with many ties.
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "tests"))

import sc_droplets as scd  # noqa: E402
from datasets import filtered_barcodes, raw_path  # noqa: E402
from sc_droplets.counts import Droplets  # noqa: E402

KB = REPO / "tests" / "oracles" / "splicing" / "pbmc_1k_v3_3_0" / "kb-0.28.2" / "counts.tsv.gz"


def write(path: Path, barcodes, nf, umi) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as handle:
        for b, f, u in zip(barcodes, nf, umi):
            handle.write(f"{b}\t{float(f)!r}\t{int(u)}\n")
    print(f"{path.name}: {len(nf)} barcodes")


def pbmc(out: Path) -> None:
    adata = scd.read(raw_path("pbmc_1k_v3_3_0"))
    umi = Droplets.from_anndata(adata).totals
    kb = pd.read_csv(KB, sep="\t", header=None, names=["barcode", "mature", "nascent", "ambiguous"], index_col=0)
    kb.index = kb.index + "-1"
    kb = kb.reindex(adata.obs_names).fillna(0)
    whole = (kb["mature"] + kb["ambiguous"] + kb["nascent"]).to_numpy(dtype=float)
    nf = np.divide(kb["nascent"].to_numpy(dtype=float), whole, out=np.zeros(whole.size), where=whole > 0)
    names = np.asarray(adata.obs_names)
    filtered = np.isin(names, list(filtered_barcodes("pbmc_1k_v3_3_0")))
    for case, keep in (("pbmc_1k_v3_filtered", filtered), ("pbmc_1k_v3_500", umi >= 500)):
        write(out / f"{case}.tsv.gz", names[keep], nf[keep], umi[keep])


def simulated(out: Path) -> None:
    rng = np.random.default_rng(20261128)
    names = lambda n: [f"S{i:06d}" for i in range(n)]  # noqa: E731
    # Two peaks: empty droplets low, cells higher.
    for seed, (n_empty, n_cell) in enumerate(((3000, 1500), (800, 2500))):
        r = np.random.default_rng(seed)
        nf = np.concatenate([r.beta(2, 20, n_empty), r.beta(12, 10, n_cell)])
        umi = np.concatenate([r.lognormal(6.3, 0.5, n_empty), r.lognormal(8.0, 0.6, n_cell)]).round()
        write(out / f"sim_bimodal_{seed + 1}.tsv.gz", names(nf.size), nf, umi)
    nf = rng.beta(12, 10, 2000)
    write(out / "sim_unimodal.tsv.gz", names(nf.size), nf, rng.lognormal(8, 0.6, nf.size).round())
    nf = np.concatenate([rng.beta(2, 20, 15), rng.beta(12, 10, 10)])
    write(out / "sim_small.tsv.gz", names(nf.size), nf, rng.lognormal(7, 1, nf.size).round())
    reads = rng.poisson(25, 3000) + 1
    nf = np.concatenate([rng.binomial(reads[:1800], 0.1), rng.binomial(reads[1800:], 0.5)]) / reads
    write(out / "sim_ties.tsv.gz", names(nf.size), nf, rng.lognormal(7, 1, nf.size).round())


if __name__ == "__main__":
    target = Path(sys.argv[1])
    pbmc(target)
    simulated(target)
