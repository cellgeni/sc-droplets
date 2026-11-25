"""Simulated raw droplet matrices with known cells, for tests and examples.

The shape follows a real raw matrix: a few thousand cells of a few types,
each with some ambient RNA; a shoulder of empty droplets holding only
ambient RNA; and a long tail of barcodes with a handful of counts (sequencing
errors, free-floating barcodes). `obs["truth"]` says which is which.
"""

from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd
import scipy.sparse as sp


def _barcodes(rng: np.random.Generator, n: int, length: int = 16) -> np.ndarray:
    letters = np.array(list("ACGT"))
    seen = set()
    out = []
    while len(out) < n:
        batch = rng.integers(0, 4, size=(2 * (n - len(out)) + 10, length))
        for row in letters[batch]:
            code = "".join(row) + "-1"
            if code not in seen:
                seen.add(code)
                out.append(code)
                if len(out) == n:
                    break
    return np.array(out, dtype=object)


def simulate_droplets(
    n_cells: int = 1500,
    n_empty: int = 8000,
    n_low: int = 30000,
    n_genes: int = 1000,
    n_types: int = 4,
    cell_umis: float = 3000.0,
    cell_sigma: float = 0.6,
    empty_umis: float = 40.0,
    empty_sigma: float = 0.5,
    low_umis: float = 2.0,
    low_sigma: float = 0.8,
    contamination: float = 0.05,
    seed: int = 0,
    splicing: bool = False,
    ensure_ambient_range: Optional[int] = None,
):
    """An AnnData of raw counts with `obs["truth"]` in {cell, empty, low}.

    Sizes are log-normal: median `cell_umis` (σ `cell_sigma`) for cells,
    `empty_umis` for empty droplets, `low_umis` for the tail. A cell's counts
    are a share `contamination` of ambient RNA and the rest from its type's
    profile; empty and low barcodes hold ambient RNA only. The ambient profile
    is the average of the type profiles, as RNA leaking from lysed cells.

    With `splicing`, layers `spliced` and `unspliced` split each count: cells
    have about 25 % unspliced, ambient RNA (cytoplasmic) about 10 %.

    `ensure_ambient_range` adds low barcodes until there are at least that
    many in total, for methods that take their ambient pool from fixed ranks
    (Cell Ranger and STARsolo use ranks 45,000 to 90,000).
    """
    import anndata as ad

    rng = np.random.default_rng(seed)
    if ensure_ambient_range is not None:
        n_low = max(n_low, ensure_ambient_range - n_cells - n_empty)
    profiles = rng.dirichlet(np.full(n_genes, 0.3), size=n_types)
    ambient = profiles.mean(axis=0)

    def sizes(n: int, median: float, sigma: float, minimum: int) -> np.ndarray:
        return np.maximum(minimum, np.round(rng.lognormal(np.log(median), sigma, size=n))).astype(np.int64)

    cell_size = sizes(n_cells, cell_umis, cell_sigma, 200)
    empty_size = sizes(n_empty, empty_umis, empty_sigma, 1)
    low_size = sizes(n_low, low_umis, low_sigma, 1)
    cell_type = rng.integers(0, n_types, size=n_cells)

    n_total = n_cells + n_empty + n_low
    rows, cols = [], []

    def draw(indices: np.ndarray, counts: np.ndarray, profile: np.ndarray) -> None:
        total = int(counts.sum())
        if total == 0:
            return
        rows.append(np.repeat(indices, counts))
        cols.append(rng.choice(n_genes, size=total, p=profile))

    cell_ambient = rng.binomial(cell_size, contamination)
    cell_own = cell_size - cell_ambient
    cell_index = np.arange(n_cells)
    for t in range(n_types):
        members = cell_type == t
        draw(cell_index[members], cell_own[members], profiles[t])
    draw(cell_index, cell_ambient, ambient)
    draw(np.arange(n_cells, n_cells + n_empty), empty_size, ambient)
    draw(np.arange(n_cells + n_empty, n_total), low_size, ambient)

    row = np.concatenate(rows)
    col = np.concatenate(cols)
    ones = np.ones(row.shape[0], dtype=np.int32)
    matrix = sp.coo_matrix((ones, (row, col)), shape=(n_total, n_genes)).tocsr()
    matrix.sum_duplicates()

    truth = np.array(["cell"] * n_cells + ["empty"] * n_empty + ["low"] * n_low, dtype=object)
    kind = np.concatenate([cell_type.astype(str), np.full(n_empty + n_low, "", dtype=object)])
    order = rng.permutation(n_total)
    matrix = matrix[order]
    obs = pd.DataFrame(
        {"truth": pd.Categorical(truth[order], categories=["cell", "empty", "low"]), "cell_type": kind[order]},
        index=pd.Index(_barcodes(rng, n_total)),
    )
    var = pd.DataFrame(
        {"gene_ids": [f"SIMG{i:05d}" for i in range(n_genes)], "feature_types": "Gene Expression", "genome": "sim"},
        index=pd.Index([f"Gene{i}" for i in range(n_genes)]),
    )
    adata = ad.AnnData(X=matrix, obs=obs, var=var)

    if splicing:
        is_cell = (obs["truth"] == "cell").to_numpy()
        share = np.where(is_cell, rng.beta(5, 15, size=n_total), rng.beta(2, 18, size=n_total))
        coo = matrix.tocoo()
        unspliced_counts = rng.binomial(coo.data, share[coo.row])
        unspliced = sp.csr_matrix((unspliced_counts, (coo.row, coo.col)), shape=matrix.shape)
        unspliced.eliminate_zeros()
        spliced = (matrix - unspliced).tocsr()
        spliced.eliminate_zeros()
        adata.layers["spliced"] = spliced
        adata.layers["unspliced"] = unspliced
    return adata
