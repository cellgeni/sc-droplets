"""What a method sees: the barcodes, their counts, and their totals.

`Droplets.from_anndata` selects the features a method counts (Gene Expression
by default, as Cell Ranger and most tools do), the layer to read, and an
optional per-barcode read count. Every method receives one of these, never the
AnnData itself, so a method cannot depend on anything the caller did not pass.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Optional, Sequence, Union

import numpy as np
import pandas as pd
import scipy.sparse as sp

#: Columns of `var` that hold a feature's type, in the order they are tried:
#: scanpy's name, then Cell Ranger's own.
FEATURE_TYPE_COLUMNS = ("feature_types", "feature_type")
GENOME_COLUMNS = ("genome",)

ALL_FEATURES = "all"


def feature_type_column(var: pd.DataFrame) -> Optional[str]:
    for column in FEATURE_TYPE_COLUMNS:
        if column in var.columns:
            return column
    return None


@dataclass
class Droplets:
    """Barcodes × selected features, with per-barcode totals.

    `counts` is CSR with one row per barcode of the input, in input order.
    `totals` are its row sums as float64. `reads` is a per-barcode read count
    when one was given. `obs` keeps per-barcode columns a method may use
    (genes detected, unspliced fraction).
    """

    barcodes: np.ndarray
    counts: sp.csr_matrix
    features: pd.DataFrame
    totals: np.ndarray
    reads: Optional[np.ndarray] = None
    obs: pd.DataFrame = field(default_factory=pd.DataFrame)
    feature_types: Optional[tuple] = None
    layer: Optional[str] = None
    genome: Optional[str] = None

    @property
    def n_barcodes(self) -> int:
        return int(self.counts.shape[0])

    @property
    def n_features(self) -> int:
        return int(self.counts.shape[1])

    def values(self, unit: str = "umis") -> np.ndarray:
        """Per-barcode totals in `unit`: UMIs, reads or genes detected."""
        if unit == "umis":
            return self.totals
        if unit == "reads":
            if self.reads is None:
                raise ValueError(
                    "this needs read counts per barcode: pass reads= (an obs column) "
                    "or --reads on the command line"
                )
            return self.reads
        if unit == "genes":
            return np.diff(self.counts.indptr).astype(np.float64)
        raise ValueError(f"unknown count unit {unit!r}")

    def order(self, unit: str = "umis") -> np.ndarray:
        """Barcode indices by decreasing total; ties keep input order (stable)."""
        values = self.values(unit)
        return np.argsort(-values, kind="stable")

    @classmethod
    def from_anndata(
        cls,
        adata: Any,
        *,
        layer: Optional[str] = None,
        feature_types: Union[None, str, Sequence[str]] = None,
        genome: Optional[str] = None,
        reads: Union[None, str, np.ndarray] = None,
    ) -> "Droplets":
        """Select what to count from an AnnData of raw (unfiltered) counts.

        `feature_types`: a type or list of types from `var["feature_types"]`;
        "all" or None counts every feature. Where `var` has no feature-type
        column, every feature is counted whatever was asked: a plain matrix
        holds one kind of feature.
        """
        matrix = adata.layers[layer] if layer is not None else adata.X
        if matrix is None:
            raise ValueError("the AnnData has no X; pass layer=")
        matrix = sp.csr_matrix(matrix) if not sp.issparse(matrix) else matrix.tocsr()
        if matrix.shape != (adata.n_obs, adata.n_vars):
            raise ValueError(f"the matrix is {matrix.shape}, the AnnData {adata.shape}")

        var = adata.var
        keep = np.ones(adata.n_vars, dtype=bool)
        chosen_types: Optional[tuple] = None
        if isinstance(feature_types, str):
            feature_types = None if feature_types == ALL_FEATURES else (feature_types,)
        column = feature_type_column(var)
        if feature_types is not None and column is not None:
            present = list(pd.unique(var[column].astype(str)))
            missing = [t for t in feature_types if t not in present]
            if len(missing) == len(feature_types):
                raise ValueError(
                    f"no features of type {', '.join(map(repr, feature_types))}; "
                    f"the input has {', '.join(map(repr, present))}"
                )
            keep &= var[column].astype(str).isin(list(feature_types)).to_numpy()
            chosen_types = tuple(t for t in feature_types if t in present)
        if genome is not None:
            gcol = next((c for c in GENOME_COLUMNS if c in var.columns), None)
            if gcol is None:
                raise ValueError("genome= was given but var has no 'genome' column")
            genomes = list(pd.unique(var[gcol].astype(str)))
            if genome not in genomes:
                raise ValueError(f"no genome {genome!r}; the input has {', '.join(genomes)}")
            keep &= (var[gcol].astype(str) == genome).to_numpy()

        counts = matrix[:, keep] if not keep.all() else matrix
        counts = counts.tocsr(copy=True)
        counts.sum_duplicates()
        # Genes detected are the stored entries per row, so none may be zero.
        counts.eliminate_zeros()
        totals = np.asarray(counts.sum(axis=1), dtype=np.float64).ravel()

        read_values: Optional[np.ndarray] = None
        if isinstance(reads, str):
            if reads not in adata.obs.columns:
                raise ValueError(f"no obs column {reads!r} for read counts")
            read_values = adata.obs[reads].to_numpy(dtype=np.float64)
        elif reads is not None:
            read_values = np.asarray(reads, dtype=np.float64)
            if read_values.shape != (adata.n_obs,):
                raise ValueError(f"reads has shape {read_values.shape}, expected ({adata.n_obs},)")

        return cls(
            barcodes=np.asarray(adata.obs_names, dtype=object),
            counts=counts,
            features=var.loc[keep].copy(),
            totals=totals,
            reads=read_values,
            obs=adata.obs,
            feature_types=chosen_types,
            layer=layer,
            genome=genome,
        )


def looks_filtered(totals: np.ndarray) -> Optional[str]:
    """A reason to think a matrix holds only called cells, or None.

    Knee methods need the empty droplets. A raw matrix has many barcodes with
    a handful of counts; a filtered one starts in the hundreds.
    """
    nonzero = totals[totals > 0]
    if nonzero.size == 0:
        return "every barcode has zero counts"
    low = int((nonzero < 10).sum())
    if nonzero.size < 20000 and low < 0.01 * nonzero.size:
        return (
            f"only {nonzero.size} barcodes, and {low} with fewer than 10 counts: "
            "this looks like a filtered matrix, and knee methods need the empty droplets too"
        )
    return None
