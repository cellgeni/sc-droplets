"""Per-barcode measures beside the counts, and how well each call separates them.

The unspliced fraction U / (S + U) is the share of a barcode's molecules from unspliced
(intron-containing) RNA: DropletQC's nuclear fraction, DropSift's `pct_intronic`. Ambient RNA is
mostly cytoplasmic and spliced, so in whole-cell data empty droplets sit low, cells higher, and
damaged cells or stripped nuclei high; in nuclei data the nuclei sit high. The mitochondrial
fraction is the share of counts in genes named `MT-` (`mt-` in mouse).

`separation` scores a call by how the measure differs between its cells and the barcodes it
left out: unsupervised indicators, not ground truth.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, List, Optional, Sequence, Union

import numpy as np
import pandas as pd
import scipy.sparse as sp

from sc_droplets.counts import Droplets, feature_type_column
from sc_droplets.result import UNS_KEY

UNSPLICED = "unspliced_fraction"
MITO = "mito_fraction"
MITO_PREFIXES = ("MT-", "mt-")
#: Barcodes with fewer counts are too few molecules to estimate a fraction from.
MIN_COUNTS = 10


class MetricError(ValueError):
    """Where a measure cannot be had."""


def barcode_key(name: str) -> str:
    """A barcode as the sequence alone: velocyto's `sample:` prefix and `x` suffix, and a GEM-group
    suffix such as Cell Ranger's `-1`, dropped."""
    core = str(name).rsplit(":", 1)[-1]
    if re.fullmatch(r"[ACGTN]+x", core):
        core = core[:-1]
    return re.sub(r"-\d+$", "", core)


def _gene_expression(adata: Any) -> np.ndarray:
    """Which features are gene expression (every feature when there is no type column)."""
    column = feature_type_column(adata.var)
    if column is None:
        return np.ones(adata.n_vars, dtype=bool)
    return (adata.var[column].astype(str) == "Gene Expression").to_numpy()


def _row_sums(matrix, keep: np.ndarray) -> np.ndarray:
    matrix = sp.csr_matrix(matrix) if not sp.issparse(matrix) else matrix.tocsr()
    if not keep.all():
        matrix = matrix[:, keep]
    return np.asarray(matrix.sum(axis=1), dtype=np.float64).ravel()


def _fraction(part: np.ndarray, whole: np.ndarray) -> np.ndarray:
    out = np.full(whole.shape, np.nan)
    np.divide(part, whole, out=out, where=whole > 0)
    return out


def has_splicing(adata: Any) -> bool:
    return {"spliced", "unspliced"} <= {k for k in adata.layers.keys() if k is not None}


def _from_layers(adata: Any) -> pd.Series:
    keep = _gene_expression(adata)
    spliced = _row_sums(adata.layers["spliced"], keep)
    unspliced = _row_sums(adata.layers["unspliced"], keep)
    return pd.Series(_fraction(unspliced, spliced + unspliced), index=adata.obs_names.copy(), name=UNSPLICED)


def _solo_base(path: Path) -> Optional[Path]:
    for base in (path, path / "Solo.out"):
        if any((base / name / "raw").is_dir() for name in ("Velocyto", "Gene", "GeneFull")):
            return base
    return None


def _gene_full_estimate(base: Path) -> pd.Series:
    """STARsolo without Velocyto: (GeneFull − Gene) / GeneFull per barcode. GeneFull counts reads on
    exons and introns, Gene those on exons, so the difference estimates the intronic share."""
    from sc_droplets.io.mtx import read_mtx_folder

    gene = read_mtx_folder(base / "Gene" / "raw")
    full = read_mtx_folder(base / "GeneFull" / "raw")
    exonic = pd.Series(_row_sums(gene.X, _gene_expression(gene)), index=gene.obs_names)
    everything = pd.Series(_row_sums(full.X, _gene_expression(full)), index=full.obs_names)
    exonic = exonic.reindex(everything.index, fill_value=0.0).to_numpy()
    whole = everything.to_numpy()
    fraction = _fraction(np.clip(whole - exonic, 0, None), whole)
    return pd.Series(fraction, index=everything.index.copy(), name=UNSPLICED)


def unspliced_fraction(source: Union[str, Path, Any]) -> pd.Series:
    """Each barcode's unspliced fraction U / (S + U), over gene-expression features; NaN where it
    has neither.

    `source` is an AnnData with `spliced` and `unspliced` layers, or a path to anything that holds
    them: a velocyto `.loom`, STARsolo's `Velocyto/raw/`, kb `--workflow nac`, alevin-fry USA mode,
    an `.h5ad`. A STARsolo `Solo.out/` (or its parent) without Velocyto but with `Gene` and
    `GeneFull` gives the estimate (GeneFull − Gene) / GeneFull instead.
    """
    if not isinstance(source, (str, Path)):
        if not has_splicing(source):
            raise MetricError("the AnnData has no 'spliced' and 'unspliced' layers")
        return _from_layers(source)
    from sc_droplets.io import read

    path = Path(source)
    base = _solo_base(path) if path.is_dir() else None
    if base is not None:
        if (base / "Velocyto" / "raw").is_dir():
            return _from_layers(read(base, solo_feature="Velocyto"))
        if (base / "Gene" / "raw").is_dir() and (base / "GeneFull" / "raw").is_dir():
            return _gene_full_estimate(base)
        raise MetricError(f"{path}: STARsolo output needs Velocyto/raw, or Gene/raw and GeneFull/raw, for splicing")
    adata = read(path)
    if not has_splicing(adata):
        raise MetricError(
            f"{path} has no spliced and unspliced counts: give a velocyto .loom, STARsolo's Velocyto/raw or "
            "Solo.out with Gene and GeneFull, kb --workflow nac, or alevin-fry in USA mode"
        )
    return _from_layers(adata)


def mito_fraction(adata: Any) -> Optional[np.ndarray]:
    """Each barcode's share of gene-expression counts in genes named `MT-` or `mt-`; None when no
    gene is so named."""
    expression = _gene_expression(adata)
    names = pd.Index(adata.var_names.astype(str))
    mito = expression & np.asarray(names.str.startswith(MITO_PREFIXES), dtype=bool)
    if not mito.any():
        return None
    matrix = adata.X
    return _fraction(_row_sums(matrix, mito), _row_sums(matrix, expression))


def _matched(target: pd.Index, values: pd.Series) -> np.ndarray:
    """`values` (indexed by barcode) on the barcodes of `target`, by name, or by `barcode_key` when
    the names do not match; NaN where a barcode is missing."""
    if values.index.isin(target).any():
        return values.reindex(target).to_numpy(dtype=np.float64)
    keyed = pd.Series(values.to_numpy(), index=[barcode_key(b) for b in values.index])
    keyed = keyed[~keyed.index.duplicated()]
    keys = [barcode_key(b) for b in target]
    if not keyed.index.isin(keys).any():
        raise MetricError("no barcode of the splicing source matches a barcode of the input")
    return keyed.reindex(keys).to_numpy(dtype=np.float64)


def add_metrics(adata: Any, splicing: Union[None, str, Path, Any, pd.Series] = None) -> List[str]:
    """Store each barcode's measures in `obs`: `unspliced_fraction` (from `splicing`, matched by
    barcode, or from the AnnData's own `spliced` and `unspliced` layers) and `mito_fraction` (when
    genes are named `MT-` or `mt-`). Returns the columns written.

    `splicing` is anything `unspliced_fraction` takes, or a Series of fractions by barcode.
    """
    written = []
    if splicing is not None:
        values = splicing if isinstance(splicing, pd.Series) else unspliced_fraction(splicing)
        adata.obs[UNSPLICED] = _matched(adata.obs_names, values)
        written.append(UNSPLICED)
    elif has_splicing(adata):
        adata.obs[UNSPLICED] = _from_layers(adata).to_numpy()
        written.append(UNSPLICED)
    mito = mito_fraction(adata)
    if mito is not None:
        adata.obs[MITO] = mito
        written.append(MITO)
    return written


def _totals(adata: Any, key: str) -> np.ndarray:
    """The totals a stored call counted: its layer and feature types."""
    record = adata.uns.get(UNS_KEY, {}).get(key, {})
    counting = record.get("counting", {}) if hasattr(record, "get") else {}
    types = counting.get("feature_types")
    if types is None:
        column = feature_type_column(adata.var)
        types = "Gene Expression" if column is None or "Gene Expression" in set(adata.var[column].astype(str)) else "all"
    elif not isinstance(types, str):
        types = [str(t) for t in np.atleast_1d(types)]
    return Droplets.from_anndata(adata, layer=counting.get("layer"), feature_types=types).totals


def auroc(higher: np.ndarray, lower: np.ndarray) -> float:
    """P(a value of `higher` exceeds one of `lower`), ties counting half: the Mann-Whitney U over
    the number of pairs."""
    from scipy.stats import mannwhitneyu

    if higher.size == 0 or lower.size == 0:
        return float("nan")
    return float(mannwhitneyu(higher, lower).statistic / (higher.size * lower.size))


def separation(
    adata: Any,
    keys: Optional[Sequence[str]] = None,
    metric: str = UNSPLICED,
    min_counts: float = MIN_COUNTS,
) -> pd.DataFrame:
    """How each stored call separates a measure: its median in the call's cells and in the barcodes
    it left out that have at least `min_counts`, and the AUROC between the two (0.5: no
    separation; above: the cells' values are higher). One row per call."""
    if metric not in adata.obs.columns:
        raise MetricError(f"no obs column {metric!r}: run scd.add_metrics first (with splicing= for the unspliced fraction)")
    values = adata.obs[metric].to_numpy(dtype=np.float64)
    keys = list(keys) if keys is not None else list(adata.uns.get(UNS_KEY, {}))
    rows = []
    for key in keys:
        cells = adata.obs[key].to_numpy(dtype=bool)
        known = ~np.isnan(values)
        others = ~cells & (_totals(adata, key) >= min_counts) & known
        inside = values[cells & known]
        outside = values[others]
        rows.append({
            "key": key,
            f"{metric}_cells": float(np.median(inside)) if inside.size else float("nan"),
            f"{metric}_others": float(np.median(outside)) if outside.size else float("nan"),
            f"{metric}_auroc": auroc(inside, outside),
        })
    return pd.DataFrame(rows).set_index("key")
