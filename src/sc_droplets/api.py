"""The Python entry points: `call`, `compare`, `agreement`, `list_methods`."""

from __future__ import annotations

from typing import Any, Dict, Iterable, List, Optional, Sequence, Union

import numpy as np
import pandas as pd

from sc_droplets import registry
from sc_droplets.counts import Droplets
from sc_droplets.result import UNS_KEY, CallResult


def list_methods() -> pd.DataFrame:
    """Every registered method: id, tool, title, category, variants, docs page."""
    rows = []
    for cls in registry.all_methods():
        info = cls.info
        rows.append(
            {
                "id": info.id,
                "tool": info.tool,
                "title": info.title,
                "category": info.category,
                "counts": info.counts,
                "deterministic": info.deterministic,
                "variants": ", ".join(f"{v.name} ({v.versions})" if info.versioned else v.name for v in info.variants),
                "status": info.status,
                "docs": info.docs_url,
            }
        )
    return pd.DataFrame(rows).set_index("id")


def call(
    adata: Any,
    method: str,
    *,
    tool_version: Optional[str] = None,
    key_added: Optional[str] = None,
    layer: Optional[str] = None,
    feature_types: Union[None, str, Sequence[str]] = None,
    genome: Optional[str] = None,
    reads: Union[None, str, np.ndarray] = None,
    random_state: Optional[int] = None,
    inplace: bool = True,
    **params: Any,
) -> CallResult:
    """Run one method on an AnnData of raw counts (barcodes × features).

    `method` is an id, optionally with the version of the original tool to
    reproduce: `"dropletutils-emptydrops@1.26"` (or pass `tool_version=`).
    Without a version, the latest known behaviour is used.

    `feature_types` defaults to what the method counts (Gene Expression for
    most); "all" counts every feature. `reads` is an obs column of read counts
    per barcode, for methods that count reads. `random_state` defaults to the
    seed the original tool uses.

    With `inplace` (the default) the call is stored on `adata` under
    `key_added` (default: the method id, plus the version when one was given;
    see `sc_droplets.result`). The `CallResult` is returned either way.
    """
    spec = registry.parse_spec(method)
    method_id = spec.id
    if spec.version and tool_version and registry.ToolVersion(spec.version) != registry.ToolVersion(tool_version):
        raise ValueError(f"{method!r} names version {spec.version}, but tool_version={tool_version!r}")
    version = tool_version or spec.version
    cls = registry.get(method_id)
    info = cls.info
    variant, warnings = info.resolve(version)
    given = {**spec.params, **{k: v for k, v in params.items() if v is not None}}
    built = cls.build_params(variant, given)
    if feature_types is None:
        feature_types = info.feature_types
    droplets = Droplets.from_anndata(adata, layer=layer, feature_types=feature_types, genome=genome, reads=reads)
    seed = info.default_seed if random_state is None else random_state
    result = cls().run(droplets, built, variant, seed)
    result.tool_version = version
    result.notes = warnings + list(result.notes)
    result.counting = {
        "unit": getattr(built, "unit", info.counts),
        "layer": layer,
        "feature_types": list(droplets.feature_types) if droplets.feature_types else "all",
        "genome": genome,
        "n_features": droplets.n_features,
        "reads": reads if isinstance(reads, str) else (None if reads is None else "array"),
        "random_state": seed,
    }
    if inplace:
        result.store(adata, key_added or registry.default_key(method_id, version, spec.params))
    return result


def compare(
    adata: Any,
    methods: Iterable[str],
    **common: Any,
) -> pd.DataFrame:
    """Run several methods and summarise them, one row each.

    Each entry of `methods` is a spec, `id[@tool_version][:name=value,...]`,
    e.g. `"top-n:n_cells=3000"` or `"dropletutils-emptydrops@1.26"`. `common`
    holds options for every call (`layer`, `feature_types`, `reads`...).
    Every call is stored on `adata`.
    """
    rows = []
    keys = []
    for spec in methods:
        result = call(adata, spec, **common)
        keys.append(result.key)
        rows.append(
            {
                "key": result.key,
                "method": result.method,
                "tool_version": result.tool_version or "",
                "variant": result.variant,
                "n_cells": result.n_cells,
                "min_cell_total": result.min_cell_total,
            }
        )
    return pd.DataFrame(rows).set_index("key")


def agreement(adata: Any, keys: Optional[Sequence[str]] = None) -> pd.DataFrame:
    """Pairwise Jaccard index of the cell sets of stored calls."""
    stored = list(adata.uns.get(UNS_KEY, {}))
    keys = list(keys) if keys is not None else stored
    missing = [k for k in keys if k not in adata.obs.columns]
    if missing:
        raise KeyError(f"no stored call {', '.join(missing)}; stored: {', '.join(stored) or 'none'}")
    calls = {k: adata.obs[k].to_numpy(dtype=bool) for k in keys}
    matrix: Dict[str, Dict[str, float]] = {}
    for a in keys:
        matrix[a] = {}
        for b in keys:
            union = np.logical_or(calls[a], calls[b]).sum()
            matrix[a][b] = float(np.logical_and(calls[a], calls[b]).sum() / union) if union else 1.0
    return pd.DataFrame(matrix).loc[keys, keys]


def stored(adata: Any) -> List[str]:
    """Keys of the calls stored on an AnnData."""
    return list(adata.uns.get(UNS_KEY, {}))
