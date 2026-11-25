"""Writing calls: an annotated .h5ad, a barcode list, a per-barcode table.

Each write goes to a temporary file beside the target and is renamed into
place only once it is complete, so an interrupted run never leaves a
truncated output where a good one was.
"""

from __future__ import annotations

import gzip
import os
from contextlib import contextmanager
from pathlib import Path
from typing import Iterable, Iterator, Optional, Sequence

import pandas as pd


@contextmanager
def _staged(path: Path) -> Iterator[Path]:
    path = Path(path)
    if path.parent and not path.parent.exists():
        raise FileNotFoundError(f"the folder {path.parent} does not exist")
    temporary = path.with_name(f".{path.name}.sc-droplets-tmp{path.suffix}")
    try:
        yield temporary
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def check_outputs(outputs: Iterable[Optional[Path]], inputs: Sequence[Path]) -> None:
    """Refuse an output that is an input, or two outputs that are one file."""
    seen = {}
    resolved_inputs = {Path(p).resolve() for p in inputs}
    for output in outputs:
        if output is None:
            continue
        resolved = Path(output).resolve()
        if resolved in resolved_inputs:
            raise ValueError(f"{output} is also the input; write somewhere else")
        if resolved in seen:
            raise ValueError(f"{output} is given as two outputs")
        seen[resolved] = output


def write_h5ad(adata, path: Path) -> None:
    with _staged(Path(path)) as temporary:
        adata.write_h5ad(temporary)


def write_barcodes(barcodes: Iterable[str], path: Path) -> None:
    """One barcode per line, gzipped when the name ends in .gz (as Cell Ranger's)."""
    text = "".join(f"{b}\n" for b in barcodes).encode("utf-8")
    path = Path(path)
    with _staged(path) as temporary:
        if path.suffix == ".gz":
            with gzip.open(temporary, "wb") as handle:
                handle.write(text)
        else:
            temporary.write_bytes(text)


def write_table(frame: pd.DataFrame, path: Path) -> None:
    """CSV, or TSV when the name ends in .tsv / .tsv.gz; gzipped for .gz."""
    path = Path(path)
    name = path.name[:-3] if path.name.endswith(".gz") else path.name
    sep = "\t" if name.endswith(".tsv") or name.endswith(".txt") else ","
    with _staged(path) as temporary:
        frame.to_csv(temporary, sep=sep, compression="gzip" if path.suffix == ".gz" else None)
