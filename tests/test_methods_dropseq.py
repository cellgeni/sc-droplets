"""Drop-seq (dropbead) on small and simulated matrices."""

from __future__ import annotations

import math
import re

import numpy as np
import pytest
from typer.testing import CliRunner

import sc_droplets as scd
from sc_droplets.cli import app
from sc_droplets.methods.dropseq import DropbeadError, estimate_cell_number

runner = CliRunner()


def _out(result) -> str:
    return re.sub(r"\s+", " ", re.sub(r"\x1b\[[0-9;]*m", "", result.output))


def _sorted(sim) -> np.ndarray:
    totals = np.asarray(sim.X.sum(axis=1)).ravel()
    return np.sort(totals)[::-1]


def test_dropbead_takes_the_first_barcode_farthest_from_the_diagonal(sim):
    counts = _sorted(sim)
    n = 2000
    cs = np.cumsum(counts[:n] / counts[:n].sum())
    x = np.arange(1, n + 1) / n
    m = (cs[-1] - cs[0]) / (n - x[0])
    expected = int(np.argmax(math.sin(math.atan(m)) * np.abs(cs - x))) + 1
    assert estimate_cell_number(counts, n) == expected
    result = scd.call(sim, "dropseq-dropbead", max_cells=n, inplace=False)
    assert result.n_cells == expected == result.thresholds["cells"]


def test_dropbead_needs_max_cells_barcodes(sim):
    with pytest.raises(DropbeadError, match="needs at least max.cells"):
        scd.call(sim, "dropseq-dropbead", inplace=False)  # 15,000 > the barcodes of `sim`


def test_cli_calls_dropseq_dropbead(sim, tmp_path):
    path = tmp_path / "sim.h5ad"
    sim.write_h5ad(path)
    result = runner.invoke(app, ["call", "dropseq-dropbead", str(path), "--json", "--max-cells", "2000"])
    assert result.exit_code == 0, _out(result)
    assert '"n_cells":' in _out(result)
