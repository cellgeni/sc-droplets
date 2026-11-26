"""bustools allowlist on small and simulated matrices."""

from __future__ import annotations

import math
import re

import anndata as ad
import numpy as np
import pytest
import scipy.sparse as sp
from typer.testing import CliRunner

import sc_droplets as scd
from sc_droplets.cli import app
from sc_droplets.methods.bustools import BLOCK, BustoolsError, sample

runner = CliRunner()


def _out(result) -> str:
    return re.sub(r"\s+", " ", re.sub(r"\x1b\[[0-9;]*m", "", result.output))


def _barcode(i: int, length: int = 16) -> str:
    return "".join("ACGT"[(i >> (2 * k)) & 3] for k in reversed(range(length)))


def _matrix(totals, barcodes=None) -> ad.AnnData:
    adata = ad.AnnData(sp.csr_matrix(np.asarray(totals, dtype=np.float64).reshape(-1, 1)))
    adata.obs_names = barcodes if barcodes is not None else [_barcode(i + 1) + "-1" for i in range(len(totals))]
    adata.var_names = ["G0"]
    adata.var["feature_types"] = "Gene Expression"
    return adata


def test_the_threshold_is_the_error_fraction_of_the_top_ten_mean():
    totals = [1000 + i for i in range(10)] + [10] * 300 + [200, 151, 150]
    result = scd.call(_matrix(totals), "bustools-allowlist", inplace=False)
    mean = sum(1000 + i for i in range(10)) // 10
    cut = int(mean * (1 - math.pow(0.99, 16)))
    assert result.thresholds["top10_mean"] == mean and result.thresholds["umis"] == cut
    assert np.array_equal(result.is_cell, np.asarray(totals) >= cut)


def test_the_sample_is_every_barcode_in_the_first_blocks_of_records():
    # 300 barcodes of 1000 records: a block holds 100, so the sample needs two blocks
    found = sample(np.full(300, 1000), np.full(300, 1000), first_is_zero=False)
    assert found.blocks == 2
    # barcode 99 ends with the first block: it enters at the block's end, and again when the
    # next block opens barcode 100
    assert len(found.entries) == 201
    assert [e for e in found.entries if e[0] == 99] == [(99, 1000), (99, 1000)]


def test_a_barcode_open_at_a_blocks_end_enters_with_its_records_so_far():
    records = np.array([BLOCK - 10, 30] + [1] * 300)
    found = sample(records, records.copy(), first_is_zero=False)
    assert found.entries[1] == (1, 10)  # 10 of its 30 records were in the first block
    assert found.entries[2] == (1, 30)  # and it enters again, whole, in the second


def test_reads_spread_evenly_over_a_barcodes_records():
    found = sample(np.array([BLOCK + 4]), np.array([3 * (BLOCK + 4) + 2]), first_is_zero=False)
    q, r = divmod(3 * (BLOCK + 4) + 2, BLOCK + 4)
    assert found.entries[0] == (0, q * BLOCK + min(BLOCK, r))


def test_a_first_barcode_of_all_as_counts_one_short():
    """bustools 0.42.0+ starts with the current barcode 0, the code of all A's, so that barcode's
    first record adds to the empty count (-1) instead of opening it."""
    barcodes = ["A" * 16] + [_barcode(i + 1) for i in range(20)]
    totals = [500] + [1000] * 10 + [10] * 10
    found = scd.call(_matrix(totals, barcodes), "bustools-allowlist", inplace=False)
    cut = found.thresholds["umis"]
    assert cut == int(1000 * (1 - math.pow(0.99, 16)))
    forced = scd.call(_matrix([cut + 1] + totals[1:], barcodes), "bustools-allowlist", inplace=False)
    edge = scd.call(_matrix([cut] + totals[1:], barcodes), "bustools-allowlist", inplace=False)
    assert forced.is_cell[0] and not edge.is_cell[0]  # cut - 1 < cut
    given = scd.call(_matrix([cut] + totals[1:], barcodes), "bustools-allowlist", threshold=cut, inplace=False)
    assert given.is_cell[0]  # with --threshold the sample is not read, and the count is whole


def test_a_given_threshold_skips_the_sample(sim):
    result = scd.call(sim, "bustools-allowlist", threshold=100, inplace=False)
    assert "sample_size" not in result.diagnostics
    assert np.array_equal(result.is_cell, np.asarray(result.totals) >= 100)


@pytest.mark.parametrize("barcodes, message", [
    (["ACGT", "ACGTA"], "barcodes of one length"),
    (["ACGN", "ACGT"], "only A, C, G and T"),
    (["ACGT-1", "ACGT-2"], "repeat once their suffix is dropped"),
])
def test_barcodes_bustools_cannot_code_are_refused(barcodes, message):
    with pytest.raises(BustoolsError, match=message):
        scd.call(_matrix([5] * len(barcodes), barcodes), "bustools-allowlist", inplace=False)


def test_fewer_than_ten_barcodes_have_no_top_ten():
    with pytest.raises(BustoolsError, match="reads past its end"):
        scd.call(_matrix([5] * 9), "bustools-allowlist", inplace=False)


def test_cli_calls_bustools_allowlist(tmp_path):
    adata = _matrix([1000 + i for i in range(10)] + [10] * 50)
    path = tmp_path / "small.h5ad"
    adata.write_h5ad(path)
    result = runner.invoke(app, ["call", "bustools-allowlist", str(path), "--json"])
    assert result.exit_code == 0, _out(result)
    assert '"n_cells": 10' in _out(result)
