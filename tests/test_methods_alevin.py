"""alevin's initial whitelist on small and simulated matrices."""

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
from sc_droplets.methods.alevin import AlevinError, expected_cells, forced_cells, left_boundary, whitelist

runner = CliRunner()


def _out(result) -> str:
    return re.sub(r"\s+", " ", re.sub(r"\x1b\[[0-9;]*m", "", result.output))


def _matrix(totals) -> ad.AnnData:
    adata = ad.AnnData(sp.csr_matrix(np.asarray(totals, dtype=np.float64).reshape(-1, 1)))
    adata.obs_names = [f"B{i}" for i in range(len(totals))]
    adata.var_names = ["G0"]
    adata.var["feature_types"] = "Gene Expression"
    return adata


def _sorted(sim) -> np.ndarray:
    totals = np.asarray(sim.X.sum(axis=1)).ravel().astype(np.int64)
    return np.sort(totals[totals > 0])[::-1]


def _get_left_boundary(f, top, shifted=False):
    """alevin's getLeftBoundary, transcribed loop for loop."""
    cum, freqs = 0.0, []
    for i in range(top):
        cum += float(f[top - i] if shifted else f[top - i - 1])
        freqs.append(math.log(cum))
    for j in range(top):
        y = freqs[j]
        if y == freqs[0]:
            continue
        slope = y / j
        if all(not (i * slope > freqs[i]) for i in range(j + 1, top)):
            cutoff = top - j
            frequency = f[cutoff]
            while f[cutoff] == frequency:
                cutoff -= 1
            return cutoff
    return 0


@pytest.mark.parametrize("shifted", [False, True])
def test_the_knee_is_alevins_tangent_from_the_origin(sim, shifted):
    f = _sorted(sim)
    top = min(3000, f.size - 1)
    assert left_boundary(f, top, shifted=shifted) == _get_left_boundary(f, top, shifted=shifted)


def test_0_11_refuses_a_top_that_reaches_the_last_barcode(sim):
    with pytest.raises(AlevinError, match="past the end"):
        left_boundary(_sorted(sim), _sorted(sim).size, shifted=True)


def test_the_density_correction_moves_the_knee_only_within_a_tenth_of_it(sim):
    result = scd.call(sim, "alevin-initial-whitelist", inplace=False)
    knee = result.thresholds["knee_rank"]
    if "gauss_rank" in result.thresholds and result.thresholds["cells"] != knee:
        assert 0.1 * knee <= result.thresholds["cells"] <= knee
    else:
        assert result.thresholds["cells"] == knee


def test_expected_cells_take_one_more_from_1_2():
    f = np.array([5000, 4000, 3000, 2900, 2000, 300, 100, 50, 20, 10, 5, 3, 2, 1, 1, 1])
    # r = max(1, int(5 * 0.01)) = 1; cut = 0.1 * 4000 = 400; the first count below it is at rank 5
    assert expected_cells(f, 5, plus_one=False) == 5
    assert expected_cells(f, 5, plus_one=True) == 6
    with pytest.raises(AlevinError, match="right Boundary"):
        expected_cells(f, 2, plus_one=True)  # nothing below the cut before rank 2 * 2


@pytest.mark.parametrize("variant, cells", [("v0.12", 6), ("v0.14", 6), ("v0.14.1", 5), ("v1.2", 6)])
def test_forced_cells_by_release(variant, cells):
    f = np.array([900, 800, 700, 600, 500, 400, 300, 200, 100, 5, 4, 3])
    assert forced_cells(f, 6, variant, 10) == cells


def test_forced_cells_below_the_read_floor_are_dropped_from_0_14():
    f = np.array([900, 800, 700, 9, 8, 7, 6, 5, 4, 3])
    assert forced_cells(f, 6, "v0.12", 10) == 6
    # f[6] = 6 < 10: lowered to the last barcode with 10, rank 2 (and one more from 1.2)
    assert forced_cells(f, 6, "v0.14", 10) == 2
    assert forced_cells(f, 6, "v1.2", 10) == 3


def test_the_whitelist_stops_above_the_counts_tied_with_its_end():
    f = np.array([100] * 10 + [50] * 5 + [7] * 300 + [3] * 100)
    # k = 10: the band is 200 wide (the floor), ending at rank 210, a 7; every 7 goes
    for variant in ("v0.12", "v1.2"):
        size, _ = whitelist(f, 10, variant, 200)
        assert size == 15


def test_keep_cb_fraction_forces_a_fraction_of_the_barcodes_seen(sim):
    n = _sorted(sim).size
    result = scd.call(sim, "alevin-initial-whitelist", keep_cb_fraction=0.02, inplace=False)
    assert result.thresholds["cells"] == int(0.02 * n)
    with pytest.raises(AlevinError, match="before 0.13.0 has no --keepCBFraction"):
        scd.call(sim, "alevin-initial-whitelist@0.12.0", keep_cb_fraction=0.02, inplace=False)
    with pytest.raises(AlevinError, match="before 0.12.0 has only the knee"):
        scd.call(sim, "alevin-initial-whitelist@0.11.3", force_cells=100, inplace=False)


def test_classes_are_the_cells_then_the_band(sim):
    result = scd.call(sim, "alevin-initial-whitelist", force_cells=150, inplace=False)
    totals, classes = np.asarray(result.totals), np.asarray(result.classes.astype(str))
    assert (classes == "cell").sum() == 150
    assert totals[classes == "cell"].min() >= totals[classes == "low_confidence"].max()
    assert (classes != "background").sum() == result.thresholds["whitelist"]


def test_diagnostics_draw_alevins_density(adata):
    import matplotlib

    matplotlib.use("Agg")
    scd.call(adata, "alevin-initial-whitelist", key_added="al")
    fig = scd.pl.diagnostics(adata, "al")
    assert fig.axes[1].get_title() == "Density of log10 counts (alevin's correction)"


def test_cli_calls_alevin_initial_whitelist(sim, tmp_path):
    path = tmp_path / "sim.h5ad"
    sim.write_h5ad(path)
    result = runner.invoke(app, ["call", "alevin-initial-whitelist", str(path), "--tool-version", "1.1.0",
                                 "--force-cells", "120", "--json"])
    assert result.exit_code == 0, _out(result)
    assert '"n_cells": 119' in _out(result)  # 0.14.1 - 1.1.0 keep F - 1
