"""DropletQC's empty droplets, checked against DropletQC itself (`oracles/dropletqc/run.sh`)."""

from __future__ import annotations

import gzip
import re
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import pytest
from typer.testing import CliRunner

import sc_droplets as scd
from sc_droplets.cli import app
from sc_droplets.methods.dropletqc import DropletQCError, identify_empty_drops, ks_behaviour

matplotlib.use("Agg")
runner = CliRunner()
ORACLE = Path(__file__).parent / "oracles" / "dropletqc"
KS_RELEASES = {"1.14.2": True, "1.15.3": False}


def _out(result) -> str:
    return re.sub(r"\s+", " ", re.sub(r"\x1b\[[0-9;]*m", "", result.output))


# ---------------------------------------------------------------------------
# Against DropletQC itself
# ---------------------------------------------------------------------------

CASES = [(v, case) for v in KS_RELEASES
         for case in pd.read_csv(ORACLE / f"ks-{v}" / "summary.tsv", sep="\t", index_col=0).index]


@pytest.mark.parametrize("version, case", CASES, ids=[f"ks{v}-{c}" for v, c in CASES])
def test_identify_empty_drops_matches_dropletqc(version, case):
    """oracles/dropletqc/run.sh: DropletQC's own identify_empty_drops, with ks 1.14.2 and 1.15.3,
    on (nuclear fraction, UMI) tables of PBMC 1k and simulated cases."""
    row = pd.read_csv(ORACLE / f"ks-{version}" / "summary.tsv", sep="\t", index_col=0).loc[case]
    table = pd.read_csv(ORACLE / "inputs" / f"{case}.tsv.gz", sep="\t", header=None, names=["barcode", "nf", "umi"])
    status, cut, derivative, runs = identify_empty_drops(table["nf"].to_numpy(), table["umi"].to_numpy(), 0.05, 1000,
                                                         KS_RELEASES[version])
    theirs = gzip.open(ORACLE / f"ks-{version}" / f"{case}.status.gz", "rt").read().split()
    assert ["cell" if s == "rescued" else s for s in status] == theirs
    assert runs == row["runs"]
    assert derivative.h == pytest.approx(row["h"], rel=1e-13)
    assert cut == pytest.approx(row["cut"], rel=1e-13)


def test_the_two_ks_releases_cut_differently():
    """ks 1.14.3 moved the counts in its convolution: the same bandwidth, another cut."""
    cuts = {v: pd.read_csv(ORACLE / f"ks-{v}" / "summary.tsv", sep="\t", index_col=0) for v in KS_RELEASES}
    assert cuts["1.14.2"].loc["sim_bimodal_1", "h"] == cuts["1.15.3"].loc["sim_bimodal_1", "h"]
    assert cuts["1.14.2"].loc["sim_bimodal_1", "cut"] != cuts["1.15.3"].loc["sim_bimodal_1", "cut"]


def test_ks_releases_are_read_by_their_change():
    assert ks_behaviour("1.13.2") and ks_behaviour("1.14.2")
    assert not ks_behaviour("1.14.3") and not ks_behaviour("1.15.3")
    with pytest.raises(DropletQCError, match="older than any release checked"):
        ks_behaviour("1.12.0")


# ---------------------------------------------------------------------------
# The method
# ---------------------------------------------------------------------------


@pytest.fixture
def bimodal():
    """Empty droplets with a low nuclear fraction and cells with a higher one, all above 500 UMIs."""
    rng = np.random.default_rng(3)
    adata = scd.simulate_droplets(n_cells=400, n_empty=1500, n_low=500, n_genes=50, seed=3)
    n = adata.n_obs
    totals = np.asarray(adata.X.sum(axis=1)).ravel()
    order = np.argsort(-totals)
    nf = np.full(n, 0.05)
    nf[order[:400]] = rng.beta(12, 10, 400)
    nf[order[400:1600]] = rng.beta(2, 20, 1200)
    adata.obs["nuclear_fraction"] = nf
    return adata


def test_the_barcodes_given_are_a_floor_or_a_call(bimodal):
    result = scd.call(bimodal, "dropletqc", min_counts=0, inplace=False)
    assert result.thresholds["given"] == bimodal.n_obs
    bimodal.obs["chosen"] = np.arange(bimodal.n_obs) < 1000
    result = scd.call(bimodal, "dropletqc", among="chosen", inplace=False)
    assert result.thresholds["given"] == 1000
    assert set(np.asarray(result.classes)[1000:]) == {"not given"}


def test_cells_are_above_the_cut_or_rescued(bimodal):
    result = scd.call(bimodal, "dropletqc", min_counts=0, umi_rescue=1e9, inplace=False)
    nf = bimodal.obs["nuclear_fraction"].to_numpy()
    cut = result.thresholds["nuclear_fraction"]
    assert np.array_equal(result.is_cell, nf >= cut)
    rescued = scd.call(bimodal, "dropletqc", min_counts=0, umi_rescue=0, nf_rescue=0.2, inplace=False)
    classes = np.asarray(rescued.classes)
    assert np.array_equal(classes == "rescued", (nf < cut) & (nf > 0.2))


def test_the_unspliced_fraction_stands_in_with_a_note(spliced):
    adata = spliced.copy()
    scd.add_metrics(adata)
    result = scd.call(adata, "dropletqc", inplace=False)
    assert any("approximation" in note for note in result.notes)
    with pytest.raises(DropletQCError, match="nuclear fraction"):
        scd.call(spliced.copy(), "dropletqc")


def test_a_fraction_table_is_matched_by_barcode(bimodal, tmp_path):
    table = tmp_path / "nf.csv"
    pd.DataFrame({"barcode": bimodal.obs_names, "nf": bimodal.obs["nuclear_fraction"]}).to_csv(table, index=False)
    adata = bimodal.copy()
    del adata.obs["nuclear_fraction"]
    from_table = scd.call(adata, "dropletqc", fraction=str(table), min_counts=0, inplace=False)
    from_obs = scd.call(bimodal, "dropletqc", min_counts=0, inplace=False)
    assert np.array_equal(from_table.is_cell, from_obs.is_cell)


def test_a_barcode_without_a_fraction_counts_as_zero(bimodal):
    adata = bimodal.copy()
    adata.obs.loc[adata.obs_names[:3], "nuclear_fraction"] = np.nan
    result = scd.call(adata, "dropletqc", min_counts=0, inplace=False)
    assert any("set to 0" in note for note in result.notes)
    assert np.all(result.columns["nuclear_fraction"][:3] == 0)


def test_diagnostics_and_the_report_show_the_derivative(bimodal):
    scd.call(bimodal, "dropletqc", min_counts=0)
    figure = scd.pl.diagnostics(bimodal, "dropletqc")
    derivative, plane = figure.axes
    assert derivative.get_ylabel() == "density derivative" and plane.get_xlabel() == "nuclear fraction"
    assert 'id="call-dropletqc"' in scd.report(bimodal)


def test_cli_takes_the_fraction_from_a_table(bimodal, tmp_path):
    adata = bimodal.copy()
    table = tmp_path / "nf.tsv"
    pd.DataFrame({"barcode": adata.obs_names, "nf": adata.obs["nuclear_fraction"]}).to_csv(table, sep="\t", index=False)
    del adata.obs["nuclear_fraction"]
    path = tmp_path / "raw.h5ad"
    adata.write_h5ad(path)
    result = runner.invoke(app, ["call", "dropletqc", str(path), "--fraction", str(table), "--min-counts", "0", "--ks-version", "1.14.2", "--json"])
    assert result.exit_code == 0, _out(result)
    assert '"variant": "v1.0"' in _out(result)
