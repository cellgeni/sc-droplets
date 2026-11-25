"""The command line: every command end to end on simulated files."""

from __future__ import annotations

import gzip
import json
import re

import anndata as ad
import pytest
from typer.testing import CliRunner

from sc_droplets import __version__
from sc_droplets.cli import app

from conftest import write_10x_h5_v3

runner = CliRunner()


def _out(result) -> str:
    """Output with ANSI codes stripped and whitespace collapsed."""
    text = re.sub(r"\x1b\[[0-9;]*m", "", result.output)
    return re.sub(r"\s+", " ", text)


@pytest.fixture
def h5(sim, tmp_path):
    return write_10x_h5_v3(sim, tmp_path / "raw_feature_bc_matrix.h5")


def test_version_prints_the_installed_version():
    result = runner.invoke(app, ["--version"])
    assert result.exit_code == 0 and result.output.strip() == __version__


def test_methods_lists_every_method_and_describes_one():
    result = runner.invoke(app, ["methods"])
    assert result.exit_code == 0 and "top-n" in result.output and "threshold" in result.output
    result = runner.invoke(app, ["methods", "top-n"])
    assert result.exit_code == 0
    assert "--force-cells" in _out(result) and "required" in _out(result)
    result = runner.invoke(app, ["methods", "nope"])
    assert result.exit_code == 1 and "Error: no method 'nope'" in _out(result)


def test_inspect_summarises_and_draws_the_knee(h5):
    result = runner.invoke(app, ["inspect", str(h5), "--width", "30"])
    assert result.exit_code == 0, result.output
    text = _out(result)
    assert "7,300" in text and "Gene Expression" in text
    assert "┤" in result.output


def test_inspect_warns_about_a_filtered_matrix(sim, tmp_path):
    path = write_10x_h5_v3(sim[sim.obs["truth"] == "cell"].copy(), tmp_path / "filtered.h5")
    result = runner.invoke(app, ["inspect", str(path)])
    assert result.exit_code == 0 and "looks like a filtered matrix" in _out(result)


def test_call_writes_every_output_it_is_asked_for(h5, tmp_path):
    result = runner.invoke(
        app,
        [
            "call", "top-n", str(h5), "--n-cells", "300",
            "-o", str(tmp_path / "calls.h5ad"),
            "--barcodes", str(tmp_path / "cells.tsv.gz"),
            "--table", str(tmp_path / "calls.tsv"),
            "--json",
        ],
    )
    assert result.exit_code == 0, result.output
    summary = json.loads(result.stdout[result.stdout.index("{"):])
    assert summary["n_cells"] == 300
    assert len(gzip.open(tmp_path / "cells.tsv.gz").read().splitlines()) == 300
    assert (tmp_path / "calls.tsv").read_text().splitlines()[0] == "barcode\ttotal\tis_cell"
    stored = ad.read_h5ad(tmp_path / "calls.h5ad")
    assert stored.n_obs == 7300 and stored.obs["top_n"].sum() == 300


def test_cells_only_keeps_the_called_barcodes(h5, tmp_path):
    result = runner.invoke(
        app, ["call", "threshold", str(h5), "--min-counts", "150", "-o", str(tmp_path / "c.h5ad"), "--cells-only"]
    )
    assert result.exit_code == 0, result.output
    stored = ad.read_h5ad(tmp_path / "c.h5ad")
    assert stored.obs["threshold"].all()


def test_call_without_outputs_says_nothing_was_written(h5):
    result = runner.invoke(app, ["call", "threshold", str(h5), "--min-counts", "150", "--strict"])
    assert result.exit_code == 0 and "Nothing written" in _out(result)


def test_choice_options_are_checked_by_the_command_line(h5):
    result = runner.invoke(app, ["call", "top-n", str(h5), "--n-cells", "5", "--ties", "random"])
    assert result.exit_code == 2


def test_errors_are_one_line_and_exit_1(h5, tmp_path):
    result = runner.invoke(app, ["call", "top-n", str(h5)])
    assert result.exit_code == 1 and "Error: top-n needs --n-cells" in _out(result)
    result = runner.invoke(app, ["call", "top-n", str(h5), "--n-cells", "5", "-o", str(h5)])
    assert result.exit_code == 1 and "also the input" in _out(result)


def test_reads_may_come_from_a_table(h5, sim, tmp_path):
    table = tmp_path / "reads.csv"
    table.write_text("barcode,reads\n" + "".join(f"{b},{i}\n" for i, b in enumerate(sim.obs_names)))
    result = runner.invoke(
        app,
        ["call", "top-n", str(h5), "--n-cells", "3", "--unit", "reads", "--reads", str(table), "--barcodes", str(tmp_path / "b.txt")],
    )
    assert result.exit_code == 0, result.output
    # the three with the most reads, written in input order
    assert (tmp_path / "b.txt").read_text().split() == list(sim.obs_names[-3:])
    result = runner.invoke(app, ["call", "top-n", str(h5), "--n-cells", "3", "--reads", "nowhere"])
    assert result.exit_code == 1 and "neither an obs column nor a file" in _out(result)


def test_compare_summarises_writes_and_plots(h5, tmp_path):
    result = runner.invoke(
        app,
        [
            "compare", str(h5), "-m", "top-n:n_cells=300", "-m", "threshold:min_counts=150",
            "--table", str(tmp_path / "t.csv"), "--plot", str(tmp_path / "knee.png"), "-o", str(tmp_path / "all.h5ad"),
        ],
    )
    assert result.exit_code == 0, result.output
    assert "top_n_n_cells_300" in _out(result)
    header = (tmp_path / "t.csv").read_text().splitlines()[0]
    assert header == "barcode,total,top_n_n_cells_300,threshold_min_counts_150"
    assert (tmp_path / "knee.png").stat().st_size > 1000
    assert set(ad.read_h5ad(tmp_path / "all.h5ad").uns["sc_droplets"]) == {"top_n_n_cells_300", "threshold_min_counts_150"}


def test_compare_json_holds_calls_and_jaccard(h5):
    result = runner.invoke(app, ["compare", str(h5), "-m", "top-n:n_cells=300", "-m", "top-n:n_cells=200", "--json"])
    assert result.exit_code == 0, result.output
    data = json.loads(result.stdout[result.stdout.index("{"):])
    assert data["calls"]["top_n_n_cells_200"]["n_cells"] == 200
    assert data["jaccard"]["top_n_n_cells_200"]["top_n_n_cells_300"] == pytest.approx(200 / 300)


def test_plot_knee_runs_methods_or_draws_stored_calls(h5, tmp_path):
    result = runner.invoke(app, ["plot", "knee", str(h5), "-m", "top-n:n_cells=300", "-o", str(tmp_path / "k.svg")])
    assert result.exit_code == 0, result.output
    assert (tmp_path / "k.svg").read_text().startswith("<?xml")
    runner.invoke(app, ["call", "top-n", str(h5), "--n-cells", "300", "-o", str(tmp_path / "c.h5ad")])
    result = runner.invoke(app, ["plot", "knee", str(tmp_path / "c.h5ad"), "-o", str(tmp_path / "k.png"), "--cumulative"])
    assert result.exit_code == 0, result.output
    result = runner.invoke(app, ["plot", "knee", str(tmp_path / "c.h5ad"), "--key", "nope", "-o", str(tmp_path / "x.png")])
    assert result.exit_code == 1 and "no stored call" in _out(result)


def test_plot_diagnostics_runs_a_method_or_draws_a_stored_call(h5, tmp_path):
    result = runner.invoke(app, ["plot", "diagnostics", str(h5), "-m", "starsolo-cellranger2.2", "-o", str(tmp_path / "d.png")])
    assert result.exit_code == 0, result.output
    assert (tmp_path / "d.png").stat().st_size > 1000
    runner.invoke(app, ["call", "top-n", str(h5), "--n-cells", "300", "-o", str(tmp_path / "c.h5ad")])
    result = runner.invoke(app, ["plot", "diagnostics", str(tmp_path / "c.h5ad"), "--key", "top_n", "-o", str(tmp_path / "e.png")])
    assert result.exit_code == 0, result.output
    result = runner.invoke(app, ["plot", "diagnostics", str(h5), "-o", str(tmp_path / "f.png")])
    assert result.exit_code == 1 and "exactly one of --method and --key" in _out(result)
    result = runner.invoke(app, ["plot", "diagnostics", str(tmp_path / "c.h5ad"), "--key", "nope", "-o", str(tmp_path / "g.png")])
    assert result.exit_code == 1 and "no stored call nope" in _out(result)
