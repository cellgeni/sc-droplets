"""Validation against the original tools on public data (marker `integration`).

Cell Ranger: our call on a dataset's raw matrix must give exactly the barcodes
of 10x's own filtered matrix, made by the release named. 10x does not publish
the `--expect-cells` of its runs; the values here are the ones that reproduce
the calls, found by trying the usual values (1,000 for the "1k" datasets, 5,000
for pbmc4k).

Cell Ranger 10.0 is the exception: it matches up to the order of tied totals at
the ambient boundary (see its test).

STARsolo: our call must give exactly the barcodes `STAR --runMode
soloCellFiltering` kept, saved in tests/oracles/starsolo by
oracles/starsolo/run.sh (STAR's static Linux binaries).

The datasets are downloaded once (tests/datasets.py) into $SC_DROPLETS_DATA.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pytest

import sc_droplets as scd

from datasets import filtered_barcodes, raw_path

pytestmark = pytest.mark.integration

ORACLES = Path(__file__).parent / "oracles"

_READ_CACHE = {}


def _read(name):
    if name not in _READ_CACHE:
        _READ_CACHE.clear()  # one raw matrix in memory at a time
        _READ_CACHE[name] = scd.read(raw_path(name))
    return _READ_CACHE[name]


def _called(result) -> set:
    return set(np.asarray(result.barcodes)[result.is_cell])


CELLRANGER_CASES = [
    # dataset, method spec, parameters, genome of a Cell Ranger 2 MEX
    ("pbmc4k_2_1", "cellranger@2.1.0", {"expect_cells": 5000}, "GRCh38"),
    ("pbmc4k_2_1", "cellranger-ordmag@2.1.0", {"expect_cells": 5000}, "GRCh38"),
    ("pbmc_1k_v2_3_0", "cellranger@3.0.0", {"expect_cells": 1000}, None),
    ("pbmc_1k_v3_3_0", "cellranger@3.0.0", {"expect_cells": 1000}, None),
    # 7.2+ estimates the expected cells itself; these runs gave none.
    ("pbmc5k_nextgem_8_0", "cellranger@8.0.0", {}, None),  # with antibody aggregates removed
    ("pbmc5k_gemx_9_0", "cellranger@9.0.0", {}, None),
]


@pytest.mark.parametrize("name, spec, params, genome", CELLRANGER_CASES, ids=[f"{c[0]}-{c[1]}" for c in CELLRANGER_CASES])
def test_cell_ranger_calls_the_cells_10x_called(name, spec, params, genome):
    result = scd.call(_read(name), spec, inplace=False, **params)
    expected = set(filtered_barcodes(name, genome))
    ours = _called(result)
    assert ours == expected, f"missing {sorted(expected - ours)[:10]}, extra {sorted(ours - expected)[:10]}"


STAR_PARAMS = {"starsolo-topcells": {"n_cells": 3000}}
STAR_ORACLES = sorted((ORACLES / "starsolo").glob("*/*.txt.gz"))


def _star_case(path: Path):
    method, _, version = path.name[: -len(".txt.gz")].partition("@")
    return path.parent.name, method, version


@pytest.mark.parametrize("path", STAR_ORACLES, ids=[f"{p.parent.name}-{p.name[:-7]}" for p in STAR_ORACLES])
def test_starsolo_calls_the_cells_star_kept(path):
    name, method, version = _star_case(path)
    result = scd.call(_read(name), f"{method}@{version}", inplace=False, **STAR_PARAMS.get(method, {}))
    expected = set(gzip.open(path).read().decode().split())
    ours = _called(result)
    assert ours == expected, f"missing {sorted(expected - ours)[:10]}, extra {sorted(ours - expected)[:10]}"


def test_there_are_star_oracles_to_check():
    assert len(STAR_ORACLES) >= 12


def test_cell_ranger_10_0_differs_only_at_the_fdr_boundary():
    """Cell Ranger 10.0 picks its ambient droplets by numpy's default sort, and
    1,509 barcodes tie at the boundary rank; which of them fall inside depends on
    the numpy build and processor 10x ran on. The OrdMag cells must match
    exactly, and every disagreement must sit at the FDR cut."""
    name = "spinal_cord_10_0"
    result = scd.call(_read(name), "cellranger@10.0.0", inplace=False)
    expected = set(filtered_barcodes(name))
    barcodes = np.asarray(result.barcodes)
    classes = np.asarray(result.classes)
    assert set(barcodes[classes == "initial"]) <= expected
    different = _called(result) ^ expected
    assert len(different) <= 0.002 * len(expected)
    fdr = result.columns["fdr"][np.isin(barcodes, list(different))]
    assert np.all((fdr > 0.0005) & (fdr < 0.002)), np.sort(fdr)


FLEX_SAMPLES = ["flex_lymphnode1_7_2", "flex_lymphnode2_7_2", "flex_spleen1_7_2", "flex_spleen2_7_2"]


def test_cell_ranger_7_2_flex_calls_every_sample_after_the_high_occupancy_filter():
    """Each sample of a 4-plex Flex run is called on its own raw matrix; the
    high-occupancy GEM filter then runs on the union, as Cell Ranger merges them."""
    from sc_droplets.methods.cellranger import high_occupancy_gems

    calls = {}
    for name in FLEX_SAMPLES:
        calls[name] = _called(scd.call(_read(name), "cellranger@7.2.0", inplace=False))
    removed, info = high_occupancy_gems(set().union(*calls.values()), "7.2.0")
    assert info["threshold"] == 3 and len(removed) == 28
    for name in FLEX_SAMPLES:
        ours = calls[name] - set(removed)
        expected = set(filtered_barcodes(name))
        assert ours == expected, f"{name}: missing {sorted(expected - ours)[:5]}, extra {sorted(ours - expected)[:5]}"
