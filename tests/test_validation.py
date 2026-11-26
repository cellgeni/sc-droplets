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

DropletUtils: knee, inflection, defaultDrops' cells and every emptyDrops and
emptyDropsCellRanger p-value, FDR and call must be those of DropletUtils itself,
release by release, saved in tests/oracles/dropletutils by
oracles/dropletutils/run.sh (bioconda's builds, R's set.seed(100)).

UMI-tools: the whitelist of UMI-tools' own functions, release by release (0.5.5,
1.0.0, 1.1.6, each with the numpy and scipy of its bioconda build), given the
matrix's per-barcode UMI totals; saved in tests/oracles/umitools by
oracles/umitools/run.sh.

Kneedle: the knee and sensitivity of kneed's own KneeLocator (0.6.0, 0.8.5, 0.8.6,
0.9.0, one per behaviour, each with its numpy, scipy and scikit-learn) on the
rank curves of the matrix's UMI totals, with yascp's search and with fixed
settings on log axes; saved in tests/oracles/kneedle by oracles/kneedle/run.sh.

cellsweep: the cutoff of cellsweep's own infer_empty_droplets (1.0.1, and 1.0.2
with each kneed it may run with) on a float32 AnnData of the matrix's UMI totals;
saved in tests/oracles/cellsweep by oracles/cellsweep/run.sh.

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


# ---------------------------------------------------------------------------
# DropletUtils
# ---------------------------------------------------------------------------

DU_ORACLES = sorted(p for p in (ORACLES / "dropletutils").glob("*/*") if p.is_dir() and p.parent.name != "fixtures")


def _du_case(path: Path):
    """(dataset, --tool-version, extra emptyDrops arguments) for one oracle folder."""
    release = path.name
    header = (path / "barcoderanks.txt").read_text().split("\n", 1)[0].split()
    extra = header[2:]
    return path.parent.name, release, extra


def _du_test(path: Path):
    lines = gzip.open(path, "rt").read().splitlines()
    meta = {line.split()[1]: " ".join(line.split()[2:]) for line in lines if line.startswith("#")}
    rows = np.array([line.split("\t") for line in lines if not line.startswith("#")])
    return meta, rows


def _same_c_library(path: Path) -> bool:
    """Whether R made this oracle on this platform's C library (glibc, or Apple's libm)."""
    import sys

    record = path / "platform.txt"
    platform = record.read_text().split("\n", 1)[0] if record.exists() else ""
    return ("linux" in platform and sys.platform.startswith("linux")) or ("darwin" in platform and sys.platform == "darwin")


def _assert_r_value(path: Path, ours: float, theirs: float) -> None:
    """Exact on R's own C library. 1.30's knee is interpolated from log10 values, so elsewhere it
    can differ in the last bits (5199.842223575369 on Linux, ...5353 on macOS, for pbmc_1k_v3)."""
    if _same_c_library(path):
        assert ours == theirs
    else:
        assert ours == pytest.approx(theirs, rel=1e-13, abs=0)


def test_there_are_dropletutils_oracles_to_check():
    assert len(DU_ORACLES) >= 13


@pytest.mark.parametrize("path", DU_ORACLES, ids=[f"{p.parent.name}-{p.name}" for p in DU_ORACLES])
def test_dropletutils_knees_and_defaultdrops_are_rs(path):
    name, release, extra = _du_case(path)
    adata = _read(name)
    lines = dict(line.split(" ", 1) for line in (path / "barcoderanks.txt").read_text().splitlines()[1:])
    ranks = scd.call(adata, f"dropletutils-barcoderanks@{release}", inplace=False)
    _assert_r_value(path, ranks.thresholds["knee"], float(lines["knee"]))
    _assert_r_value(path, ranks.thresholds["inflection"], float(lines["inflection"]))
    cells = np.array(gzip.open(path / "defaultdrops.txt.gz", "rt").read().split(), dtype=int)
    default = scd.call(adata, f"dropletutils-defaultdrops@{release}", inplace=False)
    assert np.array_equal(np.flatnonzero(default.is_cell), np.sort(cells))


DU_TESTS = [
    (p, method) for p in DU_ORACLES
    for method, file in (("dropletutils-emptydrops", "emptydrops.tsv.gz"), ("dropletutils-emptydropscellranger", "emptydropscellranger.tsv.gz"))
    if (p / file).exists()
]


@pytest.mark.parametrize("path, method", DU_TESTS, ids=[f"{p.parent.name}-{m.split('-')[1]}-{p.name}" for p, m in DU_TESTS])
def test_dropletutils_tests_give_rs_p_values_and_cells(path, method):
    name, release, extra = _du_case(path)
    meta, rows = _du_test(path / ("emptydrops.tsv.gz" if method.endswith("emptydrops") else "emptydropscellranger.tsv.gz"))
    params = {}
    if method.endswith("emptydrops") and "alphainf" in extra:
        params["alpha"] = "inf"  # 1.2's Dirichlet draws use R's rgamma, which is not reproduced
    result = scd.call(_read(name), f"{method}@{release}", inplace=False, **params)
    index = rows[:, 0].astype(int)
    pvalue = result.columns["pvalue"]
    assert np.array_equal(np.flatnonzero(~np.isnan(pvalue)), index)
    n_above = np.rint(pvalue[index] * 10001 - 1).astype(int)
    assert np.array_equal(n_above, rows[:, 2].astype(int))
    assert np.array_equal(result.columns["fdr"][index], rows[:, 3].astype(float))
    if meta.get("alpha") not in (None, "NA", "Inf"):
        assert result.diagnostics["alpha"] == float(meta["alpha"])
    if meta.get("retain") in (None, "NA"):
        # 1.0 keeps no metadata; its retain is barcodeRanks' knee, which the oracle records
        knee = dict(line.split(" ", 1) for line in (path / "barcoderanks.txt").read_text().splitlines()[1:])["knee"]
        retain = float(knee)
    else:
        retain = float(meta["retain"])
    _assert_r_value(path, result.diagnostics["retain"], retain)
    fdr = rows[:, 3].astype(float)
    cut = 0.001 if method.endswith("emptydrops") else 0.01
    expected = set(index[fdr <= cut]) | set(np.flatnonzero(np.rint(np.asarray(result.totals)) >= retain))
    assert set(np.flatnonzero(result.is_cell)) == expected


# ---------------------------------------------------------------------------
# UMI-tools
# ---------------------------------------------------------------------------

UT_ORACLES = sorted(p for p in (ORACLES / "umitools").glob("*/*/*.txt.gz"))


def _ut_case(path: Path):
    """Dataset, release, method spec and parameters of one UMI-tools oracle file."""
    case = path.name[: -len(".txt.gz")]
    rule, _, option = case.partition("_")
    params = {}
    if option.startswith("n"):
        params["n_cells"] = int(option[1:])
    elif option.startswith("expect"):
        params["expect_cells"] = int(option[len("expect"):])
    return path.parent.parent.name, path.parent.name, f"umitools-{rule}", params


def test_there_are_umitools_oracles_to_check():
    assert len(UT_ORACLES) >= 30


@pytest.mark.parametrize("path", UT_ORACLES, ids=[f"{p.parent.parent.name}-{p.parent.name}-{p.name[:-7]}" for p in UT_ORACLES])
def test_umitools_whitelist_calls_umitools_barcodes(path):
    name, release, method, params = _ut_case(path)
    expected = gzip.open(path).read().decode().split()
    if expected == ["none"]:
        with pytest.raises(ValueError, match="No local minima was accepted"):
            scd.call(_read(name), f"{method}@{release}", inplace=False, **params)
        result = scd.call(_read(name), f"{method}@{release}", inplace=False, allow_threshold_error=True, **params)
        assert result.n_cells == 0
        return
    result = scd.call(_read(name), f"{method}@{release}", inplace=False, **params)
    ours = _called(result)
    assert ours == set(expected), f"missing {sorted(set(expected) - ours)[:10]}, extra {sorted(ours - set(expected))[:10]}"


# ---------------------------------------------------------------------------
# Kneedle (kneed)
# ---------------------------------------------------------------------------

# The cases of oracles/kneedle/oracle.py, as parameters of `kneedle`.
KNEEDLE_CASES = {
    "yascp_interp1d": {"lower": 1000},
    "yascp_polynomial": {"lower": 1000, "interp_method": "polynomial"},
    "loglog_concave": {"lower": 10, "log_x": True, "log_y": True, "curve": "concave", "sensitivity": 1.0},
    "loglog_concave_online": {"lower": 100, "log_x": True, "log_y": True, "curve": "concave", "sensitivity": 1.0, "online": True},
    "loglog_convex": {"lower": 100, "log_x": True, "log_y": True, "curve": "convex", "sensitivity": 1.0},
}


def _kneedle_rows():
    rows = []
    for path in sorted((ORACLES / "kneedle").glob("*/*/knees.tsv")):
        for line in path.read_text().splitlines():
            case, rank, sensitivity = line.split("\t")
            rows.append((path.parent.parent.name, path.parent.name, case, rank, float(sensitivity)))
    return rows


KNEEDLE_ROWS = _kneedle_rows()


def test_there_are_kneedle_oracles_to_check():
    assert len(KNEEDLE_ROWS) >= 40


@pytest.mark.parametrize("name, release, case, rank, sensitivity", KNEEDLE_ROWS, ids=[f"{r[0]}-{r[1]}-{r[2]}" for r in KNEEDLE_ROWS])
def test_kneedle_finds_kneeds_knee(name, release, case, rank, sensitivity):
    params = KNEEDLE_CASES[case]
    if rank == "none":
        with pytest.raises(ValueError, match="no knee"):
            scd.call(_read(name), f"kneedle@{release}", inplace=False, **params)
        return
    result = scd.call(_read(name), f"kneedle@{release}", inplace=False, **params)
    assert result.thresholds["knee_rank"] == int(rank)
    assert result.thresholds["sensitivity"] == sensitivity
    assert result.n_cells == int(rank)


# ---------------------------------------------------------------------------
# cellsweep
# ---------------------------------------------------------------------------

# The cases of oracles/cellsweep/oracle.py: parameters, and the smallest total of the barcodes given.
CELLSWEEP_CASES = {
    "default": ({}, 0),
    "mx_filter": ({"empty_method": "mx_filter"}, 0),
    "expected": ({}, 0),  # expected_cells added from the dataset
    "mx_filter_detected": ({"empty_method": "mx_filter"}, 1),
    "mx_filter_over100": ({"empty_method": "mx_filter"}, 100),
}
CELLSWEEP_EXPECTED = {"pbmc_1k_v3_3_0": 1000, "pbmc4k_2_1": 5000}


def _cellsweep_rows():
    rows = []
    for path in sorted((ORACLES / "cellsweep").glob("*/*/cutoffs.tsv")):
        release, _, kneed_version = path.parent.name.partition("-kneed")
        for line in path.read_text().splitlines():
            case, cutoff, cells = line.split("\t")
            rows.append((path.parent.parent.name, release, kneed_version or None, case, cutoff, cells))
    return rows


CELLSWEEP_ROWS = _cellsweep_rows()


def test_there_are_cellsweep_oracles_to_check():
    assert len(CELLSWEEP_ROWS) >= 40


@pytest.mark.parametrize(
    "name, release, kneed_version, case, cutoff, cells",
    CELLSWEEP_ROWS,
    ids=[f"{r[0]}-{r[1]}{'-kneed' + r[2] if r[2] else ''}-{r[3]}" for r in CELLSWEEP_ROWS],
)
def test_cellsweep_finds_cellsweeps_cutoff(name, release, kneed_version, case, cutoff, cells):
    params, lowest = CELLSWEEP_CASES[case]
    params = dict(params)
    if case == "expected":
        params["expected_cells"] = CELLSWEEP_EXPECTED[name]
    if kneed_version:
        params["kneed_version"] = kneed_version
    adata = _read(name)
    if lowest:
        totals = np.asarray(adata.X.sum(axis=1)).ravel()
        adata = adata[totals >= lowest]
    if cutoff == "error":
        with pytest.raises(ValueError):
            scd.call(adata, f"cellsweep-empty@{release}", inplace=False, **params)
        return
    result = scd.call(adata, f"cellsweep-empty@{release}", inplace=False, **params)
    assert result.thresholds["umis"] == float(cutoff)
    assert result.n_cells == int(cells)
