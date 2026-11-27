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

CellBender: the priors and the barcodes of every class (probable cells, other
analysed droplets, surely-empty droplets; as a count and a SHA-256) of CellBender's own
SingleCellRNACountsDataset (0.1.0, 0.2.2, 0.3.0, 0.3.2, 0.4.0, each from its
tagged source with torch, pyro and numpy before 1.25) on the raw matrix; saved in
tests/oracles/cellbender by oracles/cellbender/run.sh. Counts and classes are exact;
log and exp results to the last bits, which depend on the processor (see
_assert_cellbender_value).

bustools: the list of bioconda's bustools whitelist (0.39.3, 0.41.0, 0.43.2,
0.45.1) on a BUS file written from the matrix, one record per UMI; saved in
tests/oracles/bustools by oracles/bustools/run.sh.

alevin-fry: the permit lists of programs built from alevin-fry's own source at
0.1.0, 0.11.2, 0.17.1 and 0.19.0 (its knee and rules copied out of each tag by
oracles/alevinfry/harness.py), on the matrix's per-barcode UMI totals; saved in
tests/oracles/alevinfry by oracles/alevinfry/run.sh.

alevin: the initial whitelist of programs built from salmon's own Alevin.cpp at
0.11.3, 0.13.1, 0.14.0, 1.1.0 and 1.10.3, on the matrix's per-barcode UMI totals;
saved in tests/oracles/alevin by oracles/alevin/run.sh. The cells are checked up to
the barcodes tied at the knee, whose order in alevin comes from its hash map.

Drop-seq (dropbead), dropEst and zUMIs: their own R functions (dropbead's
estimateCellNumber, dropestr's EstimateCellsNumber, zUMIs' cellBC with inflection's
uik), sourced from each tool's repository at a commit or tag and run in R from
conda-forge on the matrix's UMI totals; saved in tests/oracles/{dropseq,dropest,zumis}
by oracles/{dropseq,dropest,zumis}/run.sh.

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


# ---------------------------------------------------------------------------
# CellBender
# ---------------------------------------------------------------------------

# The cases of oracles/cellbender/oracle.py as sc-droplets parameters (None: the dataset's expected cells).
CELLBENDER_CASES = {
    "default": {},
    "expected": {"expected_cells": None},
    "total": {"total_droplets_included": 15000},
    "low30": {"low_count_threshold": 30},
    "forced": {"force_cell_umi_prior": 3000.0, "force_empty_umi_prior": 30.0},
    "expected_total": {"expected_cells": None, "total_droplets_included": 15000},
}
CELLBENDER_EXPECTED = {"pbmc_1k_v3_3_0": 1000, "pbmc4k_2_1": 5000}
# CellBender's names for what sc-droplets stores under another name; None: checked through the classes,
# or not kept (0.1-0.2's uncapped n_cells and cell_prob, 0.3's transition_point before user values).
CELLBENDER_NAMES = {"empty_UMI_threshold": "empty_umi_threshold", "n_cells": None, "cell_prob": None,
                    "transition_point": None, "n_analysed": None, "n_empty": None}


def _cellbender_rows():
    rows = []
    for path in sorted((ORACLES / "cellbender").glob("*/*/priors.tsv")):
        cases = {}
        for line in path.read_text().splitlines():
            case, name, value = line.split("\t")
            cases.setdefault(case, {})[name] = value
        rows += [(path.parent.parent.name, path.parent.name, case, values) for case, values in cases.items()]
    return rows


CELLBENDER_ROWS = _cellbender_rows()


def test_there_are_cellbender_oracles_to_check():
    assert len(CELLBENDER_ROWS) >= 52


@pytest.mark.parametrize("name, release, case, values", CELLBENDER_ROWS,
                         ids=[f"{r[0]}-{r[1]}-{r[2]}" for r in CELLBENDER_ROWS])
def test_cellbender_prior_reproduces_cellbenders_dataset(name, release, case, values):
    params = {k: CELLBENDER_EXPECTED[name] if v is None else v for k, v in CELLBENDER_CASES[case].items()}
    adata = _read(name)
    if "error" in values:
        with pytest.raises(ValueError):
            scd.call(adata, f"cellbender-prior@{release}", inplace=False, **params)
        return
    result = scd.call(adata, f"cellbender-prior@{release}", inplace=False, **params)
    for theirs, value in values.items():
        ours = CELLBENDER_NAMES.get(theirs, theirs)
        if ours is not None and not theirs.startswith("class:"):
            _assert_cellbender_value(result.thresholds[ours], float(value), theirs)
    classes = np.asarray(result.classes.astype(str))
    barcodes = np.asarray(result.barcodes)
    for kind in ("cell", "transition", "empty"):
        assert _digest(barcodes[classes == kind]) == values[f"class:{kind}"], kind


def _assert_cellbender_value(ours: float, theirs: float, name: str) -> None:
    """Counts exactly; other values to 4 units in the last place. numpy computes float64 log and
    exp with SVML on Linux processors with AVX-512, which can differ from the C library in the
    last bit, so CellBender's own priors depend on the machine (2 ulps between a Linux runner and
    macOS for pbmc4k's empty prior, 36.234075926476564 against ...56)."""
    import math

    if theirs.is_integer() or ours == theirs:
        assert ours == theirs, name
    else:
        assert abs(ours - theirs) <= 4 * math.ulp(max(abs(ours), abs(theirs))), (name, ours, theirs)


def _digest(barcodes) -> str:
    """oracles/cellbender/oracle.py's `count:sha256` of the barcodes, sorted, one per line."""
    import hashlib

    text = "\n".join(sorted(str(b) for b in barcodes)).encode()
    return f"{len(barcodes)}:{hashlib.sha256(text).hexdigest()}"


# ---------------------------------------------------------------------------
# bustools
# ---------------------------------------------------------------------------

BUSTOOLS_CASES = {"default": {}, "threshold500": {"threshold": 500}}


def _bustools_rows():
    rows = []
    for path in sorted((ORACLES / "bustools").glob("*/*/cases.tsv")):
        for line in path.read_text().splitlines():
            case, threshold, _lines, unique = line.split("\t")
            rows.append((path.parent.parent.name, path.parent.name, case, int(threshold), int(unique)))
    return rows


BUSTOOLS_ROWS = _bustools_rows()


def test_there_are_bustools_oracles_to_check():
    assert len(BUSTOOLS_ROWS) >= 16


@pytest.mark.parametrize("name, release, case, threshold, cells", BUSTOOLS_ROWS,
                         ids=[f"{r[0]}-{r[1]}-{r[2]}" for r in BUSTOOLS_ROWS])
def test_bustools_allowlist_lists_bustools_barcodes(name, release, case, threshold, cells):
    result = scd.call(_read(name), f"bustools-allowlist@{release}", inplace=False, **BUSTOOLS_CASES[case])
    assert result.thresholds["umis"] == threshold
    with gzip.open(ORACLES / "bustools" / name / release / f"{case}.txt.gz", "rt") as handle:
        expected = set(handle.read().split())
    ours = {str(b).split("-", 1)[0] for b in np.asarray(result.barcodes)[result.is_cell]}
    assert len(expected) == cells
    assert ours == expected


# ---------------------------------------------------------------------------
# alevin-fry
# ---------------------------------------------------------------------------

# The cases of oracles/alevinfry/oracle.py as sc-droplets parameters (None: the dataset's expected cells).
ALEVINFRY_CASES = {
    "knee": {},
    "expect": {"expect_cells": None},
    "force": {"force_cells": None},
    "unfiltered": {"unfiltered": True},
    "flex_knee": {"flex": True},
    "flex_expect": {"flex": True, "expect_cells": None},
    "flex_force": {"flex": True, "force_cells": None},
}
ALEVINFRY_EXPECTED = {"pbmc_1k_v3_3_0": 1000, "pbmc4k_2_1": 5000}


def _alevinfry_rows():
    rows = []
    for path in sorted((ORACLES / "alevinfry").glob("*/*/cases.tsv")):
        for line in path.read_text().splitlines():
            rows.append((path.parent.parent.name, path.parent.name, *line.split("\t")))
    return rows


ALEVINFRY_ROWS = _alevinfry_rows()


def test_there_are_alevinfry_oracles_to_check():
    assert len(ALEVINFRY_ROWS) >= 32


@pytest.mark.parametrize("name, release, case, cells, smallest, threshold", ALEVINFRY_ROWS,
                         ids=[f"{r[0]}-{r[1]}-{r[2]}" for r in ALEVINFRY_ROWS])
def test_alevinfry_permit_list_keeps_alevinfrys_barcodes(name, release, case, cells, smallest, threshold):
    params = {k: ALEVINFRY_EXPECTED[name] if v is None else v for k, v in ALEVINFRY_CASES[case].items()}
    result = scd.call(_read(name), f"alevinfry-permit-list@{release}", inplace=False, **params)
    # every rule keeps the barcodes with count >= a threshold: the number kept and the smallest
    # count kept fix them
    assert result.n_cells == int(cells)
    assert np.asarray(result.totals)[result.is_cell].min() == int(smallest)
    if threshold != "-":
        assert result.thresholds["umis"] == int(threshold)


# ---------------------------------------------------------------------------
# alevin
# ---------------------------------------------------------------------------

# The cases of oracles/alevin/oracle.py as sc-droplets parameters (None: the dataset's expected cells).
ALEVIN_CASES = {
    "knee": {},
    "expect": {"expect_cells": None},
    "force": {"force_cells": None},
    "keep": {"keep_cb_fraction": 0.005},
}
ALEVIN_EXPECTED = {"pbmc_1k_v3_3_0": 1000, "pbmc4k_2_1": 5000}


def _alevin_rows():
    rows = []
    for path in sorted((ORACLES / "alevin").glob("*/*/cases.tsv")):
        cases = {}
        for line in path.read_text().splitlines():
            case, name, value = line.split("\t")
            cases.setdefault(case, {})[name] = value
        rows += [(path.parent.parent.name, path.parent.name, case, values) for case, values in cases.items()]
    return rows


ALEVIN_ROWS = _alevin_rows()


def test_there_are_alevin_oracles_to_check():
    assert len(ALEVIN_ROWS) >= 34


@pytest.mark.parametrize("name, release, case, values", ALEVIN_ROWS, ids=[f"{r[0]}-{r[1]}-{r[2]}" for r in ALEVIN_ROWS])
def test_alevin_initial_whitelist_is_alevins(name, release, case, values):
    import hashlib

    params = {k: ALEVIN_EXPECTED[name] if v is None else v for k, v in ALEVIN_CASES[case].items()}
    adata = _read(name)
    if "error" in values:
        with pytest.raises(ValueError):
            scd.call(adata, f"alevin-initial-whitelist@{release}", inplace=False, **params)
        return
    result = scd.call(adata, f"alevin-initial-whitelist@{release}", inplace=False, **params)
    if "knee" in values:
        assert result.thresholds["knee_rank"] == int(values["knee"])
        assert result.thresholds.get("gauss_rank", 0) == int(values["gauss"])
    assert result.thresholds["cells"] == int(values["cells"])
    assert result.thresholds["low_confidence"] == int(values["low"])
    totals = np.asarray(result.totals)
    classes = np.asarray(result.classes.astype(str))
    # the cells, up to the barcodes tied at the knee: as many, the same smallest count, and every
    # barcode above it
    cells = totals[classes == "cell"]
    assert cells.size == int(values["cells"])
    assert cells.min() == int(values["cells_min"])
    assert np.sum(totals > cells.min()) == int(values["cells_above"])
    white = np.asarray(result.barcodes)[classes != "background"]
    assert white.size == int(values["whitelist"])
    assert hashlib.sha256("\n".join(sorted(map(str, white))).encode()).hexdigest() == values["whitelist_digest"]


# ---------------------------------------------------------------------------
# Drop-seq (dropbead), dropEst, zUMIs
# ---------------------------------------------------------------------------


def _r_rows(tool: str, filename: str):
    rows = []
    for path in sorted((ORACLES / tool).glob(f"*/*/{filename}")):
        rows += [(path.parent.parent.name, path.parent.name, line.split("\t")) for line in path.read_text().splitlines()]
    return rows


DROPBEAD_ROWS = _r_rows("dropseq", "cells.tsv")
DROPEST_ROWS = _r_rows("dropest", "bounds.tsv")
ZUMIS_ROWS = _r_rows("zumis", "cells.tsv")


def test_there_are_r_tool_oracles_to_check():
    assert len(DROPBEAD_ROWS) >= 4 and len(DROPEST_ROWS) >= 4 and len(ZUMIS_ROWS) >= 16


@pytest.mark.parametrize("name, release, fields", DROPBEAD_ROWS, ids=[f"{r[0]}-{r[2][0]}" for r in DROPBEAD_ROWS])
def test_dropbead_counts_dropbeads_cells(name, release, fields):
    max_cells, cells = (int(v) for v in fields)
    result = scd.call(_read(name), f"dropseq-dropbead@{release}", max_cells=max_cells, inplace=False)
    assert result.thresholds["cells"] == cells == result.n_cells


@pytest.mark.parametrize("name, release, fields", DROPEST_ROWS, ids=[f"{r[0]}-lag{r[2][0]}" for r in DROPEST_ROWS])
def test_dropest_bounds_are_dropestrs(name, release, fields):
    lag, expected, low, high = float(fields[0]), *(int(v) for v in fields[1:])
    result = scd.call(_read(name), f"dropest-bounds@{release}", lag=lag, inplace=False)
    assert (result.thresholds["expected"], result.thresholds["min"], result.thresholds["max"]) == (expected, low, high)
    assert result.n_cells == expected


@pytest.mark.parametrize("name, release, fields", ZUMIS_ROWS, ids=[f"{r[0]}-{r[1]}-{r[2][0]}" for r in ZUMIS_ROWS])
def test_zumis_cell_barcodes_are_zumis(name, release, fields):
    floor, cells, smallest = (int(v) for v in fields)
    version = release.split("-inflection")[0]
    result = scd.call(_read(name), f"zumis-cell-barcodes@{version}", n_reads_per_cell=floor, inplace=False)
    with gzip.open(ORACLES / "zumis" / name / release / f"cells_{floor}.txt.gz", "rt") as handle:
        expected = set(handle.read().split())
    ours = set(np.asarray(result.barcodes)[result.is_cell])
    assert len(expected) == cells and ours == expected
    assert result.thresholds["umis"] == smallest


# ---------------------------------------------------------------------------
# Generic rules: scuttle's isOutlier (mad), scikit-image's Otsu and multi-Otsu
# ---------------------------------------------------------------------------


def _generic_rows(prefix: str, filename: str):
    rows = []
    for path in sorted((ORACLES / "generic").glob(f"*/{prefix}-*/{filename}")):
        rows += [(path.parent.parent.name, path.parent.name, line.split("\t")) for line in path.read_text().splitlines()]
    return rows


MAD_ROWS = _generic_rows("scuttle", "mad.tsv")
OTSU_ROWS = _generic_rows("skimage", "otsu.tsv")
MULTIOTSU_ROWS = _generic_rows("skimage", "multiotsu.tsv")


def test_there_are_generic_rule_oracles_to_check():
    assert len(MAD_ROWS) >= 24 and len(OTSU_ROWS) >= 32 and len(MULTIOTSU_ROWS) >= 64


@pytest.mark.parametrize("name, release, fields", MAD_ROWS, ids=[f"{r[0]}-{r[1]}-{r[2][0]}-{r[2][1]}" for r in MAD_ROWS])
def test_mad_limit_is_scuttles_isoutlier(name, release, fields):
    floor, nmads, limit, cells = float(fields[0]), float(fields[1]), float(fields[2]), int(fields[3])
    result = scd.call(_read(name), "mad", min_counts=floor, nmads=nmads, inplace=False)
    assert result.thresholds["umis"] == limit
    assert result.n_cells == cells


@pytest.mark.parametrize("name, release, fields", OTSU_ROWS, ids=[f"{r[0]}-{r[1]}-{r[2][0]}-{r[2][1]}" for r in OTSU_ROWS])
def test_otsu_threshold_is_scikit_images(name, release, fields):
    floor, how, threshold = float(fields[0]), fields[1], float(fields[2])
    result = scd.call(_read(name), "otsu", min_counts=floor, transform=how, inplace=False)
    assert result.thresholds["transformed"] == threshold


@pytest.mark.parametrize("name, release, fields", MULTIOTSU_ROWS,
                         ids=[f"{r[0]}-{r[1]}-{r[2][0]}-{r[2][1]}-{r[2][2]}" for r in MULTIOTSU_ROWS])
def test_multiotsu_thresholds_are_scikit_images(name, release, fields):
    floor, how, classes = float(fields[0]), fields[1], int(fields[2])
    expected = [float(v) for v in fields[3].split(",")]
    result = scd.call(_read(name), "multiotsu", min_counts=floor, transform=how, classes=classes, inplace=False)
    assert [result.thresholds[f"threshold_{i + 1}"] for i in range(classes - 1)] == expected
