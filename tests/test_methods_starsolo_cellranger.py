"""STARsolo and Cell Ranger methods on small matrices: every branch, by hand.

The public-data comparisons with the original tools are in test_validation.py
(marker `integration`); these cover the cases the public data does not reach.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

import sc_droplets as scd
from sc_droplets.methods.cellranger import gem_groups, genomes, ordmag_v2


def _totals(values, barcodes=None, genome=None):
    """An AnnData whose barcodes have these totals (one gene per count, cycling over 4 genes)."""
    values = np.asarray(values, dtype=np.int64)
    rows = np.repeat(np.arange(len(values)), values)
    cols = np.arange(rows.size) % 4
    matrix = sp.csr_matrix((np.ones(rows.size, dtype=np.int64), (rows, cols)), shape=(len(values), 4))
    adata = ad.AnnData(X=matrix)
    adata.obs_names = barcodes if barcodes is not None else [f"B{i:04d}-1" for i in range(len(values))]
    adata.var["feature_types"] = "Gene Expression"
    if genome is not None:
        adata.var["genome"] = genome
    return adata


# --- STARsolo CellRanger2.2 ---------------------------------------------------


def test_cellranger22_cuts_at_the_31st_largest_over_ten():
    # 40 barcodes: the 31st largest is 1000, so the cut is 100.
    values = [5000] * 30 + [1000] + [150, 100, 99, 50, 10, 1, 1, 0, 0]
    result = scd.call(_totals(values), "starsolo-cellranger2.2", inplace=False)
    assert result.thresholds == {"robust_max": 1000.0, "umis": 100.0}
    assert result.n_cells == 33
    for version in ("2.7.3a", "2.7.5a", "2.7.8a"):
        assert scd.call(_totals(values), f"starsolo-cellranger2.2@{version}", inplace=False).n_cells == 33


def test_cellranger22_rounds_the_cut_half_up():
    values = [5000] * 30 + [1005] + [101, 100, 99]
    assert scd.call(_totals(values), "starsolo-cellranger2.2", inplace=False).thresholds["umis"] == 101.0


def test_cellranger22_with_fewer_barcodes_than_the_index_uses_the_smallest():
    result = scd.call(_totals([50, 40, 30, 0]), "starsolo-cellranger2.2", inplace=False)
    assert result.thresholds["robust_max"] == 30.0 and result.n_cells == 3


def test_starsolo_methods_with_no_counts_call_nothing():
    for spec in ("starsolo-cellranger2.2", "starsolo-topcells:n_cells=2", "starsolo-emptydrops-cr"):
        result = scd.call(_totals([0, 0, 0]), spec, inplace=False)
        assert result.n_cells == 0 and result.notes == ["no barcode has counts"]


# --- STARsolo TopCells ----------------------------------------------------------


def test_topcells_keeps_n_plus_one_and_ties():
    values = [90, 80, 70, 70, 60, 0]
    assert scd.call(_totals(values), "starsolo-topcells", n_cells=2, inplace=False).n_cells == 4


def test_topcells_before_2_7_6a_keeps_every_detected_barcode():
    values = [90, 80, 70, 60, 5, 0]
    assert scd.call(_totals(values), "starsolo-topcells@2.7.5c", n_cells=2, inplace=False).n_cells == 5
    with pytest.raises(ValueError, match="past the end"):
        scd.call(_totals(values), "starsolo-topcells@2.7.4a", n_cells=10, inplace=False)


# --- soloBasicCellFilter.awk --------------------------------------------------------


def test_basic_awk_cut_is_unrounded():
    values = [5000] * 30 + [1005] + [101, 100, 99, 0]
    result = scd.call(_totals(values), "starsolo-basic-awk", inplace=False)
    assert result.thresholds["umis"] == pytest.approx(100.5)
    assert result.n_cells == 32


def test_basic_awk_exact_cells_keeps_one_more_and_a_cut_of_zero_keeps_everything():
    values = [90, 80, 70, 60, 50, 0]
    assert scd.call(_totals(values), "starsolo-basic-awk", exact_cells=2, inplace=False).n_cells == 3
    everything = scd.call(_totals(values), "starsolo-basic-awk", exact_cells=10, inplace=False)
    assert everything.thresholds["umis"] == 50.0 and everything.n_cells == 5
    # Fewer barcodes than the robust-max index: awk reads a missing element as 0.
    assert scd.call(_totals([3, 2, 0]), "starsolo-basic-awk", inplace=False).n_cells == 3


# --- STARsolo EmptyDrops_CR --------------------------------------------------------


def test_emptydrops_cr_rescues_cells_below_the_cut_and_no_empties():
    # Cells with widely spread sizes, so some fall below the CellRanger2.2 cut.
    adata = scd.simulate_droplets(n_cells=300, n_empty=3000, n_low=3000, n_genes=200, cell_umis=1500, cell_sigma=1.0, seed=7)
    # umi_min=150 keeps the candidates to the small cells and a few empties, so 500
    # simulations can reach the FDR.
    result = scd.call(adata, "starsolo-emptydrops-cr", ind_min=2000, ind_max=5000, umi_min=150, sim_n=500)
    truth = (adata.obs["truth"] == "cell").to_numpy()
    classes = adata.obs["starsolo_emptydrops_cr_class"]
    assert (classes == "rescued").sum() > 0
    assert not (result.is_cell & ~truth).any()
    tested = classes.isin(["rescued", "candidate"]).to_numpy()
    assert np.isfinite(adata.obs["starsolo_emptydrops_cr_pvalue"].to_numpy()[tested]).all()
    assert result.diagnostics["ambient_ranks"] == [2000, 5000]


def test_emptydrops_cr_with_too_few_barcodes_keeps_the_initial_cells(adata):
    result = scd.call(adata, "starsolo-emptydrops-cr", inplace=False)
    assert any("skips EmptyDrops_CR" in note for note in result.notes)
    simple = scd.call(adata, "starsolo-cellranger2.2", inplace=False)
    assert np.array_equal(result.is_cell, simple.is_cell)


def test_emptydrops_cr_stops_without_candidates(adata):
    result = scd.call(adata, "starsolo-emptydrops-cr", ind_min=2000, ind_max=5000, umi_min=10**6, sim_n=10, inplace=False)
    assert "no candidate barcodes" in result.notes


# --- Cell Ranger -------------------------------------------------------------------


def test_gem_groups_come_from_the_barcode_suffix():
    assert gem_groups(np.array(["AAA-1", "CCC-2", "GGG", "TTT-x"])).tolist() == [1, 2, 1, 1]


def test_genomes_are_listed_in_order_of_first_appearance():
    features = pd.DataFrame({"genome": ["mm10", "hg19", "mm10"]})
    assert [name for name, _ in genomes(features)] == ["mm10", "hg19"]
    assert [name for name, _ in genomes(pd.DataFrame(index=range(2)))] == [""]


def test_ordmag_with_no_counts_calls_nothing():
    chosen, info = ordmag_v2(np.zeros(5, dtype=np.int64), 3000, np.random.RandomState(0))
    assert chosen.size == 0 and info["n"] == 0


def test_cellranger_finds_the_simulated_cells(adata):
    result = scd.call(adata, "cellranger@3.0", expect_cells=300)
    truth = (adata.obs["truth"] == "cell").to_numpy()
    assert not (result.is_cell & ~truth).any()
    assert (result.is_cell & truth).sum() > 0.9 * truth.sum()
    assert any("no barcodes to test" in note for note in result.notes), "7,300 barcodes leave no ambient ranks"


def test_cellranger_is_reproducible_and_seeded(adata):
    a = scd.call(adata, "cellranger-ordmag@2.1", inplace=False)
    b = scd.call(adata, "cellranger-ordmag@2.1", inplace=False)
    c = scd.call(adata, "cellranger-ordmag@2.1", random_state=7, inplace=False)
    assert np.array_equal(a.is_cell, b.is_cell)
    assert len(a.diagnostics["groups"]["sim-1"]["bootstrap"]) == 100
    assert c.counting["random_state"] == 7


def test_force_cells_takes_the_top_n_per_gem_group():
    barcodes = [f"B{i:03d}-{1 + i % 2}" for i in range(40)]
    adata = _totals(np.arange(40, 0, -1), barcodes=barcodes)
    result = scd.call(adata, "cellranger@3.0", force_cells=10, inplace=False)
    assert result.n_cells == 10
    groups = gem_groups(np.asarray(result.barcodes)[result.is_cell])
    assert sorted(groups.tolist()) == [1] * 5 + [2] * 5


def test_each_genome_gets_its_own_ordmag_and_the_cells_are_their_union():
    values = [800] * 50 + [5] * 400
    human = _totals(values, genome="hg")
    # 500 expected cells: every bootstrap's 5th largest is 800, so the cut is 80.
    result = scd.call(human, "cellranger-ordmag@2.1", expect_cells=500, inplace=False)
    # n is the bootstraps' mean count above the cut, so close to, not exactly, 50
    assert 45 <= result.n_cells <= 50
    assert (result.totals[result.is_cell] == 800).all()
    assert set(result.diagnostics["groups"]) == {"hg-1"}


def test_releases_without_public_source_are_refused(adata):
    with pytest.raises(scd.UnsupportedVersion, match="cannot be reproduced"):
        scd.call(adata, "cellranger@5.0")
    with pytest.raises(scd.UnsupportedVersion, match="cannot be reproduced"):
        scd.call(adata, "cellranger@1.3")


# --- Cell Ranger 7.2+ ---------------------------------------------------------------

from sc_droplets.methods import cellranger as cr  # noqa: E402
from sc_droplets.primitives import simulate as sim  # noqa: E402


@pytest.mark.parametrize(
    "release, chemistry, tags, expected",
    [
        ("7.2", "SC3Pv3", None, ((45000, 90000), 0.01, "D7")),
        ("7.2", "SC3Pv3LT", None, ((4500, 9000), 0.01, "LT")),
        ("7.2", "Fixed RNA Profiling", 2, ((45000, 90000), 0.01, "D7")),
        ("7.2", "MFRP", 3, ((67500, 135000), 0.01, "D7")),
        ("8.0.0", "SC3Pv4", None, ((80000, 160000), 0.001, "V8")),
        ("8.0.0", "SC3Pv3HT", None, ((45000, 90000), 0.01, "D7")),
        ("8.0.1", "SC3Pv3HT", None, ((80000, 160000), 0.01, "HT8")),
        ("9.0", "Single Cell 3' v3 (polyA)", None, ((45000, 90000), 0.001, "D9")),
        ("9.0", "SC3Pv3", None, ((45000, 90000), 0.001, "D9")),
        ("9.0", "SC3Pv4-polyA-OCM", 1, ((20000, 40000), 0.001, "R9")),
        ("9.0", "MFRP-RNA", 2, ((90000, 180000), 0.01, "D9")),
        ("10.0", "Flex-v2-R1", None, ((80000, 160000), 0.01, "F2")),
        ("10.0", "SC3Pv3LT", None, ((4500, 9000), 0.01, "LT")),
        ("10.1", "SC3Pv3LT", None, ((45000, 90000), 0.01, "D9")),
        ("10.1", "custom: my-chemistry", None, ((45000, 90000), 0.01, "D9")),
        ("10.1", None, None, ((45000, 90000), 0.01, "D9")),
    ],
)
def test_the_ambient_range_and_fdr_follow_each_releases_chemistry_tables(release, chemistry, tags, expected):
    assert cr.ambient_range(release, cr.chemistry_family(chemistry), tags) == expected


@pytest.mark.parametrize(
    "values, n0, rule, expected",
    [
        ([9, 8, 7, 6, 5, 4, 3], 5, "revert", 6),  # no ties: one past, as written
        ([9, 8, 7, 6, 5, 4, 3], 5, "truncate", 6),
        ([9, 8, 7, 7, 7, 7, 7, 3], 3, "revert", 3),  # run too long: 7.2 reverts
        ([9, 8, 7, 7, 7, 7, 7, 3], 3, "truncate", 3),  # 8.0+ keeps what passed (none here)
        ([9] * 10 + [5] * 3 + [1], 10, "revert", 11),  # the tie run ends at n0: one past it
        ([9, 8] + [5] * 10 + [4, 1], 10, "truncate", 12),  # truncated at 20 %
        ([9, 8] + [5] * 10 + [4, 1], 10, "revert", 10),
        ([5, 5, 5], 2, "revert", 2),
        ([9, 8, 7, 6, 5], 5, "revert", 5),  # already at the end
    ],
)
def test_the_tie_extension_takes_one_past_the_run_and_its_20_percent_rule(values, n0, rule, expected):
    assert cr._extend_ties(np.array(values), n0, rule) == expected


def test_expected_cells_are_divided_by_gem_groups_twice_as_written():
    assert cr._expected_v7(None, 1) is None and cr._expected_v7(1000, 1) == 1000
    assert cr._expected_v7(1000, 2) == 250


def test_ordmag_7_2_estimates_the_expected_cells_when_not_given(adata):
    result = scd.call(adata, "cellranger-ordmag@8.0", inplace=False)
    group = next(iter(result.diagnostics["groups"].values()))
    assert 250 <= group["estimated_cells"] <= 400
    truth = (adata.obs["truth"] == "cell").to_numpy()
    assert (result.is_cell & truth).sum() > 0.9 * truth.sum()
    given = scd.call(adata, "cellranger-ordmag@7.2", expect_cells=300, inplace=False)
    assert "estimated_cells" not in next(iter(given.diagnostics["groups"].values()))


def test_cellranger_7_2_without_enough_barcodes_for_the_ambient_range_keeps_ordmag(adata):
    result = scd.call(adata, "cellranger@9.0", chemistry="SC3Pv4", inplace=False)
    assert any("no barcodes to test" in note for note in result.notes)
    assert result.diagnostics["range_rule"] == "R9"
    assert result.thresholds.get("fdr") == 0.001


def test_cellranger_7_2_options_are_refused_by_older_releases(adata):
    with pytest.raises(ValueError, match="before 7.2 takes no chemistry"):
        scd.call(adata, "cellranger@3.0", chemistry="SC3Pv3")
    with pytest.raises(ValueError, match="before 7.2 takes no chemistry"):
        scd.call(adata, "cellranger-ordmag@2.1", chemistry="SC3Pv3")


def test_the_emptydrops_step_aborts_on_a_flat_profile_before_9_0_and_clamps_after():
    """A 200-gene simulation has a flat Good-Turing slope (-0.6): Cell Ranger 7.2
    and 8.x give up on the EmptyDrops step, 9.0+ clamps the slope to -1 and goes on."""
    from sc_droplets.counts import Droplets

    adata = scd.simulate_droplets(n_cells=300, n_empty=3000, n_low=3000, n_genes=200, cell_umis=1500, cell_sigma=1.0, seed=7)
    droplets = Droplets.from_anndata(adata)
    truth = (adata.obs["truth"] == "cell").to_numpy()
    originals = truth & (droplets.totals > 1500)
    ids = droplets.features["gene_ids"].to_numpy()

    def step(release):
        settings = dict(cr._v7_settings(release), num_sims=300)
        return cr.nonambient_v7(droplets.counts, ids, originals, settings, (2000, 5000), 0.01, 150, np.random.default_rng(42))

    assert "slope" in step("7.2")["stopped"]
    for release in ("9.0", "10.1"):
        found = step(release)
        rescued = found["candidates"][found["nonambient"]]
        assert rescued.size > 0 and truth[rescued].all(), release


def test_antibody_aggregates_need_five_signal_antibodies_and_many_top_k_hits():
    rng = np.random.default_rng(0)
    counts = rng.poisson(2, size=(2000, 6))
    counts[:, :] += (rng.random((2000, 6)) < 0.02) * 100  # a few bright barcodes per antibody, not shared
    counts[7] = 5000  # one barcode bright in every antibody: an aggregate
    aggregates, info = cr.antibody_aggregates(sp.csr_matrix(counts), None)
    assert info["n_signal_antibodies"] == 6 and aggregates[7] and aggregates.sum() == 1
    few, info = cr.antibody_aggregates(sp.csr_matrix(counts[:, :4]), None)
    assert not few.any() and info["n_signal_antibodies"] == 4


def test_an_aggregate_list_from_cell_ranger_is_removed_before_calling(adata, tmp_path):
    top = adata.obs_names[np.argsort(-np.asarray(adata.X.sum(axis=1)).ravel())[:3]]
    (tmp_path / "aggregate_barcodes.csv").write_text("barcode,library_type\n" + "".join(f"{b},Antibody Capture\n" for b in top))
    result = scd.call(adata, "cellranger@8.0", aggregate_barcodes=str(tmp_path / "aggregate_barcodes.csv"), inplace=False)
    assert result.diagnostics["aggregates_removed"] == 3
    removed = np.isin(np.asarray(result.barcodes), list(top))
    assert not result.is_cell[removed].any()
    assert (np.asarray(result.classes)[removed] == "removed").all()


def test_the_8_0_filters_remove_low_umi_and_high_mito_cells(adata):
    adata = adata.copy()
    adata.var["gene_ids"] = adata.var["gene_ids"].astype(str)
    adata.var.iloc[0, adata.var.columns.get_loc("gene_ids")] = "GRCh38_ENSG00000198888"  # a mito gene, barnyard-prefixed
    plain = scd.call(adata, "cellranger@8.0", inplace=False)
    strict = scd.call(adata, "cellranger@8.0", global_minimum_umis=2000, max_mito_percent=1.0, inplace=False)
    assert strict.n_cells < plain.n_cells
    assert strict.diagnostics["removed_global_minimum"] > 0
    kept_totals = strict.totals[strict.is_cell]
    assert (kept_totals >= 2000).all()


def test_force_cells_in_7_2_takes_the_top_n_with_ties_to_the_later_barcode():
    adata = _totals([50, 40, 40, 40, 10, 0])
    result = scd.call(adata, "cellranger@8.0", force_cells=3, inplace=False)
    assert result.is_cell.tolist() == [True, False, True, True, False, False]


def test_the_high_occupancy_filter_removes_crowded_gems():
    rng = np.random.default_rng(3)
    letters = np.array(list("ACGT"))
    gems = ["".join(letters[rng.integers(0, 4, 16)]) for _ in range(3000)]
    probes = ["ACTTTAGG", "AACGGGAA", "AGTAGGCT", "ATGTTGAC"]
    barcodes = [g + probes[i % 4] + "-1" for i, g in enumerate(gems)]
    crowded = gems[0]
    barcodes += [crowded + p + "-1" for p in probes[1:]]  # one GEM with all four samples
    removed, info = cr.high_occupancy_gems(barcodes, "7.2")
    assert set(removed) == {crowded + p + "-1" for p in probes}
    removed10, info10 = cr.high_occupancy_gems(barcodes, "10.0")
    assert info10["threshold"] >= 1
    assert cr.high_occupancy_gems([], "10.0") == ([], {"lambda": 0.0, "n_gems": 0, "threshold": 0, "n_removed": 0})


def test_batched_refills_follow_cell_ranger_7_2():
    """Gaps of 20 or more refill the block before they would reach its end,
    keeping its unused tail; single steps refill only once the block is used up."""
    profile = np.array([0.5, 0.3, 0.2])
    for batch in (None, 20):
        state = np.random.RandomState(0)
        totals = np.array([3, 25, 49, 50, 51])
        sim.cellranger_v1(profile, totals, np.full(5, -1e9), 2, state, block_size=30, batch_small=batch)
    # Reference: replay the documented schedule by hand.
    def schedule(batch):
        state = np.random.RandomState(0)
        draws = [("block", 30)]
        state.choice(3, size=30, p=profile)
        position = 0
        for _ in range(2):
            state.multinomial(3, profile)
            for step in (22, 24, 1, 1):
                if batch and step >= batch:
                    if position + step >= 30:
                        draws.append(("tail", position))
                        state.choice(3, size=position, p=profile)
                        position = 0
                    position += step
                else:
                    for _ in range(step):
                        position += 1
                        if position >= 30:
                            draws.append(("block", 30))
                            state.choice(3, size=30, p=profile)
                            position = 0
        return draws, state.random_sample()

    batched, after_batched = schedule(20)
    single, after_single = schedule(None)
    assert batched != single, "the two schedules consume the stream differently"
    for batch, after in ((20, after_batched), (None, after_single)):
        state = np.random.RandomState(0)
        sim.cellranger_v1(profile, np.array([3, 25, 49, 50, 51]), np.full(5, -1e9), 2, state, block_size=30, batch_small=batch)
        assert state.random_sample() == after
