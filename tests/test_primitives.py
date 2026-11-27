"""The shared building blocks, each against an independent reference.

References: Gale & Sampson's published prosody example (Good-Turing), STAR's
own `sgt.h` and GCC's libstdc++ compiled into a harness (their outputs are
embedded below), and the C++ standard's check value for mt19937.
"""

from __future__ import annotations

import itertools
import math

import numpy as np
import pytest
import scipy.sparse as sp
from scipy.stats import multinomial

from sc_droplets.primitives import goodturing, multinomial as mn, pvalues, rounding, simulate, stdrandom
from sc_droplets.primitives.otsu import (
    OtsuError,
    _between_class_table,
    histogram,
    multiotsu_indices,
    threshold_multiotsu,
    threshold_otsu,
)
from sc_droplets.primitives.rmath import clog2
from sc_droplets.primitives.rstats import r_mad, r_median

# Gale & Sampson (1995), Table 2: frequencies r and their counts n_r (prosody data).
PROSODY = [(1, 120), (2, 40), (3, 24), (4, 13), (5, 15), (6, 5), (7, 11), (8, 2), (9, 2), (10, 1), (12, 3), (14, 2), (15, 1), (16, 1), (17, 3), (19, 1), (20, 3), (21, 2), (23, 3), (24, 3), (25, 3), (26, 2), (27, 2), (28, 1), (31, 2), (32, 2), (33, 1), (34, 2), (36, 2), (41, 3), (43, 1), (45, 3), (46, 1), (47, 1), (50, 1), (71, 1), (84, 1), (101, 1), (105, 1), (121, 1), (124, 1), (146, 1), (162, 1), (193, 1), (199, 1), (224, 1), (226, 1), (254, 1), (257, 1), (339, 1), (421, 1), (456, 1), (481, 1), (483, 1), (1140, 1), (1256, 1), (1322, 1), (1530, 1), (2131, 1), (2395, 1), (6925, 1), (7846, 1)]  # noqa: E501
# Their published r* for the same rows.
PROSODY_RSTAR = [0.7628079, 1.706448, 2.679796, 3.663988, 4.653366, 5.645628, 6.63966, 7.634856, 8.63086, 9.627446, 11.62182, 13.61725, 14.61524, 15.61336, 16.6116, 18.60836, 19.60685, 20.6054, 22.60264, 23.60133, 24.60005, 25.5988, 26.59759, 27.59639, 30.59294, 31.59183, 32.59073, 33.58964, 35.58751, 40.58235, 42.58035, 44.57836, 45.57738, 46.57641, 49.57351, 70.55399, 83.54229, 100.5272, 104.5237, 120.5097, 123.507, 145.4879, 161.474, 192.4472, 198.4421, 223.4205, 225.4188, 253.3947, 256.3922, 338.3218, 420.2514, 455.2215, 480.2, 482.1983, 1138.636, 1254.537, 1320.48, 1528.302, 2128.788, 2392.562, 6918.687, 7838.899]  # noqa: E501
# STAR 2.7.11b sgt.h on PROSODY, compiled with g++: P0, then each row's estimate.
STAR_SGT = [0.0038832438029900978, 2.4684741251171716e-05, 5.5221274647292958e-05, 8.671917360082353e-05, 0.00011856800528718321, 0.00015058462833122393, 0.00018269457960756893, 0.00021486181093836214, 0.00024706673184407091, 0.00027929777516933813, 0.00031154766955937205, 0.0003760863107763635, 0.0004406592327428192, 0.00047295447666784515, 0.00050525418053383084, 0.00053755759207718733, 0.00060217328962278627, 0.00063448472179798807, 0.00066679810433153996, 0.00073142973663196421, 0.00076374759370774589, 0.00079606660036645436, 0.0008283866278748351, 0.00086070756603016398, 0.00089302931994167519, 0.00098999870640653501, 0.0010223230004906168, 0.0010546477909955426, 0.0010869730351159518, 0.0011516247362516514, 0.0013132597099410948, 0.0013779155174666801, 0.0014425721596966995, 0.0014749007603790162, 0.001507229531840712, 0.0016042167700342593, 0.0022831527263585173, 0.0027034589141850416, 0.0032530972801729182, 0.0033824247589868527, 0.0038997368078661108, 0.0039967331270382806, 0.0047080416234606881, 0.0052253586765531853, 0.0062276632183242073, 0.0064216579537809109, 0.0072299701302320366, 0.0072946351507120217, 0.0081999460029666769, 0.0082969436488775136, 0.010948215368236683, 0.013599490270747545, 0.014731132544454321, 0.015539448595563665, 0.015604113884087879, 0.036846675404433898, 0.040597265783653982, 0.042731222434119434, 0.049456419346593429, 0.068888359487410997, 0.077424187026014488, 0.22389123311466458, 0.25366963288845695]  # noqa: E501
# g++ 15 (libstdc++): std::discrete_distribution<uint32_t> over WEIGHTS below, drawn
# with std::mt19937(19760110 * (s + 1)), first 60 draws for s = 0, 1, 2.
DISCRETE = [
    [1451, 1210, 0, 785, 4973, 4985, 562, 4855, 3848, 3988, 2621, 857, 1714, 4177, 4285, 962, 4095, 664, 4103, 2361, 0, 3344, 2991, 1714, 3564, 1959, 2587, 2676, 4941, 4208, 857, 4565, 1573, 4904, 3829, 52, 2291, 4062, 4296, 3980, 682, 612, 2420, 2901, 3918, 3620, 2347, 2933, 997, 892, 4365, 1206, 997, 2396, 4371, 4915, 742, 1227, 2991, 2086],  # noqa: E501
    [3095, 4985, 1994, 4005, 2816, 454, 2308, 927, 472, 962, 418, 3442, 4017, 3647, 1714, 4495, 2719, 4137, 2956, 4475, 1119, 2464, 2374, 4595, 76, 1538, 3528, 262, 2074, 2165, 2628, 2991, 962, 3745, 3690, 2886, 2886, 646, 295, 927, 2200, 4855, 2104, 2536, 1852, 2956, 3953, 3609, 1313, 4232, 2991, 0, 4610, 507, 3784, 0, 831, 1158, 3274, 4810],  # noqa: E501
    [2215, 52, 927, 1924, 0, 1732, 3638, 2991, 3485, 3481, 4792, 1713, 3638, 2571, 3848, 2886, 927, 749, 1714, 52, 3883, 4985, 3480, 1038, 4985, 3841, 752, 1311, 1994, 997, 2401, 281, 4915, 1994, 4208, 2921, 2486, 0, 716, 3320, 4985, 3918, 4302, 4950, 681, 997, 2991, 4705, 3931, 2991, 892, 4058, 4506, 1749, 997, 2972, 4355, 892, 4418, 4250],  # noqa: E501
]
WEIGHTS = np.array([1.0 / (1 + (i * 7919) % 997) + 1e-6 * i for i in range(5000)])


# --- rounding ---------------------------------------------------------------


@pytest.mark.parametrize(
    "x, away, even",
    [(2.5, 3, 2), (3.5, 4, 4), (-2.5, -3, -2), (0.49999999999999994, 0, 0), (29.5, 30, 30), (30.000000000000025, 30, 30)],
)
def test_rounding_halves_away_from_zero_or_to_even(x, away, even):
    assert rounding.round_half_away(x) == away
    assert rounding.round_half_even(x) == even


# --- Good-Turing ------------------------------------------------------------


def test_cellranger_good_turing_reproduces_gale_and_sampson():
    r, nr = np.array(PROSODY).T
    rstar, p0 = goodturing.simple_good_turing_cellranger(r, nr)
    assert p0 == pytest.approx(0.003883244, abs=1e-9)
    assert np.allclose(rstar, PROSODY_RSTAR, rtol=1e-4)


def test_cellranger_proportions_sum_to_one_with_p0():
    r, nr = np.array(PROSODY).T
    frequencies = np.repeat(r, nr)
    pstar, p0 = goodturing.sgt_cellranger(frequencies)
    assert p0 + pstar.sum() == pytest.approx(1.0)
    assert pstar[0] == pytest.approx(PROSODY_RSTAR[0] / np.sum(nr * np.array(PROSODY_RSTAR)) * (1 - p0))


def test_cellranger_good_turing_needs_ten_distinct_frequencies_and_a_falling_slope():
    with pytest.raises(goodturing.GoodTuringError, match="needs 10"):
        goodturing.sgt_cellranger(np.array([1, 1, 2, 3, 4, 5, 6, 7, 8]))
    with pytest.raises(goodturing.GoodTuringError, match="slope"):
        goodturing.simple_good_turing_cellranger(np.arange(1, 12), np.arange(1, 12))
    with pytest.raises(ValueError, match="greater than zero"):
        goodturing.sgt_cellranger(np.array([0, 1]))


def test_star_good_turing_matches_stars_own_sgt_h():
    sgt = goodturing.SimpleGoodTuringStar()
    for observation, count in PROSODY:
        sgt.add(observation, count)
    assert sgt.analyse()
    ours = [sgt.estimate(0)] + [sgt.estimate(o) for o, _ in PROSODY]
    assert np.allclose(ours, STAR_SGT, rtol=1e-12, atol=0)
    assert sgt.estimate(11) == 0.0, "a frequency not seen has no estimate"


def test_star_good_turing_gives_up_below_five_frequencies():
    sgt = goodturing.SimpleGoodTuringStar()
    for observation in (1, 2, 3, 4):
        sgt.add(observation, 3)
    assert not sgt.analyse()


# --- multinomial ------------------------------------------------------------


def test_sparse_log_likelihoods_match_scipy():
    rng = np.random.default_rng(0)
    p = rng.dirichlet(np.ones(30))
    rows = rng.multinomial(50, p, size=6)
    ours = mn.multinomial_loglik_rows(sp.csr_matrix(rows), np.log(p))
    theirs = multinomial.logpmf(rows, 50, p)
    assert np.allclose(ours, theirs, rtol=1e-12)
    assert mn.multinomial_logpmf_dense(rows[0], 50, p) == pytest.approx(theirs[0], rel=1e-12)


def test_the_scipy_adjusted_profile_replaces_the_last_probability():
    p = np.array([0.2, 0.3, 0.5 + 1e-9])
    adjusted = mn.scipy_adjusted_profile(p)
    assert adjusted[-1] == pytest.approx(0.5) and p[-1] != adjusted[-1]


# --- p-values ---------------------------------------------------------------


def test_monte_carlo_p_values():
    assert pvalues.monte_carlo_pvalues(np.array([0, 9]), 9).tolist() == [0.1, 1.0]


def test_benjamini_hochberg_variants_agree_with_the_definition():
    p = np.array([0.01, 0.04, 0.03, 0.2, 0.001])
    expected = np.array([0.025, 0.05, 0.05, 0.2, 0.005])  # p·n/rank, then the running minimum from the top
    assert np.allclose(pvalues.bh_cellranger(p), expected)
    assert np.allclose(pvalues.bh_starsolo(p), expected)
    # STAR does not cap at 1; Cell Ranger does.
    assert pvalues.bh_starsolo(np.array([0.9, 0.95]))[0] == pytest.approx(0.95)
    assert pvalues.bh_cellranger(np.array([1.0, 1.0])).max() == 1.0
    assert pvalues.bh_cellranger(np.array([])).size == 0


def test_tied_p_values_get_one_adjusted_value_whatever_their_order():
    p = np.array([0.02, 0.01, 0.02, 0.02])
    adjusted = pvalues.bh_starsolo(p)
    assert adjusted[0] == adjusted[2] == adjusted[3]


# --- C++ random numbers -----------------------------------------------------


def test_mt19937_meets_the_standards_check_value():
    """[rand.predef]: the 10000th output of a default-constructed mt19937 is 4123659995."""
    assert stdrandom.mt19937_raw(5489, 10000)[-1] == 4123659995


def test_seeds_above_32_bits_wrap_as_std_mt19937_does():
    assert np.array_equal(stdrandom.mt19937_raw(19760110 * 300, 5), stdrandom.mt19937_raw((19760110 * 300) % 2**32, 5))


def test_discrete_draws_match_libstdcxx():
    distribution = stdrandom.LibstdcxxDiscrete(WEIGHTS)
    for s, expected in enumerate(DISCRETE):
        draws = distribution.sample_from_raw(stdrandom.mt19937_raw(19760110 * (s + 1), 2 * len(expected)))
        assert draws.tolist() == expected


def test_generate_canonical_never_returns_one():
    top = np.array([2**32 - 1, 2**32 - 1], dtype=np.uint64)
    assert stdrandom.generate_canonical(top)[0] < 1.0


# --- simulation helpers -----------------------------------------------------


def test_occurrence_rank_counts_earlier_equal_values():
    values = np.array([3, 1, 3, 3, 1, 2])
    assert simulate.occurrence_rank(values).tolist() == [1, 1, 2, 3, 2, 1]
    offset = np.array([0, 5, 0, 1])
    assert simulate.occurrence_rank(values, offset).tolist() == [2, 6, 3, 4, 7, 1]
    assert simulate.occurrence_rank(np.array([], dtype=int)).size == 0


def _profile(n=40, seed=1):
    return np.random.default_rng(seed).dirichlet(np.full(n, 0.5))


@pytest.mark.parametrize("engine", ["v1", "v1_batched", "v9", "starsolo"])
def test_simulated_p_values_are_uniform_for_ambient_droplets_and_small_for_cells(engine):
    """Droplets drawn from the profile itself get spread-out p-values; ones
    drawn from another profile get the smallest possible."""
    p = _profile()
    rng = np.random.default_rng(2)
    totals = rng.integers(40, 120, size=200)
    ambient = np.array([rng.multinomial(t, p) for t in totals])
    other = np.array([rng.multinomial(t, _profile(seed=9)) for t in totals])
    rows = sp.csr_matrix(np.vstack([ambient, other]))
    all_totals = np.concatenate([totals, totals])
    observed = mn.multinomial_loglik_rows(rows, np.log(p))
    num_sims = 199
    if engine == "v1":
        n_lower = simulate.cellranger_v1(p, all_totals, observed, num_sims, np.random.RandomState(0), block_size=5000)
    elif engine == "v1_batched":
        n_lower = simulate.cellranger_v1(p, all_totals, observed, num_sims, np.random.RandomState(0), block_size=5000, batch_small=20)
    elif engine == "v9":
        n_lower = simulate.cellranger_v9(p, all_totals, observed, num_sims, np.random.default_rng(42), chunk=50)
    else:
        n_lower = simulate.starsolo(p, all_totals, observed, num_sims, chunk=50)
    pv = pvalues.monte_carlo_pvalues(n_lower, num_sims)
    assert 0.3 < np.median(pv[:200]) < 0.7
    assert np.all(pv[200:] == 1.0 / (num_sims + 1))


def test_large_gaps_are_drawn_with_one_multinomial_in_cell_ranger_v1():
    p = _profile()
    totals = np.array([10, 2000, 2005])
    observed = np.full(3, -1e9)
    n_lower = simulate.cellranger_v1(p, totals, observed, 5, np.random.RandomState(0), block_size=50)
    assert n_lower.tolist() == [0, 0, 0]


def test_the_starsolo_engine_is_deterministic_and_seeded_per_simulation():
    p = _profile()
    totals = np.array([30, 60])
    observed = np.array([-80.0, -150.0])
    a = simulate.starsolo(p, totals, observed, 20, chunk=7)
    b = simulate.starsolo(p, totals, observed, 20, chunk=20)
    assert np.array_equal(a, b), "chunking must not change the streams"


# ---------------------------------------------------------------------------
# R's median and mad
# ---------------------------------------------------------------------------


def test_r_median_takes_the_mean_of_the_middle_two_and_mad_scales_by_1_4826():
    assert r_median(np.array([4.0, 1.0, 3.0, 2.0])) == 2.5
    assert r_median(np.array([5.0, 1.0, 3.0])) == 3.0
    assert math.isnan(r_median(np.array([])))
    assert r_mad(np.array([1.0, 2.0, 3.0, 4.0, 10.0]), 3.0) == 1.4826 * 1.0


def test_clog2_is_the_c_librarys_log2_with_c_semantics_at_zero():
    values = clog2(np.array([1.0, 3.0, 1e-300, 0.0, -1.0]))
    assert values[:3].tolist() == [math.log2(1.0), math.log2(3.0), math.log2(1e-300)]
    assert values[3] == -math.inf and math.isnan(values[4])


# ---------------------------------------------------------------------------
# Otsu and multi-Otsu (scikit-image)
# ---------------------------------------------------------------------------


def test_otsu_of_identical_values_is_that_value():
    assert threshold_otsu(np.full(5, 2.5)) == 2.5


def test_multiotsu_search_is_the_first_maximum_in_lexicographic_order():
    """The vectorised search must pick what scikit-image's recursive loop picks."""
    rng = np.random.default_rng(0)
    for _ in range(20):
        n = int(rng.integers(5, 24))
        counts = rng.integers(0, 6, n).astype(float)
        counts[[0, -1]] += 1
        prob = (counts / counts.sum()).astype(np.float32)
        table = _between_class_table(prob)
        for k in (1, 2, 3):
            best, chosen = np.float32(0), None
            for combo in itertools.combinations(range(n - 1), k):
                sigma = table[0, combo[0]] + table[combo[-1] + 1, n - 1]
                for a, b in zip(combo[:-1], combo[1:]):
                    sigma = np.float32(sigma + table[a + 1, b])
                if sigma > best:
                    best, chosen = sigma, list(combo)
            assert multiotsu_indices(prob, k).tolist() == chosen


def test_multiotsu_first_moment_starts_at_the_first_probability():
    """scikit-image starts its first moments at prob[0], not 0 * prob[0]; so must the table."""
    prob = np.array([0.5, 0.25, 0.25], dtype=np.float32)
    table = _between_class_table(prob)
    assert table[0, 0] == 0  # never set
    assert table[0, 1] == np.float32((0.5 + 0.25) ** 2 / 0.75)


def test_multiotsu_with_as_many_values_as_classes_cuts_at_the_first_ones():
    x = np.array([1.0] * 4 + [2.0] * 3 + [3.0] * 5)
    _, centers = histogram(x, 256)
    cuts = threshold_multiotsu(x, 3)
    assert cuts.tolist() == [centers[0], centers[np.searchsorted(centers, 2.0)]]
    with pytest.raises(OtsuError, match="only 2 different values"):
        threshold_multiotsu(np.array([1.0, 1.0, 2.0]), 3)
