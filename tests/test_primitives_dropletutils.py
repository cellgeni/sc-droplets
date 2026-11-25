"""DropletUtils' arithmetic, piece by piece, against R and DropletUtils themselves.

The fixtures in tests/oracles/dropletutils/fixtures come from
oracles/dropletutils/fixtures.R (DropletUtils 1.30, R 4.5).
"""

from __future__ import annotations

import gzip
import math
from pathlib import Path

import numpy as np
import pytest

from sc_droplets.methods.dropletutils import _ld_cumsum, r_pow10
from sc_droplets.primitives import boostrandom as br
from sc_droplets.primitives.goodturing import good_turing_proportions_edger
from sc_droplets.primitives.rrandom import RMersenneTwister, generate_seed_vectors
from sc_droplets.primitives.rstats import optimize_max, brent_fmin, p_adjust_bh, quantile7, seq_length_out, all_equal
from sc_droplets.primitives.simulate import dropletutils as simulate_dropletutils
from sc_droplets.primitives.smoothspline import nknots_smspl, smooth_spline

FIXTURES = Path(__file__).parent / "oracles" / "dropletutils" / "fixtures"


def _read(name: str) -> list:
    path = FIXTURES / name
    opener = gzip.open if name.endswith(".gz") else open
    with opener(path, "rt") as handle:
        return handle.read().splitlines()


# ---------------------------------------------------------------------------
# R's random numbers and dqrng
# ---------------------------------------------------------------------------


def test_set_seed_100_gives_rs_first_uniforms():
    u = RMersenneTwister(100).unif_rand(3)
    assert u.tolist() == [0.30776611086912453, 0.25767250102944672, 0.55232243333011866]


def test_seed_vectors_are_dqrngs_joined_into_64_bits():
    # R: set.seed(100); dqrng::generateSeedVectors(2, nwords = 2)
    seeds = generate_seed_vectors(RMersenneTwister(100), 2)
    words = [(1106677093, 2012380120), (1590529003, -1610620459 & 0xFFFFFFFF)]
    assert seeds.tolist() == [(a << 32) | b for a, b in words]


def test_pcg32_matches_the_reference_implementations_demo():
    # pcg32-demo: pcg32_srandom_r(&rng, 42u, 54u), the first six outputs.
    state, inc = br.pcg32_start(np.array([42]), np.array([54]))
    expected = [0xA15C02B7, 0x7B47F409, 0xBA1D3330, 0x83D2F293, 0xBFA4784B, 0xCBED606E]
    assert br.Pcg32Jumps().raw(state[0], inc[0], 6).tolist() == expected
    engine = br.Lockstep(state, inc)
    assert [int(engine.draw(np.array([0]))[0]) for _ in range(6)] == expected


def test_jumping_ahead_lands_where_stepping_does():
    state, inc = br.pcg32_start(np.array([2**63 + 12345]), np.array([7]))
    jumps = br.Pcg32Jumps(16)
    raw = jumps.raw(state[0], inc[0], 5000)
    after = jumps.advance(state[0], inc[0], 4000)
    assert jumps.raw(after, inc[0], 1000).tolist() == raw[4000:].tolist()


def test_int_bucket_follows_boost_for_a_32_bit_engine():
    assert br.int_bucket(2) == 2**31  # (2**32 - 1) % 2 == 1 == n - 1, so one more
    assert br.int_bucket(3) == (2**32 - 1) // 3
    assert br.int_bucket(36601) == (2**32 - 1) // 36601


@pytest.mark.parametrize("seed", range(6))
def test_the_fast_alias_table_pairs_columns_as_boosts_loop_does(seed):
    rng = np.random.default_rng(seed)
    if seed % 2:
        weights = rng.gamma(rng.uniform(0.01, 3, 400), 1.0)  # Dirichlet draws, as 1.2-1.26 build tables from
    else:
        weights = rng.exponential(1.0, 300) ** 3
    exact_threshold, exact_alias = br.alias_table(weights)
    fast_threshold, fast_alias = br.alias_table_fast(weights)
    assert np.array_equal(exact_alias, fast_alias)
    np.testing.assert_allclose(fast_threshold, exact_threshold, rtol=0, atol=1e-12)


def test_an_alias_table_of_one_weight_never_draws_a_column():
    state, inc = br.pcg32_start(np.array([1]), np.array([1]))
    jumps = br.Pcg32Jumps()
    threshold, alias = br.alias_table(np.array([3.0]))
    draws, after = br.sample_alias(jumps, state[0], inc[0], 10, threshold, alias)
    assert draws.tolist() == [0] * 10
    assert after == jumps.advance(state[0], inc[0], 10)


@pytest.mark.parametrize("alpha", ["inf", "0.5", "50", "2000"])
def test_monte_carlo_counts_are_dropletutils_own(alpha):
    """Each observed value sits between two simulated ones, so one different path changes a count."""
    ambient = np.array([float(v) for v in _read("mc_ambient.txt")])
    grid = [line.split("\t") for line in _read(f"mc_grid_{alpha}.tsv.gz")]
    totals = np.array([int(t) for t, _ in grid])
    observed = np.array([float(v) for _, v in grid])
    expected = np.array([int(v) for v in _read(f"mc_counts_{alpha}.txt")])
    got = simulate_dropletutils(ambient, totals, observed, 200, 100, alpha=float(alpha), chunk=64)
    assert np.array_equal(got, expected)


def test_the_inverse_sampler_and_rs_stream_draw_from_the_profile():
    ambient = np.array([0.7, 0.2, 0.1])
    for sampler in ("inverse", "r"):
        counts = simulate_dropletutils(ambient, np.array([2000]), np.array([0.0]), 3, 1, sampler=sampler)
        assert counts.shape == (1,)


# ---------------------------------------------------------------------------
# R's smooth.spline and statistics
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n", [30, 120, 700, 2500])
def test_smooth_spline_is_rs_to_the_last_bit(n):
    rows = np.array([[float(v) for v in line.split("\t")] for line in _read(f"spline_{n}.tsv.gz")])
    spar, df = (float(v) for v in _read(f"spline_{n}_spar.txt"))
    fit = smooth_spline(rows[:, 0], rows[:, 1], df=20)
    assert fit.spar == spar
    assert np.array_equal(fit.fitted, rows[:, 2])
    assert np.array_equal(fit.predict(deriv=1), rows[:, 3])
    assert np.array_equal(fit.predict(deriv=2), rows[:, 4])


def test_nknots_follows_rs_rule():
    expected = [int(v) for v in _read("nknots.txt")]
    assert [nknots_smspl(n) for n in (10, 49, 50, 199, 200, 799, 800, 3199, 3200, 10000)] == expected


def test_smooth_spline_needs_four_distinct_x():
    with pytest.raises(ValueError, match="four unique"):
        smooth_spline(np.array([1.0, 2.0, 3.0, 3.0]), np.array([1.0, 2.0, 3.0, 4.0]), df=20)


@pytest.mark.parametrize("case", [1, 2, 3, 4])
def test_edger_good_turing_proportions_match_r(case):
    rows = [line.split("\t") for line in _read(f"goodturing_{case}.tsv.gz")]
    counts = np.array([float(c) for c, _ in rows])
    expected = np.array([float(p) for _, p in rows])
    assert np.array_equal(good_turing_proportions_edger(counts), expected)


def test_optimize_finds_rs_optimum():
    maximum, minimum = (float(v) for v in _read("optimize.txt"))
    assert optimize_max(lambda a: -((a - 53.2) * (a - 53.2)) / (1 + a / 1000), 0.01, 10000.0) == maximum
    assert brent_fmin(lambda a: (a - 2.5) ** 2, 0.0, 10.0, np.finfo(float).eps ** 0.25) == minimum


def test_quantile_is_rs_type_7():
    expected = [float(v) for v in _read("quantile.txt")]
    q = np.concatenate(((np.arange(1, 3001) * 7919) % 10007, np.full(40, 5000))).astype(np.float64)
    assert quantile7(q, 0.99) == expected[0]
    assert quantile7(np.array([1.0, 7.0] * 50), 0.5) == expected[1]
    assert quantile7(np.arange(1, 18) ** 2 / 3, 0.37) == expected[2]
    # equal order statistics: R returns the value, not (1 - h)·9 + h·9
    assert quantile7(np.array([3.0] * 2969 + [9.0] * 31), 0.99) == expected[3]


def test_p_adjust_bh_leaves_na_out_of_n():
    p = np.array([0.01, np.nan, 0.04, 0.04, 0.2, 0, np.nan, 0.9, 0.03, 1])
    expected = np.array([float(v) if v != "NA" else np.nan for v in _read("padjust.txt")])
    np.testing.assert_array_equal(p_adjust_bh(p), expected)


def test_seq_length_out_steps_forward_from_the_start():
    assert seq_length_out(1, 10, 4).tolist() == [1.0, 4.0, 7.0, 10.0]
    assert seq_length_out(1, 2500, 184)[-1] == 2500.0


def test_all_equal_uses_the_mean_relative_difference():
    assert all_equal(np.array([1.0, 2.0]), np.array([1.0, 2.0 + 1e-9]))
    assert not all_equal(np.array([1.0, 2.0]), np.array([1.0, 2.5]))


def test_cumsum_accumulates_in_long_double_like_r():
    """R's cumsum keeps 64 bits of mantissa between steps; a double running sum loses the 1s."""
    out = _ld_cumsum(np.array([1e16, 1.0, 1.0, -1e16]))
    assert out.tolist() == [1e16, 1e16, 1.0000000000000002e16, 2.0]


def test_ten_to_a_log_is_not_always_the_total_again():
    # R on x86-64: 10^log10(1852) is 1851.9999999999993, so a barcode with 1852 UMIs is below that knee.
    assert r_pow10(math.log10(1852)) == 1851.9999999999993
    assert r_pow10(2.0) == 100.0


# ---------------------------------------------------------------------------
# R's lgamma and long-double sums
# ---------------------------------------------------------------------------


def test_lgammafn_is_rs_own_log_gamma():
    from sc_droplets.primitives.rmath import lgammafn

    rows = [line.split("\t") for line in _read("lgamma.txt")]
    x = np.array([float(a) for a, _ in rows])
    expected = np.array([float(b) for _, b in rows])
    assert np.array_equal(lgammafn(x), expected)


def test_sums_accumulate_in_long_double_as_r_does_on_x86_64():
    from sc_droplets.primitives.rmath import _ld_add, _ld_to_float, ld_cumsum, ld_sum

    values = np.array([float(v) for v in _read("sum_values.txt.gz")])
    total, head, head_cumsum = (float(v) for v in _read("sum.txt"))
    assert head == 10.0  # double arithmetic would give 0
    assert ld_sum(values[:12]) == head
    assert ld_cumsum(values[:12])[-1] == head_cumsum
    assert ld_sum(values) == total
    # the vectorised path agrees with adding one value at a time
    num, exp = 0, 0
    for v in values[:30000].tolist():
        num, exp = _ld_add(num, exp, v)
    assert ld_sum(values[:30000]) == _ld_to_float(num, exp)
