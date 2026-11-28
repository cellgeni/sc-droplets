"""R's ks estimators (`primitives/ks.py`), as DropletQC calls them; checked against ks itself in
tests/test_methods_dropletqc.py."""

from __future__ import annotations

import math

import numpy as np
import pytest

from sc_droplets.primitives import ks




def test_r_dnorm_is_the_normal_density_on_both_of_its_branches():
    from scipy.stats import norm

    for x in (0.0, 0.3, -1.7, 4.99, 5.0, 7.3, -11.2, 30.0):
        assert ks.r_dnorm(x, 0.7) == pytest.approx(norm.pdf(x, scale=0.7), rel=1e-14, abs=1e-300)
    assert ks.r_dnorm(50.0, 1.0) == 0.0


def test_dnorm_deriv_is_the_derivative_of_the_density():
    x = np.linspace(-3, 3, 13)
    step = 1e-5
    for r in (1, 2, 3):
        numeric = (ks.dnorm_deriv(x + step, 0.8, r - 1) - ks.dnorm_deriv(x - step, 0.8, r - 1)) / (2 * step)
        assert np.allclose(ks.dnorm_deriv(x, 0.8, r), numeric, rtol=1e-6, atol=1e-8)


def test_psins_is_the_normal_scale_functional():
    # psi_r = (-1)^(r/2) r! / ((2 sigma)^(r+1) (r/2)! sqrt(pi)); psi_4 of N(0, 1) is 3 / (8 sqrt(pi)).
    assert ks.psins_1d(4, 1.0) == pytest.approx(3 / (8 * math.sqrt(math.pi)))
    assert ks.psins_1d(3, 1.0) == 0.0


def test_linear_binning_keeps_the_mass():
    rng = np.random.default_rng(0)
    x = rng.normal(size=500)
    binned = ks.binning(x, h=0.3, supp=4)
    assert binned.counts.sum() == pytest.approx(500)
    assert binned.eval_points[0] == x.min() - 4 * 0.3 and binned.eval_points[-1] == x.max() + 4 * 0.3


@pytest.mark.parametrize("legacy", [True, False])
def test_the_binned_derivative_is_close_to_the_exact_kernel_estimate(legacy):
    """Binning and the FFT convolution approximate the kernel estimate at the same bandwidth,
    (1/n) sum_i phi_h'(t - x_i), to well within its smoothing; ks 1.14.3 places it a cell later."""
    rng = np.random.default_rng(1)
    x = rng.normal(size=5000)
    estimate = ks.kdde_deriv(x, 1, legacy=legacy)
    points = estimate.eval_points[150:250]
    exact = np.array([ks.dnorm_deriv(t - x, estimate.h, 1).mean() for t in points])
    step = estimate.eval_points[1] - estimate.eval_points[0]
    if not legacy and round((2 * math.ceil(5 * estimate.h / step) - 1) / 2) + 1 != math.ceil(5 * estimate.h / step):
        exact = np.array([ks.dnorm_deriv(t - step - x, estimate.h, 1).mean() for t in points])
    assert np.max(np.abs(estimate.estimate[150:250] - exact)) < 1e-3 * np.max(np.abs(exact))
