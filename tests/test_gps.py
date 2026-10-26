import numpy as np
import pytest
from scipy import stats

from pvsignals.gps import GPSPrior, fit_prior, posterior


def simulate(m: int, seed: int = 0):
    """Pairs whose true ratios come from a known two-gamma mixture (80% near 1, 20% elevated)."""
    rng = np.random.default_rng(seed)
    e = np.exp(rng.uniform(np.log(0.1), np.log(100), m))
    bulk = rng.random(m) < 0.8
    lam = np.where(bulk, rng.gamma(20, 1 / 20, m), rng.gamma(2, 1 / 0.5, m))
    return rng.poisson(lam * e), e, lam


def test_truncated_fit_recovers_the_mixture():
    n, e, _ = simulate(6000)
    keep = n >= 2  # count tables never list small counts: condition on n >= 2
    prior = fit_prior(n[keep], e[keep], np.full(keep.sum(), 2.0))
    assert prior.w == pytest.approx(0.8, abs=0.08)
    assert prior.a1 / prior.b1 == pytest.approx(1.0, abs=0.15)  # mean of the bulk component
    assert prior.a2 / prior.b2 == pytest.approx(4.0, rel=0.3)  # mean of the elevated component


def test_posterior_matches_a_single_gamma_conjugate_update():
    prior = GPSPrior(a1=2.0, b1=1.5, a2=2.0, b2=1.5, w=0.5)  # both components equal = one gamma prior
    n, e = np.array([7.0]), np.array([3.0])
    ebgm, lo, hi, _ = posterior(n, e, prior)
    post = stats.gamma(2.0 + 7, scale=1 / (1.5 + 3))
    assert lo[0] == pytest.approx(post.ppf(0.05), rel=1e-6)
    assert hi[0] == pytest.approx(post.ppf(0.95), rel=1e-6)
    from scipy.special import digamma
    assert ebgm[0] == pytest.approx(np.exp(digamma(9.0) - np.log(4.5)))


def test_small_counts_are_shrunk_more_than_large_ones():
    prior = GPSPrior(a1=20.0, b1=20.0, a2=2.0, b2=0.5, w=0.8)
    # the same observed-to-expected ratio of 5, from 5 reports and from 500 reports
    ebgm, lo, hi, _ = posterior(np.array([5.0, 500.0]), np.array([1.0, 100.0]), prior)
    assert ebgm[0] < ebgm[1] < 5.0
    assert ebgm[1] == pytest.approx(5.0, rel=0.05)
    assert np.all(lo < ebgm) and np.all(ebgm < hi)
