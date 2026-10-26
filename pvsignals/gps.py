"""Gamma-Poisson shrinker (DuMouchel 1999): empirical Bayes geometric mean (EBGM) and its 5th percentile (EB05).

Each drug-event pair has an observed count n and an expected count E. The true ratio lambda = mu / E gets a prior
that is a mixture of two gamma distributions, ``w * Gamma(a1, b1) + (1 - w) * Gamma(a2, b2)`` (shape, rate). Given
lambda, n ~ Poisson(lambda * E), so marginally n follows a mixture of negative binomials. The five prior parameters
are fitted by maximum likelihood over many pairs; the posterior of lambda is again a two-gamma mixture, and

    EBGM = exp(E[log lambda | n]),   EB05 = 5th percentile of lambda | n.

Count tables from openFDA list only the terms above a cut-off (at most 999 terms per query; zero counts are never
listed), so each pair's likelihood is conditioned on n >= t, its table's threshold: the truncated likelihood
DuMouchel used with t = 1, generalised to a per-pair t.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np
from scipy import optimize, special, stats

DUMOUCHEL_START = (0.2, 0.1, 2.0, 4.0, 1 / 3)  # (a1, b1, a2, b2, w) starting values from DuMouchel (1999)


@dataclass(frozen=True)
class GPSPrior:
    """Fitted two-gamma mixture prior."""

    a1: float
    b1: float
    a2: float
    b2: float
    w: float
    loglik: float = float("nan")
    n_pairs: int = 0

    def as_dict(self) -> dict:
        return asdict(self)


def _log_nb(n: np.ndarray, shape: float, rate: float, e: np.ndarray) -> np.ndarray:
    """log P(n) for n ~ NegBin(shape, p = rate / (rate + E)), the gamma-Poisson marginal."""
    p = rate / (rate + e)
    return special.gammaln(shape + n) - special.gammaln(shape) - special.gammaln(n + 1) + shape * np.log(p) \
        + n * np.log1p(-p)


def _log_sf(t: np.ndarray, shape: float, rate: float, e: np.ndarray) -> np.ndarray:
    """log P(n >= t) under the same negative binomial: P(n >= t) = I_{1-p}(t, shape) (regularised beta)."""
    q = e / (rate + e)
    with np.errstate(divide="ignore"):
        out = np.log(special.betainc(t, shape, q))
    bad = ~np.isfinite(out)
    if bad.any():  # underflow far in the tail: fall back to scipy's log-survival function
        out[bad] = stats.nbinom.logsf(t[bad] - 1, shape, rate / (rate + e[bad]))
    return out


def _unpack(x: np.ndarray) -> tuple[float, float, float, float, float]:
    a1, b1, a2, b2 = np.exp(x[:4])
    w = special.expit(x[4])
    return a1, b1, a2, b2, w


def negative_loglik(x: np.ndarray, n: np.ndarray, e: np.ndarray, t: np.ndarray) -> float:
    """Negative truncated log-likelihood of parameters on the log / logit scale."""
    a1, b1, a2, b2, w = _unpack(x)
    with np.errstate(all="ignore"):  # extreme trial parameters give -inf terms; they are rejected below
        return _negative_loglik(a1, b1, a2, b2, w, n, e, t)


def _negative_loglik(a1: float, b1: float, a2: float, b2: float, w: float, n: np.ndarray, e: np.ndarray,
                     t: np.ndarray) -> float:
    l1 = np.log(w) + _log_nb(n, a1, b1, e)
    l2 = np.log1p(-w) + _log_nb(n, a2, b2, e)
    ll = np.logaddexp(l1, l2)
    trunc = t > 0
    if trunc.any():
        s1 = np.log(w) + _log_sf(t[trunc], a1, b1, e[trunc])
        s2 = np.log1p(-w) + _log_sf(t[trunc], a2, b2, e[trunc])
        ll[trunc] -= np.logaddexp(s1, s2)
    val = -float(np.sum(ll))
    return val if np.isfinite(val) else 1e300


def fit_prior(n: np.ndarray, e: np.ndarray, t: np.ndarray | None = None,
              starts: tuple[tuple[float, ...], ...] | None = None) -> GPSPrior:
    """Maximum-likelihood two-gamma prior from observed counts ``n``, expected counts ``e`` and truncation ``t``.

    ``t[i]`` is the smallest count pair i could have been observed with (0 = untruncated). Several starting points
    are tried and the best optimum is kept.
    """
    n = np.asarray(n, dtype=float)
    e = np.asarray(e, dtype=float)
    t = np.zeros_like(n) if t is None else np.asarray(t, dtype=float)
    starts = starts or (DUMOUCHEL_START, (1.0, 1.0, 3.0, 1.0, 0.5), (0.5, 0.5, 5.0, 2.0, 0.8))
    best = None
    for s in starts:
        x0 = np.r_[np.log(s[:4]), special.logit(s[4])]
        res = optimize.minimize(negative_loglik, x0, args=(n, e, t), method="Nelder-Mead",
                                options={"maxiter": 4000, "xatol": 1e-6, "fatol": 1e-6})
        res = optimize.minimize(negative_loglik, res.x, args=(n, e, t), method="L-BFGS-B")
        if best is None or res.fun < best.fun:
            best = res
    a1, b1, a2, b2, w = _unpack(best.x)
    # Order components so that component 1 has the lower prior mean (identifiability for reporting only).
    if a1 / b1 > a2 / b2:
        a1, b1, a2, b2, w = a2, b2, a1, b1, 1 - w
    return GPSPrior(float(a1), float(b1), float(a2), float(b2), float(w), loglik=-float(best.fun), n_pairs=int(n.size))


def posterior(n: np.ndarray, e: np.ndarray, prior: GPSPrior, q: float = 0.05) -> tuple[np.ndarray, ...]:
    """EBGM, the q-th and (1 - q)-th posterior quantiles of lambda, and the weight of component 1, per pair.

    Pairs with no expected count (E = 0: the drug or event has no reports) get NaN.
    """
    n = np.asarray(n, dtype=float)
    e = np.asarray(e, dtype=float)
    ok = e > 0
    if not ok.all():
        out = [np.full(n.shape, np.nan) for _ in range(4)]
        for arr, val in zip(out, posterior(n[ok], e[ok], prior, q)):
            arr[ok] = val
        return tuple(out)
    l1 = np.log(prior.w) + _log_nb(n, prior.a1, prior.b1, e)
    l2 = np.log1p(-prior.w) + _log_nb(n, prior.a2, prior.b2, e)
    q1 = np.exp(l1 - np.logaddexp(l1, l2))
    s1, r1 = prior.a1 + n, prior.b1 + e
    s2, r2 = prior.a2 + n, prior.b2 + e
    ebgm = np.exp(q1 * (special.digamma(s1) - np.log(r1)) + (1 - q1) * (special.digamma(s2) - np.log(r2)))
    lo = np.array([_mixture_quantile(q, *args) for args in zip(q1, s1, r1, s2, r2)])
    hi = np.array([_mixture_quantile(1 - q, *args) for args in zip(q1, s1, r1, s2, r2)])
    return ebgm, lo, hi, q1


def _mixture_quantile(p: float, q1: float, s1: float, r1: float, s2: float, r2: float) -> float:
    """Quantile of q1 * Gamma(s1, r1) + (1 - q1) * Gamma(s2, r2) by root finding between the component quantiles."""
    x1 = stats.gamma.ppf(p, s1, scale=1 / r1)
    x2 = stats.gamma.ppf(p, s2, scale=1 / r2)
    lo, hi = min(x1, x2), max(x1, x2)
    if hi - lo < 1e-12 * max(hi, 1e-300):
        return float(lo)

    def f(x: float) -> float:
        return q1 * stats.gamma.cdf(x, s1, scale=1 / r1) + (1 - q1) * stats.gamma.cdf(x, s2, scale=1 / r2) - p

    if f(lo) >= 0:
        return float(lo)
    if f(hi) <= 0:
        return float(hi)
    return float(optimize.brentq(f, lo, hi, xtol=1e-12 * max(hi, 1e-12), maxiter=200))
