"""Frequentist and Bayesian disproportionality statistics on 2x2 tables of report counts.

For a drug D and an event E, with N reports in the database::

                    E        not E
        D           a          b        a + b = reports that mention D
        not D       c          d
                  a + c                 N = a + b + c + d

- PRR = [a / (a + b)] / [c / (c + d)], 95% CI from the log-normal approximation (Evans et al. 2001).
- ROR = (a d) / (b c), 95% CI from Woolf's log-normal approximation.
- IC = log2[(a + 0.5) / (E + 0.5)] with E = (a + b)(a + c) / N, the shrunk information component of the Uppsala
  Monitoring Centre's BCPNN (Noren et al. 2013); its 95% credibility interval comes from the Gamma(a + 0.5, E + 0.5)
  posterior of the observed-to-expected ratio.
- EBGM (Gamma-Poisson shrinker, DuMouchel 1999) lives in ``gps.py``.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

Z95 = stats.norm.ppf(0.975)


def two_by_two(a: np.ndarray, n_drug: np.ndarray, n_event: np.ndarray, n_total: float) -> tuple[np.ndarray, ...]:
    """Cells (a, b, c, d) from the co-report count and the margins."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(n_drug, dtype=float) - a
    c = np.asarray(n_event, dtype=float) - a
    d = float(n_total) - a - b - c
    if np.any(b < 0) or np.any(c < 0) or np.any(d < 0):
        raise ValueError("inconsistent margins: a exceeds a margin or margins exceed N")
    return a, b, c, d


def expected(n_drug: np.ndarray, n_event: np.ndarray, n_total: float) -> np.ndarray:
    """Expected co-reports under independence, E = (a + b)(a + c) / N."""
    return np.asarray(n_drug, dtype=float) * np.asarray(n_event, dtype=float) / float(n_total)


def prr(a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Proportional reporting ratio with its 95% CI; NaN where a = 0."""
    a, b, c, d = (np.asarray(x, dtype=float) for x in (a, b, c, d))
    with np.errstate(divide="ignore", invalid="ignore"):
        est = (a / (a + b)) / (c / (c + d))
        se = np.sqrt(1 / a - 1 / (a + b) + 1 / c - 1 / (c + d))
        lo, hi = est * np.exp(-Z95 * se), est * np.exp(Z95 * se)
    bad = a == 0
    return np.where(bad, np.nan, est), np.where(bad, np.nan, lo), np.where(bad, np.nan, hi)


def ror(a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Reporting odds ratio with its 95% CI; NaN where a = 0."""
    a, b, c, d = (np.asarray(x, dtype=float) for x in (a, b, c, d))
    with np.errstate(divide="ignore", invalid="ignore"):
        est = (a * d) / (b * c)
        se = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
        lo, hi = est * np.exp(-Z95 * se), est * np.exp(Z95 * se)
    bad = a == 0
    return np.where(bad, np.nan, est), np.where(bad, np.nan, lo), np.where(bad, np.nan, hi)


def chi2_yates(a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray) -> np.ndarray:
    """Pearson chi-square with Yates' continuity correction (1 df), as used in the Evans PRR criteria."""
    a, b, c, d = (np.asarray(x, dtype=float) for x in (a, b, c, d))
    n = a + b + c + d
    num = n * np.maximum(np.abs(a * d - b * c) - n / 2, 0) ** 2
    den = (a + b) * (c + d) * (a + c) * (b + d)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(den > 0, num / den, 0.0)


def information_component(a: np.ndarray, e: np.ndarray, level: float = 0.95) -> tuple[np.ndarray, ...]:
    """Shrunk IC = log2((a + 0.5) / (E + 0.5)) and its credibility interval from Gamma(a + 0.5, rate E + 0.5)."""
    a = np.asarray(a, dtype=float)
    e = np.asarray(e, dtype=float)
    shape, rate = a + 0.5, e + 0.5
    tail = (1 - level) / 2
    ic = np.log2(shape / rate)
    lo = np.log2(stats.gamma.ppf(tail, shape, scale=1 / rate))
    hi = np.log2(stats.gamma.ppf(1 - tail, shape, scale=1 / rate))
    return ic, lo, hi


def score_table(df: pd.DataFrame, n_total: float) -> pd.DataFrame:
    """Add PRR, ROR, IC (with 95% intervals), chi-square and E to a frame with columns a, n_drug, n_event.

    Pairs with no co-reports (a = 0) have no PRR or ROR (NaN); ``rank_score`` places them below every defined value
    when pairs are ranked. IC and EBGM are defined for a = 0.
    """
    out = df.copy()
    a, b, c, d = two_by_two(out["a"], out["n_drug"], out["n_event"], n_total)
    out["expected"] = expected(out["n_drug"], out["n_event"], n_total)
    out["prr"], out["prr025"], out["prr975"] = prr(a, b, c, d)
    out["ror"], out["ror025"], out["ror975"] = ror(a, b, c, d)
    out["chi2"] = chi2_yates(a, b, c, d)
    out["ic"], out["ic025"], out["ic975"] = information_component(a, out["expected"])
    return out


def rank_score(values: pd.Series) -> pd.Series:
    """Scores for ranking: undefined ratio statistics (a = 0) rank below every defined value."""
    return values.fillna(0.0)
