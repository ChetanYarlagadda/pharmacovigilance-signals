import numpy as np
import pandas as pd
import pytest
from scipy import stats

from pvsignals.disproportionality import (chi2_yates, information_component, prr, rank_score, ror, score_table,
                                          two_by_two)

# a = 20 co-reports, 1,000 reports of the drug, 500 of the event, 100,000 in total
A, ND, NE, N = 20, 1_000, 500, 100_000


def test_two_by_two_cells_and_ratios_match_hand_calculation():
    a, b, c, d = two_by_two(np.array([A]), np.array([ND]), np.array([NE]), N)
    assert (a[0], b[0], c[0], d[0]) == (20, 980, 480, 98_520)
    p, lo, hi = prr(a, b, c, d)
    assert p[0] == pytest.approx((20 / 1000) / (480 / 99_000))  # 4.125
    o, olo, ohi = ror(a, b, c, d)
    assert o[0] == pytest.approx(20 * 98_520 / (980 * 480))  # 4.189
    se = np.sqrt(1 / 20 + 1 / 980 + 1 / 480 + 1 / 98_520)
    assert olo[0] == pytest.approx(o[0] * np.exp(-1.959964 * se), rel=1e-5)
    assert lo[0] < p[0] < hi[0]


def test_chi2_matches_scipy_with_yates_correction():
    a, b, c, d = 12, 300, 40, 9_000
    expected = stats.chi2_contingency([[a, b], [c, d]], correction=True)[0]
    assert chi2_yates(np.array([a]), np.array([b]), np.array([c]), np.array([d]))[0] == pytest.approx(expected)


def test_information_component_is_zero_at_expectation_and_interval_brackets_it():
    ic, lo, hi = information_component(np.array([0, 10, 1000]), np.array([0.0, 10.0, 1000.0]))
    assert ic == pytest.approx([0, 0, 0])
    assert np.all(lo < ic) and np.all(ic < hi)
    # more data -> narrower credibility interval
    assert (hi - lo)[2] < (hi - lo)[1] < (hi - lo)[0]


def test_zero_co_reports_have_no_ratio_but_rank_lowest():
    a, b, c, d = two_by_two(np.array([0, 3]), np.array([100, 100]), np.array([50, 50]), 10_000)
    p, lo, _ = prr(a, b, c, d)
    assert np.isnan(p[0]) and np.isfinite(p[1])
    s = rank_score(pd.Series(lo))
    assert s.iloc[0] == 0 and s.iloc[1] > 0


def test_inconsistent_margins_are_rejected():
    with pytest.raises(ValueError):
        two_by_two(np.array([60]), np.array([50]), np.array([100]), 1_000)


def test_score_table_adds_every_statistic():
    df = pd.DataFrame({"a": [A, 0], "n_drug": [ND, ND], "n_event": [NE, NE]})
    out = score_table(df, N)
    for col in ["expected", "prr", "prr025", "ror", "ror025", "chi2", "ic", "ic025", "ic975"]:
        assert col in out
    assert out.loc[0, "expected"] == pytest.approx(ND * NE / N)
    assert out.loc[1, "ic"] < 0 < out.loc[0, "ic"]
