"""Cached openFDA client, drug-name matching rules and evaluation metrics."""

import json

import numpy as np
import pandas as pd
import pytest
from sklearn.metrics import roc_auc_score

from pvsignals.evaluation import auc, holm, paired_auc_tests, threshold_table
from pvsignals.openfda import BudgetExceeded, OpenFDA
from pvsignals.reference import canonical_drug, drug_spec, term_matches


class FakeResponse:
    def __init__(self, status: int, payload: dict):
        self.status_code, self._payload = status, payload

    def json(self):
        return self._payload

    def raise_for_status(self):
        if self.status_code >= 400 and self.status_code != 404:
            raise RuntimeError(self.status_code)


def test_client_caches_logs_and_respects_offline_and_budget(tmp_path, monkeypatch):
    client = OpenFDA(tmp_path, daily_budget=2, min_interval_s=0)
    calls = []

    def fake_get(url, timeout):
        calls.append(url)
        if "NOTHING" in url:
            return FakeResponse(404, {"error": {"code": "NOT_FOUND"}})
        return FakeResponse(200, {"meta": {"results": {"total": 42}}, "results": []})

    monkeypatch.setattr(client.session, "get", fake_get)
    assert client.total("x:1")[0] == 42
    assert client.total("x:1")[0] == 42  # second call served from the cache
    assert len(calls) == 1 and client.cache_hits == 1
    assert client.total("x:NOTHING")[0] == 0  # 404 = no matching reports
    rec = json.loads(next((tmp_path / "openfda").glob("*.json")).read_text())
    assert rec["url"].startswith("https://api.fda.gov/drug/event.json?search=") and rec["fetched_utc"]
    with pytest.raises(BudgetExceeded):  # 2 network calls are logged, the budget is 2
        client.total("x:2")
    offline = OpenFDA(tmp_path, offline=True)
    assert offline.total("x:1")[0] == 42
    with pytest.raises(FileNotFoundError):
        offline.total("x:3")


@pytest.mark.parametrize("term,drug,expected", [
    ("ESTRADIOL VALERATE", "estradiol", True),
    ("ETHINYL ESTRADIOL\\LEVONORGESTREL", "estradiol", False),
    ("ETHINYL ESTRADIOL\\LEVONORGESTREL", "levonorgestrel", True),
    ("NIACINAMIDE", "niacin", False),
    ("AMLODIPINE BESYLATE\\ATORVASTATIN CALCIUM", "atorvastatin", True),
    ("METHYLPREDNISOLONE", "prednisolone", False),
    ("DIVALPROEX SODIUM", "valproate", True),
    ("INSULIN HUMAN", "regular insulin, human", True),
    ("INSULIN HUMAN ISOPHANE", "regular insulin, human", False),
    ("LEVOTHYROXINE SODIUM", "thyroxine", True),
    ("LITHIUM CARBONATE", "lithium citrate", False),
])
def test_term_matching(term, drug, expected):
    assert term_matches(term, drug_spec(drug)) is expected


def test_pancreatic_enzymes_map_to_one_drug():
    assert {canonical_drug(x) for x in ["Amylases", "Lipase", "Endopeptidases"]} == {"pancrelipase"}
    assert drug_spec("Lipase").name == "pancrelipase"


def test_auc_matches_sklearn_with_ties():
    rng = np.random.default_rng(3)
    y = rng.random(300) < 0.4
    s = np.round(rng.normal(size=300) + y, 1)  # rounding creates ties
    assert auc(y, s) == pytest.approx(roc_auc_score(y, s))


def test_holm_and_paired_test():
    assert holm(np.array([0.01, 0.04, 0.03])) == pytest.approx([0.03, 0.06, 0.06])
    rng = np.random.default_rng(1)
    y = np.r_[np.ones(60), np.zeros(60)].astype(int)
    df = pd.DataFrame({"ground_truth": y, "good": y + rng.normal(0, 0.5, 120), "noise": rng.normal(size=120)})
    t = paired_auc_tests(df, [("good", "noise")], n_boot=300)
    assert t.loc[0, "auc_diff"] > 0.3 and t.loc[0, "p_holm"] < 0.01


def test_threshold_rules():
    df = pd.DataFrame({"ground_truth": [1, 1, 0, 0], "a": [10, 2, 10, 10], "prr": [3, 5, 1, 2.5],
                       "chi2": [9, 9, 0, 1], "ror025": [2, 3, 0.5, 1.5], "ic025": [0.5, 0.1, -1, 0.2],
                       "eb05": [2.5, 1, 0.5, 1.9]})
    t = threshold_table(df).set_index("rule")
    evans = t.loc["PRR (Evans: PRR >= 2, chi2 >= 4, a >= 3)"]
    assert (evans.true_pos, evans.false_pos) == (1, 0)  # pair 2 has a < 3, pair 4 has chi2 < 4
    assert t.loc["IC lower 95% > 0", "false_pos"] == 1
