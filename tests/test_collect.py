"""Count assembly against a fake openFDA: term sums, the 999-term cut-off, subset fill-ins."""

import pandas as pd

from pvsignals.collect import (CountTable, Fetcher, any_of, and_, drug_marginals, pair_counts, plan_batches,
                               window_clause)
from pvsignals.reference import drug_spec


class FakeClient:
    """Answers count / total queries from a dict keyed by (search, field); records every query."""

    def __init__(self, tables: dict, totals: dict, limit: int = 4):
        self.tables, self.totals, self.limit, self.queries = tables, totals, limit, []

    def count(self, search, field):
        self.queries.append(search)
        counts = dict(sorted(self.tables[(search, field)].items(), key=lambda kv: -kv[1])[: self.limit])
        return counts, {"url": f"{search}|{field}", "fetched_utc": "2026-01-01T00:00:00+00:00"}

    def total(self, search):
        self.queries.append(search)
        return self.totals[search], {"url": search, "fetched_utc": "2026-01-01T00:00:00+00:00"}


def test_query_builders():
    assert window_clause("20150101", "20151231") == "receivedate:[20150101 TO 20151231]"
    assert any_of("f", ["A B", "C"]) == 'f:("A B" "C")'
    assert and_("x", "", "y") == "x AND y"


def test_count_table_floor():
    assert CountTable({"A": 5, "B": 2}, "u", "t", limit=4).floor == 0
    full = CountTable({"A": 9, "B": 7, "C": 6, "D": 3}, "u", "t", limit=4)
    assert not full.complete and full.floor == 3


def test_marginals_sum_salt_and_combination_terms_from_the_subset_query():
    base, field = "W", "patient.drug.activesubstance.activesubstancename"
    overall = CountTable({"ASPIRIN": 900, "ATORVASTATIN": 500, "METFORMIN": 400, "ATORVASTATIN CALCIUM": 300},
                         "u", "t", limit=4)
    sub_search = and_(base, any_of(field, ["atorvastatin"]))
    tables = {(sub_search, field + ".exact"): {"ATORVASTATIN": 500, "ATORVASTATIN CALCIUM": 300, "ASPIRIN": 200,
                                               "AMLODIPINE BESYLATE\\ATORVASTATIN CALCIUM": 40, "EZETIMIBE": 10}}
    f = Fetcher(FakeClient(tables, {}, limit=10), base)
    spec = drug_spec("atorvastatin")
    batches = plan_batches([spec], overall)
    marg, terms = drug_marginals(f, [spec], overall, batches)
    assert marg.loc[0, "n_drug"] == 840 and marg.loc[0, "n_terms"] == 3
    assert "ASPIRIN" not in terms["atorvastatin"]


def test_pair_counts_fill_in_terms_below_the_event_table_cut_off():
    base, dfield = "W", "patient.drug.activesubstance.activesubstancename"
    pts = ["RHABDOMYOLYSIS"]
    ev_clause = any_of("patient.reaction.reactionmeddrapt.exact", pts)
    # the event table lists 4 terms; atorvastatin's combination term (20 co-reports) falls below its cut-off (25)
    event_table = CountTable({"ATORVASTATIN": 100, "SIMVASTATIN": 80, "ROSUVASTATIN": 60, "ASPIRIN": 25}, "u", "t",
                             limit=4)
    events = {"rhabdo": {"pts": pts, "n_event": 400, "drugs": event_table}}
    terms_of = {"atorvastatin": {"ATORVASTATIN": 5000, "AMLODIPINE BESYLATE\\ATORVASTATIN CALCIUM": 900},
                "simvastatin": {"SIMVASTATIN": 4000}}
    fill_search = and_(base, ev_clause, any_of(dfield, ["atorvastatin"]))
    tables = {(fill_search, dfield + ".exact"): {"ATORVASTATIN": 100, "AMLODIPINE BESYLATE\\ATORVASTATIN CALCIUM": 20,
                                                 "ASPIRIN": 9}}
    client = FakeClient(tables, {}, limit=4)
    f = Fetcher(client, base)
    specs = {d: drug_spec(d) for d in terms_of}
    pairs = pd.DataFrame({"drug": ["atorvastatin", "simvastatin"], "event": ["rhabdo", "rhabdo"]})
    out = pair_counts(f, pairs, specs, terms_of, events).set_index("drug")
    assert out.loc["simvastatin", "a"] == 80 and out.loc["simvastatin", "a_method"] == "table"
    assert out.loc["atorvastatin", "a"] == 120 and out.loc["atorvastatin", "a_method"] == "fill-in"
    assert out.loc["atorvastatin", "a_shortfall_max"] == 0
    assert client.queries == [fill_search]  # only the pair that needed it was re-queried
