"""Fetch openFDA counts, build 2x2 tables, score PRR / ROR / IC / EBGM, evaluate on reference controls, write results.

    python scripts/run_pipeline.py            # full run: about 300 openFDA requests with an empty cache, ~1 min cached
    python scripts/run_pipeline.py --quick    # smoke test: 2 outcomes, ~25 requests; outputs go to data/cache/smoke/
    python scripts/run_pipeline.py --offline  # re-run from the cache only (fails on any uncached query)

Every openFDA response is cached in data/raw/openfda/ with its query URL and fetch time, so a re-run makes no
network calls. The request log (data/raw/openfda_request_log.jsonl) records every call that was made.
"""

from __future__ import annotations

import argparse
import json
import platform
import re
import sys
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from pvsignals import DISCLAIMER  # noqa: E402
from pvsignals.collect import (DRUG_EXACT, DRUG_FIELD, PT_FIELD, WINDOW, Fetcher, any_of, drug_marginals,  # noqa: E402
                               outcome_tables, pair_counts, plan_batches, window_clause)
from pvsignals.disproportionality import score_table  # noqa: E402
from pvsignals.evaluation import auc_table, paired_auc_tests, roc_curve, threshold_table  # noqa: E402
from pvsignals.figures import make_all  # noqa: E402
from pvsignals.gps import fit_prior, posterior  # noqa: E402
from pvsignals.openfda import COUNT_LIMIT, OpenFDA  # noqa: E402
from pvsignals.reference import (EVENT_PTS, OUTCOME_LABEL, OUTCOME_PTS, drug_spec, load_reference,  # noqa: E402
                                 term_matches)

RAW = ROOT / "data" / "raw"
RESULTS = ROOT / "reports" / "results"
FIGS = ROOT / "reports" / "figures"

# Textbook signals named in the project brief, each with comparators from the same therapeutic area.
TEXTBOOK = {
    "tendon_rupture": {"fluoroquinolone": ["ciprofloxacin", "levofloxacin", "moxifloxacin", "ofloxacin"],
                       "other antibiotic": ["amoxicillin", "azithromycin", "doxycycline", "cephalexin",
                                            "nitrofurantoin"]},
    "rhabdomyolysis": {"statin": ["atorvastatin", "simvastatin", "rosuvastatin", "pravastatin", "lovastatin",
                                  "fluvastatin"],
                       "other lipid-lowering": ["ezetimibe"]},
    "depression": {"isotretinoin": ["isotretinoin"], "other acne drug": ["doxycycline", "minocycline"]},
}
TEXTBOOK_INDEX = [("tendon_rupture", "fluoroquinolone", "other antibiotic"),
                  ("rhabdomyolysis", "statin", "other lipid-lowering"),
                  ("depression", "isotretinoin", "other acne drug")]
SCREEN_DRUGS = ["levofloxacin", "atorvastatin", "isotretinoin"]
SCORES = ["prr", "prr025", "ror", "ror025", "ic", "ic025", "ebgm", "eb05"]
LOWER = ["prr025", "ror025", "ic025", "eb05"]
# Drugs approved inside the analysis window (FDA approval dates), for the time-since-approval view.
NEW_DRUGS = {"sacubitril": "2015-07-07", "dupilumab": "2017-03-28", "erenumab": "2018-05-17",
             "tirzepatide": "2022-05-13"}
AORTIC_PTS = ["AORTIC ANEURYSM", "AORTIC DISSECTION", "AORTIC ANEURYSM RUPTURE"]
REPORTER_CHECKS = [("levofloxacin", AORTIC_PTS), ("atorvastatin", ["TYPE 2 DIABETES MELLITUS"]),
                   ("atorvastatin", ["CHRONIC KIDNEY DISEASE"]), ("levofloxacin", ["CHRONIC KIDNEY DISEASE"]),
                   ("isotretinoin", ["COLITIS ULCERATIVE"])]
QUALIFICATION = {"1": "Physician", "2": "Pharmacist", "3": "Other health professional", "4": "Lawyer",
                 "5": "Consumer"}
SEED = 0


def collect(f: Fetcher, specs: dict, needed: set[str], events: dict[str, list[str]], pairs: pd.DataFrame) -> dict:
    """All counts under the fetcher's base filter: N, drug margins, event tables, co-report counts."""
    n_total = f.total()
    if n_total == 0:
        raise RuntimeError(f"the base filter matches no reports: {f.base}")
    overall = f.table(DRUG_EXACT)
    all_specs = sorted(specs.values(), key=lambda s: s.name)
    batches = [b for b in plan_batches(all_specs, overall) if any(s.name in needed for s in b)]
    marg, terms_of = drug_marginals(f, all_specs, overall, batches)
    marg = marg[marg.drug.isin(needed)].reset_index(drop=True)
    ev = outcome_tables(f, events)
    a = pair_counts(f, pairs, specs, terms_of, ev)
    return {"N": n_total, "overall": overall, "marginals": marg, "terms_of": terms_of, "events": ev, "a": a}


def gps_pool(data: dict) -> pd.DataFrame:
    """Term-level (active-substance term, event) pairs for fitting the GPS prior, with truncation thresholds."""
    overall = data["overall"].counts
    rows = []
    for ev, rec in data["events"].items():
        tab = rec["drugs"]
        t = tab.floor + 1  # listed terms all have more than `floor` reports; untruncated tables still omit zeros
        for term, n in tab.counts.items():
            if term in overall and n >= t:
                rows.append({"term": term, "event": ev, "n": n, "n_drug": overall[term], "n_event": rec["n_event"],
                             "t": t})
    pool = pd.DataFrame(rows)
    pool["expected"] = pool.n_drug * pool.n_event / data["N"]
    return pool


def score_pairs(pairs: pd.DataFrame, data: dict, prior) -> pd.DataFrame:
    """Join counts and margins, then add PRR / ROR / IC and the EBGM posterior."""
    m = data["marginals"][["drug", "n_drug"]]
    ev = pd.DataFrame([{"event": k, "n_event": v["n_event"]} for k, v in data["events"].items()])
    df = pairs.merge(data["a"], on=["drug", "event"], how="left").merge(m, on="drug").merge(ev, on="event")
    df = score_table(df, data["N"])
    df["ebgm"], df["eb05"], df["eb95"], _ = posterior(df["a"], df["expected"], prior)
    return df


def date_series(f: Fetcher, *clauses: str) -> pd.Series:
    """Reports per calendar quarter (from openFDA's per-day counts of receivedate)."""
    tab = f.table("receivedate", *clauses)
    s = pd.Series(tab.counts, dtype=float)
    s.index = pd.to_datetime(s.index, format="%Y%m%d")
    return s.groupby(s.index.to_period("Q")).sum()


def validation(f: Fetcher, df: pd.DataFrame, marg: pd.DataFrame, specs: dict, seen_terms: set[str],
               n_each: int = 10) -> pd.DataFrame:
    """Compare term-summed counts with exact openFDA search totals for a random sample of drugs and pairs.

    A token search also matches substances that merely contain the drug's name (``estradiol`` in ``ETHINYL
    ESTRADIOL``); such terms, among all terms seen in any count table, are listed in ``token_also_matches``.
    """
    rng = np.random.default_rng(SEED)
    n_drug = marg.set_index("drug").n_drug
    simple = [d for d in sorted(marg.drug) if specs[d].search == (d,) and " " not in d]
    big = [d for d in simple if n_drug[d] >= 100]
    rows = []
    for d in rng.choice(big, size=min(n_each, len(big)), replace=False):
        rows.append({"kind": "drug reports", "drug": d, "event": "", "term_sum": int(n_drug[d]),
                     "exact_search": f.total(any_of(DRUG_FIELD, [d]))})
    cand = df[(df.a > 0) & df.drug.isin(simple)].drop_duplicates(["drug", "event"])
    for _, r in cand.sample(n=min(n_each, len(cand)), random_state=SEED).iterrows():
        exact = f.total(any_of(DRUG_FIELD, [r.drug]), any_of(PT_FIELD, OUTCOME_PTS[r.event]))
        rows.append({"kind": "co-reports (a)", "drug": r.drug, "event": r.event, "term_sum": int(r.a),
                     "exact_search": exact})
    out = pd.DataFrame(rows)
    out["rel_diff"] = (out.term_sum - out.exact_search) / out.exact_search
    out["token_also_matches"] = [
        "; ".join(sorted(t for t in seen_terms if re.search(rf"\b{re.escape(d)}\b", t.lower())
                         and not term_matches(t, specs[d]))[:5]) for d in out.drug]
    return out


def evaluate(scored: pd.DataFrame, stratum: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """AUCs, paired tests and threshold metrics for OMOP, EU-ADR and the pooled unique pairs."""
    sets = {"OMOP": scored[scored.set == "omop"], "EU-ADR": scored[scored.set == "euadr"],
            "Pooled": scored.drop_duplicates(["drug", "outcome"])}
    aucs, tests, thr, rocs = [], [], [], []
    comps = [(LOWER[i], LOWER[j]) for i in range(4) for j in range(i + 1, 4)] + \
            [(p, lo) for p, lo in zip(["prr", "ror", "ic", "ebgm"], LOWER)]
    for name, d in sets.items():
        d = d[d.evaluable]
        x = d.copy()
        x[["prr", "prr025", "ror", "ror025"]] = x[["prr", "prr025", "ror", "ror025"]].fillna(0.0)
        aucs.append(auc_table(x, SCORES, seed=SEED).assign(set=name, stratum=stratum))
        tests.append(paired_auc_tests(x, comps, seed=SEED + 1).assign(set=name, stratum=stratum))
        thr.append(threshold_table(d).assign(set=name, stratum=stratum))
        for s in LOWER + ["prr", "ror", "ic", "ebgm"]:
            rocs.append(roc_curve(x.ground_truth, x[s]).assign(set=name, stratum=stratum, score=s))
    return pd.concat(aucs), pd.concat(tests), pd.concat(thr), pd.concat(rocs)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--quick", action="store_true", help="smoke test on two outcomes; writes to data/cache/smoke/")
    ap.add_argument("--offline", action="store_true", help="fail instead of calling openFDA for uncached queries")
    args = ap.parse_args()
    t0 = time.time()
    out_dir = ROOT / "data" / "cache" / "smoke" if args.quick else RESULTS
    out_dir.mkdir(parents=True, exist_ok=True)
    client = OpenFDA(RAW, daily_budget=600, offline=args.offline)

    ref = load_reference(RAW / "reference")
    tb_rows = [{"event": ev, "group": group, "drug": d}
               for ev, groups in TEXTBOOK.items() for group, drugs in groups.items() for d in drugs]
    textbook = pd.DataFrame(tb_rows, columns=["event", "group", "drug"])
    # Specs always cover the full drug list, so the batch plan (and its cached queries) is the same in both modes.
    specs = {s.name: s for s in map(drug_spec, sorted(set(ref.drug_name.str.lower()) | set(textbook.drug)))}
    outcomes = ["rhabdomyolysis", "anaphylaxis"] if args.quick else list(OUTCOME_PTS)
    ref = ref[ref.outcome.isin(outcomes)].reset_index(drop=True)
    events = {k: OUTCOME_PTS[k] for k in outcomes}
    if args.quick:
        textbook = textbook.iloc[0:0]
    else:
        events.update(EVENT_PTS)
    needed = set(ref.drug) | set(textbook.drug)
    pairs = pd.concat([ref[["drug", "outcome"]].rename(columns={"outcome": "event"}), textbook[["drug", "event"]]])
    pairs = pairs.drop_duplicates().reset_index(drop=True)

    # ---- every report received 2015-2025
    f = Fetcher(client, window_clause())
    data = collect(f, specs, needed, events, pairs)
    pool = gps_pool(data)
    prior = fit_prior(pool.n, pool.expected, pool.t)
    pool["ebgm"], pool["eb05"], pool["eb95"], _ = posterior(pool.n, pool.expected, prior)
    print(f"N={data['N']:,}  pool={len(pool):,}  prior={prior}")

    ref_pairs = ref.rename(columns={"outcome": "event"})
    scored = score_pairs(ref_pairs, data, prior).rename(columns={"event": "outcome"})
    scored["evaluable"] = scored.n_drug > 0
    aucs, tests, thr, rocs = evaluate(scored, "all reports")
    # Sensitivity: every term left unlisted by a count table contributes its maximum possible count.
    upper = dict(data, a=data["a"].assign(a=data["a"].a + data["a"].a_shortfall_max))
    scored_ub = score_pairs(ref_pairs, upper, prior).rename(columns={"event": "outcome"})
    scored_ub["evaluable"] = scored_ub.n_drug > 0
    aucs = pd.concat([aucs, evaluate(scored_ub, "all reports, unlisted terms at upper bound")[0]])

    # ---- outcome definitions: PT composition and union vs sum
    od_rows = []
    for ev, rec in data["events"].items():
        pc = rec["pt_counts"]
        for pt in rec["pts"]:
            od_rows.append({"event": ev, "label": OUTCOME_LABEL.get(ev, ev.replace("_", " ")), "pt": pt,
                            "pt_reports": pc.counts.get(pt, 0), "event_reports": rec["n_event"],
                            "drug_terms_listed": len(rec["drugs"].counts), "drug_table_floor": rec["drugs"].floor})
    outcome_defs = pd.DataFrame(od_rows)

    meta_extra: dict = {}
    if not args.quick:
        # ---- textbook signals with comparators
        tb = score_pairs(textbook, data, prior)
        tb_cmp = []
        for ev, target, comparator in TEXTBOOK_INDEX:
            t = tb[(tb.event == ev) & (tb.group == target)]
            for _, r in t.iterrows():
                c = tb[(tb.event == ev) & (tb.group == comparator)]
                a, b = r.a, r.n_drug - r.a
                cc, dd = c.a.sum(), (c.n_drug - c.a).sum()
                est = (a * dd) / (b * cc) if a > 0 and cc > 0 else np.nan
                se = np.sqrt(1 / a + 1 / b + 1 / cc + 1 / dd) if a > 0 and cc > 0 else np.nan
                tb_cmp.append({"event": ev, "drug": r.drug, "comparator": f"{comparator} ({', '.join(c.drug)})",
                               "a": int(a), "n_drug": int(r.n_drug), "comparator_a": int(cc),
                               "comparator_reports": int(c.n_drug.sum()), "ror_vs_all": r.ror,
                               "ror_vs_comparator": est, "ci_low": est * np.exp(-1.96 * se),
                               "ci_high": est * np.exp(1.96 * se)})
        tb_cmp = pd.DataFrame(tb_cmp)

        # ---- screening view: every PT reported with three drugs (PTs with a known overall count only: the 999 most
        # reported PTs plus every PT of the outcome definitions, whose counts come from the outcome tables)
        pt_overall = f.table(PT_FIELD)
        pt_marg = dict(pt_overall.counts)
        for rec in data["events"].values():
            if len(rec["pts"]) == 1:
                pt_marg[rec["pts"][0]] = rec["n_event"]
            else:
                pt_marg.update({p: rec["pt_counts"].counts[p] for p in rec["pts"] if p in rec["pt_counts"].counts})
        screen = []
        for d in SCREEN_DRUGS:
            tab = f.table(PT_FIELD, any_of(DRUG_FIELD, [d]))
            n_d = int(data["marginals"].set_index("drug").loc[d, "n_drug"])
            for pt, a in tab.counts.items():
                if pt in pt_marg:
                    screen.append({"drug": d, "pt": pt, "a": a, "n_drug": n_d, "n_event": pt_marg[pt]})
        screen = score_table(pd.DataFrame(screen), data["N"])
        screen["ebgm"], screen["eb05"], screen["eb95"], _ = posterior(screen.a, screen.expected, prior)
        screen["rank_eb05"] = screen.groupby("drug").eb05.rank(ascending=False, method="first")
        meta_extra["screen_pt_floor"] = pt_overall.floor

        # ---- validation of term-summed counts against exact searches
        seen = set(data["overall"].counts).union(*(rec["drugs"].counts for rec in data["events"].values()))
        val = validation(f, scored[scored.evaluable].rename(columns={"outcome": "event"}), data["marginals"], specs,
                         seen)

        # ---- reporting-bias views
        base_q = date_series(f)  # every report in the window, per quarter
        series = {"all reports": base_q,
                  "ranitidine": date_series(f, any_of(DRUG_FIELD, ["ranitidine"])),
                  "levofloxacin": date_series(f, any_of(DRUG_FIELD, ["levofloxacin"])),
                  "aortic aneurysm / dissection": date_series(f, any_of(PT_FIELD, AORTIC_PTS)),
                  "levofloxacin + aortic aneurysm / dissection": date_series(
                      f, any_of(DRUG_FIELD, ["levofloxacin"]), any_of(PT_FIELD, AORTIC_PTS))}
        for d in NEW_DRUGS:
            series[d] = date_series(f, any_of(DRUG_FIELD, [d]))
        ts = pd.DataFrame(series).fillna(0).astype(int)
        ts.index = ts.index.astype(str)
        ts.index.name = "quarter"
        mix_all = f.table("primarysource.qualification").counts
        mix_ran = f.table("primarysource.qualification", any_of(DRUG_FIELD, ["ranitidine"])).counts
        mix = pd.DataFrame([{"reporter": QUALIFICATION.get(k, k), "all_reports": mix_all.get(k, 0),
                             "ranitidine_reports": mix_ran.get(k, 0)} for k in QUALIFICATION])
        # Who reported the drug-event pairs that stand out in the time series and the screen (chosen after
        # reading those results).
        checks = []
        for drug, pts in REPORTER_CHECKS:
            q = f.table("primarysource.qualification", any_of(DRUG_FIELD, [drug]), any_of(PT_FIELD, pts)).counts
            checks.append({"drug": drug, "event": " / ".join(pts).lower(), "reports": sum(q.values()),
                           **{QUALIFICATION[k]: q.get(k, 0) for k in QUALIFICATION}})
        pd.DataFrame(checks).to_csv(out_dir / "reporter_mix_checks.csv", index=False)
        meta_extra["reports_with_reportduplicate_block"] = f.total("_exists_:reportduplicate")
        # Can a drug's role be tied to the drug in a search? If it could, "aspirin AND a suspect drug" would match
        # only the reports where aspirin itself is suspect, a fraction of all aspirin reports.
        aspirin = any_of(DRUG_FIELD, ["aspirin"])
        meta_extra["role_check_aspirin"] = {
            "aspirin_reports": f.total(aspirin),
            "aspirin_and_any_suspect_drug": f.total(aspirin, "patient.drug.drugcharacterization:1")}
        meta_extra["openfda_reports_all_years"] = client.get({"limit": 1})["data"]["meta"]["results"]["total"]
        # Why the window starts in 2015: share of reports per year with the active-ingredient field filled in.
        has_ai, _ = client.count(f"_exists_:{DRUG_FIELD}", "receivedate")
        every, _ = client.count(None, "receivedate")
        cov = pd.DataFrame({"with_active_ingredient": pd.Series(has_ai), "all_reports": pd.Series(every)}).fillna(0)
        cov = cov.groupby(cov.index.str[:4]).sum().astype(int).loc["2004":]
        cov["share"] = cov.with_active_ingredient / cov.all_reports
        cov.rename_axis("year").to_csv(out_dir / "active_ingredient_coverage.csv")

        # ---- per-outcome AUCs: OMOP pairs for its 4 outcomes, pooled unique pairs for all 10
        by_out = []
        ev_pairs = scored[scored.evaluable].copy()
        ev_pairs[["prr", "prr025", "ror", "ror025"]] = ev_pairs[["prr", "prr025", "ror", "ror025"]].fillna(0.0)
        for name, d0 in [("OMOP", ev_pairs[ev_pairs.set == "omop"]),
                         ("Pooled", ev_pairs.drop_duplicates(["drug", "outcome"]))]:
            for o, d in d0.groupby("outcome"):
                if d.ground_truth.nunique() == 2:
                    by_out.append(auc_table(d, LOWER, seed=SEED).assign(set=name, outcome=o))
        pd.concat(by_out).to_csv(out_dir / "auc_by_outcome.csv", index=False)
        # ---- sparse data: pooled AUC by how many reports the drug has
        pooled = ev_pairs.drop_duplicates(["drug", "outcome"])
        subsets = [("all evaluable pairs", pooled), ("drug has >= 1,000 reports", pooled[pooled.n_drug >= 1000]),
                   ("drug has < 1,000 reports", pooled[pooled.n_drug < 1000]),
                   ("at least 5 co-reports", pooled[pooled.a >= 5])]
        by_size = pd.concat([auc_table(d, LOWER, seed=SEED).assign(subset=name, pairs=len(d)) for name, d in subsets])
        by_size.to_csv(out_dir / "auc_by_drug_size.csv", index=False)

        tb.to_csv(out_dir / "textbook_pairs.csv", index=False)
        tb_cmp.to_csv(out_dir / "textbook_vs_comparators.csv", index=False)
        screen.sort_values(["drug", "rank_eb05"]).to_csv(out_dir / "screening_pts.csv", index=False)
        val.to_csv(out_dir / "validation_counts.csv", index=False)
        ts.to_csv(out_dir / "time_series_quarterly.csv")
        mix.to_csv(out_dir / "reporter_mix.csv", index=False)
        pd.Series(NEW_DRUGS, name="approval_date").rename_axis("drug").to_csv(out_dir / "new_drug_approvals.csv")

    # ---- write the core tables
    keep = ["set", "drug_name", "drug", "outcome", "ground_truth", "indication", "in_both", "evaluable", "a",
            "a_method", "a_shortfall_max", "n_drug", "n_event", "expected", "prr", "prr025", "prr975", "ror",
            "ror025", "ror975", "chi2", "ic", "ic025", "ic975", "ebgm", "eb05", "eb95"]
    scored[keep].to_csv(out_dir / "reference_pairs.csv", index=False)
    data["marginals"].to_csv(out_dir / "drug_marginals.csv", index=False)
    terms = pd.DataFrame([{"drug": d, "term": t, "reports": c} for d, tt in data["terms_of"].items() if d in needed
                          for t, c in sorted(tt.items(), key=lambda x: -x[1])])
    terms.to_csv(out_dir / "drug_terms.csv", index=False)
    outcome_defs.to_csv(out_dir / "outcome_definitions.csv", index=False)
    pool.to_csv(out_dir / "gps_pool.csv", index=False)
    aucs.to_csv(out_dir / "auc.csv", index=False)
    tests.to_csv(out_dir / "auc_tests.csv", index=False)
    thr.to_csv(out_dir / "signal_thresholds.csv", index=False)
    rocs.to_csv(out_dir / "roc_points.csv", index=False)
    corr = scored.loc[scored.evaluable, SCORES].fillna(0.0).corr(method="spearman")
    corr.to_csv(out_dir / "score_rank_correlation.csv")

    log = [json.loads(x) for x in (RAW / "openfda_request_log.jsonl").read_text(encoding="utf-8").splitlines()]
    fetched = sorted({t[:10] for t in client.urls_used.values()})
    meta = {
        "disclaimer": DISCLAIMER,
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "quick": args.quick,
        "source": "openFDA drug adverse event endpoint (FAERS), https://api.fda.gov/drug/event.json",
        "openfda_last_updated": client.get({"search": window_clause(), "limit": 1})["data"]["meta"]["last_updated"],
        "window_receivedate": list(WINDOW),
        "counting_unit": "report (openFDA record = latest version of a FAERS safety report ID); not patients, "
                         "not drug-event pairs; duplicate reports of one case from different senders are not merged",
        "drug_field": DRUG_FIELD, "event_field": PT_FIELD, "count_limit": COUNT_LIMIT,
        "N_reports": data["N"],
        "n_reference_pairs": {k: int(v) for k, v in ref.groupby("set").size().items()},
        "n_evaluable_pairs": {k: int(v) for k, v in scored[scored.evaluable].groupby("set").size().items()},
        "a_method_counts": dict(Counter(scored.a_method)),
        "gps_prior": prior.as_dict(), "gps_pool_pairs": len(pool),
        "requests": {"queries_used_by_this_run": len(client.urls_used),
                     "network_calls_this_run": client.network_calls, "cache_hits_this_run": client.cache_hits,
                     "network_calls_logged_total": len(log),
                     "log_days": sorted({x["fetched_utc"][:10] for x in log})},
        "fetch_dates": fetched,
        "runtime_s": round(time.time() - t0, 1),
        "python": platform.python_version(),
        **meta_extra,
    }
    (out_dir / "run_metadata.json").write_text(json.dumps(meta, indent=2, default=float), encoding="utf-8")
    prov = pd.DataFrame({"url": list(client.urls_used), "fetched_utc": list(client.urls_used.values())})
    detail = pd.DataFrame(f.provenance).drop_duplicates("url")[["url", "kind", "value", "floor"]]
    prov.merge(detail, on="url", how="left").to_csv(out_dir / "query_provenance.csv", index=False)
    print(json.dumps({k: meta[k] for k in ["N_reports", "n_evaluable_pairs", "a_method_counts", "requests"]}, indent=1))
    print(aucs[aucs.score.isin(LOWER + ["prr", "ebgm"])].round(3).to_string(index=False))
    if not args.quick:
        for p in make_all(RESULTS, FIGS):
            print("wrote", p.relative_to(ROOT))
    print(f"done in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
