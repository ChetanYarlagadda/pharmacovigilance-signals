"""All report figures, drawn from the result tables in reports/results/."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pvsignals import plotting
from pvsignals.plotting import (AXIS, INK, INK_2, METHOD_COLOR, MUTED, NEGATIVE, POSITIVE, SCORE_LABEL, SCORE_METHOD,
                                SERIES)
from pvsignals.reference import OUTCOME_LABEL

LOWER = ["prr025", "ror025", "ic025", "eb05"]
SETS = ["OMOP", "EU-ADR"]
STRATA = {"all reports": "Counts as fetched",
          "all reports, unlisted terms at upper bound": "Unlisted terms at their maximum (sensitivity)"}
EVENT_LABEL = {**OUTCOME_LABEL, "tendon_rupture": "Tendon rupture", "depression": "Depression"}

# Takeaway titles, written after reading the results.
TITLES = {
    "roc": "PRR, ROR, BCPNN IC and EBGM rank the reference controls almost identically:\n"
           "AUC about 0.70 on the OMOP set and 0.82-0.83 on EU-ADR",
    "auc": "Point estimates and lower bounds discriminate equally well (paired differences < 0.005, Holm p = 1);\n"
           "counting every unlisted term at its maximum lowers AUC by at most 0.01",
    "textbook": "The three textbook signals stand out against all reports, less so against their own therapeutic area:\n"
                "levofloxacin-tendon rupture ROR 41 (26 vs other antibiotics), isotretinoin-depression 9.4 "
                "(2.6 vs other acne drugs)",
    "shrinkage": "Shrinkage matters only for small counts: EBGM pulls pairs with a handful of co-reports toward 1,\n"
                 "while ratios backed by hundreds of reports stay at their observed-to-expected value",
    "stimulated": "Reporting waves follow news and litigation: ranitidine reports rose 18-fold after the 2019 NDMA alert "
                  "(41% from lawyers),\nand a burst of consumer reports in early 2017 pushed the levofloxacin-aortic "
                  "aneurysm PRR from about 10 to 448",
    "weber": "No common post-launch peak: only erenumab's reports peak early and then fall,\n"
             "while dupilumab and tirzepatide keep rising with uptake (dupilumab is on 8.8% of all 2025 reports)",
    "outcomes": "Detectability depends on the outcome: acute myocardial infarction is barely separable (OMOP AUC 0.59),\n"
                "kidney injury and upper GI bleeding are the easiest of the OMOP outcomes (0.81)",
}


def _suptitle(fig: plt.Figure, text: str, y: float = 0.99) -> None:
    """Left-aligned takeaway title placed above the axes (it grows upward, so two lines never cover the plot)."""
    fig.suptitle(text, x=0.01, y=y, ha="left", va="bottom", fontsize=12, fontweight="bold", color=INK)


def _subtitle(ax: plt.Axes, text: str) -> None:
    ax.set_title(text, loc="left", fontsize=9, color=MUTED, fontweight="normal")


def fig_roc(roc: pd.DataFrame, auc: pd.DataFrame, path: Path, title: str) -> Path:
    """ROC curves of the four lower-bound scores on each reference set (all reports)."""
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 5.2), sharey=True)
    for ax, name in zip(axes, SETS):
        a = auc[(auc.set == name) & (auc.stratum == "all reports")].set_index("score")
        for s in LOWER:
            r = roc[(roc.set == name) & (roc.stratum == "all reports") & (roc.score == s)]
            ax.step(r.fpr, r.tpr, where="post", color=METHOD_COLOR[SCORE_METHOD[s]], lw=1.8,
                    label=f"{SCORE_LABEL[s]}: {a.loc[s, 'auc']:.2f} [{a.loc[s, 'ci_low']:.2f}, {a.loc[s, 'ci_high']:.2f}]")
        ax.plot([0, 1], [0, 1], color=AXIS, lw=1, ls="--")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1.01)
        ax.set_aspect("equal")
        ax.set_xlabel("False-positive rate\n(share of negative controls flagged)")
        n_pos, n_neg = int(a.iloc[0].n_pos), int(a.iloc[0].n_neg)
        _subtitle(ax, f"{name}: {n_pos} positive, {n_neg} negative controls")
        ax.legend(loc="lower right", fontsize=8, handlelength=1.4, title="AUC [95% CI]", title_fontsize=8)
    axes[0].set_ylabel("True-positive rate\n(share of positive controls flagged)")
    _suptitle(fig, title)
    fig.tight_layout()
    plotting.save(fig, path)
    return path


def fig_auc(auc: pd.DataFrame, path: Path, title: str) -> Path:
    """AUC with 95% bootstrap CIs for every score, reference set and stratum."""
    scores = ["prr", "prr025", "ror", "ror025", "ic", "ic025", "ebgm", "eb05"]
    sets = ["OMOP", "EU-ADR", "Pooled"]
    fig, axes = plt.subplots(1, 3, figsize=(11, 4.6), sharey=True)
    y = np.arange(len(scores))[::-1]
    for ax, name in zip(axes, sets):
        for k, stratum in enumerate(STRATA):
            a = auc[(auc.set == name) & (auc.stratum == stratum)].set_index("score").reindex(scores)
            off = 0.17 if k == 0 else -0.17
            for yi, s in zip(y, scores):
                c = METHOD_COLOR[SCORE_METHOD[s]]
                ax.errorbar(a.loc[s, "auc"], yi + off, xerr=[[a.loc[s, "auc"] - a.loc[s, "ci_low"]],
                                                               [a.loc[s, "ci_high"] - a.loc[s, "auc"]]],
                            fmt="o" if k == 0 else "D", color=c, mfc=c if k == 0 else "white", ms=5, lw=1.2)
        ax.axvline(0.5, color=AXIS, lw=1, ls="--")
        ax.set_xlim(0.4, 1.0)
        n = auc[(auc.set == name) & (auc.stratum == "all reports")].iloc[0]
        _subtitle(ax, f"{name} ({int(n.n_pos)} positive / {int(n.n_neg)} negative)")
        ax.set_xlabel("ROC AUC (95% bootstrap CI)")
    axes[0].set_yticks(y, [SCORE_LABEL[s] for s in scores])
    labels = list(STRATA.values())
    handles = [plt.Line2D([], [], marker="o", color=INK_2, ls="none", ms=5, label=labels[0]),
               plt.Line2D([], [], marker="D", color=INK_2, mfc="white", ls="none", ms=5, label=labels[1])]
    fig.legend(handles=handles, loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.05), fontsize=9)
    _suptitle(fig, title)
    fig.tight_layout()
    plotting.save(fig, path)
    return path


def fig_textbook(tb: pd.DataFrame, path: Path, title: str) -> Path:
    """ROR (95% CI) and EBGM (EB05-EB95) for the textbook signals and their therapeutic-area comparators."""
    events = ["tendon_rupture", "rhabdomyolysis", "depression"]
    heights = [tb[tb.event == e].shape[0] for e in events]
    fig, axes = plt.subplots(len(events), 1, figsize=(9.5, 1.2 + 0.36 * sum(heights) + 0.5 * len(events)),
                             gridspec_kw={"height_ratios": heights}, sharex=True)
    for ax, ev in zip(axes, events):
        d = tb[tb.event == ev].reset_index(drop=True)
        groups = list(dict.fromkeys(d.group))
        d = pd.concat([d[d.group == groups[0]], d[d.group != groups[0]]]).reset_index(drop=True)
        y = np.arange(len(d))[::-1]
        for yi, (_, r) in zip(y, d.iterrows()):
            c = POSITIVE if r.group == groups[0] else NEGATIVE
            if r.a > 0:
                ax.plot([r.ror025, r.ror975], [yi + 0.15] * 2, color=c, lw=1.4)
                ax.plot(r.ror, yi + 0.15, "o", color=c, ms=5)
            ax.plot([r.eb05, r.eb95], [yi - 0.15] * 2, color=c, lw=1.4, alpha=0.6)
            ax.plot(r.ebgm, yi - 0.15, "D", mfc="white", mec=c, ms=5)
            ax.text(1.01, yi, f"{int(r.a):,} / {int(r.n_drug):,}", transform=ax.get_yaxis_transform(), va="center",
                    fontsize=7.5, color=MUTED)
        ax.set_yticks(y, [r.drug if r.drug == r.group else f"{r.drug} ({r.group})" for _, r in d.iterrows()],
                      fontsize=8.5)
        ax.axvline(1, color=AXIS, lw=1, ls="--")
        ax.set_xscale("log")
        ticks = [0.5, 1, 2, 5, 10, 20, 50]
        ax.set_xticks(ticks, [f"{t:g}" for t in ticks])
        ax.minorticks_off()
        _subtitle(ax, f"{EVENT_LABEL[ev]}: {groups[0]} (red) vs {groups[1]} (blue)")
    axes[-1].set_xlabel("Disproportionality against all other reports (log scale):\n"
                        "ROR with 95% CI (dots), EBGM with EB05-EB95 (diamonds)")
    axes[0].text(1.01, 1.0, "reports with\nevent / drug", transform=axes[0].transAxes, fontsize=7.5, color=MUTED,
                 va="bottom")
    _suptitle(fig, title)
    fig.tight_layout()
    plotting.save(fig, path)
    return path


def fig_shrinkage(pool: pd.DataFrame, pairs: pd.DataFrame, path: Path, title: str) -> Path:
    """Observed-to-expected ratio vs EBGM: the pairs that fitted the GPS prior and the reference pairs."""
    fig, ax = plt.subplots(figsize=(7.5, 5.8))
    oe = pool.n / pool.expected
    ax.scatter(oe, pool.ebgm, s=5, color=AXIS, alpha=0.5, lw=0,
               label=f"prior-fitting pairs ({len(pool):,} term x event pairs)")
    ref = pairs[pairs.evaluable & (pairs.a > 0)].drop_duplicates(["drug", "outcome"])
    roe = ref.a / ref.expected
    bins = [(1, 5, "1-4 co-reports", SERIES[7]), (5, 20, "5-19 co-reports", SERIES[3]),
            (20, 200, "20-199 co-reports", SERIES[2]), (200, np.inf, "200 or more co-reports", SERIES[0])]
    for lo, hi, label, c in bins:
        m = (ref.a >= lo) & (ref.a < hi)
        ax.scatter(roe[m], ref.ebgm[m], s=16, alpha=0.85, color=c, lw=0,
                   label=f"reference pairs, {label} ({m.sum()})")
    lim = [min(oe.min(), roe.min(), ref.ebgm.min()) * 0.8, max(oe.max(), roe.max(), ref.ebgm.max()) * 1.2]
    ax.plot(lim, lim, color=AXIS, lw=1, ls="--", label="no shrinkage (EBGM = O/E)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel("Observed / expected co-reports (O/E)")
    ax.set_ylabel("EBGM (empirical Bayes geometric mean)")
    _subtitle(ax, "Reference drug-outcome pairs with at least one co-report (unique pairs)")
    ax.legend(loc="upper left", fontsize=7.5, markerscale=1.5)
    _suptitle(fig, title)
    fig.tight_layout()
    plotting.save(fig, path)
    return path


def _quarter_x(index: pd.Index) -> np.ndarray:
    p = pd.PeriodIndex(index, freq="Q")
    return np.asarray(p.year + (p.quarter - 1) / 4 + 0.125, dtype=float)


def fig_stimulated(ts: pd.DataFrame, mix: pd.DataFrame, path: Path, title: str) -> Path:
    """Reporting after regulatory news: ranitidine volume and reporter mix, levofloxacin and aortic aneurysm."""
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4), gridspec_kw={"width_ratios": [1.3, 0.85, 1.3]})
    x = _quarter_x(ts.index)
    ax = axes[0]
    ax.bar(x, ts["ranitidine"] / 1e3, width=0.22, color=SERIES[1])
    top = ax.get_ylim()[1]
    for when in [2019 + 8.4 / 12, 2020 + 3 / 12]:
        ax.axvline(when, color=INK_2, lw=0.9, ls=":")
    ax.text(2019 + 8.4 / 12 - 0.15, top * 0.97, "Sep 2019: FDA\nNDMA alert", fontsize=7.5, color=INK_2, ha="right",
            va="top")
    ax.text(2020 + 3 / 12 + 0.15, top * 0.97, "Apr 2020: market\nwithdrawal request", fontsize=7.5, color=INK_2,
            ha="left", va="top")
    ax.set_ylabel("Reports per quarter (thousands)")
    ax.set_xlabel("Quarter received by FDA")
    _subtitle(ax, "Reports listing ranitidine")
    ax = axes[1]
    m = mix.set_index("reporter")[["all_reports", "ranitidine_reports"]]
    shares = m / m.sum()
    order = ["Physician", "Pharmacist", "Other health professional", "Consumer", "Lawyer"]
    colors = [SERIES[0], SERIES[2], SERIES[6], SERIES[3], SERIES[7]]
    left = np.zeros(2)
    for r, c in zip(order, colors):
        v = shares.loc[r].to_numpy()
        ax.barh([1, 0], v, left=left, color=c, height=0.6, label=r)
        for yi, vi, li in zip([1, 0], v, left):
            if vi > 0.07:
                ax.text(li + vi / 2, yi, f"{vi:.0%}", ha="center", va="center", fontsize=7.5, color="white")
        left += v
    ax.set_yticks([1, 0], ["All reports", "Ranitidine"])
    ax.set_xlim(0, 1)
    ax.set_xlabel("Share of reports by reporter type")
    ax.grid(False)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2, fontsize=7.5)
    _subtitle(ax, "Who sent the reports (2015-2025)")
    ax = axes[2]
    yr = ts.groupby(pd.PeriodIndex(ts.index, freq="Q").year).sum()
    a = yr["levofloxacin + aortic aneurysm / dissection"]
    nd, ne, n = yr["levofloxacin"], yr["aortic aneurysm / dissection"], yr["all reports"]
    c = ne - a
    d = n - nd - c
    prr = (a / nd) / (c / (c + d))
    se = np.sqrt(1 / a - 1 / nd + 1 / c - 1 / (c + d))
    ax.fill_between(yr.index, prr * np.exp(-1.96 * se), prr * np.exp(1.96 * se), color=SERIES[0], alpha=0.18, lw=0)
    ax.plot(yr.index, prr, "o-", color=SERIES[0], ms=4)
    ax.set_yscale("log")
    ax.set_ylim(0.3, 3000)
    for xi, ai in zip(yr.index, a):
        ax.text(xi, 0.4, f"{int(ai)}", ha="center", va="bottom", fontsize=6.5, color=MUTED)
    ax.vlines(2018 + 11.6 / 12, 0.7, 3000, color=INK_2, lw=0.9, ls=":")
    ax.text(2018 + 11.6 / 12 + 0.15, 2000, "Dec 2018:\nFDA warning", fontsize=7.5, color=INK_2, ha="left", va="top")
    ax.axhline(1, color=AXIS, lw=1, ls="--")
    ax.set_xlabel("Year received (small numbers: co-reports)")
    ax.set_ylabel("PRR (95% CI, log scale)")
    _subtitle(ax, "Levofloxacin and aortic aneurysm / dissection")
    _suptitle(fig, title)
    fig.tight_layout()
    plotting.save(fig, path)
    return path


def fig_weber(ts: pd.DataFrame, approvals: pd.DataFrame, path: Path, title: str) -> Path:
    """Reports per quarter for drugs approved inside the window, against years since approval."""
    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    x = _quarter_x(ts.index)
    for (drug, date), c in zip(approvals.itertuples(index=False), [SERIES[0], SERIES[1], SERIES[2], SERIES[6]]):
        t0 = pd.Timestamp(date)
        since = x - (t0.year + (t0.dayofyear - 1) / 365.25)
        share = (ts[drug] / ts["all reports"] * 1000).to_numpy()
        m = since >= -0.25
        ax.plot(since[m], share[m], color=c, lw=1.8, label=f"{drug} (approved {t0:%b %Y})")
    ax.set_xlabel("Years since FDA approval")
    ax.set_ylabel("Reports listing the drug\nper 1,000 FAERS reports that quarter")
    ax.axvline(2, color=AXIS, lw=1, ls="--")
    ax.legend(fontsize=8.5, loc="upper left")
    _subtitle(ax, "Quarterly reports 2015-2025; without prescription volumes, uptake and reporting cannot be separated")
    _suptitle(fig, title)
    fig.tight_layout()
    plotting.save(fig, path)
    return path


def fig_outcomes(pairs: pd.DataFrame, by_outcome: pd.DataFrame, path: Path, title: str) -> Path:
    """EB05 of positive and negative controls per outcome, with the OMOP per-outcome AUC of EB05."""
    d = pairs[pairs.evaluable].drop_duplicates(["drug", "outcome"])
    pooled = by_outcome[(by_outcome.score == "eb05") & (by_outcome.set == "Pooled")].set_index("outcome")
    auc = by_outcome[(by_outcome.score == "eb05") & (by_outcome.set == "OMOP")].set_index("outcome")
    order = pooled.auc.reindex(d.outcome.unique()).fillna(-1).sort_values().index.tolist()
    fig, ax = plt.subplots(figsize=(10, 6.0))
    rng = np.random.default_rng(0)
    for i, o in enumerate(order):
        for truth, c, off in [(1, POSITIVE, 0.18), (0, NEGATIVE, -0.18)]:
            v = d[(d.outcome == o) & (d.ground_truth == truth)].eb05
            ax.scatter(v, i + off + rng.uniform(-0.08, 0.08, len(v)), s=14, color=c, alpha=0.75, lw=0,
                       label=("positive control" if truth else "negative control") if i == 0 else None)
        n_pos = int((d[d.outcome == o].ground_truth == 1).sum())
        n_neg = int((d[d.outcome == o].ground_truth == 0).sum())
        lab = f"{OUTCOME_LABEL[o]} ({n_pos} + / {n_neg} -)"
        if o in pooled.index:
            lab += f"\nAUC {pooled.loc[o, 'auc']:.2f}"
        if o in auc.index:
            lab += f", OMOP only {auc.loc[o, 'auc']:.2f} [{auc.loc[o, 'ci_low']:.2f}, {auc.loc[o, 'ci_high']:.2f}]"
        ax.text(-0.01, i, lab, transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=8.5, color=INK_2)
    ax.set_yticks([])
    ax.set_xscale("log")
    ax.axvline(2, color=AXIS, lw=1, ls="--")
    ax.set_xlabel("EB05 (lower 5% of the EBGM posterior, log scale); dashed line: EB05 = 2")
    ax.legend(loc="lower right", fontsize=8.5)
    _subtitle(ax, "Each dot is one reference drug-outcome pair (both sets, unique pairs); AUC of EB05 per outcome")
    _suptitle(fig, title)
    fig.tight_layout()
    plotting.save(fig, path)
    return path


def make_all(results: Path, figs: Path) -> list[Path]:
    """Draw every figure from the CSVs in ``results``."""
    plotting.apply_style()
    auc = pd.read_csv(results / "auc.csv")
    roc = pd.read_csv(results / "roc_points.csv")
    ts = pd.read_csv(results / "time_series_quarterly.csv", index_col=0)
    return [
        fig_roc(roc, auc, figs / "01_roc.png", TITLES["roc"]),
        fig_auc(auc, figs / "02_auc.png", TITLES["auc"]),
        fig_textbook(pd.read_csv(results / "textbook_pairs.csv"), figs / "03_textbook_signals.png", TITLES["textbook"]),
        fig_shrinkage(pd.read_csv(results / "gps_pool.csv"), pd.read_csv(results / "reference_pairs.csv"),
                      figs / "04_shrinkage.png", TITLES["shrinkage"]),
        fig_stimulated(ts, pd.read_csv(results / "reporter_mix.csv"), figs / "05_stimulated_reporting.png",
                       TITLES["stimulated"]),
        fig_weber(ts, pd.read_csv(results / "new_drug_approvals.csv"), figs / "06_new_drugs.png", TITLES["weber"]),
        fig_outcomes(pd.read_csv(results / "reference_pairs.csv"), pd.read_csv(results / "auc_by_outcome.csv"),
                     figs / "07_outcomes.png", TITLES["outcomes"]),
    ]
