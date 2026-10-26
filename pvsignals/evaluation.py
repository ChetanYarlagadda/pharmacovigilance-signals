"""Discrimination of positive from negative controls: ROC AUC with bootstrap intervals, paired tests, thresholds."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

# Signal rules as commonly applied in practice.
SIGNAL_RULES = {
    "PRR (Evans: PRR >= 2, chi2 >= 4, a >= 3)": lambda d: (d["prr"] >= 2) & (d["chi2"] >= 4) & (d["a"] >= 3),
    "ROR lower 95% > 1 (a >= 3)": lambda d: (d["ror025"] > 1) & (d["a"] >= 3),
    "IC lower 95% > 0": lambda d: d["ic025"] > 0,
    "EB05 >= 2": lambda d: d["eb05"] >= 2,
}


def auc(y: np.ndarray, score: np.ndarray) -> float:
    """ROC AUC as the Mann-Whitney probability that a positive outscores a negative (ties count one half)."""
    y = np.asarray(y).astype(bool)
    s = np.asarray(score, dtype=float)
    pos, neg = s[y], s[~y]
    if pos.size == 0 or neg.size == 0:
        return float("nan")
    ranks = stats.rankdata(np.r_[pos, neg])
    return float((ranks[: pos.size].sum() - pos.size * (pos.size + 1) / 2) / (pos.size * neg.size))


def roc_curve(y: np.ndarray, score: np.ndarray) -> pd.DataFrame:
    """False- and true-positive rates at every distinct threshold (ties handled as one step)."""
    y = np.asarray(y).astype(bool)
    s = np.asarray(score, dtype=float)
    thr = np.r_[np.inf, np.unique(s)[::-1]]
    tpr = [(s[y] >= t).mean() for t in thr]
    fpr = [(s[~y] >= t).mean() for t in thr]
    return pd.DataFrame({"threshold": thr, "fpr": fpr, "tpr": tpr})


def _boot_indices(y: np.ndarray, n_boot: int, seed: int) -> list[np.ndarray]:
    """Stratified bootstrap resamples (positives and negatives resampled separately)."""
    rng = np.random.default_rng(seed)
    pos, neg = np.flatnonzero(y), np.flatnonzero(~y)
    return [np.r_[rng.choice(pos, pos.size), rng.choice(neg, neg.size)] for _ in range(n_boot)]


def auc_table(df: pd.DataFrame, scores: list[str], label: str = "ground_truth", n_boot: int = 2000,
              seed: int = 0) -> pd.DataFrame:
    """AUC of each score with a percentile bootstrap 95% CI."""
    y = df[label].to_numpy().astype(bool)
    idx = _boot_indices(y, n_boot, seed)
    rows = []
    for s in scores:
        x = df[s].to_numpy(dtype=float)
        boots = np.array([auc(y[i], x[i]) for i in idx])
        rows.append({"score": s, "auc": auc(y, x), "ci_low": np.percentile(boots, 2.5),
                     "ci_high": np.percentile(boots, 97.5), "n_pos": int(y.sum()), "n_neg": int((~y).sum())})
    return pd.DataFrame(rows)


def holm(p: np.ndarray) -> np.ndarray:
    """Holm step-down adjusted p-values."""
    p = np.asarray(p, dtype=float)
    order = np.argsort(p)
    adj = np.empty_like(p)
    running = 0.0
    for k, i in enumerate(order):
        running = max(running, (p.size - k) * p[i])
        adj[i] = min(1.0, running)
    return adj


def paired_auc_tests(df: pd.DataFrame, comparisons: list[tuple[str, str]], label: str = "ground_truth",
                     n_boot: int = 2000, seed: int = 1) -> pd.DataFrame:
    """Paired bootstrap of AUC differences (same resamples for both scores), two-sided p, Holm-adjusted."""
    y = df[label].to_numpy().astype(bool)
    idx = _boot_indices(y, n_boot, seed)
    rows = []
    for s1, s2 in comparisons:
        x1, x2 = df[s1].to_numpy(dtype=float), df[s2].to_numpy(dtype=float)
        diff = auc(y, x1) - auc(y, x2)
        boots = np.array([auc(y[i], x1[i]) - auc(y[i], x2[i]) for i in idx])
        p = 2 * min((boots <= 0).mean(), (boots >= 0).mean())
        rows.append({"score_1": s1, "score_2": s2, "auc_diff": diff, "ci_low": np.percentile(boots, 2.5),
                     "ci_high": np.percentile(boots, 97.5), "p": max(p, 1 / n_boot)})
    out = pd.DataFrame(rows)
    out["p_holm"] = holm(out["p"].to_numpy())
    return out


def threshold_table(df: pd.DataFrame, label: str = "ground_truth") -> pd.DataFrame:
    """Sensitivity, specificity and PPV of each conventional signal rule."""
    y = df[label].to_numpy().astype(bool)
    rows = []
    for name, rule in SIGNAL_RULES.items():
        flag = rule(df).fillna(False).to_numpy().astype(bool)
        tp, fp = int((flag & y).sum()), int((flag & ~y).sum())
        rows.append({"rule": name, "flagged": int(flag.sum()), "true_pos": tp, "false_pos": fp,
                     "sensitivity": tp / max(y.sum(), 1), "specificity": 1 - fp / max((~y).sum(), 1),
                     "ppv": tp / max(flag.sum(), 1)})
    return pd.DataFrame(rows)
