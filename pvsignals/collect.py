"""Build 2x2 report counts from openFDA count queries, handling the 999-term cap.

Counting unit: one FAERS report as indexed by openFDA (the latest version of each safety report ID). Drugs are
identified by FAERS's active-ingredient field (``patient.drug.activesubstance.activesubstancename``), on any role
(suspect or concomitant): openFDA cannot tie a role to a particular drug entry inside a search. Events are MedDRA
preferred terms in ``patient.reaction.reactionmeddrapt``.

Counts come from ``count=`` queries, which return the number of matching reports per exact term, at most 999 terms
per query without an API key. A drug usually has several exact terms (salts, combinations), so its count is the sum
over its terms; a report that lists two of a drug's terms is counted twice (measured in the validation step).
When a needed term falls below a table's cut-off, the same count query is re-run on the subset of reports that
mention the drugs in question; within that subset their terms rank high enough to be listed, and their counts are
exact because every report with the term is in the subset.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd

from pvsignals.openfda import COUNT_LIMIT, OpenFDA
from pvsignals.reference import DrugSpec, term_matches

DRUG_FIELD = "patient.drug.activesubstance.activesubstancename"
DRUG_EXACT = DRUG_FIELD + ".exact"
PT_FIELD = "patient.reaction.reactionmeddrapt.exact"
WINDOW = ("20150101", "20251231")


def window_clause(start: str = WINDOW[0], end: str = WINDOW[1]) -> str:
    """Reports received by FDA between two dates (inclusive, YYYYMMDD)."""
    return f"receivedate:[{start} TO {end}]"


def any_of(field_name: str, values: list[str] | tuple[str, ...]) -> str:
    """``field:("v1" "v2" ...)``: reports where the field matches any of the quoted values."""
    return f"{field_name}:(" + " ".join(f'"{v}"' for v in values) + ")"


def and_(*clauses: str) -> str:
    """Join search clauses with AND (empty clauses are skipped)."""
    return " AND ".join(c for c in clauses if c)


@dataclass
class CountTable:
    """One count query's result: report counts per exact term, with provenance."""

    counts: dict[str, int]
    url: str
    fetched_utc: str
    limit: int = COUNT_LIMIT

    @property
    def complete(self) -> bool:
        """True if openFDA returned every term (fewer than ``limit`` rows)."""
        return len(self.counts) < self.limit

    @property
    def floor(self) -> int:
        """Every term with more reports than this is listed; unlisted terms have at most this many (0 = complete)."""
        return 0 if self.complete else min(self.counts.values())


@dataclass
class Fetcher:
    """Count and total queries on top of a base filter (date window, optionally a reporter filter)."""

    client: OpenFDA
    base: str
    provenance: list[dict] = field(default_factory=list)

    def total(self, *clauses: str) -> int:
        n, rec = self.client.total(and_(self.base, *clauses))
        self.provenance.append({"kind": "total", "url": rec["url"], "fetched_utc": rec["fetched_utc"], "value": n})
        return n

    def table(self, count_field: str, *clauses: str) -> CountTable:
        counts, rec = self.client.count(and_(self.base, *clauses), count_field)
        tab = CountTable(counts, rec["url"], rec["fetched_utc"])
        self.provenance.append({"kind": "count", "url": rec["url"], "fetched_utc": rec["fetched_utc"],
                                "value": len(counts), "floor": tab.floor})
        return tab


def drug_terms(table: CountTable, spec: DrugSpec) -> dict[str, int]:
    """Terms of a count table that belong to a drug, with their counts."""
    return {t: c for t, c in table.counts.items() if term_matches(t, spec)}


def plan_batches(specs: list[DrugSpec], overall: CountTable, max_reports: int = 600_000,
                 max_drugs: int = 15) -> list[list[DrugSpec]]:
    """Group drugs for the term-discovery queries so that each subset of reports stays small.

    A drug's size is estimated from the overall table (or the table's floor when none of its terms is listed).
    The plan depends only on the drug list and the overall table, so re-runs ask the same (cached) queries.
    """
    sizes = []
    for s in specs:
        listed = sum(drug_terms(overall, s).values())
        sizes.append((max(listed, overall.floor), s.name, s))
    sizes.sort(key=lambda x: (-x[0], x[1]))
    batches: list[list[DrugSpec]] = []
    current: list[DrugSpec] = []
    load = 0
    for size, _, s in sizes:
        if current and (load + size > max_reports or len(current) >= max_drugs):
            batches.append(current)
            current, load = [], 0
        current.append(s)
        load += size
    if current:
        batches.append(current)
    return batches


def search_items(batch: list[DrugSpec]) -> list[str]:
    """De-duplicated search tokens of a batch of drugs, in a stable order."""
    items: list[str] = []
    for s in batch:
        for tok in s.search:
            if tok not in items:
                items.append(tok)
    return items


def drug_marginals(f: Fetcher, specs: list[DrugSpec], overall: CountTable,
                   batches: list[list[DrugSpec]]) -> tuple[pd.DataFrame, dict[str, dict[str, int]]]:
    """Reports per drug (sum over its terms) and each drug's term list, from batched discovery queries."""
    rows, terms_of = [], {}
    for i, batch in enumerate(batches):
        sub = f.table(DRUG_EXACT, any_of(DRUG_FIELD, search_items(batch)))
        for s in batch:
            terms = drug_terms(overall, s)
            terms.update(drug_terms(sub, s))  # exact counts; the subset lists rarer terms too
            terms_of[s.name] = terms
            rows.append({"drug": s.name, "n_drug": int(sum(terms.values())), "n_terms": len(terms),
                         "top_term": max(terms, key=terms.get) if terms else "", "batch": i,
                         "unlisted_term_max": sub.floor})
    return pd.DataFrame(rows), terms_of


def outcome_tables(f: Fetcher, events: dict[str, list[str]]) -> dict[str, dict]:
    """Per event: reports with any of its PTs (exact union), drug-term counts and PT counts."""
    out = {}
    for key, pts in events.items():
        clause = any_of(PT_FIELD, pts)
        rec = {"pts": pts, "n_event": f.total(clause), "drugs": f.table(DRUG_EXACT, clause),
               "pt_counts": f.table(PT_FIELD, clause)}
        out[key] = rec
    return out


def _shortfall(terms: dict[str, int], unlisted: list[str], floor: int) -> int:
    """Most reports the unlisted terms could add: each has at most ``floor`` co-reports and at most its own total."""
    return int(sum(min(terms[t], floor) for t in unlisted))


def _acceptable(shortfall: int, a: int) -> bool:
    """A table count is used as is when unlisted terms could add at most 3 reports or 1%; otherwise it is re-counted."""
    return shortfall <= max(3, 0.01 * a)


def pair_counts(f: Fetcher, pairs: pd.DataFrame, specs: dict[str, DrugSpec], terms_of: dict[str, dict[str, int]],
                events: dict[str, dict], max_tokens: int = 30) -> pd.DataFrame:
    """Co-report count ``a`` for each (drug, event) pair, filling in terms that fall below an event table's cut-off.

    ``pairs`` needs columns ``drug`` and ``event``. Returns one row per pair with ``a``, how it was obtained
    (``table``: from the event's count table; ``fill-in``: re-counted in a subset query) and ``a_shortfall_max``,
    the most reports that terms still unlisted could add (0 = exact up to the term summation).
    """
    counts: dict[tuple[str, str], int] = {}
    method: dict[tuple[str, str], str] = {}
    short: dict[tuple[str, str], int] = {}
    pending: dict[str, list[str]] = {}
    for drug, ev in pairs[["drug", "event"]].drop_duplicates().itertuples(index=False):
        tab: CountTable = events[ev]["drugs"]
        terms = terms_of.get(drug, {})
        missing = [t for t in terms if t not in tab.counts and not tab.complete]
        key = (drug, ev)
        counts[key] = sum(tab.counts.get(t, 0) for t in terms)
        method[key] = "table"
        short[key] = _shortfall(terms, missing, tab.floor)
        if not _acceptable(short[key], counts[key]):
            pending.setdefault(ev, []).append(drug)
    for ev, drugs in pending.items():
        still = _fill(f, ev, sorted(set(drugs)), specs, terms_of, events, counts, method, short, max_tokens)
        if still:  # the batched subset was still too large to list these drugs' terms: one drug per query
            _fill(f, ev, still, specs, terms_of, events, counts, method, short, max_tokens=1)
    rows = [{"drug": d, "event": e, "a": counts[(d, e)], "a_method": method[(d, e)],
             "a_shortfall_max": short[(d, e)]} for (d, e) in counts]
    return pd.DataFrame(rows)


def _fill(f: Fetcher, ev: str, drugs: list[str], specs: dict[str, DrugSpec], terms_of: dict[str, dict[str, int]],
          events: dict[str, dict], counts: dict, method: dict, short: dict, max_tokens: int) -> list[str]:
    """Re-count the terms of ``drugs`` within the reports of event ``ev`` that mention them, in groups of at most
    ``max_tokens`` search tokens (a single drug per query when ``max_tokens`` is 1).

    Returns the drugs whose count could still be short by more than the tolerance; terms left unlisted keep their
    possible contribution in ``short``, so the analysis can be repeated with the upper bound."""
    clause = any_of(PT_FIELD, events[ev]["pts"])
    groups, cur = [], []
    for d in drugs:
        if cur and len(search_items(cur + [specs[d]])) > max_tokens:
            groups.append(cur)
            cur = []
        cur.append(specs[d])
    if cur:
        groups.append(cur)
    still = []
    for g in groups:
        sub = f.table(DRUG_EXACT, clause, any_of(DRUG_FIELD, search_items(g)))
        for s in g:
            terms = terms_of.get(s.name, {})
            unlisted = [t for t in terms if t not in sub.counts and not sub.complete]
            key = (s.name, ev)
            counts[key] = int(sum(sub.counts.get(t, 0) for t in terms))
            method[key] = "fill-in" if len(g) > 1 or max_tokens > 1 else "fill-in (single drug)"
            short[key] = _shortfall(terms, unlisted, sub.floor)
            if not _acceptable(short[key], counts[key]):
                still.append(s.name)
    return still
