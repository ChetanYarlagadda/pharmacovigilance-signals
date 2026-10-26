"""A cached, rate-limited client for the openFDA drug adverse event endpoint (FAERS).

openFDA allows about 240 requests per minute and 1,000 per day per IP address without an API key. Every response is
stored in ``data/raw/openfda/`` together with the exact query URL and the UTC fetch time, so the analysis can be
re-run without touching the network, and every network call is appended to a request log that the pipeline reads
to report how many requests a run actually made.
"""

from __future__ import annotations

import hashlib
import json
import time
from datetime import datetime, timedelta, timezone
from pathlib import Path

import requests

API = "https://api.fda.gov/drug/event.json"
COUNT_LIMIT = 999  # without an API key, openFDA refuses count queries with limit=1000 (API_KEY_MISSING); 999 works


class BudgetExceeded(RuntimeError):
    """Raised before a network call that would exceed the self-imposed daily request budget."""


def build_url(params: dict[str, str | int]) -> str:
    """The fully encoded request URL for a parameter dict (used as the cache key and stored with the response)."""
    return requests.Request("GET", API, params=params).prepare().url


def cache_key(url: str) -> str:
    """File-name-safe key of a query URL."""
    return hashlib.sha1(url.encode("utf-8")).hexdigest()[:20]


class OpenFDA:
    """Fetch openFDA ``drug/event`` queries with an on-disk cache, a minimum interval and a daily request budget."""

    def __init__(self, raw_dir: Path, daily_budget: int = 600, min_interval_s: float = 0.4,
                 offline: bool = False) -> None:
        self.dir = Path(raw_dir) / "openfda"
        self.dir.mkdir(parents=True, exist_ok=True)
        self.log_path = Path(raw_dir) / "openfda_request_log.jsonl"
        self.daily_budget = daily_budget
        self.min_interval_s = min_interval_s
        self.offline = offline
        self.session = requests.Session()
        self._last = 0.0
        self.network_calls = 0
        self.cache_hits = 0
        self.urls_used: dict[str, str] = {}  # every distinct query asked for -> fetch time of its response

    def _requests_last_24h(self) -> int:
        if not self.log_path.exists():
            return 0
        since = datetime.now(timezone.utc) - timedelta(hours=24)
        n = 0
        for line in self.log_path.read_text(encoding="utf-8").splitlines():
            rec = json.loads(line)
            if datetime.fromisoformat(rec["fetched_utc"]) >= since:
                n += 1
        return n

    def _log(self, url: str, status: int, fetched: str) -> None:
        with self.log_path.open("a", encoding="utf-8") as f:
            f.write(json.dumps({"url": url, "status": status, "fetched_utc": fetched}) + "\n")

    def get(self, params: dict[str, str | int]) -> dict:
        """Return ``{"url", "fetched_utc", "status", "data"}`` for a query, from the cache when possible.

        A 404 from openFDA means "no matching records" and is cached as an empty result with total 0.
        """
        url = build_url(params)
        path = self.dir / f"{cache_key(url)}.json"
        if path.exists():
            self.cache_hits += 1
            rec = json.loads(path.read_text(encoding="utf-8"))
            self.urls_used[url] = rec["fetched_utc"]
            return rec
        if self.offline:
            raise FileNotFoundError(f"not cached (offline mode): {url}")
        if self._requests_last_24h() >= self.daily_budget:
            raise BudgetExceeded(f"daily budget of {self.daily_budget} openFDA requests reached")
        server_errors = 0
        for attempt in range(5):
            wait = self.min_interval_s - (time.monotonic() - self._last)
            if wait > 0:
                time.sleep(wait)
            self._last = time.monotonic()
            r = self.session.get(url, timeout=90)
            fetched = datetime.now(timezone.utc).isoformat(timespec="seconds")
            self._log(url, r.status_code, fetched)
            self.network_calls += 1
            if r.status_code >= 500:  # usually a query openFDA cannot run; retry once only
                server_errors += 1
                if server_errors >= 2:
                    break
                time.sleep(3)
                continue
            if r.status_code == 429:
                time.sleep(min(60, 2 ** attempt * 5))
                continue
            if r.status_code == 404:
                data = {"meta": {"results": {"total": 0}}, "results": []}
            else:
                r.raise_for_status()
                data = r.json()
            rec = {"url": url, "fetched_utc": fetched, "status": r.status_code, "data": data}
            path.write_text(json.dumps(rec), encoding="utf-8")
            self.urls_used[url] = fetched
            return rec
        r.raise_for_status()
        raise RuntimeError(f"openFDA kept failing: {url}")

    def total(self, search: str) -> tuple[int, dict]:
        """Number of reports matching ``search`` (``meta.results.total`` of a one-record query) and the cache record."""
        rec = self.get({"search": search, "limit": 1})
        return int(rec["data"]["meta"]["results"]["total"]), rec

    def count(self, search: str | None, field: str, limit: int = COUNT_LIMIT) -> tuple[dict[str, int], dict]:
        """Term counts of ``field`` over the reports matching ``search`` (top ``limit`` terms) and the cache record.

        Each count is a number of reports: a report that lists a term twice is counted once.
        """
        params: dict[str, str | int] = {"count": field, "limit": limit}
        if search:
            params = {"search": search, **params}
        rec = self.get(params)
        rows = rec["data"].get("results", [])
        return {str(r.get("term", r.get("time"))): int(r["count"]) for r in rows}, rec
