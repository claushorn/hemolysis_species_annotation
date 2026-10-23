"""Retrieve DBAASP peptide records, with an on-disk cache.

Read-only GETs against a public academic database, serialised with a delay and
cached per peptide so a rebuild costs nothing. `https://dbaasp.org/peptides/<id>`
answers with JSON for both the DBAASPS_ and DBAASPR_ id namespaces -- they share
the same numeric endpoint (verified over the full AMPBench reference set).
"""

from __future__ import annotations

import json
import sys
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ENDPOINT = "https://dbaasp.org/peptides/{pid}"
DELAY = 0.15      # be a polite client
WORKERS = 8       # DBAASP answers in ~5 s/request; 8 in flight is ~1.6 req/s
ATTEMPTS = 3


def fetch(pid: int, cache_dir: Path | str) -> dict | None:
    """One peptide record, from cache when present. None if DBAASP will not serve it."""
    cache = Path(cache_dir)
    cache.mkdir(parents=True, exist_ok=True)
    f = cache / f"{pid}.json"
    if f.exists():
        try:
            return json.loads(f.read_text())
        except json.JSONDecodeError:
            f.unlink()        # a truncated write; re-fetch rather than guess
    for attempt in range(ATTEMPTS):
        try:
            req = urllib.request.Request(
                ENDPOINT.format(pid=pid),
                headers={"Accept": "application/json", "User-Agent": "amp-research/0.1"},
            )
            with urllib.request.urlopen(req, timeout=30) as r:
                d = json.loads(r.read().decode())
            f.write_text(json.dumps(d))
            time.sleep(DELAY)
            return d
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as e:
            if attempt == ATTEMPTS - 1:
                print(f"    peptide {pid}: giving up ({e})", file=sys.stderr)
                return None
            time.sleep(1.5 * (attempt + 1))
    return None


def fetch_many(pids: list[int], cache_dir: Path | str,
               progress_every: int = 200) -> dict[int, dict]:
    """Fetch a list of peptides concurrently. Returns only those that resolved."""
    store: dict[int, dict] = {}
    done = 0
    with ThreadPoolExecutor(max_workers=WORKERS) as pool:
        for pid, d in zip(pids, pool.map(lambda p: fetch(p, cache_dir), pids)):
            done += 1
            if d:
                store[pid] = d
            if progress_every and done % progress_every == 0:
                print(f"  {done}/{len(pids)}", flush=True)
    return store
