"""Build the record-level target-cell annotation.

Both the AMPBench-MT and QMAP HC50 targets are aggregates over DBAASP
hemolytic/cytotoxic activity records, and both drop the identity of the cell the
assay was performed on. AMPBench keeps a two-valued `cell_category` and a count
(`n_target_cells`) but not the names; QMAP keeps `[min, max, consensus]` with no
cell field at all, while retaining full species names for its bacterial MIC
targets. The names are still in DBAASP, under
`hemoliticCytotoxicActivities[].targetCell.name`.

This module recovers them and nothing else. The output deliberately carries only
DBAASP-derived facts plus DBAASP's own identifiers as join keys: no sequences, no
HC50 values, no splits, no cluster assignments. That keeps the artifact an
ANNOTATION over records that are already openly redistributable, rather than a
republication of either benchmark, and lets each benchmark join it to whatever
it already holds.
"""

from __future__ import annotations

import re

import pandas as pd

# AMPBench encodes its provenance as DBAASP{S,R}_<peptideId>:toxicity:<recordId>.
TOKEN = re.compile(r"DBAASP[SR]?_(\d+):toxicity:(\d+)")

# Taxon of each erythrocyte donor seen in the DBAASP hemolytic records.
#
# `nucleated` is the biophysically load-bearing column, not `mammal`: avian and
# fish erythrocytes are nucleated and differ in membrane lipid composition, so
# pooling them with mammalian erythrocytes mixes two different systems. Any name
# absent from this table is reported for review rather than guessed -- see
# `unreviewed_cells`.
ERYTHROCYTE_TAXON = {
    "Human erythrocytes":      ("Homo sapiens",            "mammal", False),
    "Horse erythrocytes":      ("Equus caballus",          "mammal", False),
    "Mouse erythrocytes":      ("Mus musculus",            "mammal", False),
    "Rat erythrocytes":        ("Rattus norvegicus",       "mammal", False),
    "Sheep erythrocytes":      ("Ovis aries",              "mammal", False),
    "Pig erythrocytes":        ("Sus scrofa",              "mammal", False),
    "Rabbit erythrocytes":     ("Oryctolagus cuniculus",   "mammal", False),
    "Guinea pig erythrocytes": ("Cavia porcellus",         "mammal", False),
    "Bovine erythrocytes":     ("Bos taurus",              "mammal", False),
    "Buffalo erythrocytes":    ("Bubalus bubalis",         "mammal", False),
    "Dog erythrocytes":        ("Canis lupus familiaris",  "mammal", False),
    "Goat erythrocytes":       ("Capra aegagrus hircus",   "mammal", False),
    "Chicken erythrocytes":    ("Gallus gallus",           "bird",   True),
    "Salmon erythrocytes":     ("Salmo salar",             "fish",   True),
    "Tilapia erythrocytes":    ("Oreochromis sp.",         "fish",   True),
    "Fish erythrocytes":       ("unspecified fish",        "fish",   True),
}

# A plain 50% endpoint. DBAASP also stores the measured lysis fraction verbatim,
# so "50.84+/-1.38% Hemolysis" is a 50% endpoint reported with its uncertainty,
# whereas "9.63+/-1.50% Hemolysis" is not an HC50 at all and ">50% Hemolysis" is
# a bound rather than a point.
_HC50_ENDPOINT = re.compile(r"^\s*50(\.\d+)?\s*(±|\+/-)?\s*[\d.]*\s*%?\s*hemolysis\s*$",
                            re.IGNORECASE)

_NUMBER = re.compile(r"\d+(?:\.\d+)?(?:[Ee][+-]?\d+)?")


def parse_concentration(raw: object) -> tuple[float, str]:
    """Split DBAASP's free-text concentration into a value and a censoring code.

    DBAASP stores the tested concentration as free text, and that text is the
    only place the CENSORING survives. `256` is a measured HC50; `>256` means
    the assay never reached the endpoint at the highest concentration tested,
    so the true HC50 lies somewhere above 256. The numeric `activity` field
    DBAASP also provides silently discards the `>`, and both HC50 benchmarks
    inherit that loss: a bound becomes a point.

    The orientation is easy to get backwards and matters a great deal. A
    right-censored record is evidence that the peptide is *less* hemolytic than
    the bound, so recording it as `HC50 = bound` overstates toxicity, and it
    does so selectively, for exactly the non-hemolytic peptides a safety model
    exists to identify.

    Returns (value, censoring), with censoring one of:
      ``exact``      a point measurement, possibly reported with a spread
      ``right``      the assay did not reach the endpoint: HC50 > value
      ``left``       HC50 < value
      ``interval``   a reported range; the value is its geometric centre
      ``ambiguous``  a bound whose direction the string does not fix (``200<``)
      ``missing``    no parseable number

    Only the censoring code is published (see COLUMNS); the value is returned
    for callers that hold the measurements themselves, so that the code which
    classifies a record and the code which reads its value cannot disagree.
    """
    if not isinstance(raw, str):
        return float("nan"), "missing"
    s = raw.strip()
    if not s or s.upper() == "NA":
        return float("nan"), "missing"

    nums = _NUMBER.findall(s)
    if not nums:
        return float("nan"), "missing"
    value = float(nums[0])

    # Two-character operators first, or `>=65` and `>>600` both read as `>`.
    if s.startswith((">=", ">>")) or s.lower().startswith("up to"):
        return value, "right"
    if s.startswith(">"):
        return value, "right"
    if s.startswith("<="):
        return value, "left"
    if s.startswith("<"):
        return value, "left"
    if s.endswith(("<", ">")):
        # e.g. `200<`. Four records; the intended direction is not recoverable
        # from the string, so they are flagged rather than guessed.
        return value, "ambiguous"
    if len(nums) == 2 and re.fullmatch(r"\s*[\d.]+\s*-\s*[\d.]+\s*", s):
        lo, hi = float(nums[0]), float(nums[1])
        if lo > 0 and hi > 0:
            return float((lo * hi) ** 0.5), "interval"
        return value, "interval"
    # `221.8±7.8` and friends: a point measurement reported with its spread.
    return value, "exact"


def parse_source_ids(source_ids: pd.Series) -> pd.DataFrame:
    """AMPBench `source_ids` -> one row per referenced DBAASP toxicity record."""
    refs = []
    for i, s in zip(source_ids.index, source_ids.astype(str)):
        for tok in s.split(";"):
            m = TOKEN.match(tok.strip())
            if m:
                refs.append((i, int(m.group(1)), int(m.group(2))))
    return pd.DataFrame(refs, columns=["row", "peptide_id", "record_id"])


def index_records(store: dict[int, dict]) -> pd.DataFrame:
    """Flatten fetched peptides into one row per hemolytic/cytotoxic record."""
    rows = []
    for pid, d in store.items():
        for r in d.get("hemoliticCytotoxicActivities") or []:
            rid = r.get("id")
            if rid is None:
                continue
            rows.append(dict(
                peptide_id=int(pid),
                record_id=int(rid),
                target_cell=((r.get("targetCell") or {}).get("name") or "").strip(),
                concentration_unit=((r.get("unit") or {}).get("name") or "").strip(),
                lysis_endpoint=(r.get("activityMeasureForLysisValue") or "").strip(),
                # Carried only so `annotate` can derive the censoring code; the
                # raw string holds a measured value and is NOT published.
                _concentration_raw=r.get("concentration"),
            ))
    return pd.DataFrame(rows)


def is_hc50_endpoint(endpoint: str) -> bool:
    """Is this record's endpoint a 50% lysis point, as opposed to another level?"""
    return bool(_HC50_ENDPOINT.match(endpoint or ""))


def unreviewed_cells(rec: pd.DataFrame) -> list[str]:
    """Erythrocyte target cells with no entry in ERYTHROCYTE_TAXON.

    A name reaching the published artifact without a reviewed taxon would be a
    silent gap in exactly the column the artifact exists to provide, so the
    build refuses to write until each one is classified by hand.
    """
    seen = set(rec.target_cell.dropna())
    ery = {c for c in seen if "erythrocyte" in c.lower() or "rbc" in c.lower()}
    return sorted(ery - set(ERYTHROCYTE_TAXON))


def annotate(rec: pd.DataFrame) -> pd.DataFrame:
    """Add the taxon and endpoint columns to an indexed record table."""
    out = rec.copy()
    tax = out.target_cell.map(ERYTHROCYTE_TAXON)
    out["species"] = [t[0] if isinstance(t, tuple) else "" for t in tax]
    out["taxon_class"] = [t[1] if isinstance(t, tuple) else "" for t in tax]
    out["nucleated"] = [t[2] if isinstance(t, tuple) else None for t in tax]
    out["is_erythrocyte"] = out.target_cell.str.contains(
        "erythrocyte|RBC", case=False, na=False)
    out["endpoint_is_hc50"] = out.lysis_endpoint.map(is_hc50_endpoint)
    if "_concentration_raw" in out:
        out["censoring"] = [parse_concentration(c)[1]
                            for c in out._concentration_raw]
    else:
        out["censoring"] = "missing"
    return out.sort_values(["peptide_id", "record_id"]).reset_index(drop=True)


COLUMNS = [
    "peptide_id", "record_id", "target_cell", "species", "taxon_class",
    "nucleated", "is_erythrocyte", "concentration_unit", "lysis_endpoint",
    "endpoint_is_hc50", "censoring",
]
# `_concentration_raw` is deliberately absent: it is the measured value, and
# this artifact publishes identifiers and DBAASP-derived classifications only.
