"""Build the published annotation table, and the coverage report that goes with it.

Scope is the union of the peptides whose hemolytic records the two benchmarks
consume:

  * QMAP      -- every peptide in `DBAASPDataset().with_hc50()`; its `id` is the
                 DBAASP peptide id, so the join needs nothing else.
  * AMPBench  -- every peptide referenced by `source_ids` in the HC50 table.
                 Optional: pass --ampbench <dir> if you have that corpus. The
                 annotation itself does not depend on it, and AMPBench's own
                 columns are never copied into the output.

Usage:
    uv run python build_annotation.py [--ampbench PATH] [--cache PATH]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from dbaasp_target_cell import (
    COLUMNS,
    annotate,
    fetch_many,
    index_records,
    parse_source_ids,
    unreviewed_cells,
)

HERE = Path(__file__).resolve().parent
OUT = HERE / "annotation"
CSV = "dbaasp_hemolytic_target_cells.csv"


def qmap_hc50_peptide_ids() -> list[int]:
    from qmap import DBAASPDataset
    d = DBAASPDataset().with_hc50()
    return sorted({d[i].id for i in range(len(d))})


def ampbench_refs(ampbench_dir: Path) -> pd.DataFrame:
    df = pd.read_csv(Path(ampbench_dir) / "hc50_hemolysis" / "all.csv")
    return parse_source_ids(df.source_ids)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ampbench", type=Path, default=None,
                    help="AMPBench-MT hc50 corpus directory (optional)")
    ap.add_argument("--cache", type=Path, default=HERE / "cache" / "dbaasp")
    args = ap.parse_args()

    qmap_ids = set(qmap_hc50_peptide_ids())
    pids = set(qmap_ids)
    print(f"QMAP HC50 peptides: {len(pids)}")

    refs = None
    if args.ampbench:
        refs = ampbench_refs(args.ampbench)
        extra = set(refs.peptide_id) - pids
        print(f"AMPBench HC50 record references: {len(refs)} "
              f"over {refs.peptide_id.nunique()} peptides "
              f"({len(extra)} not already in QMAP)")
        pids |= set(refs.peptide_id)

    pids = sorted(pids)
    print(f"\nFetching {len(pids)} DBAASP peptides (cached after the first run) ...",
          flush=True)
    store = fetch_many(pids, args.cache)
    print(f"  retrieved {len(store)}/{len(pids)} peptides")

    rec = index_records(store)
    print(f"  hemolytic/cytotoxic records indexed: {len(rec)}")

    missing = unreviewed_cells(rec)
    if missing:
        raise RuntimeError(
            "These erythrocyte target cells have no reviewed taxon entry:\n  "
            + "\n  ".join(missing)
            + "\n\nAdd each to ERYTHROCYTE_TAXON in dbaasp_target_cell/annotate.py. "
              "The artifact is not written with an unclassified cell in it, because "
              "that is a silent gap in the one column it exists to supply."
        )

    ann = annotate(rec)
    OUT.mkdir(parents=True, exist_ok=True)
    ann[COLUMNS].to_csv(OUT / CSV, index=False)

    cov = coverage(ann, refs, qmap_ids)
    (OUT / "coverage.json").write_text(json.dumps(cov, indent=2) + "\n")

    print(f"\nwrote {OUT / CSV}  ({len(ann)} records, {ann.peptide_id.nunique()} peptides)")
    report(ann, cov)


def pooling_stats(elig: pd.DataFrame, prefix: str) -> dict:
    """Donor pooling at peptide granularity -- the granularity QMAP publishes."""
    per = elig.groupby("peptide_id").target_cell.nunique()
    nuc = elig.groupby("peptide_id").nucleated.agg(["any", "all"])
    return {
        f"{prefix}_peptides_with_eligible_records": int(len(per)),
        f"{prefix}_peptides_pooling_multiple_target_cells": int((per > 1).sum()),
        f"{prefix}_peptides_pooling_multiple_target_cells_pct": round(
            100.0 * (per > 1).mean(), 2) if len(per) else 0.0,
        f"{prefix}_max_target_cells_in_one_peptide": int(per.max()) if len(per) else 0,
        f"{prefix}_peptides_mixing_nucleated_and_not": int(
            (nuc["any"] & ~nuc["all"]).sum()),
    }


def coverage(ann: pd.DataFrame, refs: pd.DataFrame | None,
             qmap_ids: set[int] | None = None) -> dict:
    ery = ann[ann.is_erythrocyte]
    # Records at a 50% endpoint are the ones an HC50 target can legitimately be
    # built from. The rest are hemolysis measured at some other level, which
    # DBAASP stores deliberately -- their presence here is NOT a defect, and
    # counting them as one would overstate the problem considerably.
    elig = ery[ery.endpoint_is_hc50]
    cov: dict = {
        "records_total": int(len(ann)),
        "records_erythrocyte": int(len(ery)),
        "records_erythrocyte_hc50_endpoint": int(len(elig)),
        "records_erythrocyte_other_endpoint": int(len(ery) - len(elig)),
        "peptides": int(ann.peptide_id.nunique()),
        "target_cells_distinct": int(ery.target_cell.nunique()),
        "target_cell_named_pct": round(100.0 * (ery.target_cell != "").mean(), 2),
        "non_mammalian_records": int((ery.taxon_class != "mammal").sum()),
        # Censoring, which both benchmarks drop. A `>c` record says the assay
        # never reached the endpoint, so the HC50 is ABOVE c; recorded as the
        # point c it understates safety, selectively for the peptides that
        # matter most to a safety model.
        "records_by_censoring": ann.censoring.value_counts().to_dict(),
        "records_erythrocyte_hc50_right_censored": int(
            (elig.censoring == "right").sum()),
        "records_erythrocyte_hc50_right_censored_pct": round(
            100.0 * (elig.censoring == "right").mean(), 2) if len(elig) else 0.0,
        "nucleated_records": int(ery.nucleated.fillna(False).sum()),
        "by_target_cell": ery.target_cell.value_counts().to_dict(),
    }
    # Donor pooling, over the union of both benchmarks' peptides and again over
    # QMAP's alone. A claim addressed to one benchmark must use that benchmark's
    # own count, not the union's.
    cov.update(pooling_stats(elig, "union"))
    if qmap_ids is not None:
        cov.update(pooling_stats(elig[elig.peptide_id.isin(qmap_ids)], "qmap"))
    if refs is not None:
        named = set(ann[ann.target_cell != ""].record_id)
        cov["ampbench_refs"] = int(len(refs))
        cov["ampbench_refs_resolved"] = int(refs.record_id.isin(named).sum())
        cov["ampbench_refs_resolved_pct"] = round(
            100.0 * refs.record_id.isin(named).mean(), 2)
        # The defect claim applies only to records a benchmark CONSUMES.
        consumed = ann[ann.record_id.isin(set(refs.record_id))]
        off = consumed[~consumed.endpoint_is_hc50]
        cov["ampbench_consumed_records_other_endpoint"] = int(len(off))
        cov["ampbench_consumed_other_endpoint_detail"] = (
            off[["peptide_id", "record_id", "target_cell", "lysis_endpoint"]]
            .to_dict("records"))
    return cov


def report(ann: pd.DataFrame, cov: dict) -> None:
    ery = ann[ann.is_erythrocyte]
    print("\n" + "=" * 72)
    print("Coverage")
    print("=" * 72)
    for k in ["records_total", "records_erythrocyte",
              "records_erythrocyte_hc50_endpoint",
              "records_erythrocyte_other_endpoint", "peptides",
              "target_cells_distinct", "target_cell_named_pct",
              "ampbench_refs", "ampbench_refs_resolved",
              "ampbench_refs_resolved_pct"]:
        if k in cov:
            print(f"  {k:38s} {cov[k]}")

    print("\n" + "=" * 72)
    print("Erythrocyte records by donor")
    print("=" * 72)
    vc = ery.target_cell.value_counts()
    for name, k in vc.items():
        t = ery[ery.target_cell == name].iloc[0]
        flag = "" if t.taxon_class == "mammal" else f"   <-- {t.taxon_class}, nucleated"
        print(f"  {k:6d}  {100 * k / len(ery):5.1f}%  {name:26s}{flag}")

    print("\n" + "=" * 72)
    print("Censoring: how many records are bounds rather than measurements")
    print("=" * 72)
    for k, v in cov["records_by_censoring"].items():
        print(f"  {k:12s} {v:6d}  ({100*v/cov['records_total']:5.1f}%)")
    print(f"  of the {cov['records_erythrocyte_hc50_endpoint']} eligible HC50 records, "
          f"{cov['records_erythrocyte_hc50_right_censored']} "
          f"({cov['records_erythrocyte_hc50_right_censored_pct']}%) are right-censored")

    print("\n" + "=" * 72)
    print("Why this matters for evaluation, not for accuracy")
    print("=" * 72)
    print(f"  records on non-mammalian erythrocytes   : "
          f"{cov['non_mammalian_records']} "
          f"(nucleated: {cov['nucleated_records']})")
    for scope in ["union", "qmap"]:
        if f"{scope}_peptides_with_eligible_records" not in cov:
            continue
        print(f"  [{scope}] peptides pooling >1 donor        : "
              f"{cov[f'{scope}_peptides_pooling_multiple_target_cells']} of "
              f"{cov[f'{scope}_peptides_with_eligible_records']} "
              f"({cov[f'{scope}_peptides_pooling_multiple_target_cells_pct']}%)"
              f"  max {cov[f'{scope}_max_target_cells_in_one_peptide']} cells; "
              f"{cov[f'{scope}_peptides_mixing_nucleated_and_not']} mix nucleated")
    if "ampbench_consumed_records_other_endpoint" in cov:
        print(f"  consumed records NOT at a 50% endpoint   : "
              f"{cov['ampbench_consumed_records_other_endpoint']}")
    print("\n  A peptide-level HC50 that pools several donors cannot be")
    print("  reproduced, and a pooled label leaves no variable to condition on.")
    print("  Note: the other-endpoint erythrocyte records counted above are NOT")
    print("  a defect -- DBAASP stores hemolysis at many levels on purpose, and")
    print("  a benchmark is expected to filter to the 50% ones.")


if __name__ == "__main__":
    main()
