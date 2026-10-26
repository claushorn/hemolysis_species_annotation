# Which cell was the assay run on, and did the assay finish? A recovered annotation for DBAASP hemolytic records

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23661947.svg)](https://doi.org/10.5281/zenodo.23661947)

Both the **AMPBench-MT** and **QMAP** HC50 regression targets aggregate DBAASP
hemolytic/cytotoxic activity records while discarding the identity of the cell the assay
was performed on.

- **AMPBench-MT** keeps a two-valued `cell_category` (`human_erythrocytes` /
  `mammalian_erythrocytes`) and a count, `n_target_cells` — it counts the target cells and
  then drops their names.
- **QMAP** keeps `[min, max, consensus]` with no cell field at all, while retaining full
  species names for its bacterial MIC targets. It conditions potency on species and pools
  safety.

The names were never lost, only dropped: they are in DBAASP under
`hemoliticCytotoxicActivities[].targetCell.name`. This repository recovers them for every
record either benchmark consumes and publishes them as an annotation keyed by DBAASP's own
identifiers.

**This is an evaluation-correctness artifact, not a route to better accuracy.** On a matched
protocol, conditioning a model on the target cell is worth +0.013 Pearson *r* once the
estimator is fixed — marginal. What the annotation buys is the ability to *report* and
*reproduce* per-cell metrics, and to see that a single official fold can read −0.016 where
cluster-grouped CV over the same category reads +0.278. The artifact contains no labels, so
it cannot be used to claim predictive gains.

## The file

`annotation/dbaasp_hemolytic_target_cells.csv` — one row per DBAASP hemolytic/cytotoxic
record, **5,461 records over 3,001 peptides**.

| column | meaning |
|---|---|
| `peptide_id` | DBAASP peptide id. QMAP's `Sample.id` **is** this id; AMPBench encodes it in `source_ids` as `DBAASP{S,R}_<peptide_id>:toxicity:<record_id>`. |
| `record_id` | DBAASP activity-record id — the join key at assay granularity. |
| `target_cell` | DBAASP's verbatim `targetCell.name`, e.g. `Horse erythrocytes`. **100% populated.** |
| `species` | Binomial name of the donor, assigned by hand per target cell. |
| `taxon_class` | `mammal` / `bird` / `fish`. Empty for non-erythrocyte targets. |
| `nucleated` | Whether those erythrocytes are nucleated. **The load-bearing column** — see below. |
| `is_erythrocyte` | False for the cytotoxic cell-line records DBAASP stores alongside. |
| `concentration_unit` | DBAASP's reported unit, `µM` or `µg/ml`. |
| `lysis_endpoint` | Verbatim endpoint, e.g. `50% Hemolysis`, `9.63±1.50% Hemolysis`, `>50% Hemolysis`. |
| `endpoint_is_hc50` | Whether that endpoint is a 50% lysis point rather than another level or a bound. |
| `censoring` | Whether the tested concentration was a point (`exact`), a bound the assay never reached (`right`), `left`, an `interval`, `ambiguous`, or `missing`. **The second load-bearing column** - see below. |

`annotation/coverage.json` carries the counts below in machine-readable form.

## What it shows

### One label name, sixteen donors

Of 3,975 erythrocyte records:

| donor | records | share | |
|---|---|---|---|
| Human | 2,860 | 71.9% | |
| Horse | 228 | 5.7% | |
| Mouse | 223 | 5.6% | |
| Rat | 215 | 5.4% | |
| Pig | 185 | 4.7% | |
| Sheep | 147 | 3.7% | |
| Rabbit | 69 | 1.7% | |
| Bovine | 14 | 0.4% | |
| Dog | 9 | 0.2% | |
| Guinea pig | 6 | 0.2% | |
| Chicken | 5 | 0.1% | **bird, nucleated** |
| Fish (unspecified) | 4 | 0.1% | **fish, nucleated** |
| Salmon | 3 | 0.1% | **fish, nucleated** |
| Goat | 3 | 0.1% | |
| Tilapia | 2 | 0.1% | **fish, nucleated** |
| Buffalo | 2 | 0.1% | |

In AMPBench-MT's HC50 table the `mammalian_erythrocytes` category pools **14** of these,
the largest at **28.8%** — and includes chicken, fish and salmon. The `human_erythrocytes`
category is 100% one target cell, and is also the category that predicts.

### `nucleated`, not `mammal`, is the real distinction

Calling a bucket "mammalian" when it holds chicken and salmon is a naming error, but the
reason it matters is biophysical: **avian and fish erythrocytes are nucleated** and differ
in membrane lipid composition, so pooling them with mammalian erythrocytes mixes two
different systems in a target that a model is asked to regress. 14 records sit on
nucleated erythrocytes.

### Pooled peptide-level HC50 is not reproducible

Restricted to records at a 50% endpoint — the ones an HC50 target can legitimately be built
from — peptides whose records span more than one donor:

| scope | pooling >1 donor | max donors on one peptide | mixing nucleated & not |
|---|---|---|---|
| union of both benchmarks' peptides | 47 of 2,997 (1.6%) | 5 | 4 |
| QMAP's HC50 peptides alone | 44 of 2,797 (1.6%) | 4 | 4 |

A peptide-level `[min, max, consensus]` over such a set cannot be reconstructed, and the
pooled value leaves no variable to condition on. Both scopes are in `coverage.json`: a
claim addressed to one benchmark should use that benchmark's own count.

### A quarter of the records are bounds, not measurements

DBAASP stores the tested concentration as free text, and that text is the only
place the censoring survives. `256` is a measured HC50; **`>256` means the assay
never reached 50% lysis at the highest concentration tested**, so the true HC50
lies somewhere above 256. The numeric field DBAASP also provides drops the `>`,
and both benchmarks inherit that loss.

| form | records | meaning |
|---|---|---|
| point | 3,827 | a measurement, possibly with a reported spread |
| **`>N`, `>=N`, `>>N`, `up to N`** | **1,503 (27.5%)** | **right-censored: HC50 is above the bound** |
| `N-N` | 36 | a range; `censoring = interval` |
| `<N` | 29 | left-censored |
| `N<` | 4 | direction not recoverable from the string; flagged, not guessed |
| `NA` / empty | 62 | no value |

Restricted to the erythrocyte records at a 50% endpoint that an HC50 target is
built from, **1,171 of 3,142 (37.3%) are right-censored.**

The bias has a direction, and it is the unsafe one. A right-censored record is
evidence that the peptide is *less* hemolytic than the bound, so recording it as
`HC50 = bound` overstates toxicity - selectively, for exactly the non-hemolytic
peptides a safety model exists to identify. It also manufactures a spike in the
label distribution at whatever round concentration people stopped testing at,
which no model can explain and every model is scored against.

This column is why the artifact is worth joining even if you do not care about
the target cell. Treating these records as bounds rather than points, with a
Tobit likelihood on otherwise identical features and scored only on uncensored
test rows, is worth **+0.074 Pearson r (95% CI [+0.034, +0.119])** on QMAP's own
five splits - substantially more than conditioning a model on the target cell,
and more than any change of estimator we were able to find.

### Four consumed records are not HC50 measurements

Among the 1,852 records AMPBench-MT's HC50 table references, four carry an endpoint that is
not a 50% lysis point:

| peptide | record | donor | endpoint |
|---|---|---|---|
| 4680 | 3561 | Human | `>50% Hemolysis` — a bound, collapsed to a point |
| 10237 | 8739 | Rabbit | `9.63±1.50% Hemolysis` |
| 14040 | 12760 | Horse | `72.35±17.50% Hemolysis` |
| 21835 | 20788 | Human | `0.50±0.26% Hemolysis` |

Four rows in 1,852 changes no headline number. It is listed because the `lysis_endpoint`
column makes it checkable, and because `>50%` is a censoring indicator that both benchmarks
drop — AMPBench has 30 rows piled at exactly 100 µM with no censoring column.

**Not a defect:** 833 erythrocyte records sit at some other endpoint (10% hemolysis, and so
on). DBAASP stores hemolysis at many levels deliberately and a benchmark is expected to
filter; that count is reported for transparency, not as an error.

### Units are handled correctly

541 of AMPBench's referenced records are reported in µg/ml rather than µM. AMPBench converts
them and records how, in `conversion_methods` (`direct_uM` vs
`ug_per_ml_to_uM_estimated_unmodified_mw`). **There is no unit bug.** Stated because it is
the first thing a reader of the `concentration_unit` column will suspect.

### A bonus for cytotoxicity work

1,486 of the 5,461 records are not erythrocyte assays but cytotoxic measurements on **98
distinct cell lines**, carried along because they live in the same DBAASP field. They are
labelled `is_erythrocyte = false` and left otherwise unannotated.

## Using it

```python
import pandas as pd
ann = pd.read_csv("annotation/dbaasp_hemolytic_target_cells.csv")

# QMAP: Sample.id is the DBAASP peptide id.
from qmap import DBAASPDataset
ds = DBAASPDataset().with_hc50()
cells = (ann[ann.is_erythrocyte & ann.endpoint_is_hc50]
         .groupby("peptide_id").target_cell.agg(
             lambda s: s.iloc[0] if s.nunique() == 1 else "POOLED"))
donor = [cells.get(ds[i].id, "UNKNOWN") for i in range(len(ds))]

# AMPBench-MT: parse the provenance tokens, then join on record_id.
from dbaasp_target_cell import parse_source_ids
refs = parse_source_ids(amp_hc50.source_ids)          # -> row, peptide_id, record_id
refs = refs.merge(ann, on=["peptide_id", "record_id"], how="left")
```

Report metrics per `target_cell` with homology-grouped cross-validation, not on a single
fold: fold-to-fold drift in donor composition moves the label mean by up to 0.59 log₁₀
within one category, which is enough on its own to drive a cross-fold correlation negative.

## Rebuilding

```bash
uv run python build_annotation.py --ampbench /path/to/ampbench/hc50   # --ampbench optional
```

Fetches each DBAASP peptide once (cached in `cache/`, ~3,000 records at ~1.6 req/s, so
roughly 30 minutes cold and free thereafter), indexes every hemolytic/cytotoxic record,
assigns the taxon, and writes `annotation/`.

The build **refuses to write** if it meets an erythrocyte donor with no hand-reviewed taxon
entry, rather than emitting a blank. That gate is not theoretical: extending from
AMPBench's set to QMAP's full HC50 set surfaced `Tilapia erythrocytes`, which the gate
caught. An unclassified donor would be a silent gap in the one column this artifact exists
to supply.

## Provenance and licensing

Every annotated fact comes from **DBAASP**, whose Data Access Policy permits access,
copying, adaptation and redistribution without restriction — including after modification
or combination with other data — and requires that DBAASP be acknowledged as the source.
That acknowledgement is mandatory for any downstream use; see `LICENSE`.

> Data were obtained from DBAASP (Database of Antimicrobial Activity and Structure of
> Peptides), https://dbaasp.org. Pirtskhalava M. et al., *DBAASP v3*, Nucleic Acids
> Research, 2021.

**This annotation contains no AMPBench-MT or QMAP data.** Only DBAASP record identifiers and
DBAASP-derived fields. No sequences, no HC50 values, no splits, no cluster assignments. The
`censoring` column records whether a measurement was a bound, not what the measurement was. That
is deliberate: AMPBench-MT ships under a custom research-review-use licence that does not
permit redistribution, and an annotation keyed by identifiers avoids the question entirely
while remaining joinable to either benchmark. The annotation table is released **CC BY 4.0**
and the code **MIT**.

## Citing

See `CITATION.cff`. This annotation: **doi:10.5281/zenodo.23661947** (concept DOI, always the
latest version; v1.0.0 is doi:10.5281/zenodo.23661948). Cite DBAASP in all
cases; cite QMAP (doi:10.64898/2026.02.03.703041) if you use its splits.
