"""Recovered target-cell annotation for the DBAASP hemolytic activity records.

The AMPBench-MT and QMAP HC50 regression targets both aggregate DBAASP
hemolytic/cytotoxic records while discarding which cell the assay used. This
package recovers that field from DBAASP and publishes it as an annotation keyed
by DBAASP's own identifiers, so either benchmark can condition on it, report
per-cell metrics, or exclude the non-mammalian rows.
"""

from .annotate import (
    COLUMNS,
    ERYTHROCYTE_TAXON,
    TOKEN,
    annotate,
    index_records,
    is_hc50_endpoint,
    parse_source_ids,
    unreviewed_cells,
)
from .fetch import ENDPOINT, fetch, fetch_many

__all__ = [
    "COLUMNS", "ENDPOINT", "ERYTHROCYTE_TAXON", "TOKEN", "annotate", "fetch",
    "fetch_many", "index_records", "is_hc50_endpoint", "parse_source_ids",
    "unreviewed_cells",
]
