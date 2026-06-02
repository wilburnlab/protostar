"""Synthetic-peptide reference — sequence ↔ pool ↔ set ground truth.

Curates the ProteomeTools manuscript supplements (a one-time manual Excel input)
into a compact, version-controlled parquet per dataset (shipped under
``protostar/peptides/data/``). The pool names match the acquisition ``.raw``
stems, so this is the known-answer join from peptide → pool → run used by XIC
extraction / validation and downstream identification.

``load_reference(dataset)`` reads the committed parquet (pyarrow only);
``parse_supplement`` (in ``sources``) rebuilds it from the Excel workbooks and
needs the ``peptides`` optional extra (``openpyxl``).
"""

from .procal import (
    PROCAL_SCHEMA_VERSION,
    PROCAL_SUPPLEMENT_FILE,
    PROCAL_TABLE,
    load_procal,
    parse_procal_supplement,
    procal_path,
    procal_sequences,
    write_procal,
)
from .reference import (
    PEPTIDE_REFERENCE_SCHEMA_VERSION,
    PEPTIDE_REFERENCE_TABLE,
    POOL_PREFIX_DATASET,
    load_pool_targets,
    load_reference,
    pool_prefix,
    reference_path,
    resolve_pool_dataset,
    write_reference,
)
from .sources import DATASET_SPECS, SUPPLEMENT_FILES, parse_supplement

__all__ = [
    "DATASET_SPECS",
    "PEPTIDE_REFERENCE_SCHEMA_VERSION",
    "PEPTIDE_REFERENCE_TABLE",
    "POOL_PREFIX_DATASET",
    "PROCAL_SCHEMA_VERSION",
    "PROCAL_SUPPLEMENT_FILE",
    "PROCAL_TABLE",
    "SUPPLEMENT_FILES",
    "load_pool_targets",
    "load_procal",
    "load_reference",
    "parse_procal_supplement",
    "parse_supplement",
    "pool_prefix",
    "procal_path",
    "procal_sequences",
    "reference_path",
    "resolve_pool_dataset",
    "write_procal",
    "write_reference",
]
