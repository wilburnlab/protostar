"""Reference-library acquisition (fetch + extract).

Download the published ProteomeTools ``.msp`` spectral libraries
(https://www.proteometools.org/index.php?id=53, archived as Zenodo record
``15705607``) and extract the ``.msp`` files. Ingesting them into a
Constellation ``massspec.library.Library`` (via ``massspec.io.msp``) and
associating them with the raw acquisitions is a separate, later task.
"""

from .zenodo import (
    DEFAULT_RECORD_ID,
    ZenodoFile,
    extract_msp,
    fetch_libraries,
    list_library_files,
)

__all__ = [
    "DEFAULT_RECORD_ID",
    "ZenodoFile",
    "extract_msp",
    "fetch_libraries",
    "list_library_files",
]
