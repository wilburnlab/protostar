"""Reference-library ingestion and optional re-search.

Ingest the published ProteomeTools ``.msp`` spectral libraries
(https://www.proteometools.org/index.php?id=53) into a Constellation
``massspec.library.Library`` via ``massspec.io.msp`` as the canonical
identification/reference source. Optionally re-search the local ``.raw`` with
EncyclopeDIA/Scribe through Constellation's ``thirdparty`` wrapper.
"""
