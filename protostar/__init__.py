"""protostar — ProteomeTools Statistical Research.

A reproducible, open-science re-analysis of the ProteomeTools synthetic-peptide
datasets (Zolg 2017, Gessulat 2019, Wilhelm 2021), built on the Constellation
platform (``constellation.massspec`` + ``constellation.core``).

This package is deliberately *thin*: it holds data orchestration, experiment
drivers, and curated-record tooling only. All model / likelihood / peak-shape /
scoring code — including the Counter generative quantification model — lives in
``constellation.massspec``. When a new modeling capability is needed it is
contributed upstream to Constellation, never written here. See ``CLAUDE.md``.

Subpackages
-----------
fetch          Build the expected-file manifest; fresh-download + hash-verify
               ``.raw`` files (resumable / repairable); optional seed-from-local.
convert        Drive ``.raw`` -> mzpeak Parquet + scanmeta via Constellation.
library        Ingest the published ProteomeTools ``.msp`` reference libraries;
               optional EncyclopeDIA/Scribe re-search of the local ``.raw``.
metadata       Curate the acquisition time table (datetime + instrument + order).
intermediates  Drive MS1/MS2 chromatogram extraction (Constellation calls).
experiments    Shared experiment-harness utilities (IO, plotting helpers).
"""

__version__ = "0.1.0.dev0"
