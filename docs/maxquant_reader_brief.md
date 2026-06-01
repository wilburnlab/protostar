# Planning brief — MaxQuant search-output reader for Constellation (`massspec.io.maxquant`)

> Hand-off brief for an independent planning session. Goal: produce an implementation **plan**
> (not code) for a reader that ingests MaxQuant `combined/txt/` exports into Constellation
> domain objects, including a new `PSM_TABLE` schema. Ledger item **#11** (Lane A / foundational)
> in `docs/constellation_contributions.md`.

## Why this, why now

ProToStaR re-analyzes the ProteomeTools synthetic-peptide datasets on the **Constellation**
platform (`~/projects/constellation`, editable install, conda env `constellation`). Constellation's
`.raw → parquet` converter has now processed all ~4,213 acquisitions into per-scan bundles
(`peaks.parquet` + `scan_metadata.parquet` + `acquisition_metadata.parquet`). To *interpret* those
scans we need per-scan peptide identifications. Those live in each dataset's PRIDE **SEARCH** zips,
confirmed (probe: `pipelines/probe_search_format.py`) to be **MaxQuant `combined/txt/` exports**
for all three datasets (Zolg2017/PXD004732, Gessulat2019/PXD010595, Wilhelm2021/PXD021013), one zip
per acquisition. Reading them gives per-scan peptide assignments with **no re-search compute** — the
identification path that anchors every downstream experiment (PROCAL calibration, MS1 error model,
Counter).

## Ground truth to inherit (verify before planning — don't trust verbatim)

- **Net-new, not a port.** There is no MaxQuant reader in `~/projects/cartographer` (grep found
  none). Contrast the ledger's "port from Cartographer" note, which applies to items #3–9, not #11.
- **The MSP `.msp` library reader already exists** on constellation `main`
  (`constellation/massspec/io/msp/`: `_read.py`, `_annotate.py`, `_comment.py`, `_mods.py`,
  `adapters.py`, with `tests/test_msp_reader.py` + fixtures). It is the **structural template** —
  a thin `io.<fmt>` package whose `adapters.py` self-registers a Reader against the right tier's
  registry, with `__init__` triggering import-time registration (see `massspec/io/__init__.py`).
  Out of scope to rebuild.
- **Three-tier model is the central design decision.** Constellation separates
  `massspec.library` (theoretical / sample-agnostic), `massspec.quant` (empirical abundances), and
  `massspec.search` (identification scores). Each has its own `Reader/Writer` Protocol + suffix/format
  registry (`library/io.py`, `search/io.py`, `quant/io.py`) and Arrow schemas. **A MaxQuant export
  spans all three** (PSM/peptide/protein scores → search; intensities → quant; modseqs/spectra →
  library-ish). The plan must decide how a multi-tier source maps onto single-tier registries —
  likely *multiple registered readers over one shared parser*, mirroring how the EncyclopeDIA
  `.dlib`/`.elib` adapter (`massspec/io/encyclopedia/`) deliberately straddles tiers (read its
  `io/__init__.py` rationale). Reader-only: no new search/scoring **models**.
- **`msms.txt` should drive the deferred `PSM_TABLE`.** `massspec/search/schemas.py` currently ships
  only `PEPTIDE_SCORE_TABLE` / `PROTEIN_SCORE_TABLE` and *explicitly defers* a per-spectrum
  `PSM_TABLE` until "an actual PSM-emitting reader drives the schema design." MaxQuant `msms.txt` is
  exactly that per-scan, per-PSM table — so **defining `PSM_TABLE` is in scope** and this reader is
  its motivating consumer. (This touches a shared schema file, a slightly larger surface than a pure
  `io/` add — call it out as a deliberate decision.)
- **The payoff is scan-level join-back.** `msms.txt` carries `Raw file` + `Scan number`; the
  converter's `scan_metadata.parquet` is keyed by `scan`. A PSM row joins to a converted
  acquisition's scan on `(raw_file_stem, scan)`. **Cross-validation bonus:** `msms.txt` natively
  carries `Mass analyzer` (e.g. `FTMS`) and `Fragmentation` (e.g. `ETD`) — the exact `analyzer` +
  `activation_type` columns just added to `scan_metadata` (constellation PR #57). The plan should use
  these to *validate* the join, not just assume it.
- Modseqs (`Modified sequence`, e.g. `_AAVPQIK_`, `_(ac)…M(ox)…_`) should normalize to **ProForma
  2.0** via `core.sequence.proforma` (platform standard), reusing the MSP reader's modification
  resolution (`io/msp/_mods.py`) where applicable.

## Real `txt/` headers (sampled from `Thermo_SRM_Pool_88_01_01_ETD-1h-R2-tryptic.zip`, MaxQuant 1.5.3.30)

The zip carries 17 members. Priority files and **actual columns** below (✱ = key for the reader).

### `msms.txt` — 57 cols — **one row per PSM (per MS/MS scan)** — primary input
```
✱Raw file  ✱Scan number  Scan index  ✱Sequence  Length  Missed cleavages  Modifications
✱Modified sequence  Oxidation (M) Probabilities  Oxidation (M) Score Diffs  Oxidation (M)
✱Proteins  ✱Charge  ✱Fragmentation  ✱Mass analyzer  ✱Type  Scan event number  Isotope index
✱m/z  ✱Mass  ✱Mass Error [ppm]  Simple Mass Error [ppm]  ✱Retention time  ✱PEP  ✱Score
✱Delta score  Score diff  Localization prob  Combinatorics  PIF  Fraction of total spectrum
Base peak fraction  ✱Precursor Full ScanNumber  Precursor Intensity  Precursor Apex Fraction
Precursor Apex Offset  Precursor Apex Offset Time  ✱Matches  ✱Intensities  Mass Deviations [Da]
Mass Deviations [ppm]  Masses  Number of Matches  Intensity coverage  Peak coverage
Neutral loss level  ETD identification type  ✱Reverse  All scores  All sequences
All modified sequences  ✱id  Protein group IDs  ✱Peptide ID  ✱Mod. peptide ID  ✱Evidence ID
Oxidation (M) site IDs
```
Sample PSM: `Raw file=01640c_BH11-Thermo_SRM_Pool_88_01_01-ETD-1h-R2`, `Scan number=4621`,
`Sequence=AAVPQIK`, `Modified sequence=_AAVPQIK_`, `Charge=2`, `Fragmentation=ETD`,
`Mass analyzer=FTMS`, `m/z=363.72906`, `Mass Error [ppm]=0.53712`, `Retention time=23.966`,
`PEP=0.0060943`, `Score=38.303`, `Reverse=`(empty → target), `id=0`.
Note `Matches`/`Intensities` are `;`-delimited fragment annotations + intensities (mzPAF-adjacent;
a future MS2-spectral cross-check, not required for v1).

### `evidence.txt` — 59 cols — one row per (peptide, charge, raw-file) feature — secondary
Key cols: `✱Sequence  ✱Modified sequence  ✱Raw file  ✱Charge  ✱m/z  Mass  Mass Error [ppm]
Retention time  Calibrated retention time  ✱Intensity  ✱PEP  ✱Score  ✱MS/MS Count
✱MS/MS Scan Number  Reverse  Potential contaminant  ✱id  ✱Peptide ID  ✱Mod. peptide ID  ✱MS/MS IDs
Best MS/MS`. This is the feature/quant-ish bridge (intensity per feature, links many MS/MS IDs).

### `msmsScans.txt` — 41 cols — one row per MS/MS scan (identified or not) — optional
Key cols: `✱Raw file  ✱Scan number  Retention time  Ion injection time  Total ion current
Collision energy  Base peak intensity  ✱Identified  ✱MS/MS IDs  Sequence  ✱Fragmentation
✱Mass analyzer  Precursor full scan number  Precursor intensity  Scan event number
Modified sequence  Score  RawOvFtT  AGC Fill  Scan index  MS scan number`. Superset of attempted
fragmentations (incl. `Identified=-`); useful for negative/coverage analysis. Likely defer to v2.

### `peptides.txt` — 56 cols — one row per peptide sequence — defer (peptide-level rollup)
AA-count columns + `Sequence Proteins Start/End position Charges PEP Score Intensity id
Evidence IDs MS/MS IDs ...`. Rollup; lower priority than per-PSM.

### `parameters.txt` — key/value provenance (capture into a metadata side-table)
`Version=1.5.3.30`, `Fixed modifications=Carbamidomethyl (C)`, `PSM FDR=0.01`, `Protein FDR=0.01`,
`MS/MS tol. (FTMS)=20 ppm`, `MS/MS tol. (ITMS)=0.5 Da`, `Min. peptide Length=7`,
`Decoy mode=revert`, `Fasta file=...proteomeTools.fasta`, etc. Records search provenance — keep it.

Also present (likely ignore for v1): `allPeptides.txt`, `modificationSpecificPeptides.txt`,
`proteinGroups.txt`, `Oxidation (M)Sites.txt`, `aifMsms.txt`, `ms3Scans.txt`, `msScans.txt`,
`mzRange.txt`, `matchedFeatures.txt`, `libraryMatch.txt`, `summary.txt`, `tables.pdf`.

## Drafted `PSM_TABLE` schema (starting point — the planning session should refine)

One row per PSM, sourced from `msms.txt`. Field choices favor downstream join-back + provenance;
dtypes follow the `SCAN_METADATA_TABLE` conventions (int32 `scan`, float64 m/z, nullable where
MaxQuant may leave blanks). Foreign keys connect to the converted bundles and to the other tiers.

```python
PSM_TABLE = pa.schema([
    # ── identity / acquisition join-back ───────────────────────────────
    pa.field("psm_id", pa.int64(), nullable=False),          # msms.txt "id" (per-export)
    pa.field("raw_file", pa.string(), nullable=False),       # "Raw file" — maps to acquisition stem
    pa.field("acquisition_id", pa.int64(), nullable=True),   # resolved vs Acquisitions table (nullable pre-link)
    pa.field("scan", pa.int32(), nullable=False),            # "Scan number" — joins scan_metadata.scan
    pa.field("precursor_scan", pa.int32(), nullable=True),   # "Precursor Full ScanNumber"
    # ── peptide identity ───────────────────────────────────────────────
    pa.field("sequence", pa.string(), nullable=False),       # bare AA "Sequence"
    pa.field("modified_sequence", pa.string(), nullable=True),  # ProForma 2.0 (normalized from "_..._")
    pa.field("peptide_id", pa.int64(), nullable=True),       # "Peptide ID" (→ peptides.txt)
    pa.field("mod_peptide_id", pa.int64(), nullable=True),   # "Mod. peptide ID"
    pa.field("evidence_id", pa.int64(), nullable=True),      # "Evidence ID" (→ evidence.txt)
    pa.field("proteins", pa.string(), nullable=True),        # ";"-delimited; or normalize to a list
    pa.field("charge", pa.int8(), nullable=False),           # "Charge"
    # ── measured vs theoretical ────────────────────────────────────────
    pa.field("mz", pa.float64(), nullable=True),             # "m/z"
    pa.field("mass", pa.float64(), nullable=True),           # "Mass"
    pa.field("mass_error_ppm", pa.float64(), nullable=True), # "Mass Error [ppm]"
    pa.field("retention_time", pa.float64(), nullable=True), # "Retention time" (min — unit!)
    # ── acquisition context (cross-check vs scan_metadata, PR #57) ──────
    pa.field("fragmentation", pa.string(), nullable=True),   # "Fragmentation" ~ scan_metadata.activation_type
    pa.field("mass_analyzer", pa.string(), nullable=True),   # "Mass analyzer" ~ scan_metadata.analyzer
    pa.field("psm_type", pa.string(), nullable=True),        # "Type" (MULTI-MSMS, ...)
    # ── scores / confidence ────────────────────────────────────────────
    pa.field("score", pa.float64(), nullable=True),          # Andromeda "Score"
    pa.field("delta_score", pa.float64(), nullable=True),    # "Delta score"
    pa.field("pep", pa.float64(), nullable=True),            # "PEP" (posterior error prob)
    pa.field("is_decoy", pa.bool_(), nullable=False),        # "Reverse" == "+" → True
    # ── fragment annotations (optional v1; ";"-delimited in source) ─────
    pa.field("fragment_matches", pa.string(), nullable=True),     # "Matches"  (consider list<string>)
    pa.field("fragment_intensities", pa.string(), nullable=True), # "Intensities" (consider list<float>)
    # ── engine tag (consistency w/ PEPTIDE_SCORE_TABLE.engine) ──────────
    pa.field("engine", pa.string(), nullable=False),         # "maxquant"
], metadata={b"schema_name": b"PsmTable", b"schema_version": b"1"})
```

Open questions for the planner to resolve (flag, don't guess):
- **RT units**: MaxQuant RT is **minutes**; `scan_metadata.rt` is **seconds** — normalize on read or
  record the unit? (Recommend normalize to seconds to match the bundle.)
- **`acquisition_id` linkage**: resolve `raw_file` → acquisition at read time (needs an `Acquisitions`
  table) or leave null and link in a later join stage? How to map MaxQuant `Raw file` basenames
  (no extension) onto protostar acquisition **stems** (they appear identical in the sample — verify
  across datasets, esp. Wilhelm HLA naming).
- **`proteins`/fragment lists**: keep `;`-delimited strings (lossless, simple) vs Arrow `list<...>`
  (queryable). Decide per column.
- **decoy/contaminant**: `Reverse` and `Potential contaminant` are `"+"`-or-empty flags — both → bool.
- **modseq → ProForma**: MaxQuant uses `_PEPTIDE_` with inline `(ox)`/`(ac)` and a `Modifications`
  summary col + `Fixed modifications` from `parameters.txt` (Carbamidomethyl C is fixed, not in the
  modseq). Reuse `io/msp/_mods.py`; the fixed-mod reconstruction is the subtle part.
- **PSM vs PEPTIDE_SCORE_TABLE**: does the reader also populate the existing
  `PEPTIDE_SCORE_TABLE`/`PROTEIN_SCORE_TABLE` from `peptides.txt`/`proteinGroups.txt`, or only
  `PSM_TABLE` in v1? (Recommend `PSM_TABLE` + provenance first; peptide/protein rollups in a follow-up.)

## Conventions (must hold)
PyArrow in memory + partitioned Parquet on disk; **no pandas** in-package (prefer `pyarrow.csv` for
the TSV `txt/` files; pandas only at the very edge if unavoidable); `snake_case`/`PascalCase`;
Python 3.12; new schemas via `core.io.schemas.register_schema`. Tests mandatory + DLL-free — follow
`tests/test_msp_reader.py` with small synthetic `txt/` fixtures under `tests/data/maxquant/`
(do **not** commit real ProteomeTools exports).

## Deliverables from the planning session
1. Tier-mapping decision (which registries the adapter(s) target) + rationale.
2. Final `PSM_TABLE` schema (refine the draft above; resolve the open questions).
3. Module/file layout under `massspec/io/maxquant/` + adapter self-registration plan.
4. Acquisition-association strategy (`raw_file` ↔ acquisition stem ↔ `scan`), incl. the analyzer/
   fragmentation cross-validation against `scan_metadata`.
5. Parsing approach + which `txt/` files land in v1 (`msms.txt` + `parameters.txt` priority;
   `evidence.txt` next; `peptides`/`proteinGroups`/`msmsScans` deferred) and why.
6. Test-fixture plan (synthetic `txt/` covering: target/decoy, fixed+variable mods, multi-protein,
   unidentified scan, FTMS vs ITMS rows).
7. Staged PR breakdown.

## Notes / provenance
- Format confirmed by `protostar/pipelines/probe_search_format.py`; headers sampled live from a
  Zolg2017 SEARCH zip on OSC (MaxQuant 1.5.3.30).
- A partial **local** MaxQuant copy exists on OSC ESS
  (`/fs/ess/PAS2254/RME/2026_spring/ASMS/proteome_tools_investigation/peptide_txt_files/`, ~1,392
  zips ≈ 95% of Zolg2017, 0% of the other two) — a fallback test corpus, though the canonical fetch
  pulls all SEARCH fresh from PRIDE.
- Scope guard: reader only. No re-running MaxQuant; no new scoring models. The `.msp` library reader
  already exists (`massspec.io.msp`).
```
