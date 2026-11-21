# protostar — Claude Code Guide

**ProToStaR** (**Pro**teome **To**ols **Sta**tistical **R**esearch) is a reproducible,
open-science re-analysis of the ProteomeTools synthetic-peptide datasets, built on the
**Constellation** platform. It characterizes the physiochemistry of peptide fragmentation
and the statistics of our mass-spectrometry measurements — the calibrated priors and
"confident findings" that seed the lab's integrative models, chiefly the **Counter**
generative quantification model.

The repo name is lowercase `protostar` (Python convention); `ProToStaR` is the stylistic
rendering in prose. A protostar is the collapsing precursor from which a main-sequence star
forms — the foundational stage from which the lab's mature MS tooling and priors condense.
It lives inside the Constellation.

---

## The one rule: no model code here

**ALL model / likelihood / peak-shape / scoring code lives in Constellation**
(`constellation.massspec` + `constellation.core`) — **no exceptions.** This repository holds
only:

- **data orchestration** (fetch, convert, reference-library ingest, metadata curation),
- **experiment drivers** (scripts that call Constellation and produce artifacts), and
- **curated records** (`results/` — the "textbook").

When an experiment needs a modeling capability Constellation lacks, **build it in
Constellation first**, then import it here. Data-management *glue* specific to ProteomeTools
(file naming, PRIDE accessions, injection structure, the acquisition time table) legitimately
lives here; if a piece turns out to be generic, it is a candidate to upstream.

Track every needed-but-missing fundamental in **`docs/constellation_contributions.md`** with
its status (needed / PR-open / landed). An experiment that depends on one is gated on the
upstream landing. Short-lived local prototypes are allowed **only** with a tracked promotion
entry and **never** ship inside a `results/` record.

### What lives where

| Capability | Home |
|---|---|
| Thermo `.raw` → parquet-bundle reader/writer (`convert`/`convert_batch`), scanmeta/IIT access | `constellation` (`massspec.readers.thermo`) |
| Bundle access (v3: `read_bundle_manifest`, `peaks_dataset`, `points_dataset`, `is_profile_bundle`) | `constellation.massspec.bundle` |
| MS1/MS2 chromatogram extraction (sweep `extract_xic`, `extract_xic_scan_major`) + scoring | `constellation.massspec.quant` / `massspec.sweep` |
| Peak shapes (`HyperEMGPeak`, …) | `constellation.core.stats.peaks` |
| Optimizers (`DifferentialEvolution`, `LBFGSOptimizer`) | `constellation.core.optim` |
| Counter modules (`GlobalCalibration`, `Progenitor`, `Panel`, discovery/components/attribution, `estimate_n`) | `constellation.massspec.counter` |
| Credible intervals, identification Λ | `constellation.massspec.counter` |
| EncyclopeDIA/ScribeTwo search, MSP/MaxQuant readers | `constellation.massspec.search` / `massspec.io` + `constellation.thirdparty` |
| ProteomeTools fetch/verify, accessions, injection structure | **protostar** |
| Acquisition time-table curation | **protostar** (persists a `massspec.acquisitions.Acquisitions` table) |
| Experiment scripts + curated `results/` records | **protostar** |

---

## Repository layout

```
protostar/
├── CLAUDE.md                 # this file
├── README.md                 # open-science framing; reproduce end-to-end
├── pyproject.toml            # package `protostar`; dep constellation-bio[ms]
├── environment.yml           # env (lab: reuse the `constellation` conda env)
├── config/
│   ├── datasets.toml         # per-dataset accessions, .raw lists+checksums, .msp URLs
│   └── osc.toml              # ESS/home paths, account
├── protostar/                # slim glue package over constellation — NO model code
│   ├── bundles.py            # the ONLY way to locate/read proc/ bundles (v3 layout, refuses pre-v3)
│   ├── fetch/ convert/ library/ metadata/ intermediates/ experiments/
├── pipelines/                # CLI entry points + OSC/SLURM wrappers (one stage per file)
│   ├── 00_fetch_raw.py            # fresh download + hash verify; --seed-from optional
│   ├── 10_convert_raw.py          # .raw -> parquet bundle (proc/), rebuilt from scratch
│   ├── 15_reference_library.py    # ingest published .msp (+ optional re-search)
│   ├── 20_build_metadata.py       # acquisition time table from proc/ bundle metadata
│   ├── 30_extract_intermediates.py
│   └── experiments/               # one reproducible script per solidified experiment
├── results/                  # the "textbook": curated MD records + figures + small tables
│   ├── 0X_*.md  +  figures/  # committed. Bulk parquet artifacts are gitignored / on ESS.
├── notebooks/                # scratch/exploration sandbox — NOT source of record
├── docs/                     # records + specs:
│   ├── counter_model_decisions.md       # WHY each Counter feature is in/out (the decision log)
│   ├── history.md                       # dated history across MS_Quant/Cartographer/Constellation/protostar
│   ├── constellation_contributions.md   # the ledger (needed / landed upstream capabilities)
│   └── model_specification.md           # historical MS_Quant spec (stale values — see its banner)
└── data/                     # small local dev subsets only (gitignored); canonical on ESS
```

**`results/` is of record; `notebooks/` is scratch.** Interpretation is separated from code:
experiment scripts under `pipelines/experiments/` produce parquet + figures; the curated
conclusions + the figures they cite are synthesized into `results/*.md`. Bulk parquet artifacts
are gitignored and canonical on ESS; the `.md` chapters and `figures/` are committed.

---

## Data: raw-first

The **source of truth is the `.raw` Thermo files.** This is a hard departure from the prior
analysis (which keyed off mzML files and FragPipe searches) — none of that carries over.

- **Datasets & accessions:** Zolg2017 = `PXD004732`, Gessulat2019 = `PXD010595`,
  Wilhelm2021 = `PXD021013`. **4,213 raw acquisitions, ~2.9 TB** (RAW: 1460 / 888 / 1865;
  PRIDE-manifest counts — see `config/manifests/`; these span **all four** injection types
  **including ETD**). Each acquisition has one paired MaxQuant
  SEARCH zip (~110 GB more). NB: this corrects the earlier "~9,612 / ~7.8 TB" estimate, which
  appears to have included the out-of-scope *modified*-peptide pools. Wilhelm2021 publishes no
  checksums, so its files are size-verified only.
- **Reference library:** the published ProteomeTools `.msp` libraries from
  <https://www.proteometools.org/index.php?id=53>, ingested via `massspec.io.msp` — these
  replace the old individual searches. Optionally re-searched with EncyclopeDIA/Scribe on the
  local `.raw` via a Constellation wrapper.
- **Fetch is fresh-download + verify (canonical).** `pipelines/00_fetch_raw.py` builds a
  manifest from `config/datasets.toml` by querying the ProteomeXchange/PRIDE API for each
  dataset's `.raw` list + published checksums, downloads missing files, verifies each, and
  re-fetches on mismatch (resumable + repairable; `--dry-run` reports present/missing/corrupt).
  Repo scripts reflect *what a fresh user would do*.
- **`--seed-from <dir>` (local time-saver, not the reproducible path):** matches files already
  on disk (an existing ProteomeTools `.raw` copy; path from `--seed-from`,
  `$PROTOSTAR_SEED_FROM`, or `config/osc.toml`) by name+size (SHA-1 with `--verify`) and
  relocates them instead of downloading. Default `--seed-mode move` (empties the source);
  `hardlink`/`copy` available.
- **Canonical data root:** `<data_root>/` (set via `--data-root`, `$PROTOSTAR_DATA_ROOT`, or
  `config/osc.toml` — the lab points it at project ESS storage) with `raw/{dataset}/`,
  `search/{dataset}/`, `proc/{dataset}/<stem>/` + `proc/{dataset}/<stem>.profile/` (Constellation
  **bundle v3**, written in one pass), and `libraries/<mode>/`.
  - The `proc/` bundles are **rebuilt from scratch** with the latest Constellation reader (no reuse
    of prior caches) for downstream consistency.
  - **Read them only through `protostar.bundles`.** It knows the v3 layout and refuses pre-v3
    bundles. Never open `peaks.parquet` or glob `*/manifest.json` directly.
  - The June 2026 pre-v3 trees `proc/{dataset}/{centroid,profile}/` are legacy and are deleted per
    dataset once v3 verifies.
  - Only MS1 survey scans are recorded in **profile** mode. **All MS2 (Orbitrap, ion-trap, ETD) are
    centroided at acquisition** (Zolg et al. 2017), so the `.profile` bundle carries MS1 only.
- **Fragmentation mode is a per-scan property** derived from the scan filter string /
  scanmeta — not from pre-split files. Each pool was acquired as **4 injections** (`.raw`
  files): `DDA`, `2xIT_2xHCD`, `3xHCD`, `ETD` (Wilhelm2021 omits ETD for 60 pools). The 9
  modeled modes (CID35/HCD20–35 × DDA/Targeted × IT/Orbi) are **HCD/CID only** and interleave
  within the first three injections. **ETD is converted and kept for posterity:** its MS1 is
  the same full-profile Orbitrap survey (informative for the MS1 ion model + PROCAL) and the
  reader parses its EThcD/ETciD MS2 — but ETD MS2 falls outside the 9 modeled modes and has no
  published `.msp` reference. HCD vs CID is **not** recoverable from collision energy alone
  (CID35 and HCD35 share NCE 35); the activation type is an explicit filter-string token,
  surfaced as the `activation_type` + `analyzer` columns in `scan_metadata` (Constellation
  contribution — see `docs/constellation_contributions.md`).
- **PROCAL:** 40 synthetic calibrant peptides spiked into every run — calibration anchors and
  high-N statistical targets.

---

## Counter model — the decision log

**The status of every Counter model feature lives in
[`docs/counter_model_decisions.md`](docs/counter_model_decisions.md)**: what is Adopted, Under test,
Heuristic or Open, with the evidence and its scope. Read it before fitting N(t) or scoring MS1
observations. Change a feature's status only by adding evidence and a dated History line there. Never
edit model values into this file.

This replaces the "canonical MS1 ion model — required components" checklist, which was removed on
2026-11-20. That checklist was copied at the scaffold from an April MS_Quant agent guardrail. Its m/z
values pre-date the IIT fix, its μ̂ description did not match the model code, and two of its seven
"required" items were detection heuristics. See [`docs/history.md`](docs/history.md).

Hard rules, independent of the log's open questions:

- **IIT goes in the variance, not on α.** Ion count is `N = I·iit/α(z)` (D2).
- **α(z) is fit per file.** Do not borrow gain values across files without a test (D1).
- **Optimizer:** use Constellation's `core.optim.DifferentialEvolution`, plus LBFGS polish. **Never**
  use `scipy.optimize.differential_evolution`; in nb42c it produced artifacts that were mistaken for
  findings (D16).
- **nb42c is quarantined.** It used the wrong DE, an improper censored scaling and a mis-specified
  binary mixture. Do not resurrect it.
- **No partial models presented as results.** If an experiment omits an Adopted feature, it says so in
  its `results/` record and is labelled as a component study, not a model fit.
- **Evidence scope is stated.** The reference fit (nb42b) is one file. Single-file numbers are labelled
  as such.

**Interference** is **Under test** (D9). Two framings compete: additive progenitors (built in
Constellation, never run on PT data) and a per-point binary interfered/clean state (the PI's preferred
framing for PT reprocessing). Neither is the default until Test D. The MS_Quant notebooks are archived
read-only in `wilburnlab/ms_quant_archive`; cite them, do not extend them.

---

## Conventions (inherited from Constellation)

- **PyArrow in memory, partitioned Parquet on disk.** Never `pa.concat_tables(...)` between
  pipeline stages — open shard directories as `pa.dataset.dataset(...)` and stream / filter.
  No pandas inside the package (pandas stays at the edge).
- **PyTorch for numerical work**; `Parametric` ABC unifies distributions + peak shapes.
- **HUPO-PSI** external forms: ProForma 2.0 (modseqs), mzPAF (peak annotations), USI.
- **`snake_case`** functions, **`PascalCase`** classes. Default tolerance 20 ppm,
  `tolerance_unit ∈ {'ppm','Da'}`. Float sentinel `-1.0` = "not observed".
- **Python 3.12**, conda env `constellation` (lab) or `protostar` (external).

---

## OSC operational notes

- **Project dir:** keep the data root, results, and job logs on project ESS storage — **not**
  OSC home (strict storage/file-count limits). Concrete paths/allocations live in a local
  `config/osc.toml` (gitignored; copy `config/osc.example.toml` and fill in).
- **Account/cluster:** set your compute + home allocations and cluster in `config/osc.toml`
  (mirrors `~/.osc_helper/`); default cluster Cardinal (CPU), Ascend only for heavy GPU.
- **`osc_helper`** (optional `[osc]` extra) provides SSH / SLURM / rsync:
  `load_config()`, `ssh.check_connection(host)`, `jobs.{write_job_script,submit_job,job_status}`,
  `sync.push_code(...)` (syncs to OSC **home**, not ESS — `cp` to ESS manually).
- **SLURM:** `set -eo pipefail` (NOT `-euo` — conda breaks on unbound vars);
  `module load miniconda3/24.1.2-py310` + `source "$(conda info --base)/etc/profile.d/conda.sh"`.
  Prefer extended single allocations over many small jobs. Watchers run on the compute node.
- **libstdc++ pin (required):** after `conda activate constellation`, run
  `export LD_LIBRARY_PATH="$CONDA_PREFIX/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"` — else
  `module load` shadows the env's libstdc++ and constellation import dies with a
  `CXXABI_1.3.15` mismatch. The pipeline scripts self-heal via a one-time re-exec
  (`pipelines/_common.py`), but ad-hoc `python -c` / `convert.sbatch` set it explicitly.
- **Convert is ESS-I/O-bound:** use **`--workers 16`** (`--cpus-per-task=16`), not the whole
  node — 16 captures ~71% of 32-worker throughput at half the cores, and ≥32 risks CLR-spawn
  flakiness. The Thermo DLLs are loaded only inside spawn workers (never fork a CLR-initialised
  parent). `.raw` is the source of truth; ETD is converted + kept (MS1 + EThcD/ETciD MS2).

---

## Status

*Rewritten 2026-11-20 at the project restart. See [`docs/history.md`](docs/history.md) for the dated
record.*

**The data layer was completed in June 2026, but on a bundle format current Constellation refuses.**

- `00_fetch_raw`: all **4,213 `.raw`** files (seeded from the Cartographer tree) and all MaxQuant
  SEARCH zips (1460/888/1865) are on ESS.
- `10_convert_raw`: **4,212 / 4,213** were converted in June, centroid and profile as two passes. The
  one miss is a documented corrupt PRIDE source (Wilhelm2021 `…pool_122…3xHCD-1hnoincl-R1.raw`,
  57,874 bytes; `docs/osc_runbook.md` "Known-skip").
- **These bundles predate Constellation bundle v3** (#104, 2026-10-31). Any Constellation API given
  one raises `OutdatedBundleError`, and there is no upgrader. **The corpus must be reconverted**,
  after the convert driver is migrated to `convert_batch(keep_profile=…)`. The plan is a ~6-file pilot
  to size it, then full reconversion; v1 bundles are deleted per dataset only after v3 verifies and
  the PI confirms.
- `15_reference_library` (all 8 `.msp` modes) and `20_build_metadata` (acquisition time table) were
  done in June. Stage 20 must be re-run on v3 bundles.

**Counter: the port exists but has not been validated.** Constellation `massspec.counter` (#69–#86,
June) was validated on synthetic data plus one 10-PROCAL real-data smoke (2026-06-11). The original
notebook model (nb42b) is itself a single-file fit. The two differ, most importantly in the intensity
likelihood: Constellation is Poisson-only, while the PI adopted Student-t on 2026-11-20 (decision log
D5). The upstream fixes and features are specified as change requests (CR-*) in Constellation's
`docs/plans/counter-pt-reconciliation.md`, for a separate Constellation session to implement.

**Stage 30 and the Counter-motivation figures were built but not run at scale.**
`30_extract_intermediates` (`protostar/intermediates/ms1_extract.py`) extracts per-(peptide, charge,
scan) MS1 isotope chromatograms for PROCAL and known-clear peptides. Only shard 0 of 16 was run
(06-03). The four empirical figure scripts (`pipelines/experiments/10–14`) therefore rest on 1/16 of
the data. The α̂(z) moment estimate (script 11) is confounded and is not to be published.

**Next** (restart plan):
1. Prep PR (CI pin, PSM_TABLE v2 fixture).
2. History and decision-log docs.
3. Pipeline migration to bundle v3 (`protostar/bundles.py`; convert, metadata and stage 30 on the
   bundle API).
4. Constellation change requests, handed off as a Constellation plan (not tracked here).
5. OSC pilot, reconversion and the full stage-30 run.
6. Land the June backlog.
7. Tests A–E (intensity likelihood, m/z offset level, space charge, binary interference state, MS1→MS2
   conversion).

Operational facts (in `docs/osc_runbook.md`): the libstdc++ `LD_LIBRARY_PATH` pin; convert is
**ESS-I/O-bound**, so use 16 workers on an exclusive node (32 or more risks CLR-spawn flakiness); all
MS2 is centroided at acquisition, so profile data covers MS1 only; ETD is converted and kept.
