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
| Thermo `.raw` → parquet-bundle reader/writer (`convert`/`convert_batch`), scanmeta/IIT access | `constellation` (`massspec.io.thermo`) |
| MS1/MS2 chromatogram extraction + windowed scoring | `constellation.massspec` |
| Peak shapes (`HyperEMGPeak`, …) | `constellation.core.stats.peaks` |
| Optimizers (`DifferentialEvolution`, `LBFGSOptimizer`) | `constellation.core.optim` |
| Counter modules (`GlobalCalibration`, `Mass`, `Panel`, `PanelSet`, scoring, additive-progenitor NLL) | `constellation.massspec` |
| Laplace credible intervals, identification Λ | `constellation.massspec` |
| EncyclopeDIA/Scribe search wrapper | `constellation.thirdparty` / `massspec.io` |
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
│   ├── fetch/ convert/ library/ metadata/ intermediates/ experiments/
├── pipelines/                # CLI entry points + OSC/SLURM wrappers (one stage per file)
│   ├── 00_fetch_raw.py            # fresh download + hash verify; --seed-from optional
│   ├── 10_convert_raw.py          # .raw -> parquet bundle (proc/), rebuilt from scratch
│   ├── 15_reference_library.py    # ingest published .msp (+ optional re-search)
│   ├── 20_build_metadata.py       # acquisition time table from .raw headers
│   ├── 30_extract_intermediates.py
│   └── experiments/               # one reproducible script per solidified experiment
├── results/                  # the "textbook": curated MD records + figures + small tables
│   ├── 0X_*.md  +  figures/  # committed. Bulk parquet artifacts are gitignored / on ESS.
├── notebooks/                # scratch/exploration sandbox — NOT source of record
├── docs/                     # model_specification.md (empirical values); the ledger:
│   └── constellation_contributions.md
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
  PRIDE-manifest counts — see `config/manifests/`). Each acquisition has one paired MaxQuant
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
  on ESS (e.g. `/path/to/ProteomeTools/{ds}/raw/`) by
  name+size (SHA-1 with `--verify`) and relocates them instead of downloading. Default
  `--seed-mode move` (empties the source); `hardlink`/`copy` available.
- **Canonical data root (ESS):** `/fs/ess/<allocation>/<group>/protostar/data/` with
  `raw/{dataset}/`, `search/{dataset}/`, `proc/{dataset}/{centroid,profile}/<stem>/`,
  `libraries/<mode>/`. The `proc/` parquet bundles are **rebuilt from scratch** with the
  latest Constellation reader (no reuse of prior caches) for downstream consistency.
- **Fragmentation mode is a per-scan property** derived from the scan filter string /
  scanmeta — not from pre-split files. Each pool was acquired as 3 injections (`.raw` files);
  the 9 modes (CID35/HCD20–35 × DDA/Targeted × IT/Orbi) interleave within them.
- **PROCAL:** 40 synthetic calibrant peptides spiked into every run — calibration anchors and
  high-N statistical targets.

---

## Canonical MS1 ion model — required components (guardrail)

Any experiment that fits N(t) or scores MS1 observations **must** include all of the
following. These were established across the prior analysis (nb13–nb42b) and validated
extensively. Do not simplify, omit, or approximate — partial models produce misleading
results. The implementations live in Constellation; this checklist is the acceptance gate for
any experiment here.

1. **Charge-dependent α(z):** `α(z) = softplus(α₀ + α₁·z)`, linear in z (charge-only, no IIT
   power-law). **Fit per file.** Intensity is per-time; ion count `N = I·iit/α(z)`. Reference:
   DDA 60K `α(z=2) ≈ 19–22` (MLE); targeted 120K `α₀ ≈ −6, α₁ ≈ 16.9`.
2. **IIT in the variance (not on α):** `Var[I_k] ∝ α/iit`. The missing `/iit` caused 10–100×
   α inflation (the nb42a bug). IIT from scanmeta (`*.scanmeta`, level-1 `iit` ms).
3. **Fitted isotope fractions p_k (not theoretical):** `expected_isotope_envelope()` has a
   systematic mass-dependent bias (dp_0 ≈ +0.0075, dp_2 ≈ −0.0064). Use per-peptide fitted
   fractions; renormalize to sum 1.
4. **m/z error decomposition:** (a) per-file isotope-spacing correction
   `da_err_k − k·d_mz_k/charge`; (b) per-scan systematic offset `μ̂(t)` = mean residual across
   ions in the scan (calibration drift, ~0.47 ppm between-scan σ — NOT interference);
   (c) per-ion residual `ε_k = mz_err_corr − μ̂(t)`, the interference-sensitive quantity,
   `ε_k ~ t_ν_mz(0, √(c_mz/N_k^α_mz))` with c_mz per-peptide, α_mz ≈ 0.93, ν_mz ≈ 3.06.
5. **Per-ion intensity model (Student-t):** `I_pred = N(t)·p_k·α(z)`;
   `Var[I_k] = α(z)·ΣI_pred·p_k·(1−p_k)/iit`; `I_obs ~ t_ν_I(I_pred, √Var)` with ν_I
   **per-peptide** (median ~6.5 at 60K DDA; resolution-scales as (R/R_ref)^η, η≈0.87).
6. **Sliding-window contextual features** (W=15 leave-one-out): `μ_smooth`, `μ_sd`,
   `kl_smooth`.
7. **KL divergence with intensity-weighted charge aggregation** against fitted fractions;
   missing isotopes filled with pseudocount; sigmoid + charge weights.

**Optimizer rule:** for any `nn.Module` fit, use Constellation's
`core.optim.DifferentialEvolution` — **never** `scipy.optimize.differential_evolution` (the
nb42c lesson: scipy-DE produced optimizer artifacts mistaken for structural findings).

**Counter / interference:** the current direction is the **additive-progenitor** joint
likelihood (interference as a sum of co-contributing progenitors, not a binary signal/noise
mixture; zero observations as proper Poisson log-probabilities; overlapping center-weighted
2D m/z×RT panels). The math is in `docs/` (and Constellation's
`docs/plans/cartographer-counter-port.md`). **nb42c is quarantined** — wrong DE, improper
censored scaling, binary-mixture mis-specification; do not resurrect it.

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

- **Project dir (ESS):** `/fs/ess/<allocation>/<group>/protostar/` (data root, results, job logs). Do
  **not** store project files in OSC home (strict storage/file-count limits).
- **Account/cluster:** account `<allocation>`, home account `<home-allocation>`, default cluster Cardinal
  (CPU); Ascend only for heavy GPU. Config in `config/osc.toml` (mirrors `~/.osc_helper/`).
- **`osc_helper`** (optional `[osc]` extra) provides SSH / SLURM / rsync:
  `load_config()`, `ssh.check_connection(host)`, `jobs.{write_job_script,submit_job,job_status}`,
  `sync.push_code(...)` (syncs to OSC **home**, not ESS — `cp` to ESS manually).
- **SLURM:** `set -eo pipefail` (NOT `-euo` — conda breaks on unbound vars);
  `module load miniconda3/24.1.2-py310` + `source "$(conda info --base)/etc/profile.d/conda.sh"`.
  Prefer extended single allocations over many small jobs. Watchers run on the compute node.

---

## Status

Data layer implemented. `00_fetch_raw` (PRIDE v3 manifest + resumable download + verify +
`--seed-from` relocation), `10_convert_raw` (drives `convert_batch` into `proc/`, centroid +
profile, resume + SLURM sharding), and `15_reference_library` (Zenodo `.msp` fetch + extract)
are written and committed; the per-dataset expected manifests live in `config/manifests/`. The
SEARCH files are confirmed MaxQuant (`probe_search_format.py`) → ledger #11 targets a
`massspec.io.maxquant` reader. Constellation's Thermo `.raw` reader is verified production-ready
(ledger #1/#2 landed).

Next (run on OSC): seed → download → convert (centroid all three, then the profile pass after a
1–2 file validation). Then `20_build_metadata` → `30_extract_intermediates`, driving the Counter
modules into Constellation as the experiments need them.
