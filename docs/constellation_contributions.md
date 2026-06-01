# Constellation contribution ledger

protostar holds **no model code** (see [`../CLAUDE.md`](../CLAUDE.md)). Every modeling /
likelihood / peak-shape / scoring / reader capability that an experiment needs but
Constellation does not yet provide is tracked here. An experiment that depends on an item is
**gated on that item landing upstream**. Short-lived local prototypes are allowed only with a
row here (status `prototyping`) and a promotion plan, and must never ship inside a `results/`
record.

**Status:** `needed` (not started) · `prototyping` (local, must promote) · `pr-open`
(upstream review) · `landed` (available in Constellation; cite the module).

| # | Capability | Target location in Constellation | Drives which protostar stage | Status | Notes |
|---|---|---|---|---|---|
| 1 | Thermo `.raw` reader (`ThermoReader`, `convert`) | `massspec.io.thermo` | `10_convert_raw`, `20_build_metadata` | **landed** | Verified production-ready: `constellation massspec convert` CLI + `ThermoReader`/`convert_batch`; centroid (per-peak resolution/noise/baseline) + profile (raw FT grid). |
| 2 | Convert bundle: `peaks.parquet` + `scan_metadata.parquet` (IIT/TIC/`filter_string`) + `acquisition_metadata.parquet` + `manifest.json` | `massspec.io.thermo` | `10_convert_raw` | **landed** | `convert`/`convert_batch`; RT-binned row groups; source SHA-256 in manifest; per-scan `filter_string` recovers fragmentation mode. ("mzpeak" retired — no HUPO-PSI standard exists.) |
| 3 | MS1/MS2 chromatogram extraction + **MS2 spectrum retrieval** (per-scan centroid peak lists by precursor isolation m/z + RT + scan-filter mode) + `rt_range` predicate pushdown + windowed MS1 scoring | `massspec` | `30_extract_intermediates`, exp 01–08 | needed | Port from `cartographer/data/chromatogram.py` + `data/spectra.py`. `MS1TensorResult` + `rt_range` pushdown were pruned from Cartographer (see massspec CLAUDE.md "worth surfacing"). **MS2-retrieval slice is Lane A** (Part I); MS1-extraction + windowed-scoring slice is Lane B (Part II). |
| 4 | `HyperEMGPeak` (+ `WarpedEMGPeak`, `SplinePeak`) | `core.stats.peaks` | experiments (peak fitting) | needed | Only `GaussianPeak`/`EMGPeak` shipped; `core.optim.DifferentialEvolution` is available, so these are the next peak-numerics PR. |
| 5 | `GlobalCalibration` module (α₀/α₁/d_mz/ρ_R/ν_mz/α_mz, per file) | `massspec` | `04_procal_calibration` experiment | needed | Spec: Counter port doc §Architecture. |
| 6 | `Mass` module (per-progenitor: f_z, p_k, HyperEMG, ν_I, c_mz; MS2 channel) | `massspec` | Counter experiments | needed | The core progenitor object. |
| 7 | `Panel` + `PanelSet` (2D m/z×RT panels; additive-progenitor `forward()`; center-weighted gradients; Poisson zero-cell term) | `massspec` | Counter experiments | needed | Math: [`additive_progenitor_likelihood.md`](additive_progenitor_likelihood.md). |
| 8 | Seeded/global-agnostic scoring function + identification Λ (Counter v2 §10) | `massspec` | `05_counter_validation` | needed | `Λ_q = L_P(Q) − L_P(Q\{q})`. |
| 9 | Laplace credible intervals on N_total | `massspec` | uncertainty experiments | needed | 6×6 Hessian at MAP; well-conditioned now that α is stable. |
| 10 | EncyclopeDIA / Scribe search wrapper (reads `.raw` natively) | `thirdparty` + `massspec.io` | `15_reference_library` (optional re-search) | needed | EncyclopeDIA adapter (`massspec.io.encyclopedia`) + `thirdparty` registry exist; a search-invocation wrapper may need adding. Optional re-search **cross-check to #11**; not on the critical path. |
| 11 | MaxQuant search-output reader (`txt/` export: `msms`/`evidence`/`peptides`/`parameters`/…) | `massspec.io.maxquant` | exp 01–02 (PROCAL & extend identification); search↔raw association | needed | **Foundational / Lane A** — reprioritized. Confirmed by `pipelines/probe_search_format.py`: all three datasets' PRIDE SEARCH zips are MaxQuant combined-`txt` exports, already paired with every raw file. Reading `msms.txt`/`evidence.txt` gives per-scan peptide assignments with no re-search compute — the identification path that anchors every experiment. |
| 12 | MS2 spectral-similarity suite: cosine similarity, normalized dot product, spectral-entropy similarity (Li 2021), + a unified `compare_spectra` over aligned vectors | `core.stats.losses` (extend) / new `massspec.spectra.similarity` | exp 01 (`01_ms2_spectral_scoring`) | needed | **Lane A, net-new but small.** Reuse `kld`, `spectral_angle`, `l1/l2_normalize` already in `core/stats/losses.py`; only raw cosine / dot / entropy-*similarity* + the comparator are missing (no `massspec.spectra` package exists yet). |
| 13 | Consensus / aggregated-spectrum builder (align replicate spectra to a reference fragment ladder; aggregate intensities mean/median + per-fragment dispersion) | new `massspec.spectra.consensus` | exp 01–02 | needed | **Lane A, net-new.** Reuse `match_mz`/`assign_fragments` (`massspec.peptide.match`); reference impl `cartographer/data/spectra.py`. |
| 14 | Multinomial fragment-intensity **generative** model (peptide → expected fragment propensities → `Multinomial` over channels) | `massspec` | exp 01 (extend) | needed | **Deferred (later phase).** Near-term exp 02 only *characterizes* multinomial shot noise using the existing `Multinomial` (`core/stats/distributions.py`); this predictive model is not built yet. |
| 15 | Graceful skip of empty/aborted `.raw` in `convert`/`convert_batch` (catch `SelectInstrument` `ArgumentOutOfRangeException` → emit a `skipped`/`corrupt` BatchResult, not an unhandled .NET stack trace) | `massspec.io.thermo` | `10_convert_raw` | needed | **Low-priority robustness.** One real case: Wilhelm2021 `…pool_122…3xHCD-1hnoincl-R1.raw` is a 57,874-byte aborted acquisition (valid OLE2 `Finnigan` magic, no MS device) that PRIDE hosts as-is (re-download byte-identical; no published checksum). `convert_batch` should classify it as a clean skip so a full pass reports 4,212/4,213 without a hard error. Documented in [`osc_runbook.md`](osc_runbook.md). |

## Build lanes (priority)

The MS2-first manuscript ([`manuscript_roadmap.md`](manuscript_roadmap.md)) reprioritizes these
items into three lanes by what they unblock:

- **Lane A — foundational, unblocks Part I (MS2).** #11 (MaxQuant reader — identification),
  #12 (spectral-similarity suite), #13 (consensus builder), and the **MS2-retrieval slice of
  #3**. Mostly net-new and small; the minimal path to the first runnable experiment.
- **Lane B — unblocks Part II (MS1 error structure).** #3 (MS1 extraction + windowed scoring),
  the `error_model` α(z)/Student-t likelihoods, #5 (GlobalCalibration).
- **Lane C — unblocks Part III (Counter).** #4 (HyperEMG peak shapes), #6 (`Mass`), #7
  (`Panel`/`PanelSet`), #8 (scoring + Λ), #9 (Laplace credible intervals).
- **Deferred:** #14 (multinomial generative), #10 (EncyclopeDIA/Scribe re-search).

**Items #3–#9 are a Cartographer→Constellation port, not new invention.** The predecessor
engine lives at `/home/dbwilburn/projects/cartographer/cartographer/` — `calibration.py` (#5),
`error_model.py` (the α(z) + Student-t m/z & intensity likelihoods under #5/#6), `masses.py`
(#6), `scoring.py` (#8), `data/chromatogram.py` (#3), `data/spectra.py` (MS2 retrieval) — and
the migration is scoped in Constellation's `docs/plans/cartographer-counter-port.md`. Effort is
*port-and-validate* against the nb42b empirics in
[`model_specification.md`](model_specification.md). Caveat: that port doc rates the low-level
likelihoods "exists cleanly," but that describes Cartographer/notebooks — they are **not yet
importable from Constellation**.

## Reference architecture

- Counter engine architecture + Counter v2 → Constellation gap audit:
  `~/projects/constellation/docs/plans/cartographer-counter-port.md`.
- Additive-progenitor joint NLL: [`additive_progenitor_likelihood.md`](additive_progenitor_likelihood.md).
- Empirical model values + canonical-component checklist: [`model_specification.md`](model_specification.md)
  and [`../CLAUDE.md`](../CLAUDE.md).
