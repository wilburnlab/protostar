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
| 3 | MS1/MS2 chromatogram extraction + `rt_range` predicate pushdown + windowed MS1 scoring | `massspec` | `30_extract_intermediates` | needed | `MS1TensorResult` + `rt_range` pushdown were pruned from Cartographer (see massspec CLAUDE.md "worth surfacing"). |
| 4 | `HyperEMGPeak` (+ `WarpedEMGPeak`, `SplinePeak`) | `core.stats.peaks` | experiments (peak fitting) | needed | Only `GaussianPeak`/`EMGPeak` shipped; `core.optim.DifferentialEvolution` is available, so these are the next peak-numerics PR. |
| 5 | `GlobalCalibration` module (α₀/α₁/d_mz/ρ_R/ν_mz/α_mz, per file) | `massspec` | `04_procal_calibration` experiment | needed | Spec: Counter port doc §Architecture. |
| 6 | `Mass` module (per-progenitor: f_z, p_k, HyperEMG, ν_I, c_mz; MS2 channel) | `massspec` | Counter experiments | needed | The core progenitor object. |
| 7 | `Panel` + `PanelSet` (2D m/z×RT panels; additive-progenitor `forward()`; center-weighted gradients; Poisson zero-cell term) | `massspec` | Counter experiments | needed | Math: [`additive_progenitor_likelihood.md`](additive_progenitor_likelihood.md). |
| 8 | Seeded/global-agnostic scoring function + identification Λ (Counter v2 §10) | `massspec` | `05_counter_validation` | needed | `Λ_q = L_P(Q) − L_P(Q\{q})`. |
| 9 | Laplace credible intervals on N_total | `massspec` | uncertainty experiments | needed | 6×6 Hessian at MAP; well-conditioned now that α is stable. |
| 10 | EncyclopeDIA / Scribe search wrapper (reads `.raw` natively) | `thirdparty` + `massspec.io` | `15_reference_library` (optional re-search) | needed | EncyclopeDIA adapter (`massspec.io.encyclopedia`) + `thirdparty` registry exist; a search-invocation wrapper may need adding. |
| 11 | MaxQuant search-output reader (`txt/` export: `msms`/`evidence`/`peptides`/`parameters`/…) | `massspec.io.maxquant` | search↔raw association (future) | needed | Confirmed by `pipelines/probe_search_format.py`: all three datasets' PRIDE SEARCH zips are MaxQuant combined-`txt` exports. Needed to associate identifications with acquisitions (out of scope for the data-layer session). |

## Reference architecture

- Counter engine architecture + Counter v2 → Constellation gap audit:
  `~/projects/constellation/docs/plans/cartographer-counter-port.md`.
- Additive-progenitor joint NLL: [`additive_progenitor_likelihood.md`](additive_progenitor_likelihood.md).
- Empirical model values + canonical-component checklist: [`model_specification.md`](model_specification.md)
  and [`../CLAUDE.md`](../CLAUDE.md).
