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
| 1 | Thermo `.raw` reader | `core.io` / `massspec.io` | `10_build_mzpeak`, `20_build_metadata` | needed (verify) | Env ships `pythonnet` + Thermo CommonCore DLLs; massspec CLAUDE.md lists `readers/` as TODO. **First action: confirm whether a working `.raw` reader already exists.** |
| 2 | mzpeak Parquet writer/reader + scanmeta (IIT/TIC/filter_string) | `massspec.io` | `10_build_mzpeak` | needed | Cartographer's mzpeak was an internal Parquet cache; rebuild fresh under Constellation. Per-scan `filter_string` needed for per-scan mode assignment. |
| 3 | MS1/MS2 chromatogram extraction + `rt_range` predicate pushdown + windowed MS1 scoring | `massspec` | `30_extract_intermediates` | needed | `MS1TensorResult` + `rt_range` pushdown were pruned from Cartographer (see massspec CLAUDE.md "worth surfacing"). |
| 4 | `HyperEMGPeak` (+ `WarpedEMGPeak`, `SplinePeak`) | `core.stats.peaks` | experiments (peak fitting) | needed | Only `GaussianPeak`/`EMGPeak` shipped; `core.optim.DifferentialEvolution` is available, so these are the next peak-numerics PR. |
| 5 | `GlobalCalibration` module (α₀/α₁/d_mz/ρ_R/ν_mz/α_mz, per file) | `massspec` | `04_procal_calibration` experiment | needed | Spec: Counter port doc §Architecture. |
| 6 | `Mass` module (per-progenitor: f_z, p_k, HyperEMG, ν_I, c_mz; MS2 channel) | `massspec` | Counter experiments | needed | The core progenitor object. |
| 7 | `Panel` + `PanelSet` (2D m/z×RT panels; additive-progenitor `forward()`; center-weighted gradients; Poisson zero-cell term) | `massspec` | Counter experiments | needed | Math: [`additive_progenitor_likelihood.md`](additive_progenitor_likelihood.md). |
| 8 | Seeded/global-agnostic scoring function + identification Λ (Counter v2 §10) | `massspec` | `05_counter_validation` | needed | `Λ_q = L_P(Q) − L_P(Q\{q})`. |
| 9 | Laplace credible intervals on N_total | `massspec` | uncertainty experiments | needed | 6×6 Hessian at MAP; well-conditioned now that α is stable. |
| 10 | EncyclopeDIA / Scribe search wrapper (reads `.raw` natively) | `thirdparty` + `massspec.io` | `15_reference_library` (optional re-search) | needed | EncyclopeDIA adapter (`massspec.io.encyclopedia`) + `thirdparty` registry exist; a search-invocation wrapper may need adding. |

## Reference architecture

- Counter engine architecture + Counter v2 → Constellation gap audit:
  `~/projects/constellation/docs/plans/cartographer-counter-port.md`.
- Additive-progenitor joint NLL: [`additive_progenitor_likelihood.md`](additive_progenitor_likelihood.md).
- Empirical model values + canonical-component checklist: [`model_specification.md`](model_specification.md)
  and [`../CLAUDE.md`](../CLAUDE.md).
