# Constellation contribution ledger

protostar holds **no model code** (see [`../CLAUDE.md`](../CLAUDE.md)). This ledger tracks every
modelling, likelihood, peak-shape, scoring or reader capability that an experiment needs but
Constellation does not yet provide.

- **Gating:** an experiment that depends on an item waits until the item lands upstream.
- **Local prototypes:** allowed only with a row here (status `prototyping`) and a plan for moving the
  code upstream. They never ship inside a `results/` record.

**Status:** `needed` (not started) · `prototyping` (local; must move upstream) · `pr-open` (upstream
review) · `landed` (available in Constellation; cite the module).

> **Ledger refreshed 2026-11-20** against `constellation@77f540f` (main). Each claim was checked
> against the source.
>
> This refresh corrects the 2026-09-05 refresh, which overstated the Counter port:
> - **"Ported and validated" is wrong.** The port (#69–#86) was validated on synthetic data plus **one**
>   real-data smoke (10 PROCAL, 2026-06-11).
> - PanelSet, center-weighted gradients and identification Λ were **never built**.
> - The CLI's interval path is **Laplace**, not VB.
> - The intensity term is **Poisson-only**; the Student-t code exists but is never called.
>
> It also records the bundle-v3 break. Since #104 (2026-10-31), bundles written before then raise
> `OutdatedBundleError`, so protostar's June corpus must be reconverted.
>
> Open Counter work is specified as change requests (CR-*) in Constellation's
> [`docs/plans/counter-pt-reconciliation.md`](https://github.com/wilburnlab/constellation/blob/main/docs/plans/counter-pt-reconciliation.md)
> (committed with the first PR of that set),
> and feature status lives in [`counter_model_decisions.md`](counter_model_decisions.md) (D-numbers).

| # | Capability | Target location in Constellation | Drives which protostar stage | Status | Notes |
|---|---|---|---|---|---|
| 1 | Thermo `.raw` reader (`ThermoReader`, `convert`, `convert_batch`) | `massspec.readers.thermo` | `10_convert_raw`, `20_build_metadata` | **landed** (API changed) | Revised in #103 (2026-10-27) and #104 (2026-10-31).<br>• **`profile=` was removed.** The call is now `convert_batch(..., keep_profile=True, centroid_algo=…, precursor_major=True, allow_empty_levels=False)`.<br>• **One pass** writes the centroid bundle plus a sibling `<stem>.profile/` bundle.<br>• protostar's convert driver still passes `profile=`, so it raises `TypeError`. Being migrated (restart PR "bundle-v3 pipelines"). |
| 2 | Convert bundle: peaks, `scan_metadata`, `acquisition_metadata`, `manifest.json` | `massspec.bundle` + `massspec.readers.thermo` | `10_convert_raw` and every reader of bundles | **landed** (v3) | **Bundle v3** (#104):<br>• Peaks live in `peaks/part-NNNNN.parquet`, not `peaks.parquet`.<br>• Profile data is in `<stem>.profile/profile/part-*`.<br>• An optional `precursor_major/` copy (#109).<br>• Access goes through `massspec.bundle`: `read_bundle_manifest`, `peaks_dataset`, `points_dataset`, `is_profile_bundle`, `related_bundles`.<br>• Older manifests raise `OutdatedBundleError` (no upgrader).<br>`SCAN_METADATA` is now v5, which adds `intensity_floor` (#111). |
| 3 | MS1/MS2 chromatogram extraction + MS2 spectrum retrieval + RT-window pushdown | `massspec.quant.chromatogram`, `massspec.sweep` | `30_extract_intermediates`, experiments | **landed** (XIC) · needed (A3) | `extract_xic_scan_major` (#61) still exists and now accepts a bundle path. **`extract_xic_indexed` was removed** (#109, Mode B retired).<br>New sweep `extract_xic(peaks, targets, *, levels, tolerances, rt_window, layout, claims, fragment_abundance, …)` (#105) returns a `SweepResult`. The `QueryBatch.dense` views (#110) provide `Dense.parent_slots` and `Dense.brackets` for interpolating MS1 at MS2 scans.<br>**A3 `retrieve_spectra` (full per-scan MS2 peak lists) still does not exist.** Re-evaluate whether the sweep's level-2 output covers that need. |
| 4 | `HyperEMGPeak` (+ `WarpedEMGPeak`, `SplinePeak`) | `core.stats.peaks` | experiments (peak fitting) | **landed** (HyperEMG) | #67. Area-parameterised left+right EMG mixture, vmap-safe, fit with `core.optim.DifferentialEvolution` plus LBFGS polish. 4894d07 (2026-11-07) changed "the EMG mode", so June fits need re-checking. `WarpedEMGPeak` and `SplinePeak` are still `needed`. |
| 5 | `GlobalCalibration` (α₀/α₁, d_mz, m/z offset, α_mz/ν_mz) | `massspec.counter.calibration` | PROCAL calibration | **landed**, with gaps | #72–#73, fit by `StagedCalibration`. Gaps:<br>• Defaults α₀ = 1, α₁ = 15 (α(2) = 31) are not evidence-based.<br>• The CLI fits α only with `--gain`, and there is no prior (CR-B3, decision log D1).<br>• Only **one file-level** `mz_offset_ppm` (CR-C1, D4b).<br>• `rho` is a dead parameter (CR-A5).<br>α_mz/ν_mz defaults 0.6/7.8 match nb42b. |
| 6 | Per-progenitor model (f_z, p_k, HyperEMG, ν_I, c_mz) | `massspec.counter.model` | Counter experiments | **landed** (MS1), with gaps | #70, shipped as **`Progenitor`**. p_k is theoretical plus a learned k ≥ 1 correction (D3).<br>**`log_nu_intensity` is a dead parameter:** the intensity term is a Poisson PMF on `I·τ/α`, and `channels.intensity_log_prob` (Student-t) is never called (CR-B1, D5).<br>The CLI fits **one charge per target** (CR-B4, D14). No MS2 channel. |
| 7 | `Panel` + additive-progenitor score + zero-cell term (`PanelSet`, center-weighted gradients) | `massspec.counter.panel` | Counter experiments | **landed** (Panel) · **not built** (PanelSet, center-weighting) | #77. `Panel`, `panel_log_prob`; censoring uses the Poisson void probability with **`floor_ions = 1.0`**, while the real floor is about 200 ions at z = 2 (CR-B2, D8).<br>Multi-progenitor machinery: `discover` (#78), `candidates` (#82), `component` (#83–#85).<br>**Bug:** `component` stacks members by channel *index*, while `_channels_overlap` matches any channel to any channel by m/z (CR-A2).<br>**Never run on PT data.** |
| 8 | Scoring + identification Λ | `massspec.counter` | Counter validation | **landed** (scoring) · **needed** (Λ) | `panel_log_prob`, `attribution` (#79), CLI emit (#86). Λ is deferred in `counter.discover`. |
| 9 | Credible intervals on N_total | `massspec.counter.orchestrate` | uncertainty experiments | **landed**, with gaps | `estimate_n` defaults to VB, but **the CLI calls `inference="map"`, which gives a Laplace interval** (`cli/counter.py:187`). `discover` and `component` report Laplace intervals *conditional* on the other progenitors. No coverage test on real data (CR-D1, D12). |
| 10 | EncyclopeDIA / ScribeTwo search (reads `.raw` natively) + predicted libraries | `massspec.search.encyclopedia` + `thirdparty.encyclopedia`; `massspec.library` | `15_reference_library` (optional re-search) | **landed** | History:<br>• Wrapper ported 2026-05-27 (7328487, EncyclopeDIA 6.5.15).<br>• Koina prediction backend #89 and mass-tolerance units #90 (2026-09-05).<br>• **#100 (2026-10-19)** pinned EncyclopeDIA 6.9.23, added ScribeTwo (DDA spectrum-centric with FDR) and `PSM_TABLE` v2 (`precursor_id`, `qvalue`).<br>Search-space and library prediction: #106–#107. RT alignment: `rt.align` (#102).<br>(The 2026-09-05 refresh attributed the wrapper to #89/#90.) |
| 11 | MaxQuant search-output reader | `massspec.io.maxquant` | experiments (PROCAL and extended identification) | **landed** | `read_maxquant_search` returns `Search.psms`, using `PSM_TABLE`, now **v2** (two added nullable columns). protostar's PSM test fixture was updated for v2 in restart PR #18. |
| 12 | MS2 spectral-similarity suite | `core.stats.losses` + `massspec.spectra.similarity` | exp 01–02 | **landed** | #65. 4894d07 (2026-11-07) made the similarities scale-free, so re-check the exp02 and exp03 numbers. |
| 13 | Consensus / aggregated-spectrum builder | `massspec.spectra.consensus` | exp 01–02 | **landed** | #65. `refine_consensus` is deferred. |
| 14 | Multinomial fragment-intensity **generative** model | `massspec` | MS2 experiments; MS1→MS2 conversion (D13) | needed | Should be over-dispersed (exp03: n_eff in the hundreds). Folded into CR-D3. |
| 15 | Graceful handling of empty or aborted `.raw` | `massspec.readers.thermo` | `10_convert_raw` | **landed** (verify) | beedabb (2026-10-26) checks the instrument count before `SelectInstrument` and raises "holds no mass-spectrometer data", which `convert_batch` reports as `status='error'` instead of a .NET stack trace. **Verify on the known-corrupt Wilhelm2021 file during the reconversion pilot.** |
| 16 | **Orbitrap m/z ↔ frequency reconstruction** | `massspec.readers.thermo` (coefficient capture) + a new transform / panel builder | profile frequency-grid experiment ([`orbitrap_profile_frequency.md`](orbitrap_profile_frequency.md), pending commit) | **prototyping** | The reader captures only `ft_conversion_a/b/c` (`_trailer.py:178-184`). The June prototype found a uniform Δf = 1/512 grid. The payoff would be structural, not better per-peak precision. **Deferred**: profile is a referee, not the PT path (D17). |
| 17 | **`CandidateScreenReducer`**: global MS1 screen as a sweep reducer (matches only; per-query candidate regions with screen features) | `massspec.sweep` / `massspec.search` | Global extraction, Tier 1 ([`extraction_design.md`](extraction_design.md)) | **needed** | Specified from pilot P2 in Constellation's `docs/plans/global-extraction.md`. The prototype is `pipelines/experiments/pilots/_pilot_reducers.py` (pilot only; promote, then delete). |
| 18 | **Tier-2 model-input reducer**: the dense record for N(t) fits (MS1 cells with censored zeros, MS2 with brackets and isolation) | `massspec.sweep` (+ `counter` input) | Global extraction, Tier 2; Counter fits (D18) | **needed** | Today `extract_xic(rule="nearest", drop_unmatched=False)` covers it at 55–130 KB per (sequence, z) (P4). |
| 19 | Sweep gaps found by the pilots:<br>• MS2 tolerance per analyzer (ion trap 0.5 Da vs Orbitrap 20 ppm in one run)<br>• isolation gating by the whole isotope envelope<br>• zero rows under rule `all`<br>• `rt.align` thread count<br>• scan-major does not scale with workers for global sweeps | `massspec.sweep`, `massspec.rt` | Tiers 1–2, RT transfer | **needed** | Measurements are in [`extraction_design.md`](extraction_design.md):<br>• rule `all` + `drop_unmatched=False` writes 75–1,319 of 3,230 MS1 rows per channel;<br>• a 40-point `rt.align` takes 66 s at 104 threads vs 0.04 s at 1;<br>• scan-major takes 690–745 s at 4/8/16 workers vs 30–78 s precursor-major. |

## Build lanes (priority)

From the MS2-first manuscript ([`manuscript_roadmap.md`](manuscript_roadmap.md)), restated for the
restart:

- **Lane 0 — data (blocking everything).** Migrate the protostar pipelines to the bundle-v3 API (#1,
  #2), then reconvert the corpus. A Constellation-side nicety is CR-A4: `convert_batch` should
  reconvert outdated bundles rather than report them `skipped`.
- **Lane A — Part I (MS2).** #11, #12 and #13 have landed; re-check the numbers after 4894d07. A3
  (#3) is still open but may be covered by the sweep.
- **Lane B — Part II (MS1 error structure).** Extraction (#3) has landed. **The model is not ready:**
  the intensity likelihood (CR-B1, D5), the detection floor (CR-B2, D8), the α prior (CR-B3, D1) and
  the m/z offset level (CR-C1, D4b) are open. The empirical figures in
  [`counter_figures_notes.md`](counter_figures_notes.md) are descriptive only until these land.
- **Lane C — Part III (Counter).** Correctness fixes come first: ProForma parsing in the CLI (CR-A1),
  component channel alignment (CR-A2), and v3 inputs (CR-A3). After that come cross-file validation,
  the interference choice (D9), uncertainty (CR-D1), and the MS1→MS2 conversion model (D13, CR-D3).

## Reference architecture

- **Counter engine design:** `~/projects/constellation/docs/plans/cartographer-counter-port.md` and
  `counter-real-data-deconstruction.md`.
- **Bundle v3 rationale:** `~/projects/constellation/docs/plans/profile-and-peak-lists.md`. Search
  stack: `sweep.md`, `search-space.md`, `peak-detection.md`.
- **Feature decisions:** [`counter_model_decisions.md`](counter_model_decisions.md).
- **History:** [`history.md`](history.md).
- **Historical specs** (banners list their stale values):
  [`model_specification.md`](model_specification.md) and
  [`additive_progenitor_likelihood.md`](additive_progenitor_likelihood.md).
- **Original prototyping:** the private archive `wilburnlab/ms_quant_archive`. Its
  `scripts/calibration_model.py` imports Cartographer's `feat/chromatogram-logL-scoring` branch.
