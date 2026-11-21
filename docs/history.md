# Development history and where things live

This is a reconstructed history of the Counter / ProteomeTools work across four codebases. It was
written at the project restart (2026-11-20) because the record was split between two repositories,
an un-versioned notebook tree, and a superseded library. All dates are 2026 unless stated.

- **Git repositories:** dates come from commit history.
- **MS_Quant:** it had no version control, so its dates are file modification times, preserved in the
  archive's `MANIFEST.tsv`.

The *reasons* each model feature is in, or out, are recorded in
[`counter_model_decisions.md`](counter_model_decisions.md). This file records *when and where*.

## Where things live

| Location | What it holds | Role now |
|---|---|---|
| **MS_Quant archive**: GitHub `wilburnlab/ms_quant_archive` (private, archived). Original: lab OneDrive `MS_Quant/ProteomeTools` | nb01–nb43 with outputs; 61 scripts including `scripts/calibration_model.py` (the nb41e/nb42b stepwise-EM model); notes (`CLAUDE.md`, `session_notes.md`, theory notes); result files under 5 MB; Counter manuscript drafts | **Frozen provenance.** Cite it; never extend it. The snapshot keeps original mtimes in `MANIFEST.tsv`. Raw data, mzML/mzpeak caches and result files over 5 MB were not archived. |
| **MS_Quant OSC mirror**: project ESS storage (a `ProteomeTools/` tree beside the protostar root) | A fuller copy of the above, plus the mzpeak caches the notebooks read | Authoritative copy of the large MS_Quant results |
| **Cartographer**: `wilburnlab/cartographer` (last commit 04-21) | The predecessor library the notebooks imported. `HyperEMGPeak` and the DE optimizer used by `calibration_model.py` exist **only** on the pushed, unmerged branch `feat/chromatogram-logL-scoring` | Superseded by Constellation; history only |
| **Constellation**: `wilburnlab/constellation` | All model, reader and search code: `massspec.counter`, `readers.thermo`, `bundle`, sweep/`extract_xic`, the EncyclopeDIA/ScribeTwo search stack | The only home for model code (see `CLAUDE.md`, "The one rule") |
| **protostar**: this repo | Data orchestration, experiment drivers, the curated `results/` record, and the decision record (this file, the decision log, the ledger) | Active |
| **protostar data root**: `<data_root>` on project ESS (set in `config/osc.toml`) | Raw corpus (4,213 `.raw` plus MaxQuant zips); converted bundles; stage-30 shard 0; Counter harness outputs | Data root. Bundles predate bundle v3 and must be reconverted (see 2026-10-31 below) |

## Timeline

### Before protostar: Cartographer and MS_Quant (2024-04 → 2026-04-21)

| Dates | Where | What |
|---|---|---|
| 2024-04 → 2025-09 | Cartographer | Sporadic library work: PyTorch fragment ladder, Koina |
| 02-15 → 03-13 | Cartographer | mzpeak format, chromatogram extraction, isotope envelopes and KL; `calibration.py` and `error_model.py`, a Student-t m/z error model (03-08); scanmeta |
| 03-04, 03-13 | MS_Quant | Theory notes: `fragment_indexing_explained`, `ion_calibration_theory` |
| 03-06 → 03-13 | MS_Quant | nb01–nb12: smoke test, m/z error characterization, isotope spacing, neutron mass, PROCAL isotope-ratio variance and ion counting |
| 03-14 | MS_Quant | **nb13** shot-noise calibration, **nb14** multinomial, **nb15** fitted isotope fractions |
| 03-15 → 03-16 | MS_Quant | nb16–nb24: pooled and cross-mode α, Student-t model comparison (nb20), charge-linear α (nb22–23). **nb25** per-peptide ν, **nb27** resolution scaling |
| 03-17 → 03-25 | MS_Quant | nb28–nb36k: detection and scoring (KL, LLR, the W=15 sliding window, `ms1_scoring_test.py`), scaled MS1 detection |
| 03-28 | MS_Quant | **nb37 / nb37b**: m/z error vs abundance, including the high-intensity upturn. **nb38**: joint m/z–intensity model, α_mz ≈ 0.93 by BIC, fitted without IIT |
| 03-29 → 03-31 | MS_Quant | nb38b/c; nb39a–m peak shape (EMG → warped → spline → **HyperEMG**, 03-31). Cartographer `feat/chromatogram-logL-scoring` (03-30, 04-01) |
| ~03-28 → 04-02 | MS_Quant | The **"Canonical MS1 Ion Model — Required Components"** section is written into `ProteomeTools/CLAUDE.md`, a Claude-maintained guardrail. nb40a cites it as "all 7 components" |
| 04-02 → 04-04 | MS_Quant | **nb40a** per-ion HyperEMG with the canonical model; **nb41a–e** PROCAL calibration and the interference EM. nb41e's stepwise EM is the breakthrough |
| 04-04 → 04-05 | MS_Quant | **nb42a**: the missing `/iit` in the intensity variance had inflated α 10–100×; β (the IIT power law) is unnecessary. Stubs `42b_ms2_incorporation` and `42e_laplace_uncertainty` (never executed) |
| 04-05 → 04-09 | MS_Quant | **nb42b** corrected EM: the reference fit, on **one file** (Zolg2017 `01625b_GA1-TUM_first_pool_1_01_01-DDA-1h-R2`, 29 PROCAL). Converged values: α₀ = 1.06, α₁ = 10.2 (α(2) = 21.5), α_mz = 0.604, ν_mz = 7.77, median ν_I = 6.4. The α_mz and ν_mz values were only printed (`results/calibration/nb42b/42b_global_params_printed.csv` in the archive). `scripts/calibration_model.py` (04-05); Counter manuscript outline and technical reference v1/v2 (04-08/09) |
| 04-11 → 04-16 | MS_Quant | nb42c peak-boundary refinement, **quarantined** (`docs/nb42c_state_of_play.md`, 04-15); `additive_progenitor_likelihood.md` (04-16) |
| 04-17 → 04-18 | MS_Quant | Centroid vs profile comparison (`260417`), then **nb43 profile vs centroid: the last MS_Quant analysis** |
| 04-18 → 04-21 | Cartographer | Thermo `.raw` reader, float64 mzpeak, run → acquisition rename. **Last commit 04-21** |

### Constellation and protostar (2026-04-23 → today)

| Dates | Where | What |
|---|---|---|
| 04-23 → 05-08 | Constellation | Scaffold; `core.chem`, `core.sequence`, `core.stats`; `docs/plans/cartographer-counter-port.md` and the nb42c / profile-vs-centroid research notes (05-08) |
| 05-02 → 05-27 | Constellation | Sequencing, transcriptome, viz, reference cache (#1–#35); MSP reader (#19–#20, 05-19) |
| **05-25** | protostar | **Scaffold** (66fe805). The canonical section, `model_specification.md` and `additive_progenitor_likelihood.md` are copied from MS_Quant, with the evidence range relabelled "nb13–nb42b" without re-checking the values |
| 05-27 → 06-05 | Constellation | Cartographer port, phases 1–8 (#36–#50, including the EncyclopeDIA wrapper, 05-27); Thermo reader and `convert_batch` (#51–#52); scan-metadata activation/analyzer (#56–#57); MaxQuant reader (#58); acquisition order (#60); XIC extraction (#61); MS2 similarity and consensus (#65); `HyperEMGPeak` (#67); thermo moved `io` → `readers` (#68) |
| 05-29 → 06-10 | protostar | #1–#17: CI and tooling; PRIDE fetch, seed/relocate, convert; libstdc++ fix and OSC runbook; MaxQuant; acquisition metadata; peptide/PROCAL reference; exp02 MS2 multinomial; exp30 calibrant MS2 extract; thermo import path (#17, merged 06-10 — **last push before the restart**). The full corpus is fetched and converted on ESS: 4,212/4,213 bundles, centroid and profile, in the pre-v3 bundle format |
| 06-01 → 06-17 | protostar (uncommitted) | Stage-30 MS1 extraction (`ms1_extract.py`; only shard 0 of 16 was run, 06-03); empirical figure scripts 10–14; exp02/exp03 MS2 notebooks; profile → frequency-grid prototype; Counter XIC validation harness (one 10-PROCAL smoke, 06-11) |
| 06-07 → 06-21 | Constellation | **Counter port**, #69–#86: foundations, model, VB, calibration (#72 makes the intensity term Poisson-only, motivated by synthetic data), staged calibration, gain Jacobian, real-data pivot (#75), panel discovery, same-grid, components, CLI, attribution. Design doc `counter-real-data-deconstruction.md` (06-18). **`massspec/counter/` has had only a one-line change since** (2150a53, 11-09) |
| 07-18 → 10-10 | Constellation | Viz, Koina prediction backend (#89), mass-tolerance units (#90), transcriptome / genome assembly / de novo EM |
| 09-05 | protostar (uncommitted) | Contribution-ledger refresh. It overstated the port ("ported and validated"); corrected in this PR |
| 10-11 → 10-23 | Constellation | EncyclopeDIA pinned to 6.9.23, ScribeTwo, `PSM_TABLE` v2 (#100); ms-browser (#101); `rt.align` (#102) |
| 10-26 → 11-02 | Constellation | Thermo reader revision (#103: centroiding, workers, part files). **Bundle v3** (#104, feature commit 10-31): profile written as a sibling `<stem>.profile/` in the same pass; older bundles raise `OutdatedBundleError`, and the design doc assumed "no outside users", so protostar's bundles must be reconverted |
| 11-03 → 11-11 | Constellation | Search stack: sweep `extract_xic` (#105), search space (#106), library prediction (#107), query plumbing (#108), precursor-major layout with Mode B / `extract_xic_indexed` retired (#109), `QueryBatch.dense` (#110), KL scorer and `PeakPicker`, `SCAN_METADATA` v5 with `intensity_floor` (#111) |
| **11-20** | all | **Restart.** MS_Quant archived to GitHub; protostar's canonical-7 guardrail retired in favour of the decision log; Student-t intensity adopted (PI decision); bundle-v3 pipeline migration; Counter change requests handed off as a Constellation plan (`docs/plans/counter-pt-reconciliation.md`) |
| 11-20 → 11-21 | protostar (OSC) | **Bundle-v3 reconversion.** A 5-file pilot (job 17014192: the nb42b file, Zolg 2xIT_2xHCD and ETD, Gessulat 3xHCD, Wilhelm DDA) is **bit-identical** to the June bundles (every m/z and intensity, MS1, MS2 and profile points) despite RawFileReader 1.4.5 → 8.0.42. v3 is 0.82× the June size; the precursor-major copy adds 15% of the June size and 28% of the time and is kept (PI decision). The full pass (jobs 17014460–2, 16 workers, Constellation 6e47177) converts **4,212 / 4,213** in 2 h 04 min to 4 h 47 min per dataset; the corrupt Wilhelm2021 stub fails cleanly. Stage 20 rebuilt; every acquisition row matches June (only the `source_kind` label changed) |

## What this history implies

- **Neither Counter lineage has been validated across files.** nb42b is one file. The Constellation
  port was tested on synthetic data plus one 10-peptide real-data smoke.
- **The protostar "canonical 7" was a copy of an April agent guardrail, not a reviewed specification.**
  Its m/z values pre-date the IIT fix, and its description of μ̂ does not match `calibration_model.py`.
  See decision-log entries D4b and D4c.
- **The protostar data layer was finished against a bundle format that no longer opens.** It was
  reconverted on 2026-11-20/21. A 5-file check found identical data, so June-era analyses should
  stand at the data level; their code had to move to the v3 bundle API.
