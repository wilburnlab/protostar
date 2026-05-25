# 04 — PROCAL per-file calibration

**Status:** planned (stub).

## Scope
Per-file calibration of the global parameters (`GlobalCalibration`: α₀, α₁, d_mz_1, d_mz_2,
ρ_R, ν_mz, α_mz) using the 40 PROCAL synthetic calibrants spiked into every run. Establish
cross-mode and cross-dataset consistency of the calibration, and characterize how the global
parameters track instrument state and acquisition order.

## Source
- Experiment script: `pipelines/experiments/` (TBD)
- Intermediates: stage 30 PROCAL chromatograms; stage 20 acquisition time table
- Model: `GlobalCalibration` (Constellation; ledger item #5), stepwise EM
- Empirical record: `../docs/model_specification.md` §6–7.

## Key questions
- α(z) and d_mz stability across 365/222/481 replicates per mode; cross-mode correlation.
- Calibration drift vs acquisition order (carryover / batch effects via stage 20).
- Per-peptide ν_I, c_mz reproducibility across files.

_Findings to be written once per-file calibration runs at scale._
