# 06 — Chromatographic peak shape

**Status:** planned (stub).

## Scope
Characterize MS1 chromatographic peak shape across the three datasets and nine fragmentation
modes: peak width (FWHM), asymmetry / tailing, signal-to-noise, and the HyperEMG(1,1) fit
quality. Establish reference shape priors (σ, τ_R, τ_L, η) and how they vary with mode,
instrument, and acquisition order.

## Source
- Experiment script: `pipelines/experiments/` (TBD)
- Intermediates: stage 30 MS1 extracted chromatograms (PROCAL + library peptides)
- Model: `HyperEMGPeak` (Constellation `core.stats.peaks`; ledger item #4)

## Key questions
- Distribution of HyperEMG parameters per mode; when does τ_R run to a bound (interference)?
- Reproducibility of shape across the 365/222/481 replicates per mode (CV%).

_Findings to be written once the peak-fit experiment runs._
