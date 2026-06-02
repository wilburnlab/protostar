# 05 — Resolution and the transient (feature c)

**Status:** planned (stub). _Part II — MS1: two errors, one problem._

## Scope
**Feature (c): resolution = transient length = number of FFT samples.** The FT physics by which a
longer transient simultaneously sharpens the m/z estimate (finer frequency resolution), improves
the amplitude (intensity) estimate, and resolves ¹³C/¹⁵N isotope fine structure — so **fitted**
isotope fractions `p_k` beat theoretical ones. Quantify across 60K (DDA) vs 120K (targeted): the
intensity-variance resolution scaling `ρ_R` (≈0.5, √R theory), the ν resolution scaling
(η≈0.87), and the near-Cauchy `ν_mz` tails from the FT Lorentzian peak shape. This chapter carries
the m/z-error decomposition (per-file `d_mz`, per-scan drift `μ̂(t)`, interference-sensitive
residual `ε_k ~ t_ν_mz(0, √(c_mz/N^α_mz))`) and **closes Part II by deriving Counter** from the
three aligned feature results (03 → 04 → 05).

## Source
- Experiment script: `pipelines/experiments/05_resolution_and_transient.py` (TBD)
- Intermediates: stage 30 MS1 chromatograms + scanmeta
- Model: Student-t m/z residual + resolution scaling (Constellation; ledger items #3, #5).
  Empirical record: `../docs/model_specification.md` §4–5, §8.

## Key questions
- d_mz per mode/dataset; ν_mz near-Cauchy behavior at 120K vs 60K; α_mz precision-vs-N exponent.
- ρ_R and η resolution scaling; do fitted p_k beat theoretical once fine structure resolves?
- Within-scan correlation of m/z errors (per-scan calibration state), distinguished from interference.

_Findings to be written once the resolution/transient experiment runs._
