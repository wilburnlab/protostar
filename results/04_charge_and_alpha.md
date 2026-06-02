# 04 — Charge and the gain α (feature b)

**Status:** planned (stub). _Part II — MS1: two errors, one problem._

## Scope
**Feature (b): charge rescales the ion count behind a given signal.** Fit the charge-dependent
intensity→ion-count gain `α(z) = softplus(α₀ + α₁z)` and show how charge changes the `N` implied
by a fixed observed signal, and therefore shifts **both** error channels together — the m/z
residual scale and the intensity shot-noise scale move with z through the same `N`. Characterize
per-peptide `ν_I` behavior across charge states. PROCAL's identical-given-resolution MS1 across
many replicates and charges is the clean leverage for this fit.

## Source
- Experiment script: `pipelines/experiments/04_charge_and_alpha.py` (TBD)
- Intermediates: stage 30 MS1 chromatograms + scanmeta (IIT, charge)
- Model: `α(z)` within `GlobalCalibration` (Constellation; ledger items #3, #5). Empirical record:
  `../docs/model_specification.md` §3, §8.

## Key questions
- α₀, α₁ per mode/instrument (target: targeted 120K α₀≈−6, α₁≈16.9; DDA 60K α(z=2)≈19–22).
- How do both error scales shift with z, as the shared-N picture predicts?
- ν_I across charge; reproducibility across replicates.

_Findings to be written once the charge/α experiment runs._
