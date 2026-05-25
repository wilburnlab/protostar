# 02 — MS1 intensity model (ion-count calibration)

**Status:** planned (stub).

## Scope
The intensity → ion-count relationship: charge-dependent gain `α(z) = softplus(α₀ + α₁z)`,
**IIT in the variance** (`Var ∝ α/iit`), per-peptide intensity degrees of freedom `ν_I`, and
fitted (not theoretical) isotope fractions `p_k`. Re-establish the reference values on the
complete datasets and across resolution settings (60K DDA vs 120K targeted).

## Source
- Experiment script: `pipelines/experiments/` (TBD)
- Intermediates: stage 30 MS1 chromatograms + scanmeta (IIT)
- Model: Student-t intensity likelihood (Constellation; ledger items #5–6). Empirical record:
  `../docs/model_specification.md` §3.

## Key questions
- α(z=2) per mode/instrument; charge scaling α₁; resolution scaling ρ_R.
- ν_I distribution and its correlates (RT, tailing, interference); reproducibility.
- Isotope-fraction bias dp_k vs theory; mass dependence.

_Findings to be written once the calibration experiment runs. Guardrail: the canonical
ion-model checklist in `../CLAUDE.md` (never omit /iit, per-peptide ν_I, fitted p_k)._
