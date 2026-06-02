# 08 — Counter validation (additive-progenitor quantification)

**Status:** planned (stub).

## Scope
Validate the Counter generative quantification model on ProteomeTools: the additive-progenitor
joint likelihood over overlapping center-weighted 2D (m/z × RT) panels, with zero observations
as proper Poisson log-probabilities. Produce N_total estimates with Laplace credible intervals
and the identification likelihood ratio Λ, first on PROCAL (cross-validating the prior nb42b
results) then scaling to pool peptides; add the MS2 fragment channel.

## Source
- Experiment script: `pipelines/experiments/` (TBD)
- Intermediates: stage 30 MS1 (+MS2) chromatograms; stage 04 calibration; stage 15 library
- Model: `Mass` / `Panel` / `PanelSet` / scoring + Λ + Laplace CIs (Constellation; ledger
  items #6–9)
- Math: `../docs/additive_progenitor_likelihood.md`

## Key questions
- Do α, d_mz, ν_I recover nb42b within uncertainty? Do peak shapes converge without τ_R-at-bound?
- Are discovered interferer progenitors (synthesis artifacts) a useful diagnostic?
- MS1-only vs MS1+MS2 N_total agreement; credible-interval calibration.

_Findings to be written once the Counter modules land in Constellation and the experiment runs.
This chapter is the bridge to the Counter manuscript (`../../Counter_Manuscript/`)._
