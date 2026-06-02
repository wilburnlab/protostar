# 03 — Intensity and ion count (feature a)

**Status:** planned (stub). _Part II — MS1: two errors, one problem._

## Scope
**Feature (a): the observed signal is `≈ α·N`.** Across PROCAL's wide intensity dynamic range,
show that *both* error channels move with the single ion count `N`: m/z precision
`σ_mz = √(c_mz / N^α_mz)` tightens while intensity shot-noise variance
`Var[I] ∝ α·ΣI·p(1−p)/iit` grows — two consequences of one `N`. Establish the intensity→ion-count
gain `α(z) = softplus(α₀ + α₁z)`, the **IIT-in-variance** correction (`Var ∝ α/iit`; the nb42a
lesson), per-peptide intensity degrees of freedom `ν_I`, and fitted (not theoretical) isotope
fractions `p_k`. This chapter opens Part II and receives the multinomial bridge from
[`02_ms2_multinomial.md`](02_ms2_multinomial.md): the finite precursor pool *is* this `N`.

## Source
- Experiment script: `pipelines/experiments/03_intensity_and_ion_count.py` (TBD)
- Intermediates: stage 30 MS1 chromatograms + scanmeta (IIT)
- Model: Student-t intensity likelihood + `α(z)` (Constellation; ledger items #3, #5–6).
  Empirical record: `../docs/model_specification.md` §3, §8.

## Key questions
- α(z=2) per mode/instrument (target ~19–28); charge scaling α₁; the IIT-in-variance fix.
- ν_I distribution and its correlates (RT, tailing, interference); reproducibility.
- Isotope-fraction bias dp_k vs theory (dp₀≈+0.0075, dp₂≈−0.0064); mass dependence.
- Does m/z precision and intensity variance co-move with N exactly as the shared-N physics predicts?

_Findings to be written once the intensity/ion-count experiment runs. Guardrail: the canonical
ion-model checklist in `../CLAUDE.md` (never omit /iit, per-peptide ν_I, fitted p_k)._
