# 03 — m/z error model

**Status:** planned (stub).

## Scope
The three-level m/z error decomposition: per-file isotope-spacing correction `d_mz`, per-scan
systematic offset `μ̂(t)` (calibration drift, not interference), and the interference-sensitive
per-ion residual `ε_k ~ t_ν_mz(0, √(c_mz / N_k^α_mz))`. Re-establish d_mz, ν_mz, α_mz, and the
per-peptide c_mz on the complete datasets; quantify how m/z precision scales with ion count.

## Source
- Experiment script: `pipelines/experiments/` (TBD)
- Intermediates: stage 30 MS1 chromatograms + scanmeta
- Empirical record: `../docs/model_specification.md` §4.

## Key questions
- d_mz per mode/dataset; ν_mz near-Cauchy behavior at 120K vs 60K.
- Within-scan correlation of m/z errors (per-scan calibration state).
- c_mz drivers; precision vs ion count exponent α_mz.

_Findings to be written once the m/z experiment runs._
