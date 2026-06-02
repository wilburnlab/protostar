# 02 — MS2 fragmentation is multinomial

**Status:** planned — in active development (Part I). _Part I — MS2._

## Scope
The mechanism behind [`01_ms2_spectral_scoring.md`](01_ms2_spectral_scoring.md): an MS2 spectrum
is a **`Multinomial(N, p)`** draw — a finite precursor pool of `N` ions fragments into `K`
competing channels, so the fragment channels are **mutually exclusive** (`Σp = 1`) with negative
covariance `Cov[x_i, x_j] = −N·p_i·p_j` and per-channel variance `Var[x_k] = N·p_k(1−p_k)`.
Characterize, per fragmentation mode, (i) the **stability of the precursor→fragment conversion
ratios** `p̄_k` across replicates of a peptide, and (ii) the **across-replicate variance law** —
that observed `Var[p̂_k]` tracks `p̄_k(1−p̄_k)/N` (and `Cov[p̂_i,p̂_j]` tracks `−p̄_ip̄_j/N`) with a
single fitted per-precursor slope `1/N_eff`. That slope, plotted against total intensity, is the
**bridge to Part II**: the finite pool size *is* the MS1 ion count `N`, and `signal ≈ α·N`
([`03_intensity_and_ion_count.md`](03_intensity_and_ion_count.md)). Robustness: use the **median**
consensus and *report the per-replicate deviance-from-bulk distribution* so the inlier assumption
is shown to be checked, not assumed (the airtight iterative-screening refinement is a deferred
follow-up).

This chapter *characterizes* multinomial shot noise using the existing `Multinomial`; the
predictive fragment-propensity **generative** model (ledger #14) is deferred.

## Source
- Experiment script: `pipelines/experiments/02_ms2_multinomial.py` (TBD)
- Identifications + intermediates: as in chapter 01 (MaxQuant #11; XIC level-2; consensus #13).
- Model: `core.stats.distributions.Multinomial` + `massspec.spectra.similarity.multinomial_deviance`
  (ledger #12, **prototyping**); the retained `per_replicate[R,K]` matrix from
  `massspec.spectra.consensus.build_consensus` (ledger #13, **prototyping**). Any fitted slope uses
  `core.optim` — never scipy.

## Key questions
- How stable are the per-mode conversion ratios `p̄_k` across replicates and across the intensity range?
- Does the across-replicate (co)variance follow the multinomial law with a single slope `1/N_eff`?
- Does fitted `1/N_eff` track total intensity (the `signal ≈ α·N` foreshadowing)?
- Is `2·N·KL` (the multinomial deviance) calibrated to χ²_{K−1} on confident replicates?

_Figures (planned):_ (2.1) per-mode conversion-ratio stability `p̄_k` ± across-replicate error;
(2.2) the mean–variance law `Var[p̂_k]` vs `p̄_k(1−p̄_k)` and `Cov` vs `−p̄_ip̄_j`, single fitted
slope; (2.3) fitted `1/N_eff` vs total intensity (the bridge).

_Findings to be written once exp02 runs (gated on #13 + the existing `Multinomial`)._
