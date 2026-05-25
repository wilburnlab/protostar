# MS1 Ion Model — Technical Specification

A statistical model for quantifying peptide ion counts from Orbitrap MS1 data,
developed through notebooks nb13–nb42b on the ProteomeTools synthetic peptide dataset.

> **Scope note (protostar).** This document is the **empirical record** of the canonical
> MS1 ion model — the parameter forms and reference values established in the prior analysis.
> Two reframings apply under the current architecture:
>
> 1. **Implementation lives in Constellation.** The model code (α(z), the Student-t intensity
>    and m/z likelihoods, the HyperEMG peak, the Counter modules) is being built in
>    `constellation.massspec`/`core`, not in this repo. See
>    [`constellation_contributions.md`](constellation_contributions.md). protostar experiments
>    import and validate it.
> 2. **Raw-first data construction.** Observations are now built from Thermo `.raw` via
>    Constellation readers (rebuilt mzpeak + scanmeta), not from pre-split mzML / FragPipe
>    inputs. The model is unchanged; only the ingestion path differs. References to
>    `extract_ms1_chromatogram()` below denote the equivalent Constellation extraction.
>
> The §7 Phase-C global `InterferenceMixture` and the binary signal/noise framing are
> **superseded** by the **additive-progenitor** joint likelihood
> ([`additive_progenitor_likelihood.md`](additive_progenitor_likelihood.md), mirrored in
> Constellation's `docs/plans/cartographer-counter-port.md`). nb42c is quarantined.

---

## 1. Overview

The model estimates N(t), the chromatographic ion flux profile (ions/sec) for each
peptide, by jointly fitting intensity and m/z observations across scans. It consists
of three components: (1) an intensity variance model that relates observed intensities
to ion counts via the detector gain parameter alpha, (2) an m/z error model that uses
ion count-dependent precision scaling, and (3) a chromatographic peak shape model
(HyperEMG) that parameterizes N(t) as a function of retention time.

### Physical setting

- **Instrument:** Thermo Orbitrap (FT mass analyzer)
- **Reported intensity:** Image current amplitude from the Fourier transform, proportional
  to ions x charge. Units are per-time (intensity = ion_count x alpha / iit).
- **Ion injection time (IIT):** Time the ion trap accumulates ions before Orbitrap injection.
  Controlled by AGC (automatic gain control) targeting a fixed total ion count.
  Range in DDA: 0.2–34 ms (156x ratio).
- **Resolution:** Determined by transient length. MS1 at 60K (DDA) or 120K (targeted).
  Higher resolution = narrower FT peaks = better intensity precision.

### What the model estimates

The target quantity is **N_total** = the integral of N(t) over the chromatographic peak,
representing the total ion flux for a peptide. This is proportional to the amount of
peptide present, modulo ionization efficiency.

---

## 2. Observations

Each observation is a (scan, charge, isotope) triplet where all three isotopes (M+0,
M+1, M+2) are detected within the extraction tolerance.

**Per observation:**
- I_k: observed intensity of isotope k (k = 0, 1, 2)
- mz_err_k: observed m/z error in ppm (after d_mz isotope spacing correction)
- z: precursor charge state (1-4)
- iit: ion injection time (seconds)
- t: retention time (seconds)

**Data construction:**
- Extracted from mzML via `extract_ms1_chromatogram()` at 20 ppm tolerance
- Only complete isotope series (all 3 isotopes observed in one scan at one charge)
- RT window: ±60s from Bayes3-scored apex

---

## 3. Intensity model

### 3.1 Detector gain alpha

The conversion between intensity and ion count depends on charge state:

```
alpha(z) = softplus(alpha_0 + alpha_1 * z)
```

where softplus ensures positivity. The charge dependence arises because the Orbitrap
image current is proportional to ion_count x charge; more charge per ion produces more
signal per ion.

**Reference values (DDA 60K Orbitrap):** alpha_0 ~ -3 to 0, alpha_1 ~ 10-15,
giving alpha(z=2) ~ 20.

Alpha is fit per file — it depends on resolution, AGC settings, and instrument state.

### 3.2 Ion count from intensity

Intensity is in per-time units. The actual ion count accumulated in a scan is:

```
N_ions = I * iit / alpha(z)
```

This is the number of ions of a given isotope/charge that were present in the Orbitrap
during measurement. Longer IIT = more ions accumulated for the same ion flux.

### 3.3 Intensity prediction

**Round 0 (no peak model):**
```
I_pred_k = p_k * sum(I_obs)    [fraction of observed total]
```

**Round 1+ (peak model available):**
```
I_pred_k = N(t) * p_k * alpha(z)    [from chromatographic model]
```

where N(t) is the HyperEMG-predicted ion flux and p_k are the fitted isotope fractions.

### 3.4 Intensity variance

The variance of observed intensity arises from multinomial shot noise on ion counts,
converted back to intensity units:

```
Var(I_k) = alpha(z)^2 * N_ions * p_k * (1 - p_k) / iit
         = alpha(z) * sum_I * p_k * (1 - p_k) / (iit * rho_R)
```

where:
- The `/ iit` normalizes for accumulation time: short IIT scans have fewer ions
  and therefore more relative noise (this is the key correction discovered in nb42a)
- `rho_R = (R / R_ref)^rho` is a resolution scaling factor (rho ~ 0.5, matching
  sqrt(R) theory from FT transient length). Higher resolution = longer transient =
  narrower FT peaks = better intensity precision per ion.

**Physical derivation:**
1. N ions distributed among isotopes with probabilities p_k (multinomial)
2. Variance of count in bin k: Var(n_k) = N * p_k * (1 - p_k)
3. Intensity = count * alpha / iit, so Var(I_k) = (alpha/iit)^2 * Var(n_k)
4. Substituting N = sum_I * iit / alpha:
   Var(I_k) = (alpha/iit)^2 * (sum_I * iit / alpha) * p_k * (1-p_k)
            = alpha * sum_I * p_k * (1-p_k) / iit

### 3.5 Intensity likelihood

```
I_obs_k ~ Student-t(nu_I, I_pred_k, sqrt(Var(I_k)))
```

The Student-t distribution (rather than Gaussian) accounts for the fact that the
variance estimate itself is uncertain. The "true" noise is Gaussian (Poisson
approximation for large N), but the Orbitrap FT peak shape introduces additional
variation that broadens the tails. The result is a convolution of Gaussian shot noise
with the Lorentzian/Cauchy FT peak profile, which is well-approximated by a Student-t.

**nu_I is per-peptide** (not global). Established in nb25 as highly peptide-specific:
- Median ~ 6-8 for 120K targeted, ~3-7 for 60K DDA
- Range: 2-20+ across PROCAL peptides
- Correlates with RT (r=-0.40), tailing (r=-0.31), interference level (r=0.30)
- Resolution scaling: nu scales as (R/R_ref)^eta, eta ~ 0.87
- Physical interpretation: different peptides produce different FT transient profiles,
  different space-charge environments, different isotope fine structure interference

### 3.6 Isotope fractions p_k

Per-peptide fitted fractions, not theoretical. Parameterized as softmax over logits:

```
logits = [eps_0, eps_1, 0]    [2 free parameters per peptide]
p = softmax(logits)           [ensures sum to 1, all positive]
```

Initialized from `expected_isotope_envelope()` but allowed to deviate. The theoretical
envelope has systematic mass-dependent bias (nb15): dp_0 ~ +0.0075, dp_2 ~ -0.0064.
Using theoretical fractions creates systematic residuals that dominate over interference.

---

## 4. m/z error model

### 4.1 Error decomposition

Three levels of m/z error structure:

**(a) Isotope spacing correction d_mz:**
```
mz_err_corrected_k = mz_err_raw_k - k * d_mz_k / charge
```
Corrects for the difference between pure 13C spacing and the actual centroid spacing
(influenced by 15N and other isotopes). Fit per file from high-quality scans.

**(b) Per-scan systematic offset mu_hat:**
```
mu_hat(t) = mean(mz_err_corrected_k across all isotopes in scan t)
```
Captures instrument calibration drift (~0.47 ppm between-scan sigma). This is NOT
interference — it's a shared shift across all ions in a scan.

**(c) Per-ion residual epsilon:**
```
epsilon_k(t) = mz_err_corrected_k(t) - mu_hat(t)
```
This is the interference-sensitive quantity. After removing the shared drift,
the remaining error reflects per-ion measurement precision.

### 4.2 m/z precision model

The precision of the m/z measurement depends on the number of ions contributing
to the FT peak:

```
sigma_mz_k = sqrt(c_mz / N_k^alpha_mz)
```

where:
- N_k = I_k * iit / alpha(z) is the ion count for isotope k (intensity is per-time)
- c_mz is per-peptide (captures peptide-specific m/z precision)
- alpha_mz ~ 0.93 (global) is the power-law exponent

More ions → narrower FT peak → better m/z precision. The alpha_mz exponent being
slightly less than 1.0 suggests diminishing returns at very high ion counts
(possibly due to space-charge broadening).

**After peak model is available (round 1+):** N_k should use the model-predicted
ion count N(t) * iit * p_k rather than the observed I_k-based estimate, to avoid
circular dependency between m/z scoring and intensity noise.

### 4.3 m/z likelihood

```
epsilon_k ~ Student-t(nu_mz, 0, sigma_mz_k)
```

where nu_mz ~ 3-5 (global). The near-Cauchy tails (especially at 120K) arise from
the Orbitrap FT peak shape — the frequency-domain representation of the transient
has Lorentzian tails that produce heavy-tailed m/z error distributions.

---

## 5. Chromatographic peak model (HyperEMG)

### 5.1 Functional form

The ion flux profile N(t) is modeled as a Hyper-Exponentially Modified Gaussian
with one right-tailing and one left-tailing component:

```
N(t) = N_total * [eta * EMG_right(t; mu, sigma, tau_R)
                  + (1-eta) * EMG_left(t; mu, sigma, tau_L)]
```

where EMG_right is a Gaussian convolved with a right-side exponential decay
(produces chromatographic tailing) and EMG_left is the mirror image (produces
fronting). Each component integrates to 1, so N_total is the exact integrated
ion flux.

### 5.2 Parameters (6 per peptide)

| Parameter | Description | Typical values |
|-----------|-------------|----------------|
| N_total | Integrated ion flux (ions/sec, integrated over time) | varies |
| mu | Center of underlying Gaussian (NOT the peak apex) | RT in seconds |
| sigma | Gaussian core width | 2-5 s |
| tau_R | Right exponential decay constant | 2-10 s (real), >30 s indicates interference |
| tau_L | Left exponential decay constant | 1-5 s |
| eta | Mixing weight (eta=1: pure right-tailing) | 0.7-1.0 |

### 5.3 Why HyperEMG

Chosen over alternatives in nb39 series:
- vs EMG (4p): HyperEMG wins 92% of peptides (median loss ratio 0.55). Independent L/R tails.
- vs WarpedEMG (5p): HyperEMG wins 83%. No monotonicity/fold-back artifacts.
- vs SplinePeak (14p): Similar flexibility with half the parameters. Closed-form integral.
- N_total as direct parameter simplifies uncertainty quantification (Laplace on one param).

### 5.4 Fitting

Per-peptide DE+polish optimization on log-space squared error:
```
loss = sum(w_fit * sqrt(N_obs) * (log(N_pred) - log(N_obs))^2)
```
where N_obs = sum_I * iit / alpha at the dominant charge state.

---

## 6. Parameter summary

### Global parameters (4, fit in phase A)

| # | Parameter | Role |
|---|-----------|------|
| 1 | alpha_0 | Charge intercept for detector gain |
| 2 | alpha_1 | Charge slope for detector gain |
| 3 | alpha_mz | Power-law exponent for m/z precision vs ion count |
| 4 | nu_mz | Student-t df for m/z residuals |

### Per-peptide parameters (4 per peptide, fit in phase A)

| # | Parameter | Role |
|---|-----------|------|
| 5 | eps_fracs (x2) | Isotope fraction logits |
| 6 | c_mz | m/z variance amplitude |
| 7 | nu_I | Intensity Student-t df (per-peptide tail behavior) |

### Per-peptide peak parameters (6 per peptide, fit in phase B)

| # | Parameter | Role |
|---|-----------|------|
| 8 | N_total | Integrated ion flux (**target quantity**) |
| 9 | mu | Peak center RT |
| 10 | sigma | Gaussian core width |
| 11 | tau_R | Right exponential tail |
| 12 | tau_L | Left exponential tail |
| 13 | eta | Right/left mixing weight |

### Fixed constants

| Constant | Value | Source |
|----------|-------|--------|
| rho | 0.5 | Resolution exponent for intensity variance (nb27, sqrt(R) theory) |
| R_ref | 120,000 | Reference resolution |
| d_mz_1, d_mz_2 | ~3.7e-4, ~5.2e-4 Da | Isotope spacing corrections (per file) |

### Removed parameters

| Parameter | Why removed |
|-----------|-------------|
| beta (IIT power-law on alpha) | Unnecessary with /iit in variance (nb42a). Beta absorbs ~0 when variance is correctly normalized. |
| gamma (IIT power-law on nu) | Unnecessary with per-peptide nu_I (nb42a). Per-peptide nu captures all relevant variation. |

---

## 7. Optimization architecture

### EM-style stepwise algorithm

**Phase A — Global + per-peptide non-peak parameters:**
- Differential Evolution (pop=200, gen=500, patience=40)
- Round 0 uses only ±5s apex window (bootstrap from Bayes3 candidate)
- Round 1+ uses full ±60s data with interference weights

**Phase B — Per-peptide peak shapes:**
- Independent HyperEMG DE+polish per peptide (pop=100, gen=1000, patience=30)
- Uses alpha from phase A to convert I → N_obs = sum_I * iit / alpha
- Warm-started from previous round

**Phase C — Interference detection (E-step):**
- Per-observation logL evaluated using fitted model
- Signal/noise classification via mixture model or per-peptide threshold
- Weights applied to next round's phases A and B
- Currently using global InterferenceMixture (to be replaced by per-peptide approach, nb42c)

### Convergence

Typical run: 2 warmup rounds (phases A+B, uniform weights) + 5 EM rounds (A+B+C).
Alpha converges by round 2. Peak shapes stabilize by round 3. Total runtime ~30 min
for 31 PROCAL peptides on a single CPU node.

---

## 8. Connection between intensity and m/z through N(t)

N(t) is the single link between the two likelihood components:

**Intensity side:**
- N(t) predicts how many ions are present at each scan time
- Combined with alpha and p_k, predicts the expected intensity and its variance
- The variance scales with N(t) (more ions → more absolute noise, but less relative noise)

**m/z side:**
- N(t) * iit * p_k gives the actual ion count per isotope per scan
- More ions → narrower FT peak → better m/z precision (sigma_mz scales as N_k^(-alpha_mz/2))
- After peak model is available, N_k from N(t) should replace the I_obs-based estimate

In round 0 (no peak model), both likelihoods use I_obs-based approximations:
- Intensity: I_pred_k = p_k * sum_I_obs (fraction of observed total)
- m/z: N_k = I_k * iit / alpha (ion count from observation)

In round 1+ (peak model available), both should use the model-predicted N(t):
- Intensity: I_pred_k = N(t) * p_k * alpha
- m/z: N_k = N(t) * iit * p_k

---

## 9. Key findings and open questions

### Resolved
- **IIT normalization** (nb42a): Intensity variance MUST include /iit. Without it, alpha
  inflates 10-100x to compensate for the missing variance contribution from short-IIT scans.
- **Beta unnecessary** (nb42a): The IIT power-law on alpha was compensating for the missing
  /iit in the variance. With proper normalization, beta ~ 0.
- **Per-peptide nu essential** (nb42a): Global nu causes alpha to absorb peptide-specific
  variance. Per-peptide nu_I captures real physical variation (peptide-specific FT behavior).
- **Alpha well-determined** (nb42a/b): alpha(z=2) ~ 20 for DDA 60K, consistent across
  all ablation conditions (±m/z, ±beta, ±priors).

### Open
- **Peak boundary detection** (nb42c planned): The global interference mixture is mis-specified.
  Per-peptide thresholds based on seed-window logL distribution outperform global approach
  but need improved stopping criteria to avoid bridging to distant interference.
- **MS2 integration** (nb42d planned): Fragment ion intensities and m/z errors provide
  independent evidence for peak location and shape, using the same N(t) link with
  resolution-scaled alpha and nu.
- **Uncertainty on N_total** (nb42e planned): Laplace approximation at MAP via 6x6 Hessian.
  Now that alpha is stable, this should be well-conditioned.
