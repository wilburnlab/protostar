# Counter model — decision log

This is the record of **why each feature is in the Counter model, or out of it.** It replaces the
"canonical MS1 ion model — required components" checklist that protostar's `CLAUDE.md` carried from
its scaffold until 2026-11-20. That checklist was a copy of an April guardrail written by an agent in
the MS_Quant notebooks, and some of its values were stale; see [`history.md`](history.md).

**How to use this log**
- Every feature has one entry here.
- Change an entry's status only with new evidence (a notebook, a `results/` record or a PR) and a dated
  line in its **History**.
- Never delete an entry. Superseded entries stay, marked as such.
- Experiments cite entry IDs.

**Statuses**

| Status | Meaning |
|---|---|
| **Adopted** | Required in any model fit that is reported |
| **Under test** | A named experiment decides it; until then, report results under each option |
| **Heuristic** | Allowed for detection, seeding or QC; never a likelihood term |
| **Open** | Not yet modelled; known to matter or suspected to |
| **Rejected** / **Superseded** | Kept for the record |

**Evidence scope matters.** The reference fit is nb42b: **one file**, Zolg2017
`01625b_GA1-TUM_first_pool_1_01_01-DDA-1h-R2`, DDA at 60K, 29 PROCAL peptides. No Counter feature has
yet been validated across files.

**Where the evidence lives**
- `archive:` means the private MS_Quant archive (`wilburnlab/ms_quant_archive`, commit 30d4fc5, under
  `ProteomeTools/`).
- `constellation:` means `constellation/massspec/` at 77f540f.
- `CR-*` IDs are the upstream change requests in Constellation's
  [`docs/plans/counter-pt-reconciliation.md`](https://github.com/wilburnlab/constellation/blob/main/docs/plans/counter-pt-reconciliation.md).
- Tests A–E are protostar's restart experiments: A intensity likelihood (D5), B m/z offset level
  (D4b), C space charge (D10), D binary interference state (D9), E MS1→MS2 conversion (D13).

---

## D1 — Charge-dependent gain α(z) = softplus(α₀ + α₁·z), fit per file

- **Status:** Adopted for the functional form. **Under test** for the values and whether they transfer
  across files.
- **Claim:** intensity per ion scales with charge, roughly linearly. Ion count is
  `N = I·iit/α(z)`, where Thermo's reported intensity is a rate **per second**: **iit in seconds**
  (the trailer's "Ion Injection Time (ms)" divided by 1000).
- **Evidence:**
  - nb22–23 established α as charge-linear.
  - nb42b fit α₀ = 1.06 and α₁ = 10.2, so α(2) = 21.5 (archive:
    `results/calibration/nb42b/42b_global_params_printed.csv`), with iit in seconds
    (`scripts/calibration_model.py`, `iit_s`). That value has the resolution factor ρ_R = 0.707
    divided out.
  - Independent agreement: Hsu et al. 2025 (J. Proteome Res. 24:5742) calibrate `I = α·N/t` in ions/s.
    - Their code (`uw-maccosslab/manuscript-ion-counting`, `Ion_Calibration_Table1_FigS3-5.Rmd`)
      computes `intensity * (injectionTime / 1000)`.
    - They report Orbitrap α = 9.4–15.4 by instrument (Orbitrap Fusion 15.44, Exploris 10.95).
    - Their method, the shot-noise variance of ion ratios, is nb42b's. Their ions are fragments of
      Glu-fibrinopeptide B; the charge is not stated, so read it as roughly z = 1.
  - α₀ and α₁ are **not separately identified** on DDA PROCAL: there are only about 30 z = 1
    observations.
- **Constellation:**
  - `counter/iit.py` uses iit in **milliseconds** ("ms convention"). With α of the evidenced size,
    that makes N 1000× too large and the shot-noise variance (D2) 1000× too small (CR-A6).
  - `counter/calibration.py:44-56`, `GlobalCalibration`, defaults to α₀ = 1 and α₁ = 15, giving
    α(2) = 31. That is seconds-sized, but it is applied in ms, and no evidence supports it.
  - α is fit only when the CLI is given `--gain`.
  - There is no prior; nb41 needed a log-normal prior on α₁ to avoid a degenerate large-α mode.
- **Next:** CR-A6 (seconds), then CR-B3 (a per-file fit with a prior). Fit α per file across the PROCAL
  files and check how stable it is.
- **History:**
  - 2026-11-21: the unit was settled as seconds (nb42b code; Hsu et al. 2025 code), and the
    Constellation mismatch was found (CR-A6).

## D2 — IIT enters the variance, not α

- **Status:** Adopted.
- **Claim:** `Var[I_k] = α(z)·ΣI_pred·p_k(1−p_k)/iit`. An IIT power law on α (β) is not needed.
- **Evidence:**
  - nb42a's β ablation (archive: `results/calibration/nb42a/beta_ablation7.png`) gave β ≈ −0.01 once
    `/iit` was in the variance.
  - Without `/iit`, α inflated 10–100× and the per-peptide ν collapsed to its floor.
  - DDA IIT spans about 0.2–34 ms.
- **Constellation:** the expected count is `λ = N(t)·f_z·p_k·τ`, with τ the scan's IIT, and the
  observation is `n = I·τ/α`. That is consistent in form, but τ is in ms, not seconds, so the variance
  is 1000× too small at evidenced α (D1, CR-A6).

## D3 — Fitted isotope fractions p_k, not theoretical

- **Status:** Adopted.
- **Evidence:** the theoretical envelope is biased with mass: nb15 and nb19 found dp₀ ≈ +0.0075 and
  dp₂ ≈ −0.0064. Left uncorrected, this produces log-ratio offsets of 0.3–0.7 that swamp the
  interference signal.
- **Constellation:** `counter/model.py:116-167` starts from the theoretical log-weights
  (`peptide_envelope`) and adds a learned correction for each progenitor for k ≥ 1. Consistent with
  this entry.

## D4a — Per-file isotope-spacing correction d_mz

- **Status:** Adopted.
- **Evidence:**
  - nb04, nb06 and nb07 (isotope spacing and neutron mass).
  - nb37 `fig1_dmz_by_charge.png`.
  - Without the correction, M+1/M+2 m/z errors are biased by about 0.5 ppm.
- **Constellation:** `GlobalCalibration.d_mz_da`, with M+0 pinned to 0.

## D4b — Level of the m/z offset that is removed before scoring ε

- **Status:** **Under test** (Test B).
- **The record disagrees with itself:**
  - The old checklist described μ̂(t) as the mean residual across **all ions in a scan**, i.e. drift.
  - The model code, archive `scripts/calibration_model.py:116-120`, actually uses the mean of the
    **peptide's own ≤3 isotope errors** in each (peptide, scan, charge), on raw errors with no d_mz
    correction. With one valid isotope, ε = 0.
  - Constellation uses **one `mz_offset_ppm` per file** (`counter/calibration.py:54,85`).
- **Why it matters:**
  - The peptide's own isotope mean absorbs any shift common to the whole envelope. That includes the
    pull an interferer puts on the centroid and any self-space-charge shift.
  - So it removes part of the very signal D9 needs.
  - A scan-level offset estimated from *other* ions does not have this problem. It was never built.
  - Measured between-scan σ was 0.59–0.75 ppm. The "0.47 ppm" in the old checklist is a Student-t
    scale from nb07/nb08.
- **Next:** CR-C1 makes the offset selectable: file-level, per scan from other ions, or per-peptide
  isotope mean. Test B refits α_mz and ν_mz under each option.

## D4c — Per-ion m/z residual ε ~ Student-t(ν_mz, scale √(c_mz/N^α_mz))

- **Status:** Adopted for the form. Values: α_mz = 0.60 and ν_mz = 7.8 (nb42b), **pending a per-file
  refit** under the D4b choice.
- **Evidence:**
  - nb38 chose the power-law precision model by BIC. Its α_mz = 0.929 was fit **without IIT**, and
    ν_mz = 3.06 was a stored per-peptide value that nb40a hard-coded.
  - nb41b, nb41e and nb42b refit both freely. nb42b converged to **α_mz = 0.6037, ν_mz = 7.77**.
    These values exist only in the notebook's printed output.
  - c_mz is per peptide (nb38b/c).
- **Constellation:** defaults α_mz = 0.6 and ν_mz = 7.8 (`counter/calibration.py:52-53`), matching
  nb42b. The old checklist's 0.93 / 3.06 should not be used.

## D5 — Intensity noise is Student-t with per-peptide ν_I

- **Status:** **Adopted.** PI decision, 2026-11-20.
- **Claim:** `I_obs ~ t_{ν_I}(I_pred, √Var)` with `Var` as in D2 and ν_I per peptide.
- **Evidence:**
  - nb20 compared Student-t with Gaussian.
  - nb25 found ν highly peptide-specific (`fig07_nu_predictors.png`).
  - nb27 found ν scales with resolution, η ≈ 0.87.
  - nb42b: median ν_I = 6.4, range 2.3–20.8. In that file ν_I correlates with RT (r = −0.40), tailing
    (−0.31) and interference (0.30), so part of ν_I may be absorbing interference. That should be
    revisited once D9 is settled.
  - protostar: the variance of isotope shares falls with summed intensity at slope −0.48, not the −1
    pure counting noise predicts (stage-30 shard 0, June).
  - exp03: MS2 spectra are over-dispersed relative to counting, with n_eff in the hundreds.
- **Rejected alternative: Poisson-only intensity.**
  - Constellation #72 (753dc88, June) made the intensity term a Poisson PMF on `I·τ/α`. The reason
    was a detection-edge bias seen on **synthetic data only**. It has never been compared with
    Student-t on real data.
  - In the current code (`counter/panel.py`):
    - `channels.intensity_log_prob` (the Student-t, `counter/channels.py:51`) is never called.
    - `Progenitor.log_nu_intensity` (`counter/model.py:160`) is still searched by DE
      (`counter/orchestrate.py:144`) but affects nothing.
- **Next:** CR-B1 makes the likelihood selectable, with Student-t as the default and Poisson and
  NB/Dirichlet-multinomial as comparisons. Test A compares all of them on held-out PROCAL points and
  checks calibration.

## D6 — Sliding-window context features (W = 15, leave-one-out: μ_smooth, μ_sd, kl_smooth)

- **Status:** **Heuristic.** It was listed as "required" in the old checklist, but it is a detection
  and scoring device, not a likelihood term.
- **Source:** `ms1_scoring_test.py` (archive `scripts/`, 03-23); nb36c.

## D7 — Isotope KL with intensity-weighted charge aggregation

- **Status:** **Heuristic.** It was listed as "required" in the old checklist, but it too is for
  detection and seeding.
- **Note:** Constellation's search stack has its own KL scorer (#111), which is used for detection.

## D8 — Zero and below-floor observations are censored

- **Status:** Adopted in principle; **the floor value is wrong in Constellation**.
- **Evidence:**
  - Unobserved cells carry information.
  - nb42c tried a censored term with an arbitrary `cens_weight`, which was not a proper likelihood;
    that notebook is quarantined.
  - The derivation is in [`additive_progenitor_likelihood.md`](additive_progenitor_likelihood.md).
  - The floor is an *intensity*; its ion equivalent scales with iit.
    - In nb42c (archive `notebooks/42c_peak_boundary_refinement.ipynb`, cell 14) the floor, taken as
      the 0.1th percentile of matched MS1 intensities, was 2,959.
    - At the median iit of 1.4 ms and α(2) = 21.3, the notebook printed `N_FLOOR = 0 ions`, which is
      about 0.2 ions with iit in seconds (D1). At a 30 ms iit the same floor is about 4 ions.
    - The "about 200 ions" in this entry before 2026-11-21 (from `docs/nb42c_state_of_play.md`) put iit
      in ms into seconds-convention arithmetic.
    - A sub-ion floor suggests that the lowest matched intensities are noise peaks. The floor should
      come from each scan's own record, not from a percentile.
- **Constellation:**
  - Censoring uses the Poisson void probability with `floor_ions = 1.0` (`counter/panel.py:51`), one
    fixed ion count for every scan.
  - `SCAN_METADATA` v5 (#111) now records each scan's `intensity_floor`.
- **Next:**
  - CR-B2 takes the floor from the scan, in intensity units, converted to ions per scan with iit in
    seconds (CR-A6).
  - Under D5, the censored term must be the Student-t CDF below the floor, so it is consistent with the
    data term (part of CR-B1).
- **History:**
  - 2026-11-21: the "about 200 ions" floor was corrected to about 0.2 ions at the nb42c median iit;
    the ion-equivalent floor varies with iit.

## D9 — Interference

- **Status:** **Under test** (Test D). Two framings compete.
- **(a) Additive progenitors.**
  - Interferers are explicit extra sources, with soft attribution `γ_q = I_q/ΣI`.
  - Built in Constellation: `discover`, `candidates`, `component`, `attribution` (#77–#86).
  - **Never run on PT data.**
  - It has an index-alignment bug: components are stacked by channel *index* even though channels are
    matched by m/z (CR-A2).
- **(b) A per-point binary state, interfered or not.**
  - This is the PI's preferred framing for PT reprocessing (2026-11-20). Centroiding makes it hard to
    assign one signal to several sources, and Orbitrap resolution makes true co-centroiding rare.
  - The earlier attempt (nb41c/e, nb42b) was a two-Gaussian mixture over each ion's *summed log
    likelihood*. It flagged 3.7–6% of points.
  - It was mis-specified rather than refuted. The log-likelihood level differs by about 7 nats between
    peptides and correlates with log N (|ρ| ≈ 0.54), so the classifier sorted peptides and penalised
    peak apexes.
- **Physics both framings must respect:**
  - Interference only **adds** intensity.
  - It pulls the centroid toward the interferer in proportion to the interferer's share.
  - (b) is (a) with one unidentified source averaged out, so a proper (b) is a generative mixture
    whose interfered component has a positive intensity excess coupled to an m/z pull.
- **Gates:** D5 and D10 must be settled first. Without them, bright clean apexes look interfered.
- **Validation design:** PROCAL is in every run, but the background differs by pool. Real interference
  should recur across injections of the *same* pool at the same RT and channel. Add synthetic
  spike-ins, and use profile data as a referee on a subset.

## D10 — Space charge and trap load

- **Status:** **Open**, and believed to matter.
- **Evidence:**
  - In nb37, ε² turns **up** in the highest intensity bins, including the cleanest isotope-fit stratum
    (archive: `results/figures/nb37/fig4_power_law_comparison.png`, `fig5b_power_law_by_kl.png`), and
    ν falls with intensity (`fig6_nu_vs_intensity.png`).
  - At the time this was blamed on interference and absorbed into a floor term
    (`fig7b_power_law_plus_floor.png`).
  - Trap charge (TIC × IIT) was **never tested** as a covariate.
  - On the intensity side, nb25 `fig03_z_vs_tic_aggregate.png` showed no TIC effect.
- **Constellation:** absent; the Counter never reads TIC.
- **Next:** Test C. Plot raw per-ion data first: m/z residual mean and spread, and isotope-ratio bias,
  against the ion's own N, against TIC × IIT, and against Thermo's per-scan compensation fields. Then
  CR-C2 exposes those covariates to the m/z precision term.

## D11 — Elution profile N(t) is a HyperEMG

- **Status:** Adopted.
- **Evidence:** nb39a–m went from EMG to warped EMG to spline to HyperEMG (nb39m, 03-31) and compared
  DE optimizers (nb39k).
- **Constellation:** `core.stats.peaks.HyperEMGPeak` (#67), which is area-parameterised and fit to
  **flux**, `λ = N(t)·τ`.
- **Unverified hypothesis (agent inference):** nb42b may have fit the peak to per-scan counts rather
  than flux. With DDA IIT spanning more than two orders of magnitude (about 0.2–34 ms), that could
  inflate the fitted tails. Check before comparing nb42b peak widths with Constellation's.

## D12 — Uncertainty on N_total

- **Status:** **Open.**
- **History:**
  - MS_Quant planned Laplace intervals (nb42e), but that notebook was never executed.
  - In Constellation, `estimate_n` defaults to VB, but the CLI calls `inference="map"`
    (`cli/counter.py:187`), which gives Laplace.
  - `discover` and `component` report Laplace intervals *conditional* on the other progenitors
    (`counter/discover.py:295`).
  - There has been no coverage test on real data.
- **Next:** CR-D1 adds joint uncertainty and a coverage test, after D5 lands.

## D13 — MS1 → MS2 conversion rates

- **Status:** **Open.** This is the near-term goal.
- **History:** `42b_ms2_incorporation` was never executed, and no conversion rate has been estimated
  anywhere.
- **Identifiability:** only `c_j·α_MS2/α_MS1` is identifiable, unless the MS1 gain carries over to
  MS2. That is plausible for the Orbitrap (same detector) but not for the ion trap.
- **Requires:**
  - v3 bundles, using the sweep with `levels=[1, 2]` and `Dense.brackets` to interpolate N(t) at each
    MS2 scan.
  - Each MS2 scan's own IIT and isolation transmission.
  - Masking the ion-trap low-mass cutoff.
  - An over-dispersed fragment model (see D5's exp03 evidence).
- **Start with** 3xHCD Orbitrap. Constellation work is CR-D3.

## D14 — Fit all charge states of a target jointly

- **Status:** Adopted.
- **Evidence:** charge-dependent gain (D1) is identified only across charges, and the minor charge
  states carry evidence.
- **Constellation:** the CLI calibrates on charges 1–4 (`cli/counter.py:510`) but fits **one charge
  per target** (`cli/counter.py:175,516`). CR-B4.

## D15 — Resolution scaling (ρ_R for variance, η for ν)

- **Status:** **Open.**
- **Evidence:** nb27. ρ and η were never jointly identified (nb27 free fit: ρ = 0.11, η = 1.10). The
  old checklist used ρ = 0.5 and η ≈ 0.87.
- **Constellation:** `GlobalCalibration.rho` exists (`counter/calibration.py:56,100`) but affects
  nothing (CR-A5). This matters once files at different resolutions are compared.

## D16 — Optimisation uses Constellation's differential evolution

- **Status:** Adopted, as a rule.
- **Evidence:** in nb42c, scipy's differential evolution produced optimizer artifacts that were
  mistaken for structural findings. Use `core.optim.DifferentialEvolution`, plus LBFGS polish.

## D17 — Data representation for PT reprocessing: centroid

- **Status:** Adopted for the PT reprocessing. Profile is **deferred**, as a referee on subsets.
- **Evidence:**
  - nb43 (profile vs centroid) was inconclusive, and its implementation had bugs.
  - The 04-17 comparison found centroid intensity ≈ profile *apex height*, not area.
  - The protostar profile → frequency-grid prototype (June) found the Orbitrap frequency grid uniform
    across scans (Δf = 1/512), with only sub-bin calibration drift.
- **PI view (2026-11-20):** profile data would make interference easier to model, but aligning scans in
  RT is still a matching problem. For consolidating PT, use centroids with D9.

---

## Change log

| Date | Change |
|---|---|
| 2026-11-20 | Log created from the restart review. Replaces the CLAUDE.md "canonical 7". D5 adopted by PI decision (Student-t; Poisson-only rejected). D6/D7 reclassified as heuristics. D4b and D4c corrected against the code and the nb42b output. |
