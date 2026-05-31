# Manuscript roadmap

The experiment plan for the ProToStaR re-analysis, ordered to match the companion manuscript
(the open-science narrative paired with the Counter technical work in
`../../Counter_Manuscript/`). This document is the **plan of record** for *what we run and in
what order*; the curated conclusions land in [`../results/`](../results/), and every modeling
dependency is tracked in [`constellation_contributions.md`](constellation_contributions.md).

> **Status.** Forward-looking. The data layer (fetch / convert / reference-library) is built;
> data migration + `.raw → proc/` conversion are in progress on OSC, and the
> Cartographer→Constellation port is in flight. No experiment here can run until its gate items
> land — this is the roadmap, not a record of results. The `results/` renames and experiment
> stubs it implies are **deferred** until experiments are closer to runnable.

---

## 1. Thesis and arc

The spine of the manuscript is one claim: **m/z error and intensity error are not two problems
but one.** Both are governed by the same physical/instrumental quantities, and once that is
visible the Counter model — joint inference of the ion count **N** through both error channels —
feels inevitable rather than clever.

The reading order is built to deliver that realization to a mass spectrometrist on familiar
ground first:

1. **MS2 (Part I).** Lead where ProteomeTools is best known — the data that catalyzed Prosit.
   Frame it as a *comparison of spectral scoring properties*, and show that **KL divergence**
   scores similar spectra more sensibly than cosine / dot-product measures. That observation
   motivates treating an MS2 spectrum as a **multinomial** draw: an expected precursor→fragment
   conversion, with fragments competing for a **finite precursor pool**.
2. **The bridge.** The size of that finite precursor pool *is* the MS1 ion count **N**. So the
   MS2 story hands off directly to the MS1 question: what governs N and the errors on it?
3. **MS1 (Part II).** The less-explored half. Organized **by feature** — intensity, charge,
   resolution — each shown to move **both** error channels together, in alignment with what
   physics predicts. The three results compound into an "almost *well, duh*" derivation of
   Counter.
4. **Counter (Part III).** Joint additive-progenitor inference of N with credible intervals and
   an identification statistic; interference detection as the headline demonstration.

---

## 2. Data strategy — PROCAL first, for statistical power

The 40 **PROCAL** synthetic calibrant peptides are spiked into *every* run. Across the
~4,213 raw acquisitions that yields an enormous replicate count spanning a wide **intensity
dynamic range** — exactly the leverage needed to fit variance-versus-intensity behavior. An
individual PROCAL peptide's **MS1 signature should be identical contingent on resolution**;
only the **MS2 fragmentation profile varies by acquisition mode**. PROCAL is therefore the
high-N proving ground: establish the patterns there, then **extend to the chemically diverse
but sparse single-pool peptides**.

PROCAL-first is a *statistics* choice, **not** an identification shortcut. We still need a
search engine to know which scans are which peptide. That identification comes from the
**MaxQuant SEARCH zip already paired with every raw file** — a MaxQuant reader
([ledger #11](constellation_contributions.md)) pulls per-scan assignments
(`msms.txt`/`evidence.txt`) with no re-search compute. The published MSP libraries (fetched by
stage 15) provide reference consensus spectra + RT per mode as a guide/filter; an
EncyclopeDIA/Scribe re-search ([#10](constellation_contributions.md)) is an optional
cross-check, not on the critical path.

---

## 3. Target `results/` structure (MS2-first)

Flat `0N_*.md`, file order = reading order. The existing MS1-first chapters are renumbered (via
`git mv`, history preserved) and the two MS2 chapters lead. **Deferred** — captured here as the
target, not yet materialized.

| New | File | Provenance | Scope |
|----|------|-----------|-------|
| 01 | `ms2_spectral_scoring.md` | new | KL vs cosine / dot / spectral-angle / entropy-similarity on PROCAL replicate spectra; why KL is the right comparator (cosine over-penalizes deviations in modestly-abundant ions via its squared penalty) |
| 02 | `ms2_multinomial.md` | new | MS2 as multinomial: per-mode precursor→fragment conversion ratios, finite-pool competition, Var ∝ p(1−p)·N — the bridge to N |
| 03 | `intensity_and_ion_count.md` | old `02_ms1_intensity_model` + intensity-axis of old `03` | **Feature (a):** observed signal ≈ α·N (neither defined yet). Over PROCAL's dynamic range, *both* m/z precision and intensity shot noise move with N; the IIT-in-variance correction |
| 04 | `charge_and_alpha.md` | new (+ α(z) from old `02`) | **Feature (b):** charge scales α = softplus(α₀+α₁z), hence the N behind a given signal; *both* error scales shift with z |
| 05 | `resolution_and_transient.md` | old `03_mz_error_model` + resolution scaling | **Feature (c):** resolution = transient length = number of FFT samples — the FT physics that sharpens *both* errors and resolves isotope fine structure (so fitted p_k beat theoretical); ρ_R, η, ν_mz tails |
| 06 | `chromatogram_shape.md` | old `01` | HyperEMG(1,1) → N(t), N_total; shape priors (σ, τ_R, τ_L, η); τ_R-at-bound interference diagnostic |
| 07 | `procal_calibration.md` | old `04` | per-file GlobalCalibration; cross-mode / cross-dataset consistency; instrument-state and acquisition-order tracking |
| 08 | `counter.md` | old `05` | additive-progenitor joint model: N_total + Laplace credible intervals; identification Λ; interference demos |

**On the restructure.** The old `02_ms1_intensity_model` / `03_mz_error_model` split is
**dissolved**. Per the thesis we do not characterize the two errors separately; Part II is
organized by the three features, and each chapter shows that feature's effect on **both** errors
together. Chapters 03 → 04 → 05 compound into the Counter derivation that opens Part III. The
per-component empirics in [`model_specification.md`](model_specification.md) feed these chapters
(its §8, "Connection between intensity and m/z through N(t)", is the seed of the thesis), but
the chapter *organization* is by feature, not by error channel.

---

## 4. Experiment roadmap

Each chapter is produced by a driver under `pipelines/experiments/` (named below; not yet
written). Drivers call Constellation, emit partitioned parquet + figures; the chapter
synthesizes the conclusions. Gates reference [`constellation_contributions.md`](constellation_contributions.md).

### Part I — MS2 (high-N PROCAL proving ground)

**exp01 → `01_ms2_spectral_scoring`.** Compare similarity metrics (KL divergence, cosine,
normalized dot, spectral angle, spectral-entropy similarity) on PROCAL replicate MS2 spectra,
**stratified by fragmentation mode** (MS2 profiles vary by mode). Demonstrate cosine's
squared-penalty over-sensitivity to deviations in modestly-abundant ions versus KL's
robustness, on within-peptide (replicate) vs between-peptide pairs.
- *Data:* identify PROCAL scans from the paired MaxQuant searches (#11), guided/filtered by the
  MSP reference (`massspec.io.msp`); aggregate replicates → consensus (#13); score every pair
  under each metric (#12). *Extend* to single-pool peptides via the same #11 identifications.
- *Inputs:* `proc/{dataset}/centroid` bundles; `libraries/{mode}/*.msp`; paired MaxQuant searches.
- *Outputs:* per-pair metric table (parquet) + figures (metric distributions; KL-vs-cosine
  scatter; an illustrative modestly-abundant-ion deviation case).
- *Gates:* **#11, #12, #13, #3 (MS2 retrieval).** Not gated on Counter.

**exp02 → `02_ms2_multinomial`.** Characterize fragment intensities as multinomial draws:
per-mode conversion-ratio stability and across-replicate variance versus the p(1−p)·N
prediction; introduce the finite-precursor-pool competition that ties fragments to the MS1 N.
- *Reuse:* the existing `Multinomial` distribution (`core/stats/distributions.py`). The
  predictive fragment-propensity generative model (#14) is **deferred** — this chapter
  *characterizes*, it does not build a predictor.
- *Outputs:* per-fragment variance-vs-p(1−p)N table + figures (conversion-ratio stability;
  multinomial shot-noise overlay).
- *Gates:* **#11, #13** + existing `Multinomial`.

### Part II — MS1: two errors, one problem

Each driver computes **both** the per-ion m/z residual ε and the intensity variance from the
*same* PROCAL observations, as a function of one feature, overlaying the physically-derived
expectation on the data. PROCAL's identical-given-resolution MS1 + huge dynamic range is what
makes these fits clean.

**exp03 → `03_intensity_and_ion_count` (feature a).** Across PROCAL's intensity dynamic range,
show m/z precision σ_mz = √(c_mz / N^α_mz) tightening **and** intensity shot-noise variance
Var[I] ∝ α·ΣI·p(1−p)/iit growing — two consequences of one N. Establish signal ≈ α·N and the
IIT-in-variance correction (the nb42a lesson). Targets: α(z=2) ≈ 19–28, per-peptide ν_I median
~6.5.

**exp04 → `04_charge_and_alpha` (feature b).** Fit α(z) = softplus(α₀+α₁z); show how charge
rescales the N behind a fixed signal and therefore shifts **both** error channels; per-peptide
ν_I behavior across charge.

**exp05 → `05_resolution_and_transient` (feature c).** The transient/FFT exposition: a longer
transient = more samples = finer frequency (m/z) resolution, a better amplitude (intensity)
estimate, and resolved ¹³C/¹⁵N fine structure (so fitted p_k beat theoretical; dp₀≈+0.0075,
dp₂≈−0.0064). Quantify across 60K (DDA) vs 120K (targeted): intensity-variance resolution
scaling ρ_R (ρ≈0.5, √R theory), ν resolution scaling (η≈0.87), and the near-Cauchy ν_mz tails
from the FT Lorentzian peak shape. **Closes Part II by deriving Counter** from the three aligned
results.

- *Gates for exp03–05:* **#3 (MS1 extraction), the `error_model` likelihoods + #5
  (GlobalCalibration).**

### Part III — Counter

**exp06 → `06_chromatogram_shape`.** HyperEMG(1,1) fits → N(t) and N_total; shape priors per
mode; τ_R-at-bound interference diagnostic. *Gates:* **#4, #3.**

**exp07 → `07_procal_calibration`.** Per-file GlobalCalibration (α₀/α₁/d_mz/ρ_R/ν_mz/α_mz);
cross-mode / cross-dataset consistency; acquisition-order and instrument-state tracking (needs
the stage-20 Acquisitions table). *Gates:* **#5**, stage 20.

**exp08 → `08_counter`.** Additive-progenitor joint likelihood over overlapping center-weighted
2D (m/z × RT) panels; N_total + Laplace credible intervals; identification Λ; MS1-only vs
MS1+MS2 agreement; interference detection demos. Math:
[`additive_progenitor_likelihood.md`](additive_progenitor_likelihood.md). *Gates:* **#6, #7,
#8, #9, #3.**

---

## 5. Constellation dependency map (build lanes)

All modeling/scoring/likelihood code lives upstream in Constellation (protostar's one rule). The
gap items group into three lanes by what they unblock; the full ledger is in
[`constellation_contributions.md`](constellation_contributions.md).

- **Lane A — foundational, unblocks Part I (mostly net-new, small).** #11 MaxQuant reader
  (identify PROCAL + all scans), #12 spectral-similarity suite (cosine/dot/entropy-similarity on
  top of the existing `kld`/`spectral_angle`), #13 consensus builder, and the MS2-retrieval
  slice of #3. Minimal path to the first runnable experiment.
- **Lane B — Cartographer port (error/calibration), unblocks Part II.** #3 (MS1 extraction +
  windowed scoring, from `cartographer/data/chromatogram.py`), the `error_model.py` α(z) +
  Student-t likelihoods, #5 GlobalCalibration (`cartographer/calibration.py`).
- **Lane C — Cartographer port (Counter engine), unblocks Part III.** #4 peak shapes /
  HyperEMG, #6 `Mass` (`cartographer/masses.py`), #7 `Panel`/`PanelSet` (new additive-progenitor
  architecture), #8 scoring + Λ (`cartographer/scoring.py`), #9 Laplace credible intervals.
- **Deferred:** #14 multinomial generative model; #10 EncyclopeDIA/Scribe re-search wrapper.

The MS1/Counter lanes are a **port-and-validate** effort, not invention: the predecessor engine
lives at `/home/dbwilburn/projects/cartographer/cartographer/` and the migration is scoped in
Constellation's `docs/plans/cartographer-counter-port.md`. Validation targets are the empirical
values in [`model_specification.md`](model_specification.md) (the nb42b-validated MS1 model).

---

## 6. See also

- [`constellation_contributions.md`](constellation_contributions.md) — the gap ledger (items
  referenced throughout).
- [`model_specification.md`](model_specification.md) — empirical record of the MS1 ion model
  (parameter forms + reference values) that Part II validates.
- [`additive_progenitor_likelihood.md`](additive_progenitor_likelihood.md) — the Counter joint
  NLL (Part III math).
- [`../CLAUDE.md`](../CLAUDE.md) — the "no model code here" rule and the canonical MS1-model
  guardrail.
- [`../results/`](../results/) — where the curated conclusions land.
