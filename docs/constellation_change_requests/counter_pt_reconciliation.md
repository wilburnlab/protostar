# Counter ↔ ProteomeTools reconciliation: Constellation change requests

**Status:** proposed, 2026-11-20. Written by protostar at the project restart, for a dedicated
Constellation session.

**Scope:** `constellation/massspec/counter/` and the parts of `readers` and `cli` it touches, at
`constellation@77f540f`.

**Why:** the Counter port (#69–#86, June 2026) was built against synthetic data and one real-data
smoke. Since then the PT restart review found:
- correctness bugs;
- one modelling departure the PI has now rejected: Poisson-only intensity;
- several places where the notebook evidence and the port disagree.

Each CR names the decision-log entry it feeds ([`../counter_model_decisions.md`](../counter_model_decisions.md),
"D-n"). The dated history is in [`../history.md`](../history.md).

**Evidence locations**
- **archive:** the private `wilburnlab/ms_quant_archive` (commit 30d4fc5), paths under `ProteomeTools/`.
  The reference fit is **nb42b**, on one file, Zolg2017 `01625b_GA1-TUM_first_pool_1_01_01-DDA-1h-R2`
  (60K DDA, 29 PROCAL). Its converged globals are in
  `results/calibration/nb42b/42b_global_params_printed.csv`.
- **The original model code:** `scripts/calibration_model.py`. It imports Cartographer branch
  `feat/chromatogram-logL-scoring`.

**Order and gating**
- Tier A: correctness. Do it first; it is small.
- Tier B: model features with evidence. CR-B1 is the priority.
- Tier C: switches that protostar experiments need. No model decision is baked in.
- Tier D: research. **Do not start**; it is gated on protostar Tests A–D.
- One PR per tier, or per CR where large, each with tests. Accepted CRs are copied into a Constellation
  `docs/plans/` doc before coding, per Constellation convention.

---

## Tier A — correctness

### CR-A1 · The Counter CLI must parse ProForma modified sequences
- **Problem:** `cli/counter.py:175,225,290,516` build `Peptidoform(sequence=modseq)`, treating a
  ProForma string as a bare sequence. Any modified target crashes or is mis-massed. MaxQuant and
  ScribeTwo PSMs carry ProForma `modified_sequence`.
- **Change:** parse with `core.sequence.proforma.parse_proforma`, as other massspec code does.
- **Test:** a CLI round-trip with a carbamidomethyl-C and an oxidised-M target. Theoretical m/z must
  match `massspec.peptide.mz.precursor_mz` for the parsed peptidoform.

### CR-A2 · Component co-fit must align members' channels by m/z, not by index
- **Problem:**
  - `counter/candidates.py:217` `_channels_overlap` puts two targets in one component if **any**
    channel of one is within `collide_ppm` of **any** channel of the other.
  - `counter/component.py:56-61` then checks only that members have the same *number* of channels,
    and stacks them by channel index.
  - So when B's M+0 overlaps A's M+1, or a different charge state, B's predicted signal is added to
    the wrong cells.
- **Change:** build the member → reference channel map from the per-channel m/z centres. Score each
  member's channels at their own m/z on the reference grid; channels outside the grid are scored or
  censored, never silently re-indexed. Alternatively, refuse components whose overlaps are not
  index-aligned until this is supported.
- **Test:** a synthetic two-progenitor component where B's M+0 = A's M+1 (and a cross-charge case).
  Recover both N_total within the CI, and attribute the overlap cell to both members.
- **Feeds:** D9(a).

### CR-A3 · Counter end-to-end on bundle-v3 inputs and the current XIC schema
- **Problem:**
  - The Counter froze in June. Bundle v3 (#104), `XIC_TRACE` v3 and the sweep `extract_xic` (#105)
    came later.
  - `counter estimate` loads `--trace` through `load_xic`, but `counter calibrate`
    (`cli/counter.py:503-504`) uses `pq.read_table` on `--trace` and `--scan-metadata`.
  - No test runs the Counter from a v3 bundle.
  - The columns Counter reads (`counter/observe.py:120-133,281-292`) still exist in `XIC_TRACE` v3,
    so this is about the integration path, not a known schema break.
- **Change:** accept a bundle directory, or a sweep `SweepResult` trace, for both commands. Read
  scan metadata through the bundle API. Use one loader for both commands.
- **Test:** an integration test that runs `convert`, then `extract_xic`, `counter calibrate` and
  `counter estimate` on a tiny synthetic v3 bundle. If a fixture builder exists, use the same one
  protostar's tests use.

### CR-A4 · `convert_batch` should not report an outdated bundle as "already converted"
- **Problem:** `readers/thermo/_batch.py:63-67` `_already_converted` checks only that
  `manifest.json` exists. A pre-v3 bundle in the output parent is reported as `skipped`, and the
  run's downstream users then raise `OutdatedBundleError`.
- **Change:** check the manifest version with `bundle.read_bundle_manifest` or `_check_current`.
  Treat an outdated bundle as pending, or return `status='error'` with an "outdated, pass force=True"
  detail. Do not report it as skipped.
- **Test:** a v1 or v2 manifest stub in the output parent. Without `force`, the result is not
  `skipped`.
- **Note:** protostar's convert driver handles this itself (its restart PR 2), so this CR is
  hygiene, not a blocker.

### CR-A5 · Remove or wire dead parameters
- **Problem:**
  - `GlobalCalibration.rho` (`counter/calibration.py:56,100`) is persisted but has no effect.
  - `Progenitor.log_nu_intensity` (`counter/model.py:160`) is searched by DE
    (`counter/orchestrate.py:144`) but has no effect.
  - Dead parameters inflate the DE search and mislead readers of fit tables.
- **Change:**
  - `log_nu_intensity` becomes live under CR-B1.
  - `rho`: remove it, or wire it as the resolution factor ρ_R (see D15). Removing it is preferred
    until a protostar test identifies it.
- **Test:** every parameter in `GlobalCalibration` and `Progenitor` changes the log-probability. A
  parametrised gradient-nonzero test on synthetic data.

---

## Tier B — model features with evidence

### CR-B1 · Selectable intensity likelihood, defaulting to Student-t with per-progenitor ν_I (priority)
- **Decision:** D5 was **Adopted by the PI on 2026-11-20**. Poisson-only (#72, 753dc88) was
  motivated by a detection-edge bias on synthetic data and was never compared with Student-t on real
  data.
- **Evidence:**
  - nb20 compared Student-t with Gaussian.
  - nb25 found ν highly peptide-specific (archive `results/figures/nb25/fig07_nu_predictors.png`).
  - nb42b: median ν_I 6.4, range 2.3–20.8.
  - protostar shard 0: the variance of isotope shares scales with summed intensity at slope −0.48,
    not −1.
  - exp03: MS2 is over-dispersed by orders of magnitude.
- **Current code:**
  - `counter/panel.py` scores intensity as a Poisson PMF on `n = I·τ/α`, plus a Jacobian.
  - `counter/channels.py:51` `intensity_log_prob` (Student-t) is exported and never called.
- **Change:**
  - Add `intensity_likelihood ∈ {"student_t", "poisson", "negbin", "dirichlet_multinomial"}`
    (names are the implementer's call), with **`"student_t"` the default**.
  - Student-t form, as nb42b and `calibration_model.py` define it:
    `I_obs ~ t_{ν_I}(I_pred, √Var)`, with `I_pred = N(t)·f_z·p_k·α(z)` and
    `Var[I_k] = α(z)·ΣI_pred·p_k(1−p_k)/τ`.
  - ν_I is per progenitor (`log_nu_intensity`, revived), with sensible bounds.
  - In the additive-progenitor panel, `I_pred` and `ΣI_pred` sum over progenitors per cell.
  - **The censored term must match the chosen likelihood.** For Student-t, use
    `log P(I < floor) = log F_t((floor − I_pred)/√Var; ν_I)` instead of the Poisson void probability
    (see CR-B2 for the floor).
  - The gain-Jacobian bookkeeping (#74) must stay correct for each option.
- **Acceptance:**
  1. **Simulation:** recover ν_I, α and N_total under each likelihood, including with the other
     likelihood misspecified, with CI coverage ≥ 0.9 on synthetic panels.
  2. **Reproduce nb42b** on the reference file (as a **v3 bundle**; protostar's reconversion pilot
     produces it), with the 29 PROCAL targets and Student-t intensity:
     - α(2) within about 15% of 21.5 (note that nb42b's α has ρ_R = 0.707 divided out; state the
       convention used);
     - median ν_I within [5, 8];
     - α_mz within [0.5, 0.7];
     - ν_mz within [6, 10].
     Report any systematic difference and why, for example flux vs counts (D11) or the offset level
     (D4b).
  3. **Comparison hook for protostar Test A:** the same fit can be run under each likelihood, giving
     per-point held-out log-likelihoods and PIT values.
- **Feeds:** D5, D8, and Test A.

### CR-B2 · Detection floor from the scan, not 1 ion
- **Problem:** `counter/panel.py:51` sets `_DEFAULT_FLOOR_IONS = 1.0`. The real floor is about 200
  ions at z = 2 (about 4,000 counts in the nb42b file). A floor of 1 ion makes every empty cell a very
  strong "no signal" statement.
- **Change:**
  - Take each scan's floor from `scan_metadata.intensity_floor` (`SCAN_METADATA` v5, #111; the
    lowest stored intensity) and convert it to ions per channel via α(z) and τ.
  - Let it be overridden per run.
  - Record the floor in the fit output.
- **Test:** censoring with the scan floor removes the downward bias in N_total on synthetic panels
  with realistic floors. Repeat on PROCAL across files.
- **Feeds:** D8.

### CR-B3 · α(z) fit per file, with a prior
- **Problem:**
  - `GlobalCalibration` defaults to α₀ = 1, α₁ = 15 (α(2) = 31), which no evidence supports.
  - α is fit only with `--gain`, and there is no prior.
  - nb41 found a degenerate large-α mode without a prior. α₀ and α₁ are not separately identified on
    DDA PROCAL (too few z = 1).
- **Change:**
  - Fit α per file by default in `calibrate`.
  - Add an optional prior on α(2), or on α₁, e.g. log-normal centred near the nb42b value with a
    wide scale. Document it as a weak prior, not data.
  - Report α(z) at each charge with intervals.
  - Change the default initialisation to the nb42b values.
- **Test:** on synthetic data, α is recovered with and without the prior. With the prior off, the
  degenerate mode is detected and reported, not silently returned.
- **Feeds:** D1.

### CR-B4 · Fit all observed charge states of a target
- **Problem:** the CLI calibrates on charges 1–4 (`cli/counter.py:510`) but estimates **one charge
  per target** (`cli/counter.py:175,516`, `[z]`). The minor charge states' evidence is dropped.
- **Change:** targets carry a charge set, defaulting to the charges observed in the trace, and
  `Progenitor.for_peptide` is built with all of them. `f_z` is already modelled.
- **Test:** a 2-charge synthetic target gives a narrower N_total CI than the dominant charge alone, at
  the same truth.
- **Feeds:** D14.

---

## Tier C — switches for protostar experiments

### CR-C1 · Selectable m/z offset level
- **Problem:** the record disagrees with itself (D4b).
  - The old guardrail said "per-scan offset across all ions".
  - The notebook code used the per-(peptide, scan, charge) mean over the peptide's own isotopes.
  - Constellation has one `mz_offset_ppm` per file (`counter/calibration.py:54,85`).
- **Change:** `mz_offset ∈ {"file", "scan_external", "peptide_isotope_mean"}`.
  - `scan_external` is estimated per scan from *other* confident ions in that scan, for example the
    calibration targets, leaving the scored peptide out. It is passed in as a per-scan table.
  - `file` stays the default until protostar Test B decides.
- **Test:** on synthetic data with an injected per-scan drift, `scan_external` removes it, and
  `peptide_isotope_mean` also removes an injected envelope-wide interference pull. That shows why it
  is the wrong choice for detecting interference.
- **Feeds:** D4b, D4c, Test B.

### CR-C2 · Scan covariates available to the m/z precision term
- **Problem:** the high-intensity rise in m/z error (nb37 `fig4_power_law_comparison.png`,
  `fig5b_power_law_by_kl.png`) may be space charge. The Counter never reads TIC or Thermo's
  compensation fields.
- **Change:**
  - Thread per-scan covariates from `scan_metadata` into the observation: TIC, τ, and
    space-charge / lock-mass compensation where present.
  - Add an optional precision extension, e.g.
    `scale² = c_mz/N^α_mz + c_sc·(TIC·τ)^γ` (form open), **off by default**.
  - The point is to make protostar Test C's model comparison possible, not to choose the form.
- **Test:** on synthetic data, an injected covariate-driven variance is recovered when enabled, and
  the extension has no effect when disabled.
- **Feeds:** D10, Test C.

---

## Tier D — research (do not start; gated on protostar Tests A–D)

- **CR-D1 · Joint uncertainty and coverage.**
  - Today: the CLI uses `inference="map"`, giving Laplace (`cli/counter.py:187`). `discover` and
    `component` report Laplace intervals *conditional* on the other progenitors
    (`counter/discover.py:295`).
  - Wanted: joint intervals (VB on the panel, or a joint Laplace) plus a coverage study, synthetic
    first, then PROCAL replicate injections.
  - Feeds D12.
- **CR-D2 · Per-point binary interference state as a generative mixture.**
  - The interfered component has a positive intensity excess, with a centroid pull coupled to that
    excess. Equivalently, it is the additive model with one unidentified progenitor averaged out.
  - Only after protostar settles D5 and D10.
  - Feeds D9(b).
- **CR-D3 · MS1 → MS2 conversion model.**
  - Use the precursor's N(t), interpolated at each MS2 scan through `QueryBatch.dense` /
    `Dense.brackets`, together with the MS2 scan's τ and isolation transmission.
  - Use an over-dispersed fragment model (Dirichlet-multinomial).
  - Only `c_j·α_MS2/α_MS1` is identifiable unless the gain carries over between levels.
  - Mask the ion-trap low-mass cutoff.
  - Feeds D13.

---

## Out of scope for the Constellation session

- **protostar edits of any kind.** protostar owns the experiments (Tests A–E), the decision log and
  the pipelines.
- Changing defaults beyond what each CR says. In particular, the interference framing (D9) stays
  undecided.
- Profile-mode / frequency-grid work (ledger #16, D17), which is deferred.

---

## Agent brief (paste into a fresh Claude Code session in `~/projects/constellation`)

> You are implementing a set of change requests to Constellation's Counter
> (`constellation/massspec/counter/`, plus the `cli/counter.py` and `readers/thermo/_batch.py` pieces
> named below).
>
> The spec is
> `~/projects/protostar/docs/constellation_change_requests/counter_pt_reconciliation.md`. Read it in
> full. Also read:
> - its companion decision log, `~/projects/protostar/docs/counter_model_decisions.md`;
> - Constellation's own `CLAUDE.md` files and `docs/plans/cartographer-counter-port.md` /
>   `counter-real-data-deconstruction.md`.
>
> **Rules**
> - Work **tier by tier: A, then B, then C**. Do not start Tier D.
> - Before coding a tier, write a Constellation `docs/plans/` entry for it (accepted CRs, design,
>   tests), following that repo's conventions.
> - Open one PR per tier, or per CR if large. Each PR needs tests, and must pass the full test suite
>   and the linters.
> - Do not merge your own PRs. The PI reviews each one.
> - Do **not** edit the protostar repo.
> - When a CR's acceptance result would change a decision-log entry, do not edit the log. Put the
>   proposed new status and evidence in the PR description, for the PI to carry into protostar.
>
> **Priorities**
> - **CR-B1** (Student-t intensity, adopted by the PI on 2026-11-20) is the most important change.
> - Its real-data acceptance test needs the nb42b reference file as a **bundle-v3** conversion, which
>   protostar's reconversion pilot will produce. Until it exists, complete the simulation acceptance
>   and leave the real-data check as a documented follow-up. Do not substitute a different file.
>
> **Evidence and conventions**
> - The original notebook evidence is in the private repo `wilburnlab/ms_quant_archive`. Clone it
>   read-only if needed.
> - Cite evidence by path. Do not copy notebook code; re-derive.
> - Use `core.optim.DifferentialEvolution`, never scipy's.
> - Physics-first, generative framing. Every new parameter must change the likelihood (see CR-A5).
>
> **Report back** with:
> - the PR links;
> - for each CR: done, partial or blocked, and why;
> - the CR-B1 simulation coverage numbers;
> - any place where the spec was wrong about the code.
