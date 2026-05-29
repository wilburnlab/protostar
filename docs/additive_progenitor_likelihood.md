# Additive-Progenitor Joint Likelihood

> **Scope note (protostar).** This is the mathematical foundation for the Counter
> `Panel.forward()` implementation, which lives in **Constellation** (`constellation.massspec`).
> The canonical home for this spec is Constellation's `docs/plans/cartographer-counter-port.md`;
> this copy is kept for in-repo reference and traceability. nb42c is quarantined and superseded
> by the formulation below.

**Purpose.** Extends the Counter v2 joint likelihood (§6) and interference mixture (§8)
to the additive-progenitor setting, and specifies the zero-region censored term on the
same log-probability scale as the data term. This document is the mathematical foundation
for the `Panel.forward()` implementation (architecture in Constellation's
`docs/plans/cartographer-counter-port.md`).

**Scope.** MS1 Orbitrap. MS2 is a direct generalization (different
`p_k`-analogue for fragment ions) and is noted where it enters. Non-Orbitrap analyzers
are deferred.

**Relation to Counter v2.** The per-observation marginals (intensity Student-t,
m/z Student-t, conditional independence) are unchanged. What changes is that a single
observation is now explained by a *sum* of progenitors rather than attributed to a
single peptide. The binary signal/noise mixture of v2 §8 becomes a special case
(one target progenitor + one noise progenitor with fixed background rate).

---

## 1. Notation

**Progenitors.** Indexed by `q ∈ Q_P`, the set of progenitors active in panel `P`.
Each progenitor has per-file state (learned):

- `m_q` — monoisotopic mass
- `f_{z,q}` — charge fractions (simplex over `z`)
- `p_{k,q}` — isotope fractions (simplex over `k`, per charge or shared)
- `N_q(t)` — HyperEMG chromatographic profile; integrates to `N_total,q`
- `ν_{I,q}` — intensity degrees of freedom (per-peptide)
- `c_{mz,q}` — m/z precision constant (per-peptide)
- `m_{z,k,q} = m_q / z + k · Δm_iso` — predicted m/z of the (z, k) channel
  (with per-file `d_{mz,k}` isotope-spacing correction applied)

**Global per-file parameters** (from `GlobalCalibration`): `α₀, α₁, ρ_R, d_{mz,1},
d_{mz,2}, ν_{mz}, α_{mz}`.

**Per-observation indexing.** An observation `i` is a centroided peak at
`(t_i, z_i, m_{obs,i}, I_{obs,i})`. Its IIT is `τ_i = τ_{IIT}(t_i)`.

**Per-cell indexing.** A zero-cell is a panel-local grid cell with no observation,
at `(t_c, z_c, m_{cell,c})`.

**Contribution kernel.** `φ_q(m; z)` is the instrument m/z response kernel for
progenitor `q` at channel `(z, k_q(m))`, centered on `m_{z,k,q}`. In View A (graph,
resolved isotopes) it reduces to an indicator of whether `m_{z,k,q}` falls in the
observation's m/z window. In View B (convolution, partially-resolved), it's the
smoothed-square-wave kernel (plan file §View B).

---

## 2. Additive intensity prediction

### 2.1 Per-progenitor contribution

For progenitor `q` at channel `(z, k)` at time `t`, the predicted intensity contribution
is (Counter v2 Eq 14):

$$
I_q(t; z, k) = \alpha(z) \cdot N_q(t) \cdot f_{z,q} \cdot p_{k,q}
\tag{1}
$$

### 2.2 Total predicted intensity at an observation

For observation `i` at `(t_i, z_i, m_{obs,i})`, let `A_i = \{q : φ_q(m_{obs,i}; z_i) > 0\}`
be the set of progenitors whose (z_i, k_q)-channel falls in the observation's m/z
window. The total predicted intensity is:

$$
\boxed{
I_{pred,i} = \sum_{q \in A_i} I_q(t_i;\, z_i, k_q(m_{obs,i}))
}
\tag{2}
$$

This is the single substantive change from v2 Eq 13: sum over contributors rather than
attribute to one peptide.

### 2.3 Total predicted variance

Ion counts across progenitors are independent Poisson variables (the source rates
`λ_q(t)` are independent per peptide; their sum is Poisson with summed rate). So
variances add:

$$
\mathrm{Var}[I_{obs,i}]
= \sum_{q \in A_i} \mathrm{Var}\!\left[I_q(t_i; z_i, k_q)\right]
$$

Using v2 Eq 16 for each per-progenitor term:

$$
\boxed{
\mathrm{Var}[I_{obs,i}]
= \frac{\alpha(z_i)}{\tau_i \cdot \rho_R}
\sum_{q \in A_i} I_q(t_i; z_i, k_q) \cdot (1 - p_{k_q,q})
}
\tag{3}
$$

The `(1 - p_k)` factor stays per-progenitor because it reflects the within-peptide
multinomial correction, which is a within-peptide property.

**Conditioning choice (flagged).** Eq 3 inherits Counter v2 Eq 6's conditioning on
`N_z` (ion count *after* charge partitioning), which treats `f_{z,q}` as a
deterministic partition rather than a multinomial draw. The more complete form
conditions on `N_{tot}` (before charge partitioning, v2 Eq 5), giving
`Var ∝ f_{z,q} · p_{k,q} · (1 − f_{z,q} · p_{k,q})` per progenitor. The difference
matters only for dominant charges (where `f_z` is large); for minor charges
`(1 − f_z·p_k) ≈ (1 − p_k)`. We keep the simpler form for the MVP and revisit if
residuals at dominant charges show systematic under-dispersion. Physically, the
fuller form is the right one — proton loading during ESI is stochastic per ion —
but the MVP's approximation is consistent with how v2 is currently implemented.

### 2.4 Per-observation intensity likelihood

Each per-progenitor contribution `I_q` is marginally Student-t via the v2 §4.4
scale-mixture argument, with its own `(μ_q, σ_q^2, ν_{I,q})`. The observed total
`I_{obs,i} = \sum_q I_q` is therefore a **sum of independent Student-t random
variables** with different means, scales, and degrees of freedom. This sum is *not*
itself Student-t in closed form.

**Welch–Satterthwaite approximation.** For `X_q \sim t_{\nu_q}(\mu_q, \sigma_q^2)`
independent, the sum is well approximated by a Student-t with matched first two
moments and an effective degrees-of-freedom derived by matching the variance of
the variance (the Satterthwaite formula):

$$
\sum_q X_q \;\approx\; t_{\nu_{\text{eff}}}\!\left(\sum_q \mu_q,\ \sum_q \sigma_q^2\right),\qquad
\frac{1}{\nu_{\text{eff}}} = \sum_{q \in A_i} \frac{w_q^2}{\nu_{I,q}},\quad
w_q = \frac{\sigma_q^2}{\sum_{q' \in A_i} \sigma_{q'}^2}
\tag{4}
$$

with `σ_q^2 = (α(z_i)/(τ_i ρ_R)) · I_q · (1 − p_{k_q,q})` (the per-progenitor
contribution to Eq 3).

Applied to the total observation:

$$
I_{obs,i} \sim t_{\nu_{\text{eff},i}}\!\Big(I_{pred,i},\ \sqrt{\mathrm{Var}[I_{obs,i}]}\Big)
\tag{4'}
$$

**Limits.** When one progenitor dominates (`σ_{q^\star}^2 \gg \sum_{q \ne q^\star}
\sigma_q^2`): `w_{q^\star} \to 1`, `w_{\text{rest}} \to 0`, so
`ν_{\text{eff}} \to ν_{I, q^\star}` — the dominant-progenitor rule emerges
automatically. When progenitors have equal variances: `w_q = 1/|A_i|`, so
`ν_{\text{eff}} = |A_i| \cdot \bar{\nu}_{\text{harm}}` where `\bar{\nu}_{\text{harm}}`
is the harmonic mean of the `ν_q`'s — more independent contributors push the
aggregate closer to Gaussian, which is physically sensible.

This formula replaces the hand-waved "dominant vs harmonic mean" choice from the
initial draft. It is a genuine approximation (the sum of Student-t's is not
Student-t in general), but it is the standard, principled moment-matching choice
and reduces cleanly to the single-progenitor case.

---

## 3. Additive m/z prediction

**Honest framing.** Unlike intensity, where Poisson independence gives us clean
additivity for both mean (Eq 2) and variance (Eq 3), the m/z observable is the
output of vendor centroiding — a black box that typically fits a peak shape and
reports the fitted apex. We do not have a physical model of how that algorithm
responds to an asymmetric composite peak. The "intensity-weighted mean" ansatz
below (Eq 5) is a *first-order approximation* to what centroiding does in the
limit of two well-matched Gaussian contributors; it is not a derivation from
FT physics. For contaminated observations, §3.6 describes a masking strategy
that avoids placing weight on this uncertain prediction.

### 3.1 Expected centroid (first-order ansatz)

As a first-order model, a centroided FT peak containing contributions from
multiple progenitors is taken to have expected centroid equal to the
intensity-weighted mean of the contributors' predicted m/z:

$$
\boxed{
m_{pred,i} = \frac{\sum_{q \in A_i} I_q(t_i; z_i, k_q) \cdot m_{z_i, k_q, q}}
{\sum_{q \in A_i} I_q(t_i; z_i, k_q)}
}
\tag{5}
$$

When `|A_i| = 1`, this reduces to `m_{pred,i} = m_{z,k,q}`, the v2 single-progenitor
case.

### 3.2 Per-scan systematic offset

The per-scan drift `\hat{\mu}(t)` (v2 §5.2b) is unchanged. It is computed from the
pooled set of per-ion residuals across all progenitors at scan `t`, after Eq 5 is
applied:

$$
\hat{\mu}(t) = \text{mean}_{i \in t} \big[m_{obs,i} - m_{pred,i} - d_{mz,k_q(i)}/z_i\big]
\tag{6}
$$

### 3.3 Per-observation m/z residual

$$
\varepsilon_i = m_{obs,i} - m_{pred,i} - d_{mz,k_q(i)}/z_i - \hat{\mu}(t_i)
\tag{7}
$$

### 3.4 m/z residual variance

Variance of the centroid of a mixture of FT peaks. Each per-progenitor centroid has
variance `c_{mz,q} · N_q^{-α_{mz}}` (v2 Eq 25). For the intensity-weighted mean of
independent centroid estimates with weights `w_q = I_q / \sum_q I_q`:

$$
\boxed{
\sigma^2_{\varepsilon,i}
= \sum_{q \in A_i} w_q^2 \cdot c_{mz,q} \cdot N_{z_i, k_q, q}(t_i)^{-\alpha_{mz}}
}
\tag{8}
$$

where `N_{z,k,q}(t) = N_q(t) · f_{z,q} · p_{k,q}` is the per-channel ion count.
Single-progenitor limit: `w_q → 1`, recovering v2 Eq 25.

### 3.5 Per-observation m/z likelihood

$$
\varepsilon_i \sim t_{\nu_{mz}}(0,\ \sigma_{\varepsilon,i})
\tag{9}
$$

with global `ν_{mz} ≈ 3.06` (v2 §5.4).

### 3.6 Contamination masking

Because the centroid prediction (Eq 5) is not backed by a real centroiding model,
we restrict the m/z likelihood to observations where a single progenitor
*clearly dominates* — cases in which Eq 5 reduces to the well-understood
single-progenitor prediction `m_{pred,i} = m_{z,k,q^\star}`.

Define the dominant-progenitor attribution (from Eq 11 below):

$$
\gamma_{i,q^\star} = \max_q \gamma_{i,q},\qquad q^\star = \arg\max_q \gamma_{i,q}
$$

**m/z masking rule:**

- If `γ_{i,q^\star} \ge γ_{\text{thresh}}` (e.g. 0.9): the observation is
  "cleanly attributed" to `q^\star`. Use the single-progenitor prediction
  `m_{pred,i} = m_{z_i, k_{q^\star}, q^\star}` and the single-progenitor variance
  `σ_{\varepsilon,i}^2 = c_{mz, q^\star} · N_{z_i, k_{q^\star}, q^\star}^{-α_{mz}}`.
  The m/z term contributes to `\ell_i`.
- If `γ_{i,q^\star} < γ_{\text{thresh}}`: the observation is contaminated.
  **Drop the m/z term** from `\ell_i`. Only the intensity Student-t contributes.

This is a pragmatic hybrid that blurs Views A and B: the attribution (`γ`) is
computed in View A style (single-progenitor m/z likelihoods), but the decision
to include or drop the m/z term implicitly acknowledges that View B-style
convolutional centroiding doesn't have a validated model yet.

**Notes.** The intensity contribution always remains (Eqs 2–4). The masking is
m/z-only. `γ_thresh` is a tuning parameter; 0.9 is a starting point, but the
right value depends on how sharp `γ` distributions are empirically — candidate
for sensitivity analysis on PROCAL. As the m/z kernel model matures (View B),
this masking can be replaced with a proper likelihood contribution under the
mixture centroid.

---

## 4. Joint per-observation log-likelihood

Conditional independence of intensity and m/z given ion counts (v2 Eq 26) survives the
additive-progenitor generalization because it holds per-progenitor and independent
progenitor contributions factor. With the m/z mask from §3.6:

$$
\boxed{
\ell_i = \log t_{\nu_{\text{eff},i}}(I_{obs,i};\ I_{pred,i},\ \sigma_{I,i})
\;+\; \mathbb{1}[\gamma_{i,q^\star} \ge \gamma_{\text{thresh}}] \cdot
\log t_{\nu_{mz}}(\varepsilon_i;\ 0,\ \sigma_{\varepsilon,i})
}
\tag{10}
$$

This extends Counter v2 Eq 27: the intensity term always contributes (additive
model is clean there), and the m/z term contributes only when the observation is
cleanly attributed to a dominant progenitor. The indicator function is deterministic
given `γ_{i,q^\star}` (which is a post-hoc computation from Eq 11), so the
likelihood remains a proper joint log-density for the observations in each category.

### 4.1 Posterior progenitor attribution (by-product)

If needed for diagnostics (or for the identification score below), the posterior
probability that observation `i` came from progenitor `q` is:

$$
\gamma_{i,q} = \frac{I_q(t_i; z_i, k_q) \cdot \mathcal{L}_{mz,q}(m_{obs,i})}
{\sum_{q' \in A_i} I_{q'}(t_i; z_i, k_{q'}) \cdot \mathcal{L}_{mz,q'}(m_{obs,i})}
\tag{11}
$$

where `\mathcal{L}_{mz,q}(m) = t_{\nu_{mz}}(m - m_{z,k,q}; 0, \sqrt{c_{mz,q} · N^{-α_{mz}}})`
is the single-progenitor m/z likelihood. Note: `γ_{i,q}` is a post-hoc diagnostic,
not a parameter of the loss.

---

## 5. Zero-observation likelihood

**Terminology note.** Earlier drafts called this the "censored likelihood"
following nb42c's usage. In classical statistics, *censoring* refers to the case
where a measurement is known to fall below (or above, or within) a threshold
without its exact value being observed — the likelihood contribution is the
integral of the density over the unobserved region (a CDF evaluated at the
threshold). That terminology applies correctly to the optional detection-floor
form (§5.3 below), but is a stretch for the default form (§5.2), which is an
*exact point-probability* under the Poisson model — no threshold, no integration
over a latent true value. We therefore use "zero-observation" for the overall
section, reserving "censored" for the detection-floor subcase where it applies
in its proper sense.

### 5.1 Panel discretization

A panel tiles `(t, z, m)` into cells. For each (scan `t_c`, charge `z_c`, m/z cell
`m_{cell,c}`) in the panel:

- **Observed cells:** contribute Eq 10.
- **Zero cells:** contribute the log-probability that no ion was reported given the
  predicted rate.

"Cell" here is a resolution-scaled m/z bin at a given scan; its width is the
instrument-kernel width from §1. Only cells where at least one progenitor has
non-negligible predicted contribution (`I_{pred,c} > \epsilon`) enter the zero
term — otherwise the prediction trivially matches the non-observation and the
contribution is identically zero.

### 5.2 Poisson zero-event probability (default)

Ion counts are Poisson (v2 Eq 2). The probability that the underlying count
`N_{pred,c}^{total} = \sum_q N_{q,z_c,k_q}(t_c)` takes the value zero is an exact
discrete point-probability:

$$
\boxed{
\log P(N = 0 \mid N_{pred,c}^{total}) = -N_{pred,c}^{total}
}
\tag{12}
$$

This is *not* a censored likelihood in the classical sense — it is the exact
probability of the discrete event `N = 0` under the Poisson model. In point-process
language, this is the *void probability* of the predicted intensity measure over
the cell. It has the right units (log-probability in nats), it sums directly with
Eq 10 on the same scale — no `cens_weight` fudge factor — and it has no
divide-by-zero pathology because it is a log-exp rather than a standardized CDF.

### 5.3 Left-censored likelihood at the detection floor (optional)

Eq 12 identifies "no observation" with "zero ions". Real Orbitrap data has an
intensity floor `I_{floor}` below which centroided peaks aren't reported: a cell
with a small but nonzero ion count still produces no observation if the resulting
intensity fails to clear the detection threshold. In this case the observation is
genuinely *left-censored* — we know the true intensity was below `I_{floor}` but
not by how much — and the proper likelihood contribution is the integral of the
data density over `(−∞, I_{floor})`:

$$
\log P(I_{obs,c} < I_{floor} \mid I_{pred,c})
= \log F_{t_{\nu_{\text{eff}}}}\!\left(\frac{I_{floor} - I_{pred,c}}{\sigma_{I,c}}\right)
\tag{13}
$$

Eq 13 is a proper left-censored Student-t log-likelihood. At low `N_{pred}` it
approaches Eq 12; at moderate `N_{pred}` it additionally captures cells where a
few ions are present but the intensity fails to clear the floor.

**MVP recommendation.** Use Eq 12 (Poisson zero-event). It is the exact zero-ion
probability and dominates the information content of zero cells for our data —
almost all zero cells correspond to `N_{pred}` values where the Poisson void
probability is the binding constraint, not the detection threshold. The
left-censored form (Eq 13) costs a CDF evaluation per zero-cell and should be
adopted only if Eq 12 underweights zeros empirically on PROCAL.

### 5.4 Summing over zero cells

The panel-local zero-observation contribution to the log-likelihood is:

$$
\mathcal{L}_{\text{zero}, P}
= \sum_{c \in \text{zero cells of } P} -N_{pred,c}^{total}
= -\sum_{c} \sum_{q} N_{q}(t_c; z_c, k_q(m_{cell,c}))
\tag{14}
$$

**Interpretation.** This is equivalent to subtracting the integrated predicted ion
count over the panel's unobserved volume from the log-likelihood. The more predicted
signal the model places in a cell with no observation, the larger the penalty — and
the penalty is linear in predicted intensity, which is a clean, well-behaved gradient.

### 5.5 Relation to nb42c's "censored" term

nb42c used `log F_{t_ν}((I_{floor} − I_{pred})/\sqrt{N_{pred}})` *per peptide per scan*
with `cens_weight` as a scalar fudge factor. The term was properly called censored
(it was the left-censored Student-t CDF of Eq 13), but three structural issues made
it unworkable:

1. **Grid was peptide-local, not panel-local.** nb42c had to decide what region
   belonged to each peptide — the problem §5.1 dissolves by making the grid panel-
   local.
2. **`cens_weight` was a fudge.** The "weighted squared log-residual" data term and
   the censored log-CDF term were on incomparable scales; `cens_weight` papered
   over this. In Eq 12 + Eq 10, both terms are proper log-probabilities and sum
   directly.
3. **Divide-by-zero.** nb42c's denominator `\sqrt{N_{pred}}` collapsed at small
   predictions. Eq 12 has no denominator.

The optional Eq 13 here inherits nb42c's censored form but on a proper likelihood
scale and with the panel-local grid.

---

## 6. Panel-level loss

### 6.1 Full panel log-likelihood

Panel `P` with observation set `O_P`, zero-cell set `Z_P`, and progenitor set `Q_P`:

$$
\boxed{
\mathcal{L}_P = \sum_{i \in O_P} \ell_i + \sum_{c \in Z_P} \log P(\text{no obs} \mid c)
}
\tag{15}
$$

where `\ell_i` is Eq 10 and the zero-term is Eq 12 or Eq 13.

### 6.2 Run-level loss

Panels compose additively across the run (global parameters are shared; progenitor
parameters are owned by one home panel each; see §7):

$$
\mathcal{L}_{\text{run}} = \sum_P \mathcal{L}_P
+ \log \pi(\theta_{\text{global}}) + \sum_q \log \pi(\theta_q)
\tag{16}
$$

with informative priors where Counter v2 specifies them (log-normal on `α_1` at 20,
simplex priors on `f_z` and `p_k`, etc.).

---

## 7. Center-weighted gradients for overlapping panels

### 7.1 The mechanics

Panels overlap in `(m/z, t)`. Each progenitor `q` has a single **home panel** `P^{\star}(q)`
— the one whose center `(m_P, t_P)` is nearest to `m_q` and the apex of `N_q(t)`. In
the home panel, `q`'s parameters receive full gradient. In overlapping non-home panels,
`q` still contributes to the forward prediction (Eqs 2, 5, and 14 all sum over
*every* progenitor whose kernel overlaps the panel), but its gradient from those
panels is either zero (hard gate) or attenuated (soft).

### 7.2 Hard gate (MVP default)

For panel `P` and progenitor `q`:

$$
\frac{\partial \mathcal{L}_P}{\partial \theta_q}^{\text{effective}}
= \begin{cases}
\partial \mathcal{L}_P / \partial \theta_q & \text{if } P = P^{\star}(q) \\
0 & \text{otherwise}
\end{cases}
\tag{17}
$$

Implementation: in non-home panels, wrap `q`'s predicted contribution in `.detach()`
before it enters the forward computation, so it contributes to `I_{pred}` and
`m_{pred}` numerically but carries no gradient.

### 7.3 Soft weight (optional refinement)

$$
w_{\text{center}}(q, P) = \exp\!\left(
- \frac{(m_q - m_P)^2}{2 h_m^2}
- \frac{(\mu_q - t_P)^2}{2 h_t^2}
\right)
\tag{18}
$$

with `h_m`, `h_t` half the panel's m/z and RT half-widths. Use `w_{\text{center}}` as
a multiplier on the gradient (or equivalently, as a multiplier on the panel's
contribution to q's loss). Hard gate is recovered in the limit `w \in \{0,1\}` with
Voronoi boundaries between panel centers.

### 7.4 Why this doesn't double-count

Each `q` has exactly one home panel. Its parameters are updated only from that one
panel's loss. Non-home panels contribute correct forward predictions (the additive
model requires all progenitors overlapping the panel to be summed, for correctness of
`I_{pred}` at observations in the overlap region) but do not influence `q`'s learned
parameters. Observations in overlap regions contribute to whichever home panels
contain them (typically one — panels can be set to have centers Voronoi-far apart
while still overlapping at edges).

---

## 8. Reduction to Counter v2 special cases

### 8.1 Single-progenitor, no interferers

`|Q_P| = 1`, `|A_i| = 1` for all `i`:

- Eq 2 → Eq 13 (v2): `I_{pred,i} = α(z_i) · N_q(t_i) · f_{z_i,q} · p_{k_q,q}`.
- Eq 3 → Eq 16 (v2) after substitution.
- Eq 5 → `m_{pred,i} = m_{z,k,q}`.
- Eq 8 → Eq 25 (v2) variance.
- Eq 11 → `γ_{i,q} = 1` (trivial).

### 8.2 Counter v2 binary mixture

Two progenitors: target `q_t` with HyperEMG, and "noise" `q_n` with flat
`N_{q_n}(t) = \lambda_n` and no peak structure. Eq 11 then recovers the v2 §8
two-component posterior `γ_i`. We consider this reduction a consistency check,
not a preferred model — the additive form with physically meaningful noise
progenitors (synthesis artifacts, modifications with their own learned peak shapes)
is strictly more informative.

### 8.3 Identification log-likelihood ratio

Counter v2 §10.1 `Λ` for peptide `q` is recovered as:

$$
\Lambda_q = \mathcal{L}_P(Q_P) - \mathcal{L}_P(Q_P \setminus \{q\})
\tag{19}
$$

i.e., the change in panel log-likelihood when `q` is removed from the progenitor set.
This is the marginal value of progenitor `q` in explaining the panel's data — the
natural identification score that also accounts for interference from the other
progenitors in the panel.

---

## 9. Open items and flagged approximations

1. **Welch–Satterthwaite ν_eff as approximation** (§2.4). The sum of independent
   Student-t's is not Student-t; Eq 4 is a moment-matched approximation. Adequate
   in the regimes we care about (one dominant progenitor, or many contributors
   pulling toward Gaussian), but should be sanity-checked against Monte Carlo on
   representative panels.
1b. **Charge-partition conditioning** (§2.3). MVP uses `(1 − p_k)` (conditioning on
    `N_z`, v2 Eq 6). Full multinomial on `N_tot` would give `(1 − f_z · p_k)` and
    is physically more correct (proton loading is stochastic per ion). Revisit if
    dominant-charge residuals show systematic under-dispersion.
2. **m/z kernel `φ_q` for View B** (§1). MVP uses View A (indicator-style membership);
   the smoothed-square-wave form for low-resolution / fine-structure is not yet
   parameterized. Plan: specify shoulder shape as an instrument-level prior per
   analyzer family once IT/TOF data are characterized.
3. **m/z masking threshold `γ_thresh`** (§3.6). Pragmatic hybrid that admits we
   don't have a centroiding model for composite peaks. Starting value 0.9;
   candidate for sensitivity analysis. Retires once a validated View B kernel
   model is in hand.
3b. **Intensity floor `I_{floor}` as metadata** (§5.3). Currently a per-notebook
    empirical value (`≈ 4000 counts` for DDA 60K Zolg2017). Proper home: per-file
    `proc/` bundle metadata, characterized per-file or per-mode. Affects whether the
    left-censored form (Eq 13) is needed on top of the Poisson zero-event (Eq 12).
4. **Per-cell "meaningful prediction" threshold** (§5.1). Which cells enter the zero
   term? Proposal: cells where `I_{pred,c}^{total} > ε · I_{floor}` for some small
   `ε` (e.g. 0.01), so cells with negligible predicted signal contribute negligibly.
   This needs an empirical pass on PROCAL to avoid under- or over-inclusion.
5. **`A_i` construction (which progenitors belong to observation `i`)** (§2.2). In
   View A, `A_i = \{q : |m_{obs,i} - m_{z,k,q}| < \Delta_{\text{res}}\}` where
   `\Delta_{\text{res}} = m / R` is the resolution window. This is the m/z-side
   analogue of the panel tiling. Note that membership is a function of predicted m/z
   (not observed), so the dependency is one-way and differentiable.
6. **Per-panel vs per-progenitor ν_I fit** (§2.4). Currently ν_{I,q} is a Mass-level
   parameter. When the same progenitor has a home panel and contributes to neighbors,
   its ν_I is only updated from its home panel's loss (consistent with §7).
7. **MS2 extension.** Fragment ions enter as extra `(z, k)`-like channels on the
   `Mass` module with their own `α_{MS2}(z)` (resolution-scaled), their own
   `p_{k,MS2,q}` (fragment intensity fractions, analogous to isotope fractions for
   MS1), and share the same `N_q(t)` peak profile. Eqs 2, 3, 5, 8 carry through
   directly with expanded channel indices.

---

## 10. Summary of key equations (quick reference)

**Additive intensity prediction:**
$$
I_{pred,i} = \sum_{q \in A_i} \alpha(z_i) \cdot N_q(t_i) \cdot f_{z_i,q} \cdot p_{k_q,q}
$$

**Additive intensity variance:**
$$
\sigma_{I,i}^2 = \frac{\alpha(z_i)}{\tau_i \rho_R}
\sum_{q \in A_i} I_q(t_i; z_i, k_q)(1 - p_{k_q,q})
$$

**Additive m/z centroid prediction:**
$$
m_{pred,i} = \frac{\sum_q I_q \cdot m_{z,k,q}}{\sum_q I_q}
$$

**Additive m/z residual variance (intensity-weighted mixture):**
$$
\sigma_{\varepsilon,i}^2 = \sum_q w_q^2 \cdot c_{mz,q} \cdot N_{z,k,q}^{-\alpha_{mz}}
$$

**Welch–Satterthwaite effective d.o.f.:**
$$
\frac{1}{\nu_{\text{eff},i}} = \sum_{q \in A_i} \frac{w_q^2}{\nu_{I,q}},\quad
w_q = \sigma_q^2 / \sum_{q'} \sigma_{q'}^2
$$

**Joint per-observation log-likelihood (with m/z mask):**
$$
\ell_i = \log t_{\nu_{\text{eff},i}}(I_{obs,i}; I_{pred,i}, \sigma_{I,i})
+ \mathbb{1}[\gamma_{i,q^\star} \ge \gamma_{\text{thresh}}] \cdot
\log t_{\nu_{mz}}(\varepsilon_i; 0, \sigma_{\varepsilon,i})
$$

**Zero-cell log-probability (Poisson zero-event):**
$$
\log P(N = 0 \mid c) = -N_{pred,c}^{total}
$$

**Panel log-likelihood:**
$$
\mathcal{L}_P = \sum_{i \in O_P} \ell_i - \sum_{c \in Z_P} N_{pred,c}^{total}
$$

**Identification score:**
$$
\Lambda_q = \mathcal{L}_P(Q_P) - \mathcal{L}_P(Q_P \setminus \{q\})
$$
