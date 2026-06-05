# Spectral Similarity: L1 vs L2 Geometry and Why Dot-Product Measures Mislead

> **Scope note (protostar).** Reference/theory note, not model code. The similarity
> functions themselves (`spectral_angle`, `cosine_similarity`, `kld`,
> `multinomial_deviance`, `compare_spectra`) live in **Constellation**
> (`constellation.core.stats.losses`, `constellation.massspec.spectra.similarity`).
> This document records *why* the exp02/exp03 MS2 work prefers the information-geometry
> comparators (KL / multinomial deviance) over the L2 / dot-product family (cosine,
> normalized dot, spectral angle) for fragment-ion spectra, and writes down the exact
> formulas as we use them.

**TL;DR.** A centroided MS2 spectrum is a draw of `N` ions partitioned over `K` fragment
channels — a *multinomial* sample on the probability simplex, not a vector in Euclidean
space corrupted by additive noise. The matched divergence is therefore the KL / multinomial
deviance (an L1-/simplex-native, information-geometry measure), not the cosine or spectral
angle (an L2-/sphere-native, least-squares measure). Dot-product measures are the maximum-
likelihood comparator under the *wrong* noise model: they are dominated by the few most
intense ions, over-penalize the *expected* sampling wobble of those dominant ions, stay
nearly blind to low-abundance channels, and are not comparable across peptides.

---

## 1. Notation

For two non-negative intensity vectors `a`, `b` over `K` fragment channels (e.g. an observed
per-scan spectrum and a consensus/reference):

- **L1-normalized proportions** (the simplex representation):
  `p̂_k = a_k / Σ_j a_j` and `q_k = b_k / Σ_j b_j`; both sum to 1.
- **L2-normalized directions** (the sphere representation):
  `â = a / ‖a‖₂` and `b̂ = b / ‖b‖₂`, where `‖a‖₂ = sqrt(Σ_k a_k²)`.
- **Counts and total:** `x_k = N · p̂_k` with `N = Σ_k x_k` the total ion count. In
  ProteomeTools we only observe `N_proxy = total_intensity · iit` (detector gain unresolved),
  so `N` enters as a *relative* count scale — label any `N` axis "N_proxy (a.u., gain unknown)".

---

## 2. The two geometries

| | **L1 / simplex** | **L2 / sphere** |
|---|---|---|
| Norm | `‖x‖₁ = Σ_k x_k` (intensities ≥ 0) | `‖x‖₂ = sqrt(Σ_k x_k²)` |
| Normalized object | proportions on the simplex (`Σ p_k = 1`) | unit direction on the sphere (`‖x̂‖₂ = 1`) |
| Implied generative model | **Multinomial**: `x ~ Multinomial(N, p)` | **Gaussian**: `a = signal + iid additive noise` |
| Per-channel variance | `p_k(1−p_k)/N` (coupled to the mean) | `σ²` (constant, magnitude-independent) |
| Channel correlations | negative (sum constraint) | none (independent) |
| Matched objective (MLE) | KL / relative entropy → multinomial deviance | least squares → cosine / dot product |
| Sampling-noise scale | `~1/sqrt(N)` shot noise | fixed `σ` |

The physical question is just *which row matches MS2 fragmentation*. Fragment ions are
produced by (largely without-replacement) sampling from a finite precursor pool and then
counted — the **multinomial / L1** row. exp02 established this as bedrock (multinomial beats
iid-Gaussian for fragment proportions; the empirical covariance is consistent with an
imperfect multinomial). The additive-Gaussian / L2 row is the assumption *baked into every
dot-product measure*.

---

## 3. The underlying distribution shapes

**L1 norm → proportions → Multinomial.**

- Normalizing by L1 puts a spectrum on the probability simplex — the natural object when the
  data are *counts* of a fixed total partitioned over channels.
- Under `x ~ Multinomial(N, p)`, each proportion has variance `p_k(1−p_k)/N`, so a **dominant
  ion (large `p_k`) carries large absolute sampling variance**; channels are negatively
  correlated; relative shot noise shrinks as `1/sqrt(N)`.
- The matched (MLE) divergence is the **KL divergence** / relative entropy; its count-scaled
  form `2N·KL` is the multinomial likelihood-ratio statistic `G²`, asymptotically `χ²_{K−1}`
  and therefore *flat in N* — a peptide-comparable scale.

**L2 norm → unit vector → Gaussian.**

- Normalizing by L2 puts a spectrum on a sphere — the natural object when the data are a
  *direction in Euclidean space* and the relevant noise is additive.
- The implied model is `a_k = s_k + ε_k` with `ε_k ~ iid N(0, σ²)`: a **fixed-size fluctuation
  is equally probable on a large or a tiny ion** (variance independent of magnitude), and the
  channels are independent and unconstrained.
- The matched (MLE) objective is **least squares**; the matched similarity is the **cosine /
  dot product**, and the **spectral angle** is a monotone reparametrization of the cosine.

---

## 4. Why dot-product measures are mis-specified for fragment spectra

(cosine, normalized dot product, Pearson, spectral angle — the L2 / cosine family)

- **Wrong likelihood.** They are the MLE-matched similarity under *additive Gaussian* noise,
  but fragment proportions are *multinomial*. The metric therefore mis-weights each channel's
  disagreement relative to its true sampling variance.
- **Dominated by the most intense ions.** `‖a‖₂²` is carried by the few largest channels, so
  the cosine is set mostly by the dominant ions. A small *relative* wobble in a dominant ion —
  which is *expected* multinomial noise (variance `∝ p(1−p)/N`, largest exactly where `p` is
  large) — rotates the unit vector and depresses the score despite carrying no real signal.
  (The exp03 "SA-only" failure mode: a butterfly pair that looks near-identical, spectral
  angle ≈ 0.94, KL negligible.)
- **Nearly blind to low-abundance channels.** A present/absent difference in a small ion barely
  changes `‖a‖₂`, so cosine/SA miss real differences there that forward KL/G² resolve. (The
  exp03 "KL-only" failure mode: high spectral angle, tail KL.)
- **Not comparable across peptides.** Because L2 sensitivity depends on how *concentrated* the
  spectrum is (a dominant-ion-heavy vs a flat spectrum), a fixed cosine/SA threshold means
  different things for different peptides. exp03 quantifies this: peptide identity explains
  ~5× more of the spectral-angle variance than of the KL variance (one-way ANOVA `η²` ≈ 0.11
  vs ≈ 0.02; lower = more peptide-invariant).
- **The angular transform is cosmetic.** `spectral_angle` is a *monotone* function of the
  cosine (it stretches the high-similarity region for readability); it inherits all of the L2
  mis-weighting and changes none of the underlying model.

> These are *model-specification* points, not a claim that cosine is "always wrong" — for
> direction-only matching of well-resolved, high-N library spectra it is often adequate. The
> argument is that for *per-scan* fragment counts with real shot noise, KL/G² is the matched
> comparator and the dot-product family systematically mis-weights both the dominant and the
> faint ions.

---

## 5. Formulas (exactly as implemented)

All three are dispatched through
`compare_spectra(query, reference, *, method, as_counts, pseudocount)`
(`constellation/massspec/spectra/similarity.py`), which broadcasts a `(K,)` reference against
a `(B,K)` query.

### 5.1 Spectral angle (L2 / cosine family)

The cosine of the angle between the two L2-normalized intensity vectors, mapped onto `[0,1]`:

$$
\cos\theta = \frac{\langle a, b\rangle}{\lVert a\rVert_2\,\lVert b\rVert_2}
           = \sum_{k} \hat a_k\,\hat b_k,
\qquad
\mathrm{SA}(a,b) = 1 - \frac{2}{\pi}\,\arccos\!\big(\cos\theta\big)
\tag{1}
$$

- `1` = identical direction; scale-invariant (`SA(a, c·a) = 1` for `c > 0`).
- The plain **cosine similarity** is `cosθ` itself; the **normalized dot product** is the same
  quantity. "Not observed" channels (sentinel `-1.0`) are zeroed before normalization via
  `mask = (target+1)/(target+1+eps)` in `losses.spectral_angle`.
- Call through the `compare_spectra` dispatcher, **not** `losses.spectral_angle` directly, on a
  `(B,K)`-vs-`(K,)` batch (the direct `batch_dim` path raises).

### 5.2 KL divergence (L1 / information-geometry family)

Forward direction — observed proportions against the pseudocount-smoothed reference, in nats:

$$
\mathrm{KL}(\hat p \,\Vert\, q) = \sum_{k} \hat p_k \,\log\!\frac{\hat p_k}{q_k}
\;\ge\; 0,
\qquad
q_k \leftarrow \frac{q_k + \varepsilon}{1 + K\varepsilon}
\tag{2}
$$

- `compare_spectra(query, ref, method="kld")` returns `KL(query ‖ ref)` — the forward,
  multinomial-MLE direction; it is **asymmetric**. `0` = identical, lower = more similar.
- The reference is smoothed by a flat-on-simplex pseudocount `ε` so `log q_k` stays finite when
  the reference leaves a channel empty that the observation populates; the forward direction is
  exactly the one sensitive to that case. At `ε = 0`, `kld` can return `Inf` — prefer
  `multinomial_deviance` (which smooths the reference only) when an N-aware statistic is wanted.

### 5.3 Multinomial deviance (the count-aware bridge)

With counts `x_k = N · p̂_k` and `N = Σ_k x_k`:

$$
G^2 = 2\sum_{k} x_k \,\log\!\frac{\hat p_k}{q_k}
    = 2N \cdot \mathrm{KL}(\hat p \,\Vert\, q)
    \;\xrightarrow[\text{Wilks}]{}\; \chi^2_{K-1}
\tag{3}
$$

- This is the multinomial likelihood-ratio statistic and the identity that ties §5.2 to the
  generative model: `E[G²] ≈ K−1` *independent of N* (here `K = 40` for the 11-mer PROCAL b/y
  basis, so `E[G²] ≈ 39` on a true single-`p` group at matched N — the sanity anchor).
- `compare_spectra(..., method="multinomial_deviance", as_counts=True)` with
  `counts = props · ion_proxy`. Real per-scan groups are over-dispersed relative to a perfect
  multinomial at `N_proxy` (gain unresolved); the exp03 Gamma/χ² fit of the per-spectrum
  KL-from-consensus recovers an *effective* `ν ≈ 2` and effective `N ≈ 10²`, i.e. the spectra
  behave like a much smaller multinomial than `N_proxy` would imply.

---

## 6. Practical upshot for exp02/exp03

- **Default to KL / multinomial deviance** for per-scan fragment-spectrum comparison; reserve
  cosine/spectral angle for direction-only library matching where N is high and uniform.
- **Lead with the per-pair / per-spectrum distribution**, not a single tuned scalar. A KL/G²
  threshold is peptide-comparable (`2N·KL ~ χ²_{K−1}`); a cosine/SA threshold is not.
- **Gate shot-noise-dominated scans** (exp03 uses `total_intensity ≥ 1e6`); below that the
  multinomial variance dominates and *every* comparator inflates.
- **Separate terminal ions** (y1/y2/b1/b2, known-unstable per the exp02 bedrock) from
  mid-ladder when attributing disagreement — terminal-only divergence is an artifact, not
  structure.

---

## See also

- [`model_specification.md`](model_specification.md) — the canonical MS1 ion model (the
  N-aware, multinomial-rooted intensity / m·z likelihoods).
- [`additive_progenitor_likelihood.md`](additive_progenitor_likelihood.md) — the Counter joint
  likelihood that uses these per-observation marginals.
- exp02 `SESSION_NOTES.md` / `memory/ms2-fragmentation-bedrock.md` — the empirical basis for
  the multinomial framing (multinomial > Gaussian; covariance ≈ imperfect multinomial).
- Constellation: `core/stats/losses.py`, `massspec/spectra/similarity.py`,
  `massspec/spectra/consensus.py`.
