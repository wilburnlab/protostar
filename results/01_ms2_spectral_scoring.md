# 01 — MS2 spectral scoring: why L2 measures mislead

**Status:** planned — in active development (Part I lead). _Part I — MS2._

## Scope
A comparison of MS2 **spectral-similarity scoring properties** on replicate spectra of the same
peptide. Score within-peptide (replicate) and between-peptide pairs under each of cosine
similarity, normalized dot product, Pearson correlation, spectral angle, spectral-entropy
similarity (Li 2021), KL divergence, and the **multinomial deviance** `2·N·KL(p̂ ‖ p_ref)`,
**stratified by fragmentation mode** (MS2 profiles vary by mode). Demonstrate that the L2-geometry
measures (cosine / dot / Pearson / spectral angle) **inflate apparent differences for intense
fragment ions** — their squared penalty is dominated by the shot-noise deviation of the largest
channels — and, decisively, that this inflation **scales with N** (the total ion count): cosine
dissimilarity between two replicates shrinks as N grows, a pure sampling artifact, while
`2·N·KL` is flat in N (asymptotically χ²_{K−1}). This is the empirical face of the analytical
spine carried into [`02_ms2_multinomial.md`](02_ms2_multinomial.md): KL *is* the multinomial-scale
comparator, because the multinomial NLL of `x` vs `p_ref` is `N·KL(p̂ ‖ p_ref)` up to a constant.

### Sub-analysis — what is the published library, really?
Resolve the provenance of the published ProteomeTools `.msp` reference. **Resolved (Branch B):**
a direct check of the real libraries on OSC (`FTMS_HCD_35_annotated_2019-11-12.msp`, 479,261
spectra) shows every `Comment:` carries **only** `Parent`/`Mods`/`Modstring`/`iRT` — **no
`Scan`, `RTInSeconds`, or collision-energy field** — and the peak intensities are **large
non-normalized integers** (e.g. 14792, 111970, 307036), not 0–1 / 0–10000 reals. So the library
cannot be traced to a source scan by metadata, and its count-space integer intensities are
already strong evidence it is a **summed/consensus** spectrum, not a normalized model prediction
or a single-scan dump. The forensic therefore runs by scoring, exactly as designed:
(i) tabulate the MSP intensity quantization (integer count-space ⇒ consensus); (ii) score the MSP
against the **sum-consensus of all PSM-matched spectra** and show the match improves as replicate
count R grows (the consensus fingerprint); (iii) MSP-vs-best-single-scan (absence of one perfect
match rules out a single representative). The `iRT`-only RT (predicted, not measured) corroborates
an aggregate product. This justifies using the PSM-sum consensus as the orthogonal `p_ref` the
multinomial chapter compares against.

## Source
- Experiment script: `pipelines/experiments/01_ms2_spectral_scoring.py` (TBD)
- Identifications: paired MaxQuant searches (`massspec.io.maxquant`, ledger #11 **landed**),
  PROCAL/calibrant + recurring peptides; guided/filtered by the MSP reference
  (`massspec.io.msp`).
- Intermediates: MS2 fragment intensities via XIC level-2 `assigned_scans_only`
  (`massspec.quant.chromatogram`); analyzer-appropriate tolerance (FTMS 20 ppm / ITMS 0.5 Da).
- Model: `massspec.spectra.similarity.compare_spectra` (ledger #12, **prototyping**) +
  `massspec.spectra.consensus.build_consensus` (ledger #13, **prototyping**); the existing
  `core.stats.distributions.Multinomial`.

## Key questions
- Per mode, how separable are within- vs between-peptide pairs under each metric?
- How strongly does each metric's within-pair score depend on N (the artifact magnitude)?
- Is the published `.msp` a consensus sum, a single representative scan, or a prediction?

_Figures (planned):_ (1.1) metric distributions, within vs between, faceted by mode;
(1.2) KL-vs-cosine scatter colored by N; (1.3) one illustrative peptide where cosine penalizes an
intense ion's shot noise but KL does not.

_Findings to be written once exp01 runs (gated on #12/#13 landing)._
