# results/ — the ProToStaR "textbook"

Curated records of the confident, established findings from this re-analysis. **Code is
separated from interpretation:** experiment scripts (`pipelines/experiments/`) produce
parquet + figures; each chapter here synthesizes the conclusions and embeds the figures it
cites. A chapter is a reference others (and the Counter manuscript, `../../Counter_Manuscript/`)
can point to — not a lab notebook.

**Committed:** the `.md` chapters, `figures/`, and small summary tables. **Not committed:**
bulk parquet artifacts (canonical on OSC ESS; gitignored).

Each chapter should state, for every finding: the experiment script that produced it, the
datasets/modes it covers, the figure(s), and the confidence/caveats. A finding that depends
on a not-yet-landed Constellation capability does not belong here yet (see
`../docs/constellation_contributions.md`).

## Chapters

File order = reading order (MS2-first; see `../docs/manuscript_roadmap.md` §3). Part I (MS2)
leads, Part II (MS1: two errors, one problem) is organized by feature, Part III is Counter.

| File | Topic |
|---|---|
| `01_ms2_spectral_scoring.md` | **Part I.** MS2 spectral scoring: KL / multinomial deviance vs cosine / dot / spectral angle; why L2 measures inflate differences (and depend on N); MSP-library provenance forensic |
| `02_ms2_multinomial.md` | **Part I.** MS2 as multinomial: per-mode conversion ratios, Var ∝ p(1−p)·N finite-pool competition — the bridge to the MS1 ion count N |
| `03_intensity_and_ion_count.md` | **Part II (a).** Signal ≈ α·N: both m/z precision and intensity shot noise move with N; the IIT-in-variance correction |
| `04_charge_and_alpha.md` | **Part II (b).** Charge scales α = softplus(α₀+α₁z), hence the N behind a signal; both error scales shift with z |
| `05_resolution_and_transient.md` | **Part II (c).** Resolution = transient length = FFT samples: sharpens both errors, resolves isotope fine structure (fitted p_k); m/z error decomposition; closes Part II → Counter |
| `06_chromatogram_shape.md` | **Part III.** Chromatographic peak shape: width, asymmetry, tailing, SNR; HyperEMG fits → N(t) |
| `07_procal_calibration.md` | **Part III.** Per-file PROCAL calibration (GlobalCalibration); cross-mode / cross-dataset consistency |
| `08_counter.md` | **Part III.** Counter additive-progenitor quantification: N_total + credible intervals; identification Λ |
