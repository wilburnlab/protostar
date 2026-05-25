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

| File | Topic |
|---|---|
| `01_chromatogram_shape.md` | Chromatographic peak shape: width, asymmetry, tailing, SNR; HyperEMG fits |
| `02_ms1_intensity_model.md` | MS1 intensity → ion count: α(z), IIT variance, per-peptide ν_I, fitted isotope fractions |
| `03_mz_error_model.md` | m/z error structure: d_mz spacing, per-scan drift μ̂, per-ion residual ε, precision vs ion count |
| `04_procal_calibration.md` | Per-file PROCAL calibration (GlobalCalibration); cross-mode / cross-dataset consistency |
| `05_counter_validation.md` | Counter additive-progenitor quantification: N_total + credible intervals; identification Λ |
