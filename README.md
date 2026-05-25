# protostar

**ProToStaR** — **Pro**teome **To**ols **Sta**tistical **R**esearch.

A reproducible, open-science re-analysis of the [ProteomeTools](https://www.proteometools.org)
synthetic-peptide datasets, characterizing the physiochemistry of peptide fragmentation and
the statistics of Orbitrap mass-spectrometry measurements. These are the calibrated priors and
confirmed findings that seed the lab's integrative models — chiefly the **Counter** generative
quantification model.

ProteomeTools peptides are synthetic with known sequences, so they isolate instrument behavior,
fragmentation physics, and chromatographic reproducibility from biological variation. Because
the data are public, this re-analysis is published in the spirit of open science: every figure
and finding traces back to a script you can re-run from the raw data.

> *A protostar is the gravitationally-collapsing cloud from which a main-sequence star forms —
> the foundational, formative stage. ProToStaR is that stage for the lab's MS tooling: the
> place where raw observations condense into the priors the mature tools rely on. It lives
> inside the [Constellation](https://github.com/wilburnlab/constellation).*

## Built on Constellation

protostar is a **thin** layer over [Constellation](https://github.com/wilburnlab/constellation)
(`constellation.massspec` + `constellation.core`). It contains **no model code**: all
likelihood, peak-shape, calibration, and scoring machinery — including the Counter model —
lives in Constellation. This repo holds data orchestration, experiment drivers, and curated
records. When an analysis needs a new modeling capability, it is contributed upstream to
Constellation rather than written here. See [`CLAUDE.md`](CLAUDE.md) for the architecture and
the contribution workflow.

## Datasets

| Dataset | Publication | Accession |
|---|---|---|
| Zolg 2017 | Zolg et al., *Nat. Methods* 2017 | [PXD004732](https://www.ebi.ac.uk/pride/archive/projects/PXD004732) |
| Gessulat 2019 | Gessulat et al. (Prosit), *Nat. Methods* 2019 | [PXD010595](https://www.ebi.ac.uk/pride/archive/projects/PXD010595) |
| Wilhelm 2021 | Wilhelm et al., *Nat. Commun.* 2021 | [PXD021013](https://www.ebi.ac.uk/pride/archive/projects/PXD021013) |

Reference spectral libraries: the published ProteomeTools `.msp` files from
<https://www.proteometools.org/index.php?id=53>.

## Install

**Lab (recommended):** reuse the shared `constellation` conda env and add this repo editable.

```bash
conda activate constellation
pip install -e .
# optional, for OSC job/sync helpers:
pip install -e ~/projects/osc_helper   # provides the [osc] extra
```

**External / reproducible:** create a self-contained env that pulls Constellation from PyPI.

```bash
conda env create -f environment.yml
conda activate protostar
```

## Pipeline (raw-first)

```
00_fetch_raw         download .raw from PRIDE + verify checksums (resumable/repairable)
10_build_mzpeak      .raw -> mzpeak Parquet + scanmeta (rebuilt fresh via Constellation)
15_reference_library ingest published .msp libraries (+ optional EncyclopeDIA/Scribe re-search)
20_build_metadata    acquisition time table (datetime + instrument + acquisition order)
30_extract_intermediates   common MS1/MS2 extracted chromatograms (PROCAL + library peptides)
experiments/         one reproducible script per solidified analysis -> parquet + figures
results/             curated "textbook" records synthesizing the confident findings
```

Large data and bulk artifacts are canonical on OSC ESS
(`/fs/ess/<allocation>/<group>/protostar/`); only small dev subsets live locally. The curated
`results/` records and their figures are committed.

## Status

Scaffold phase — directory architecture and conventions are in place; pipeline stages and
experiments are stubs. See [`CLAUDE.md`](CLAUDE.md) and
[`docs/constellation_contributions.md`](docs/constellation_contributions.md).

## License

[Apache-2.0](LICENSE) — matching [Constellation](https://github.com/wilburnlab/constellation).
See [`NOTICE.md`](NOTICE.md) for dataset citations and attributions.
