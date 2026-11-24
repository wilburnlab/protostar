# Global chromatogram extraction for N(t): design and pilot record

*2026-11-23. Approved design (PI, 2026-11-21); pilots P1–P4 run on OSC 2026-11-23. Drivers are in
`pipelines/experiments/pilots/`. The Constellation side is specified for a dedicated session in
Constellation's `docs/plans/global-extraction.md`.*

## Why

The goal is a latent ion-count trajectory N(t) per peptide per injection, plus the MS1→MS2
conversion. The PI chose three things:
- a **global search**: every synthesized sequence in every run, every isotopologue × charge 1–6 in
  MS1, and fragments in MS2;
- **in-pool / off-pool** labels as the confidence model, not decoys (decision log D19);
- an extraction that leaves three questions to the model, which MS_Quant never settled: window
  boundaries, which points enter the likelihood, and how many zeros to include (D18).

This replaces the June stage 30 (`30_extract_intermediates`), which extracted only identified PSMs
of seed peptides at z = 1–4 within ±60 s.

## Data flow

1. **Tier 0 (protostar).**
   - `protostar.search_space`: 824,863 sequences as ProForma, with masses and composition keys;
     `targets_for_run` gives (sequence, z = 1..6) inside the run's MS1 range.
   - `protostar.pools`: stem → pool and injection, run order per instrument across datasets, and
     per-run labels (`standard`, `in_pool`, `carryover`, `shadowed`, with `off_pool` implied).
2. **Tier 1 (Constellation `CandidateScreenReducer`, every query in every run).**
   - An MS1 sweep, precursor-major, keeping matches only. The XIC is reduced inside the reducer and
     discarded; only candidate rows are written.
3. **RT transfer (protostar over `rt.align`).**
   - Per-run maps from PROCAL plus the run's own identifications.
   - Each peptide gets a reference RT from its in-pool runs.
4. **Tier 2 (Constellation model-input reducer, a bounded set per run), written to disk.**
   - In-pool and standard peptides over the whole gradient.
   - Off-pool, carryover and shadowed controls at the transferred RT ± 6 min.
   - Zeros are kept as censored cells, and MS2 is linked to its MS1 brackets.
5. **Counter fits from the Tier-2 record**, as separate jobs. Window width, point inclusion and zero
   handling are chosen at fit time.

## Pilot results

All pilots ran on Zolg2017 runs with Constellation 6e47177 on Cardinal.

### P1: cost of the global MS1 screen (job 17098364)

Setup:
- Full screen: **2.48M targets per run** (824,863 sequences × z = 1..6 inside 360–1,300 m/z).
- Channels M−1..M+3 at 20 ppm, rule `nearest`, Constellation's `QuerySummaryReducer` (matches only).

| Run (pool 1) | MS1 scans | Matches | Precursor-major, 16 workers |
|---|---|---|---|
| DDA | 3,901 | 4.0 × 10⁹ | 87 s |
| 2xIT_2xHCD | 3,183 | 4.7 × 10⁹ | 113 s |
| 3xHCD | 3,230 | 5.0 × 10⁹ | 121 s |
| ETD | 2,303 | 2.8 × 10⁹ | 68 s |

- **Worker scaling (precursor-major, DDA):** 4 / 8 / 16 workers → 257 / 143 / 87 s.
- **Compile:** about 60 s, identical for every run with the same MS1 range, so compile once and reuse.
- **Scan-major does not scale with workers.** At 30% of the space it took 690–745 s at 4, 8 and 16
  workers alike, spilled 16.5 GB and used 27.5 GB per worker. Precursor-major took 30–78 s, spilled
  0.03 GB and used under 10 GB per worker.
- **Dense (zero-filled) screening is 25–30× dearer.** 525M cells took 31.7 s at 16 workers, which
  projects to about 48 min per run. The screen therefore keeps matches only.
- **Matches are 8–12% of MS1 cells at 20 ppm.** MS1 scans are dense: about 1,300 centroids each.
- **Corpus estimate:** about 100 s per run on 16 CPUs, about **120 node-hours (about 1,900 core-hours)**
  for 4,212 runs.

### P2: what separates in-pool from off-pool (job 17104326)

Setup:
- **Block:** 21 consecutive runs on one instrument, all DDA (Zolg2017 first_pool_1–11 and
  second_pool_33–42), one run per pool. Each block peptide is in-pool in one run and off-pool or
  carryover in the other 20.
- **Queries:** the 21,050 block-pool sequences, the standards, and a hashed 2% null (16,257
  sequences), at z = 1–6.
- **Rows:** a pilot reducer kept every (query, MS1 scan) with at least two of M0–M+2 matched: about
  19M rows per run. Features, regions and thresholds were computed offline.
- **Truth:** MaxQuant in-pool identifications at PEP ≤ 0.01. A region is "true" when it contains the
  PSM RT ±30 s.

**Recall** (top 5 regions per query, gap ≤ 2 MS1 scans):

| KL max | \|ppm\| max | min scans | In-pool recall | PROCAL recall | Null regions per query |
|---|---|---|---|---|---|
| 0.1 | 10 | 2 / 3 / 5 | 0.951 / 0.951 / 0.950 | 0.998 | 2.43 / 1.87 / 1.29 |
| 0.5 | 10 | 3 | 0.948 | 0.995 | 2.22 |
| ∞ | 10 | 3 | 0.949 | 0.995 | 2.27 |
| any | 20 | 3 | 0.928–0.934 | 0.990–0.993 | 2.38–2.69 |

A 20 ppm window lowers recall: chance matches crowd the top 5. The run's PROCAL mass offset
(0.65–2.1 ppm across the block) is subtracted before the cut.

**Separation.** AUC, true in-pool regions vs regions of clean off-pool sequences, at KL ≤ 0.1,
\|ppm\| ≤ 10, ≥ 3 scans:

| Feature | AUC |
|---|---|
| log Σ flux (Σ I·iit, iit in s) | 0.913 |
| qualifying scans | 0.899 |
| \|median ppm\| (offset-corrected) | 0.882 |
| median KL | 0.802 |
| median cosine | 0.792 |
| max cosine | 0.755 |
| median ppm spread across isotopes | 0.671 |
| ppm MAD across scans | 0.608 |
| max M−1 ratio | 0.404 (inverted) |
| count-aware deviance (α = 10·z, seconds) | 0.349 (inverted) |

- **The count-aware deviance points the wrong way.** True, bright signals score worse. This is nb42b's
  failure mode (decision log D9 history) reproduced on 21 runs. A screen must not use count-scaled
  goodness of fit.
- **The M−1 ratio also points the wrong way.** True peptides show an M−1 neighbour more often.
- **A combined score helps.** A logistic score with run-held-out cross-validation reaches **AUC
  0.965**. The dominant weights are flux, |ppm|, max cosine and KL.

The trade-off against volume, at the loosest setting (about 7.9M null regions per run before a
threshold):

| True regions kept | Off-pool regions kept | Rows per run | Corpus (≈ 60–100 B/row) |
|---|---|---|---|
| 99% | 55% | about 4.4M | 1.3–1.9 TB |
| 95% | 20% | about 1.5M | 0.4–0.6 TB |
| 90% | 8% | about 0.6M | about 0.2 TB |

**Carryover and shadowing.**
- Carryover regions sit at the peptide's own-pool RT (±60 s) in 4.5% (one pool back) and 4.3% (two
  pools back) of cases. A rough chance rate is about 3%, so carryover from the previous one or two
  pools is weak at most.
- Shadowed sequences' regions are RT-coincident (±30 s) with the in-pool sequence they shadow in
  5–7% of cases, against about 2% by chance. Shadowing is real but a minority. Labels keep
  `shadow_of`, so it can be checked.

**Limits.**
- DDA only, one instrument block.
- PSM RTs are imperfect truth (see P3).
- The other injection types must be checked before the score is fixed.

### P3: RT transfer (job 17104327)

Setup:
- 140 runs (the P2 block's pools with all their injections, 5 Zolg2017 ↔ Gessulat2019 pool pairs,
  and 2 Wilhelm2021 `noincl` reruns), MaxQuant PSMs at PEP ≤ 0.01.
- PROCAL's published 60-minute RTs as the reference clock.
- `rt.align` maps per run: PROCAL alone, or "refined" (plus the run's own pool IDs, 5-fold held out).

Residuals, observed − transferred (s), refined map:

| Regime | n | Median | SD (MAD) | p95 \|r\| | p99 \|r\| | > 300 s |
|---|---|---|---|---|---|---|
| Within pool (block DDA runs ← other injections) | 17,620 | 1.3 | 10.5 | 243 | 557 | 4.1% |
| Cross pool (multi-pool peptides) | 5,016 | 0.3 | 21.5 | 108 | 348 | 1.4% |
| Cross dataset (Wilhelm `noincl` ← Zolg) | 1,705 | 0.6 | 10.4 | 24 | 96 | 0.06% |

- **Map choice.** PROCAL-only maps leave a median offset of about +5 s; refining removes it.
- **The tails are peptide-specific.** Refining the map does not shrink them. The within-pool tail is
  heavier than the cross-dataset one, which points at the truth rather than the transfer: a DDA PSM's
  RT can be a re-sampling in the peak's tail, or a wrong identification.
- **Tier-2 window.** ±6 min covers ≥ 98.6% cross-pool and 99.9% cross-dataset. Off-pool detections
  outside the window are caught by Tier 1, not Tier 2.
- **Re-measure against screen apexes.** Replace PSM RTs with Tier-1 apexes before setting the window
  for production.
- **Constellation finding.** With torch's default thread count, a 40-point `rt.align` fit took 66 s on
  a 104-core node, against 0.04 s on one thread. The driver pins one thread.

### P4: the Tier-2 model-input record (job 17105067)

Setup:
- One pool's 1,135 sequences (pool 1 plus PROCAL) at z = 1..6, M−1..M+3, fragments b/y (c/z for ETD)
  up to min(3, z−1).
- `extract_xic(levels=[1, 2], rule="nearest", drop_unmatched=False)`, precursor-major.
- MS2 swept at 20 ppm and at 0.5 Da, each row kept by its scan's analyzer.

| Record | Targets per run | MB per run | KB per (seq, z) | MS2 scans per target | MS1 brackets present |
|---|---|---|---|---|---|
| In-pool, whole gradient (4 runs) | 2,322 | 128–302 | 55–130 | 55–260 | 0.999+ |
| Off-pool, PSM RT ± 6 min (8 runs) | 1,969 | 31–54 | 16–28 | 26–43 | 0.997+ |

- **Speed:** each run took 14–49 s for both sweeps at 16 workers.
- **Projection:** about 430 MB per run for Tier 2, **about 1.8 TB** for the corpus. Trimming M−1/M+3
  or co-isolated MS2 would cut that.
- **Co-isolation.** DDA isolation windows of *other* precursors cover a peptide's m/z in many MS2 scans
  (55–260 per target). These rows are the raw material for interference and co-isolation modelling.
- **Constellation finding.** With `rule="all"`, `drop_unmatched=False` writes **no zero rows**: 75–1,319
  of 3,230 MS1 scans per channel, in both layouts. `nearest` fills all 3,230. The first P4 submission
  used `all` and was cancelled. Counting the peaks within tolerance therefore waits on Constellation.

## Decisions taken from the pilots

1. **The Tier-1 screen is precursor-major and keeps matches only, with no zero-fill.** About 100 s per
   run.
2. **Screen features:**
   - region qualifying-scan count, flux, offset-corrected |ppm|, KL and cosine of the renormalised
     isotope shares, ppm spread and MAD;
   - the M−1 ratio, kept for diagnosis;
   - no count-scaled deviance.

   The components are written, and the combined score is learned from labels (refit per injection
   type).
3. **Thresholds by label.**
   - Queries labelled in-pool, standard or carryover are emitted without a threshold: about 3k
     sequences per run, so cheap.
   - Shadowed queries (about 240k sequences per run) and off-pool queries are thresholded at a
     combined-score cut. The PI chooses between 95% (about 0.5 TB) and 90% (about 0.2 TB) of true
     regions kept.
   - A hashed 1% sentinel set is emitted unthresholded as the unbiased null.
4. **The Tier-2 window is the transferred RT ± 6 min,** using refined maps (PROCAL plus the run's own
   IDs).
5. **The Tier-2 record uses rule `nearest`,** with zeros as acquired-but-unmatched cells. MS2 tolerance
   is per analyzer.

## Open before production

- Repeat P2 on a non-DDA block (2xIT_2xHCD, 3xHCD, ETD) and refit the combined score per injection
  type.
- Re-measure RT transfer against Tier-1 apexes instead of PSM RTs.
- **Constellation (hand-off):**
  - `CandidateScreenReducer` and the Tier-2 reducer;
  - MS2 tolerance per analyzer;
  - isolation gating by the whole isotope envelope;
  - zero rows under rule `all`;
  - `rt.align`'s thread count;
  - a note that scan-major does not scale for global sweeps;
  - CR-A6 (Counter injection time in seconds).
- **P5,** testing window width, inclusion and zero handling with N(t) fits on Tier-2 data, after CR-A6,
  CR-B1 and CR-B2 land.
