# OSC execution runbook

Run the data layer on OSC: relocate existing `.raw`, download anything missing, fetch the
reference libraries, and convert to parquet. Paths/allocations come from a local, gitignored
`config/osc.toml` (copy `config/osc.example.toml`), so the commands below need no path flags.

## 0. Get the repo
```bash
cd ~ && git clone git@github.com:wilburnlab/protostar.git && cd ~/protostar
# already cloned? the history has been rewritten, so reset rather than `git pull`:
#   git fetch origin && git reset --hard origin/main
```

## 1. Local config (do this first)
```bash
cp config/osc.example.toml config/osc.toml   # then edit: username, allocations,
                                             # [paths].data (your ESS data root), [paths].seed
```
`[paths].data` is the default data root; `[paths].seed` is the default `--seed-from` source.

## 2. Environment
```bash
module load miniconda3/24.1.2-py310
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate constellation                 # shared lab env: constellation + Thermo DLLs
# Put the env's libstdc++ first, else sqlite3 -> libicui18n -> libstdc++ fails with a
# CXXABI_1.3.15 mismatch (module load puts the older system lib ahead). The pipeline scripts
# self-heal via a one-time re-exec, but exporting here also fixes ad-hoc `python -c` / the CLI:
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd ~/protostar

# sanity:
ls config/manifests/         # Zolg2017.json Gessulat2019.json Wilhelm2021.json (committed)
python -c "from constellation.massspec.readers.thermo import is_thermo_available; print('thermo:', is_thermo_available())"
```

## 3. Move existing data + download the rest (stage 00)
```bash
python pipelines/00_fetch_raw.py --dataset all --seed-from --dry-run   # preview (no changes)

tmux new -s ptfetch                                                    # SEARCH download is ~110 GB
python pipelines/00_fetch_raw.py --dataset all --seed-from --workers 8
#  - the move runs fully before any download (Ctrl-C after "seeded N" = relocate only)
#  - resumable + idempotent

python pipelines/00_fetch_raw.py --dataset all --dry-run --verify      # integrity (Wilhelm2021: size-only)
```

## 4. Reference libraries (stage 15)
```bash
python pipelines/15_reference_library.py --modes all --extract         # ~2.3 GB -> data/libraries/<mode>/
```

## 5. Convert to parquet (stage 10, bundle v3; needs a compute node for the Thermo DLLs)

Constellation ≥ 2026-10-31 writes **bundle v3** in **one pass**: `proc/<dataset>/<stem>/` (the
centroid peaks, `peaks/part-*.parquet`, plus scan and acquisition metadata) and a sibling
`proc/<dataset>/<stem>.profile/` (the MS1 profile points, plus the stored centroids of every MS2
scan, since MS2 was centroided at acquisition). With `precursor_major` on (the default), the peak
bundle also gets a `precursor_major/` copy for MS2 sweeps. Pre-v3 bundles, i.e. the June 2026
corpus in `proc/<dataset>/{centroid,profile}/<stem>/`, raise `OutdatedBundleError` in Constellation.

**What the peaks are (`centroid_algo=thermo`).** A profile scan with a centroid stream (every
ProteomeTools MS1 survey in the pilot) gets the instrument's own centroids, with `peak_resolution`,
`peak_noise` and `peak_baseline`. A centroid-acquired scan (all MS2) gets its stored centroids,
which carry no resolution, noise or baseline: those columns are null. Thermo's `Scan.ToCentroid`
is used only for a profile scan with no centroid stream; no pilot file had one.

Stage 10 never writes into the old trees and never deletes them. `--dry-run` reports them as
"legacy". Delete a dataset's legacy trees only after its v3 conversion verifies (step 6), and only
with the PI's go-ahead.

```bash
# 5a. Pilot (sizes + fidelity). Pick ~6 files: each dataset × injection type, including the
#     nb42b reference file 01625b_GA1-TUM_first_pool_1_01_01-DDA-1h-R2. Run each with and
#     without the precursor-major copy, and record disk use + wall time.
salloc -A <allocation> -p cardinal -N 1 -n 1 -c 16 -t 2:00:00
conda activate constellation
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd ~/protostar
python pipelines/10_convert_raw.py --dataset Zolg2017 --raw-name <stem> --no-precursor-major
du -sh <data_root>/proc/Zolg2017/<stem> <data_root>/proc/Zolg2017/<stem>.profile
#   Also check: peaks agree with the June centroid bundle of the same file (scan count,
#   peaks per scan, m/z and intensity), and MS1 peak_resolution/noise/baseline are populated.
#   MS2 resolution/noise/baseline are null in v3 (June stored NaN). If they disagree -> STOP.
#   Also convert the known-corrupt Wilhelm2021 file (below): expect status=error
#   ("holds no mass-spectrometer data"), not a .NET stack trace.
exit   # leave salloc

# 5b. Full pass (one SLURM job per dataset; convert.sbatch sets LD_LIBRARY_PATH itself).
#     EXTRA_ARGS=--no-precursor-major skips the precursor-major copy (about +15% disk and
#     +28% time over the June bundles; kept for the 2026-11 pass, PI decision).
#     With precursor-major on, convert_batch reads EVERY file in the job before it builds any
#     layout or writes any peak manifest, so progress shows "0/N" for most of the run and a
#     killed job leaves the whole dataset unfinished (stage 10 clears and redoes it on resubmit).
#     Memory peaks while reading: 41-66 GiB at 16 workers, inside the default 77 GiB.
#     2026-11 pass: Gessulat 2h04, Zolg 3h27, Wilhelm 4h47. Wilhelm ends FAILED (exit 1)
#     because of the known-corrupt file below; check the log's "N ok, N skipped, N failed".
for DS in Zolg2017 Gessulat2019 Wilhelm2021; do
  sbatch --account=<allocation> --export=ALL,PROTOSTAR_DIR=$HOME/protostar,DATASET=$DS \
    pipelines/slurm/convert.sbatch
done
squeue --me

# 5c. Acquisition tables: rebuild from v3. A June table is not "present" (it lacks the
#     bundle-schema stamp), so no --force is needed. This overwrites
#     proc/<dataset>/acquisitions.parquet: copy the old table aside first to compare.
#     (2026-11: every row matched June except source_kind, thermo_raw -> thermo_acquisition.)
python pipelines/20_build_metadata.py --dataset all
```

## 6. Reconcile
```bash
python pipelines/00_fetch_raw.py  --dataset all --dry-run    # fetch completeness
python pipelines/10_convert_raw.py --dataset all --dry-run   # v3 current / outdated / pending + legacy
python pipelines/20_build_metadata.py --dataset all --dry-run
```
Expect:
- Raw counts of about 1460 / 888 / 1865 (Zolg / Gessulat / Wilhelm), and v3 current == raw present.
- The **one known-bad source file** below is the exception. Wilhelm2021 converts 1864 / 1865, so
  its dry-run permanently shows "1 pending". That is expected, not a failure.
- The acquisition counts should match the June tables.

## Known-skip: one corrupt PRIDE source file (Wilhelm2021)

`03174_GF7-TUM_second_pool_122_01_01-3xHCD-1hnoincl-R1.raw` is a **57,874-byte aborted
acquisition** — a valid Thermo OLE2 container (correct `Finnigan` magic) with no MS instrument
device, so `convert` raises `ArgumentOutOfRangeException: Instrument index not available …`
(`RawFileReader.SelectInstrument`) and the bundle is never written. **No data is lost:** the same
pool's R2 re-acquisition and its other three injections all convert fine.

Verified faithful to the source (2026-05-31): a fresh re-download is byte-identical to our copy
(`sha256 083d9238…`), and PRIDE's live `Content-Length` is 57,874 — so this is a corrupt file at
the **source archive**, not a truncated download. Wilhelm2021 publishes no checksums, so the fetch
size-check (57,874 == expected) passes and the stub is retained. Do not chase the perpetual
"1 pending"; re-fetching reproduces the identical stub. Since Constellation beedabb (2026-10-26)
the reader checks for an MS device before selecting it, so this file comes back as a per-file
`status=error` ("holds no mass-spectrometer data") rather than a .NET stack trace (confirmed in the
2026-11-20 reconversion pilot; ledger #15). No bundle is written.

## Notes
- Idempotent throughout — re-running after an interruption is the recovery path.
- `convert.sbatch` leaves `#SBATCH --account` commented; pass `--account=<allocation>` on submit
  (or uncomment it). For very large datasets, fan out: `--array=0-19 ... ,N_SHARDS=20`.
- Everything keys off the `.raw` files; `proc/` bundles are rebuilt from scratch (no cache reuse).
- Read bundles only through `protostar.bundles` (it knows the v3 layout and refuses pre-v3 bundles);
  never open `peaks.parquet` or glob `*/manifest.json` directly.

## Troubleshooting
- **`CXXABI_1.3.15 not found` on import** — the env's libstdc++ is being shadowed by the
  spack/system one. The `export LD_LIBRARY_PATH=...` in step 2 fixes it; the pipeline scripts
  also self-heal (one-time re-exec). Affects ad-hoc `python -c` only if you skipped the export.
- **`thermo: False`** — `is_thermo_available()` requires all three: (1) `pythonnet` importable,
  (2) the Thermo DLL pack discoverable, (3) `Data` + `RawFileReader` + `OpenMcdf` DLLs all present
  (the registry only markers on `RawFileReader.dll`, so a pack missing the other two still reads
  False). Pinpoint which:
  ```bash
  python - <<'PY'
  try:
      import pythonnet; print("pythonnet: OK")
  except Exception as e:
      print("pythonnet: MISSING ->", e)
  from constellation.massspec.readers.thermo._netruntime import get_dll_dir, _missing_dlls
  d = get_dll_dir(); print("dll_dir:", d)
  print("missing DLLs:", _missing_dlls(d) if d else "registry found nothing")
  PY
  ```
  - pythonnet MISSING → `conda env update -f environment.yml --prune` in the constellation repo
    (pulls in `pythonnet` + `dotnet-runtime`), or `pip install pythonnet` into the env.
  - dll_dir None → `bash scripts/install-thermo-dlls.sh` (in the constellation repo), or
    `export CONSTELLATION_THERMO_HOME=<constellation>/third_party/thermo/current`.
  - missing non-empty → the DLL pack is incomplete; re-fetch it.
