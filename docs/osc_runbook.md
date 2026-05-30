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
python -c "from constellation.massspec.io.thermo import is_thermo_available; print('thermo:', is_thermo_available())"
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

## 5. Convert to parquet (stage 10 — needs a compute node for the Thermo DLLs)
```bash
# 5a. Smoke-test one file:
salloc -A <allocation> -p cardinal -N 1 -n 1 -c 8 -t 1:00:00
conda activate constellation
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd ~/protostar
python pipelines/10_convert_raw.py --dataset Zolg2017 --limit 1

# 5b. Validate profile on 1-2 files (gating step):
python pipelines/10_convert_raw.py --dataset Zolg2017 --limit 2 --profile
#   confirm proc/Zolg2017/profile/<stem>/ exists, more peak rows/scan than centroid,
#   peak_resolution/noise/baseline null. If wrong -> STOP.
exit   # leave salloc

# 5c. Full centroid pass (one SLURM job/dataset; convert.sbatch sets LD_LIBRARY_PATH itself):
for DS in Zolg2017 Gessulat2019 Wilhelm2021; do
  sbatch --account=<allocation> --export=ALL,PROTOSTAR_DIR=$HOME/protostar,DATASET=$DS \
    pipelines/slurm/convert.sbatch
done
squeue --me

# 5d. Full profile pass (only after 5b checks out):
for DS in Zolg2017 Gessulat2019 Wilhelm2021; do
  sbatch --account=<allocation> --export=ALL,PROTOSTAR_DIR=$HOME/protostar,DATASET=$DS,MODE_FLAG=--profile \
    pipelines/slurm/convert.sbatch
done
```

## 6. Reconcile
```bash
python pipelines/00_fetch_raw.py  --dataset all --dry-run             # fetch completeness
python pipelines/10_convert_raw.py --dataset all --dry-run           # centroid converted vs planned
python pipelines/10_convert_raw.py --dataset all --dry-run --profile # profile
```
Expect raw counts ~ 1460 / 888 / 1865 (Zolg / Gessulat / Wilhelm); converted == present per mode.

## Notes
- Idempotent throughout — re-running after an interruption is the recovery path.
- `convert.sbatch` leaves `#SBATCH --account` commented; pass `--account=<allocation>` on submit
  (or uncomment it). For very large datasets, fan out: `--array=0-19 ... ,N_SHARDS=20`.
- Everything keys off the `.raw` files; `proc/` bundles are rebuilt from scratch (no cache reuse).

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
  from constellation.massspec.io.thermo._netruntime import get_dll_dir, _missing_dlls
  d = get_dll_dir(); print("dll_dir:", d)
  print("missing DLLs:", _missing_dlls(d) if d else "registry found nothing")
  PY
  ```
  - pythonnet MISSING → `conda env update -f environment.yml --prune` in the constellation repo
    (pulls in `pythonnet` + `dotnet-runtime`), or `pip install pythonnet` into the env.
  - dll_dir None → `bash scripts/install-thermo-dlls.sh` (in the constellation repo), or
    `export CONSTELLATION_THERMO_HOME=<constellation>/third_party/thermo/current`.
  - missing non-empty → the DLL pack is incomplete; re-fetch it.
