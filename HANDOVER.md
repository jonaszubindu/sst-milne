# Handover: state, history, open items (2026-10-14)

Read `README.md` first for installation and usage. This file records what
exists where, what has been run, the problems already solved, and what is
still open. Facts marked **(to confirm)** were not visible from the code
side and should be checked with the owner (Jonas Zbinden) or on the machines.

## Scientific context

ME inversions of SST/CRISP Fe I 6302 time series of forming sunspots,
first AR 13010 (2022-05-16; partial penumbra forming 08:30–13:00 UT in the
southern sector, filaments linked to MMFs; Zbinden et al. 2026, A&A 705,
A35), then 2–3 comparison datasets. Purpose: maps of B, inclination,
azimuth, v_los and intensity per scan, as the basis for tracking MMFs and
granules at SST resolution. **What exactly to track or measure has not
been decided yet** — the owner wants to see the overviews/movies first.

It complements the HMI study in
[moat-flow-onset](https://github.com/jonaszubindu/moat-flow-onset)
(results in its `VETTING.md`): at HMI resolution the moat outflow starts
~6 h and tracked MMFs ~9 h after the penumbra onset (sections J–R), and
for AR 13010 the HMI penumbra onset and first MMFs fall inside the SST
formation interval (sections M, R). Planned "test 3" there: correlation
tracking on the SST continuum, degraded to HMI resolution, to see whether
an outflow or edge inflow is visible at SST but lost at HMI.

## Where things are

| what | where |
|---|---|
| cluster login | `zbindenj_aiub@147.86.10.94` (owner's SSH alias `fhnw_server`), Slurm partition `aiub`, scratch `/unibe/scratch/zbindenj` |
| this repo | GitHub `jonaszubindu/sst-milne` (private); Mac `~/Desktop/sst-milne`; sml server `/sml/zbindenj/sst-milne`; cluster `/unibe/scratch/zbindenj/sst-milne` |
| pyMilne (upstream, unmodified) | Mac `~/Documents/GitHub/pyMilne` (built: `pyMilne.cpython-311-darwin.so`, conda env `Milne`); server `/sml/zbindenj/pyMilne` (venv `~/venvs/milne`, Python 3.8); cluster `/unibe/scratch/zbindenj/pyMilne` (conda env `milne`) -- all built |
| AR 13010 cube | sml `/sml/zbindenj/idl/SST_data/4th_batch/2022-05-16/CRISP/nb_6302_2022-05-16T08:28:21_08:28:21=0-19_stokes_corrected_export2023-02-24T09:32:31_im.fits` — 20 scans, 13 tunings 6302.198–6302.758 Å, 1315×1313 px at 0.0592″/px, NaN padding |
| 2025-06-17 cube | cluster, `20250617/crisp/nb_6302_2025-06-17T16:00:07_16:00:07=0-50_stokes_corrected_im.fits` relative to the job's submit dir — 51 scans, 12 tunings 6302.202–6302.682 Å, 2943×2945 px at 0.0441″/px |
| other comparison cubes | **(to confirm)** list them in `cubes.txt` on the cluster |
| results | sml: `/sml/zbindenj/me_AR13010` (default two-pass run, possibly stopped/partial), `/sml/zbindenj/me_all` (single-pass batch, if started); cluster: `/unibe/scratch/zbindenj/me_out/<cube tag>/` (job 110587) and `/unibe/scratch/zbindenj/me_out_17/<cube tag>/` (job 110736) -- completeness: `tools/me_status.py` |
| prepared AR 13010 maps (test data only) | Google Drive `My Drive/Daniel-Master-Arbeit/Maps_for_Daniel/Prepared_map_SIR_AR13010_nb_6302_2022-05-16T08:28:21_frame_{0,11,19}.fits` — (4, 12, 550, 600); 13-point grid 6302.4931 + (−0.28…+0.20) Å in 40 mÅ steps, slot 1 a −1 placeholder |

## Run history

1. Mac: pyMilne built (see below), verified on pyMilne's `example_CRISP`
   and for 1- vs 8-thread identity; pipeline tested on SSTRED-layout test
   cubes from the prepared maps (`tools/make_test_cube.py`), including
   added NaN/zero padding and an empty tuning.
2. sml server, AR 13010: first run crashed after scan 0 (NaN in the
   quiet-Sun median from the NaN padding) → fixed (commit 9b1cf9e). A
   default run into `me_AR13010` and a `--single-pass` folder batch into
   `me_all` were started; the owner later stopped the server run to move
   to the cluster.
3. Cluster (Slurm array, job 110736, partition `aiub`, 384 threads/node):
   most cubes ran fine; the 2025-06-17 cube skipped every scan ("only 0
   valid pixels") because one tuning had no data → fixed (977c12c): such
   tunings are now left out per scan. The skipped tasks need to be
   resubmitted (`sbatch --array=<tasks> ...`, same `cubes.txt` and output
   folder) -- done (confirmed 2026-10-14). The pyMilne builds on the
   server and the cluster are in the listed places (confirmed).
4. Earlier cluster job 110587 (output `/unibe/scratch/zbindenj/me_out`,
   2026-10-13, before the empty-tuning fix) also skipped cubes with "only
   0 valid pixels", e.g. `20250625/crisp/nb_6302_2025-06-25T10:54:22...`
   (3 scans, 12 tunings, 2805x2821 px): rerun those tasks with that job's
   cube list and `me_out` (`tools/me_status.py` lists them).

5. 2026-10-15: four 2025-06-17 cubes have **no polarimetry** in the
   reduced data -- Q/U/V NaN everywhere, only Stokes I (and that in ~36 %
   of the field, the rest padding): `6302_2025-06-17T10-26-05` (27 scans),
   `...T12-38-46` (10), `...T14-53-01` (9), `...T16-00-07` (51). Not an
   inversion problem: they need a re-reduction / re-export with Stokes
   (SSTRED polcal) before ME can be run. `...T08-03-03` (8/9 scans) and
   `...T09-20-37` (19/31) are probably partly affected the same way.
   Skips are now recorded in `skipped_scans.json` per cube (see
   `tools/me_status.py`).

## Problems already solved (do not re-debug)

- **pyMilne on Apple silicon** (`pymilne_build/setup_macos_arm64.py`):
  upstream `setup.py` uses x86-only flags; Homebrew GCC 14.2 cannot compile
  on macOS 27 (its include-fixed headers clash with the SDK); with Apple
  clang: C++17 (Eigen 3.3.7 needs `std::result_of`), `std::abs<T>` → `std::abs`
  and `schedule(dynamic)` → `schedule(static)` in a build copy (the conda
  env's llvm-openmp 16 lacks `__kmpc_dispatch_deinit`), `-DEIGEN_DONT_PARALLELIZE`
  plus `-include omp.h`, `-include iostream`. Results identical across thread counts.
- **`crisp.py`** (CRISP transmission profile) ships only in pyMilne's
  `example_CRISP/`; `me_invert.py` searches there too.
- **SSTRED padding**: NaN/zero pixels around the rotated field → excluded
  from statistics, given a quiet-Sun profile during the fit, NaN in output.
- **Empty tunings**: left out per scan; the far-wing reference moves to the
  farthest tuning with data (header `WINGWAV`, `NTUNE`).
- **Split observations** (SSTRED splits a too-large series into cubes with
  the same start time, e.g. `...T06:55:31=0-39` and `=40-79`): cube_tag
  used to ignore the scan range, so both parts shared one folder and
  overwrote / skipped each other's `me_scanNNNN.fits` (2025-06-25 06:55 in
  `me_out_25`). Now parts not starting at scan 0 get `_s<a>-<b>`, a scan
  file whose `CUBE` header names another cube is redone, and me_status
  warns about shared folders. Merging the parts into one time series is
  not done (scan numbers restart at 0 in each part; the true numbers are
  in the cube's `VAR-EXT-SCANNUM`).
- **Scans without Q/U/V** (Stokes-I-only reductions): skipped with the
  reason "no Q/U/V (Stokes I only)", recorded in `skipped_scans.json`;
  `me_status.py` counts them as done-but-skipped, not as missing.
- **FITS headers reject NaN** → `fin()` writes 'NaN' strings.
- **Background runs**: `nohup sh -c "python ..."` picked a Python without
  the packages; use the full interpreter path or `python -u` directly;
  Python buffers redirected output (`python -u`, start-up lines flushed).
- **pyMilne's progress counter** (carriage returns) floods logs; the Slurm
  script filters it, otherwise `tr '\r' '\n' < log | grep '^scan'`.
- **ffmpeg**: `/usr/local/bin/ffmpeg` on the Mac is an x86 build that does
  not run; `me_quicklook.py` picks the first ffmpeg that actually runs.

- **Pulling quicklooks**: run `tools/pull_quicklooks.sh` on the Mac (rsync
  of `quicklook/` folders only). Do not use `scp -r` on the output folders:
  it copies the ~0.4-0.5 GB `me_scan*.fits` files, and with several remote
  arguments scp copies them INTO the last folder on the cluster.

## Open items / decisions for the owner

1. **Keep the better of the two passes per pixel** on the first scan
   (pass 2 currently replaces pass 1 unconditionally; ~3 lines, like
   `repair_outliers`). Offered, not decided.
2. **Seeing quality per scan** (r0 from `VAR-EXT-ATMOS_R0`, quiet-Sun
   granulation contrast, median chi2 vs time) to flag bad scans before
   tracking. Offered, not decided.
3. **What to track in the ME maps** (MMFs in B_los, granules in I_wing):
   to be decided after looking at the overviews/movies.
4. Coherent longitudinal-solution patches where Q/U are weak: pyMilne's
   spatially regularised inversion if they matter.
5. v_los zero point (currently relative to the quiet-Sun median per scan);
   azimuth/axes in the image frame — SST images are not necessarily
   rotated to solar north.
6. No automated tests: check changes with `tools/make_test_cube.py` (+
   padding / empty tuning variants as in the history above) and a short
   `--scans 0-1 --crop ...` run.
