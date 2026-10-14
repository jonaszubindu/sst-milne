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
| this repo | GitHub `jonaszubindu/sst-milne` (private); Mac `~/Desktop/sst-milne`; sml server `/sml/zbindenj/sst-milne`; cluster `/unibe/scratch/zbindenj/sst-milne` |
| pyMilne (upstream, unmodified) | Mac `~/Documents/GitHub/pyMilne` (built: `pyMilne.cpython-311-darwin.so`, conda env `Milne`); server `/sml/zbindenj/pyMilne` (venv `~/venvs/milne`, Python 3.8); cluster `/unibe/scratch/zbindenj/pyMilne` (conda env `milne`) **(to confirm: builds on server/cluster)** |
| AR 13010 cube | sml `/sml/zbindenj/idl/SST_data/4th_batch/2022-05-16/CRISP/nb_6302_2022-05-16T08:28:21_08:28:21=0-19_stokes_corrected_export2023-02-24T09:32:31_im.fits` — 20 scans, 13 tunings 6302.198–6302.758 Å, 1315×1313 px at 0.0592″/px, NaN padding |
| 2025-06-17 cube | cluster, `20250617/crisp/nb_6302_2025-06-17T16:00:07_16:00:07=0-50_stokes_corrected_im.fits` relative to the job's submit dir — 51 scans, 12 tunings 6302.202–6302.682 Å, 2943×2945 px at 0.0441″/px |
| other comparison cubes | **(to confirm)** list them in `cubes.txt` on the cluster |
| results | sml: `/sml/zbindenj/me_AR13010` (default two-pass run, possibly stopped/partial), `/sml/zbindenj/me_all` (single-pass batch, if started); cluster: `/unibe/scratch/zbindenj/me_out_17/<cube tag>/` (job 110736) **(to confirm which are complete)** |
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
   folder) **(to confirm whether done)**.

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
- **FITS headers reject NaN** → `fin()` writes 'NaN' strings.
- **Background runs**: `nohup sh -c "python ..."` picked a Python without
  the packages; use the full interpreter path or `python -u` directly;
  Python buffers redirected output (`python -u`, start-up lines flushed).
- **pyMilne's progress counter** (carriage returns) floods logs; the Slurm
  script filters it, otherwise `tr '\r' '\n' < log | grep '^scan'`.
- **ffmpeg**: `/usr/local/bin/ffmpeg` on the Mac is an x86 build that does
  not run; `me_quicklook.py` picks the first ffmpeg that actually runs.

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
