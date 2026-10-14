# SST/CRISP Fe I 6302: Milne-Eddington inversions

Scan-by-scan ME inversions of SSTRED narrowband cubes with
[pyMilne](https://github.com/jaimedelacruz/pyMilne) (de la Cruz Rodríguez
2019, A&A 631, A153), plus overview plots and movies. Meant for AR13010
and 2–3 comparison datasets; tracking (MMFs, granules) is deliberately
not part of this yet. Split out of
[moat-flow-onset](https://github.com/jonaszubindu/moat-flow-onset), the
HMI moat-flow study it complements.

Agents and new maintainers: read `HANDOVER.md` (state, solved problems,
open decisions) and `CLAUDE.md` (working conventions).

## pyMilne on Apple silicon

pyMilne itself is not part of this repository: clone
https://github.com/jaimedelacruz/pyMilne (here `~/Documents/GitHub/pyMilne`;
elsewhere set `PYMILNE_DIR`). Its `setup.py` does not build on an
Apple-silicon Mac with a current macOS; use the build file in this repo:

```bash
conda activate Milne          # numpy, scipy, astropy, matplotlib, cython, eigen3
python pymilne_build/setup_macos_arm64.py ~/Documents/GitHub/pyMilne
```

It leaves `pyMilne.cpython-311-darwin.so` next to `MilneEddington.py`. The
build file documents every deviation from `setup.py` (Apple clang +
libomp, C++17, three `std::abs<T>` calls and three dynamic OpenMP
schedules patched in a build copy, Eigen threading off). Checked on
pyMilne's own `example_CRISP` (same maps as the reference) and for
1- vs 8-thread identity.

## On a Linux server

```bash
# code (once)
git clone git@github.com:jonaszubindu/sst-milne.git
git clone https://github.com/jaimedelacruz/pyMilne.git

# python env (once) -- either a venv ...
python3 -m venv ~/venvs/milne && source ~/venvs/milne/bin/activate
pip install "numpy<2" scipy astropy matplotlib "cython<3"
# ... or conda (also brings FFTW and Eigen 3.4 without root):
# conda create -n milne -c conda-forge python=3.11 "numpy<2" scipy astropy \
#     matplotlib "cython<3" fftw eigen && conda activate milne

# build pyMilne (once per machine type: it compiles with -march=native)
sst-milne/pymilne_build/build_linux.sh pyMilne
```

`build_linux.sh` checks GCC >= 10, finds FFTW3 and Eigen >= 3.4 (conda
env or system), links Eigen as `pyMilne/eigen3`, builds with pyMilne's
own `setup.py` and runs a smoke test (an inversion with 1 and N threads).
If something is missing it stops and prints what to install. On Linux
none of the macOS workarounds are needed (GCC + libgomp).

Then, with `PYMILNE_DIR` pointing at the checkout (a single cube; for a
whole folder see *All cubes in a folder* below):

```bash
export PYMILNE_DIR=$PWD/pyMilne
nohup python sst-milne/me_invert.py /path/to/nb_6302_..._im.fits \
    --out me_AR13010 --nthreads $(nproc) > me_AR13010.log 2>&1 &
tr '\r' '\n' < me_AR13010.log | grep '^scan'          # progress
python sst-milne/me_quicklook.py me_AR13010 --movie   # when done
```

## All cubes in a folder

```bash
python me_invert.py /path/to/SST_data/ --out me_all --single-pass --nthreads $(nproc)
python me_quicklook.py me_all --movie
```

Folders are searched recursively for `nb_6302*_im.fits` (`--pattern`);
each cube gets its own sub-folder in `--out` (e.g.
`me_all/6302_2022-05-16T08-28-21/`). Files and folders can be mixed. A
failing cube is reported at the end and the others still run; finished
scans are skipped on a rerun. `--single-pass` inverts every scan once
(first scan: one pass, no outlier re-inversion; later scans still start
from the previous one) -- faster, with more stuck pixels than the default.

## On a Slurm cluster

pyMilne parallelises with OpenMP threads inside one process (not MPI), so
the unit of parallel work across nodes is the cube: `slurm/me_invert_array.sbatch`
runs a job array with one cube per task on a whole node (all its CPUs as
pyMilne threads); scans within a cube stay sequential.

```bash
# once: environment (conda-forge) and a pyMilne build on a compute node
conda create -n milne -c conda-forge python=3.11 "numpy<2" scipy astropy \
    matplotlib "cython<3" fftw eigen ffmpeg
sbatch sst-milne/slurm/build_pymilne.sbatch          # -march=native: build on a node

# per batch of cubes
find /path/to/SST_data -name 'nb_6302*_im.fits' | sort > cubes.txt
sbatch --array=0-$(( $(wc -l < cubes.txt) - 1 )) \
    sst-milne/slurm/me_invert_array.sbatch cubes.txt /path/to/me_out --single-pass
```

Defaults: partition `aiub`, 24 h, `PYMILNE_DIR=/unibe/scratch/zbindenj/pyMilne`,
`SST_MILNE_DIR=/unibe/scratch/zbindenj/sst-milne`, `CONDA_ENV=milne`
(override with `sbatch --export=ALL,PYMILNE_DIR=...`; uncomment the
reservation line if needed). Each task writes `<out>/<cube tag>/` with
the scans and the quicklook overview + movie (`--per-cube-dirs
--quicklook`), and logs `me_<job>_<task>.out` with pyMilne's progress
counter filtered out. Throttle with `--array=0-K%N`; reruns skip finished
scans.

## Run

```bash
conda activate Milne
python me_invert.py nb_6302_2022-05-16T08:28:21_..._im.fits --out me_AR13010
python me_quicklook.py me_AR13010 --movie            # + --scan N, --all-overviews
```

Useful options of `me_invert.py`: `--scans 0-99`, `--crop x0 x1 y0 y1`,
`--bin 2`, `--nthreads 8`, `--no-profile` (skip the CRISP transmission
profile: faster, less exact), `--overwrite`. Existing scans are skipped,
so an interrupted run resumes.

What it does per scan: read (memory-mapped) → leave out tunings without
data in most of the field, mark padding pixels (NaN/zero) invalid →
normalise to the quiet-Sun intensity at the far-wing point (farthest
tuning with data) → noise in Q, U, V measured in the quiet Sun → observed
points on a regular 10 mÅ grid convolved with the CRISP dual-FPI profile
(as in pyMilne's `example_CRISP`) → inversion in row blocks (first scan:
two passes, 4 randomisations; later scans start from the previous scan;
padding pixels get a quiet-Sun profile) → re-inversion of outlier pixels
→ padding written as NaN.

Output `me_scanNNNN.fits`: 15 planes — B [G], inc, azi [rad], vlos
[km/s], vdop [Å], eta0, damp, S0, S1, chi2, I_wing (normalised far-wing
intensity), B_los, B_trans [G], inc, azi [deg] (plane names in `PAR00`…).
Header: `DATE-OBS`, `SCAN`, `MU`, `XCEN`/`YCEN`, `PIXSCALE`, `BIN`, `CROP`,
`ICQS` (quiet-Sun normalisation), `WINGWAV` (far-wing tuning used),
`NTUNE` (tunings used), `NVALID` (pixels with data), `SIGI/Q/U/V`,
`PROFILE`, `NITER`, `PASSES`, `NREPAIR` (outliers re-inverted), `VQS`
(quiet-Sun median vlos, not subtracted).

Speed (M-series, 8 threads, CRISP profile on): ~3400 px/s for scans
started from the previous one, i.e. ~5 min per scan for a 1000×1000 FOV;
the first scan takes ~4x longer.

## Caveats

- v_los has no absolute zero point; plots show it relative to the
  quiet-Sun median of each scan.
- Azimuth and image axes are in the cube's image frame.
- Where linear polarisation is weak compared with the noise, ME can fall
  into a longitudinal-field solution (inclination 0/180°, B_trans ≈ 0,
  small Doppler width) in coherent patches; the outlier pass only fixes
  isolated pixels. pyMilne's spatially regularised inversion is the
  remedy if these patches matter.
- Developed on SSTRED-layout test cubes from the prepared AR13010 maps
  (`tools/make_test_cube.py`; 2022-05-16 08:28:21, frames 0/11/19;
  11 tunings −280…+160 mÅ) and run on the original SSTRED cubes (AR13010
  2022-05-16, 2025-06-17 and others; see `HANDOVER.md`).
- With a tuning left out, the far-wing reference sits deeper in the line
  wing, so the measured noise (`SIGQ/U/V`) comes out higher.

## Files

| file | |
|---|---|
| `me_invert.py` | ME inversion of an SSTRED cube, one FITS per scan |
| `me_quicklook.py` | overview plots and movies of the results |
| `sstred.py` | SSTRED cube reader (wavelength, time, pointing, mu) |
| `pymilne_build/setup_macos_arm64.py` | builds pyMilne on Apple silicon |
| `pymilne_build/build_linux.sh` | checks dependencies and builds pyMilne on Linux |
| `tools/make_test_cube.py` | SSTRED-layout test cube from the prepared AR13010 maps |
| `slurm/me_invert_array.sbatch` | job array, one cube per node |
| `slurm/build_pymilne.sbatch` | builds pyMilne on a compute node |
