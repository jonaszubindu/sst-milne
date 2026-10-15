#!/usr/bin/env python
"""Milne-Eddington inversion of an SSTRED CRISP Fe I 6302 cube, scan by scan.

Uses pyMilne (J. de la Cruz Rodriguez 2019, A&A 631, A153) with the CRISP
dual-FPI transmission profile (crisp.py from the pyMilne repository), as in
pyMilne's example_CRISP/invert_crisp.py. Run in the conda env that holds
pyMilne (here: `Milne`).

  python me_invert.py nb_6302_..._im.fits --out me_AR13010 \
      [--scans 0-99] [--crop x0 x1 y0 y1] [--bin 2] [--nthreads 8]
  python me_invert.py /data/2022-05-16/CRISP/ /data/2022-05-17/ --out me_all
      # every nb_6302*_im.fits in the folders (recursive), one sub-folder
      # per cube; --single-pass inverts each scan exactly once

Per scan:
  1. read (memory-mapped), crop, bin;
  2. normalise to the quiet-Sun intensity at the wavelength farthest from
     the line core(s) (quiet Sun = low |V|/I, near-median intensity);
  3. noise in Q, U, V from the quiet Sun at that wavelength (robust);
  4. place the observed points on a regular fine grid (one region per line,
     unobserved points get sigma = 1e32) and convolve with the CRISP
     profile;
  5. invert in row blocks (memory). The first scan starts from a generic
     model with --nrandom-first randomisations, then re-inverts from its
     own median-filtered result (Jaime's two-pass recipe); every later
     scan starts from the previous scan's filtered result with --nrandom;
  6. re-invert outlier pixels (chi2 > 3x median or inclination > 30 deg
     off its 5x5 median) from their neighbourhood, keeping improvements.
Writes <out>/me_scanNNNN.fits: (15, ny, nx) float32 with the 9 ME
parameters, chi2, the normalised far-wing intensity, B_los, B_trans and
the inclination/azimuth in degrees; header: time, mu, pointing, plate
scale, normalisation, noise, settings. Existing scans are skipped
(resume) unless --overwrite.
"""

import argparse
import json
import os
import pathlib
import sys
import time

import numpy as np
from astropy.io import fits
from scipy import ndimage

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from sstred import SSTRedCube  # noqa: E402

PYMILNE_DIR = os.environ.get("PYMILNE_DIR",
                             str(pathlib.Path.home() / "Documents/GitHub/pyMilne"))
LINES = {6301: 6301.4995, 6302: 6302.4931}       # pyMilne's line list
PARAMS = ["B", "inc", "azi", "vlos", "vdop", "eta0", "damp", "S0", "S1"]
EXTRA = ["chi2", "I_wing", "Blos", "Btrans", "inc_deg", "azi_deg"]
M0 = np.float64([400.0, 1.5, 0.5, 0.0, 0.03, 15.0, 0.2, 0.3, 0.7])


def parse_scans(spec, n):
    if spec in (None, "all"):
        return list(range(n))
    out = []
    for part in spec.split(","):
        a, _, b = part.partition("-")
        out += list(range(int(a), int(b or a) + 1))
    return [k for k in out if 0 <= k < n]


def line_segments(wav):
    """Split the tuning points into per-line regions; returns
    [(pyMilne line id, indices into wav), ...]."""
    segs = []
    for lid, lam in LINES.items():
        idx = np.flatnonzero(np.abs(wav - lam) < 0.6)
        if idx.size:
            segs.append((lid, idx))
    if not segs:
        raise SystemExit(f"no Fe I 6301/6302 tuning points in {wav}")
    return segs


def fine_grid(w, dw, extra=5):
    """Regular grid of step dw containing every point of w (+extra margin);
    positions of w in it (Jaime's findgrid)."""
    nw = np.rint((w - w[0]) / dw).astype(int)
    if np.max(np.abs(w - (w[0] + nw * dw))) > 0.15 * dw:
        raise SystemExit(f"tuning points {w} are not on a {dw} A grid; "
                         "choose --dw")
    grid = w[0] + (np.arange(nw[-1] + 1 + 2 * extra) - extra) * dw
    return grid, nw + extra


def build_regions(wav, dw, use_profile):
    import crisp
    regions, pos, lines, off = [], [], [], 0
    for lid, idx in line_segments(wav):
        grid, p = fine_grid(wav[idx], dw)
        tr = None
        if use_profile:
            n = int(round(0.3 / dw))
            tw = (np.arange(2 * n + 1) - n) * dw
            tr = crisp.crisp(6302.0).dual_fpi(tw, erh=-0.001)
            tr = tr / tr.sum()
        regions.append([grid, tr])
        pos.append((idx, p + off))
        lines.append(lid)
        off += grid.size
    return regions, pos, lines, off


def tuning_good(d):
    """(ny, nx, ntune): data present at that pixel and tuning -- all four
    Stokes finite and I > 0. SSTRED pads the rotated/aligned field with NaN
    (or zeros) and fills a tuning that was not recorded the same way."""
    with np.errstate(invalid="ignore"):
        return np.isfinite(d).all(axis=2) & (np.nan_to_num(d[:, :, 0]) > 0)


def quiet_sun(I, V, cont, valid):
    """Mask of quiet-Sun pixels: weak polarisation, near-median intensity
    (valid pixels only)."""
    ic = I[..., cont]
    pol = np.zeros(ic.shape)
    pol[valid] = np.abs(V[valid]).max(axis=-1) / I[valid].max(axis=-1)
    med = np.median(ic[valid])
    near = valid & (np.abs(ic / med - 1) < 0.1)
    qs = near & (pol < np.percentile(pol[valid], 30))
    return qs if qs.sum() > 200 else near


def fin(x):
    """FITS headers take no NaN: non-finite floats become the string 'NaN'."""
    return float(x) if np.isfinite(x) else "NaN"


def invert_blocks(me, model, obs, sig, mu, nrandom, niter, rows):
    out_m = np.empty_like(model)
    out_chi = np.empty(obs.shape[:2], dtype=np.float64)
    for y0 in range(0, obs.shape[0], rows):
        y1 = min(y0 + rows, obs.shape[0])
        m, _, chi2 = me.invert(np.ascontiguousarray(model[y0:y1]),
                               np.ascontiguousarray(obs[y0:y1]), sig,
                               nRandom=nrandom, nIter=niter, chi2_thres=1.0,
                               mu=mu)
        out_m[y0:y1], out_chi[y0:y1] = m, chi2
    return out_m, out_chi


def repair_outliers(me, m, chi2, obs, sig, mu, nrandom, niter, valid):
    """Re-invert pixels stuck in a wrong minimum (chi2 > 3x the median, or
    inclination > 30 deg off its 5x5 median) from the median-filtered
    neighbourhood; keep the new fit only where chi2 improves."""
    inc_med = ndimage.median_filter(m[..., 1], 5)
    bad = valid & ((chi2 > 3 * np.median(chi2[valid]))
                   | (np.abs(m[..., 1] - inc_med) > np.radians(30)))
    n = int(bad.sum())
    if n == 0:
        return m, chi2, 0
    init = smooth_model(m, size=5)[bad][None]
    mb, _, cb = me.invert(np.ascontiguousarray(init),
                          np.ascontiguousarray(obs[bad][None]), sig,
                          nRandom=nrandom, nIter=niter, chi2_thres=1.0, mu=mu)
    better = cb[0] < chi2[bad]
    m, chi2 = m.copy(), chi2.copy()
    idx = np.flatnonzero(bad.ravel())[better]
    m.reshape(-1, m.shape[-1])[idx] = mb[0][better]
    chi2.ravel()[idx] = cb[0][better]
    return m, chi2, n


def smooth_model(m, size=3):
    out = m.copy()
    for k in range(m.shape[-1]):
        if k == 2:      # azimuth wraps at pi: filter cos/sin of 2*azi
            c = ndimage.median_filter(np.cos(2 * m[..., k]), size)
            s = ndimage.median_filter(np.sin(2 * m[..., k]), size)
            out[..., k] = (0.5 * np.arctan2(s, c)) % np.pi
        else:
            out[..., k] = ndimage.median_filter(m[..., k], size)
    return out


def find_cubes(inputs, pattern):
    """Cube files from a mix of file and folder arguments (folders: pattern,
    recursively), sorted, without duplicates."""
    found = []
    for x in inputs:
        p = pathlib.Path(x)
        found += sorted(p.rglob(pattern)) if p.is_dir() else [p]
    seen, out = set(), []
    for p in found:
        if p.resolve() not in seen:
            seen.add(p.resolve())
            out.append(p)
    return out


def cube_tag(path):
    """Output folder name for a cube: 'nb_6302_2022-05-16T08:28:21_...'
    -> '6302_2022-05-16T08-28-21'. A cube whose scan range does not start
    at 0 ('...=40-79_...', the later part of an observation SSTRED split
    in two with the same start time) gets the range appended,
    '..._s40-79', so the parts do not share a folder; tags of cubes
    starting at scan 0 are unchanged (existing results stay valid)."""
    import re
    parts = path.name.split("_")
    tag = "_".join(parts[1:3]) if len(parts) > 2 else path.stem
    m = re.search(r"=(\d+)-(\d+)", path.name)
    if m and int(m.group(1)) != 0:
        tag += f"_s{m.group(1)}-{m.group(2)}"
    return tag.replace(":", "-")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+",
                    help="SSTRED cube(s) and/or folders holding them")
    ap.add_argument("--pattern", default="nb_6302*_im.fits",
                    help="cube file pattern inside folders")
    ap.add_argument("--out", required=True,
                    help="output folder (one sub-folder per cube for folders / several cubes)")
    ap.add_argument("--per-cube-dirs", action="store_true",
                    help="one sub-folder per cube in --out even for a single cube "
                         "(Slurm array tasks)")
    ap.add_argument("--quicklook", action="store_true",
                    help="after each cube: overview of the middle scan and the movie "
                         "(me_quicklook.py)")
    ap.add_argument("--single-pass", action="store_true",
                    help="invert every scan once: no second pass on the first scan, "
                         "no outlier re-inversion (faster, more outlier pixels)")
    ap.add_argument("--scans", default="all", help="e.g. 0-99 or 0,5,10-20")
    ap.add_argument("--crop", type=int, nargs=4, metavar=("X0", "X1", "Y0", "Y1"))
    ap.add_argument("--bin", type=int, default=1)
    ap.add_argument("--nthreads", type=int, default=os.cpu_count())
    ap.add_argument("--nrandom-first", type=int, default=4)
    ap.add_argument("--nrandom", type=int, default=2,
                    help="randomisations for scans initialised from the previous one")
    ap.add_argument("--niter", type=int, default=25)
    ap.add_argument("--dw", type=float, default=0.01, help="fine grid step [A]")
    ap.add_argument("--no-profile", action="store_true",
                    help="skip the CRISP transmission profile (much faster, less exact)")
    ap.add_argument("--sigma-i", type=float, default=5e-3,
                    help="Stokes I noise (model error dominates; QUV are measured)")
    ap.add_argument("--rows", type=int, default=0,
                    help="row block size (0: ~150k pixels per block)")
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()

    # pyMilne ships crisp.py (CRISP transmission profile) only in example_CRISP/
    sys.path[:0] = [PYMILNE_DIR, os.path.join(PYMILNE_DIR, "example_CRISP")]
    try:
        import MilneEddington as ME
        if not args.no_profile:
            import crisp  # noqa: F401
    except ImportError as e:
        raise SystemExit(f"{e}: set PYMILNE_DIR to the pyMilne checkout "
                         f"(now {PYMILNE_DIR}); crisp.py is expected in it or "
                         "in its example_CRISP/") from e

    cubes = find_cubes(args.inputs, args.pattern)
    if not cubes:
        raise SystemExit(f"no cubes matching {args.pattern} in {args.inputs}")
    many = (args.per_cube_dirs or len(cubes) > 1
            or any(pathlib.Path(x).is_dir() for x in args.inputs))
    print(f"{len(cubes)} cube(s); results in {args.out}"
          + ("/<cube>/" if many else ""), flush=True)
    failed = []
    for n, path in enumerate(cubes, 1):
        out = pathlib.Path(args.out) / cube_tag(path) if many else pathlib.Path(args.out)
        print(f"\n=== [{n}/{len(cubes)}] {path.name} -> {out}", flush=True)
        try:
            invert_cube(path, out, args, ME)
        except Exception as e:                    # keep going with the others
            import traceback
            traceback.print_exc()
            failed.append((path.name, f"{type(e).__name__}: {e}"))
    if failed:
        print("\nFAILED:")
        for name, msg in failed:
            print(f"  {name}: {msg}")
        sys.exit(1)


def invert_cube(path, out, args, ME):
    cube = SSTRedCube(path)
    out = pathlib.Path(out)
    out.mkdir(parents=True, exist_ok=True)
    scans = parse_scans(args.scans, cube.nscan)
    wav = cube.wavelengths_A(scans[0])
    segs = line_segments(wav)
    dist_line = np.min([np.abs(wav - LINES[l]) for l, _ in segs], axis=0)
    cont = int(np.argmax(dist_line))
    regions, pos, lines, nfine = build_regions(wav, args.dw, not args.no_profile)
    me = ME.MilneEddington(regions, lines, nthreads=args.nthreads,
                           precision="float32")
    x_c, y_c, scale = cube.pointing(scans[0])
    print(f"{cube.path}\n  {cube.nscan} scans, {cube.ntune} tunings "
          f"({wav.min():.3f}-{wav.max():.3f} A), FOV {cube.nx}x{cube.ny} px at "
          f"{scale:.4f}\"/px; lines {lines}; far-wing point {wav[cont]:.3f} A; "
          f"fine grid {nfine} points; CRISP profile "
          f"{'off' if args.no_profile else 'on'}", flush=True)

    prev = None
    # scans skipped for missing data, with the reason -- tools/me_status.py
    # counts them as done instead of listing the cube for a rerun
    skip_file = out / "skipped_scans.json"
    skipped = json.loads(skip_file.read_text()) if skip_file.exists() else {}
    for k in scans:
        fout = out / f"me_scan{k:04d}.fits"
        if fout.exists() and not args.overwrite:
            src = fits.getheader(fout).get("CUBE", "")
            if src and src != os.path.basename(cube.path):
                # written by another cube that shared this folder (split
                # observations before cube_tag told the parts apart): redo
                print(f"scan {k}: {fout.name} belongs to {src} -- redone for "
                      "this cube", flush=True)
            else:
                # resume: the next scan starts from this scan's saved result
                prev = np.ascontiguousarray(
                    fits.getdata(fout)[:len(PARAMS)].transpose(1, 2, 0),
                    dtype=np.float64)
                print(f"scan {k}: exists, skipped", flush=True)
                continue
        t0 = time.time()
        d = cube.scan(k, crop=args.crop, binning=args.bin)   # (ny, nx, 4, nw)
        ny, nx = d.shape[:2]
        # tunings without data in most of the field (not recorded / empty)
        # are left out of this scan, as if not observed; pixels are judged
        # on the others
        good = tuning_good(d)
        frac = good.mean(axis=(0, 1))
        use = frac >= 0.5 * frac.max() if frac.max() > 0 else np.zeros(len(frac), bool)
        if (~use).any() and use.any():
            print(f"scan {k}: tunings {', '.join(f'{w:.3f}' for w in wav[~use])} A have "
                  "no data in most of the field -- left out", flush=True)
        valid = good[..., use].all(axis=-1) if use.any() else np.zeros((ny, nx), bool)
        if valid.sum() < 1000:
            with np.errstate(invalid="ignore"):
                fin_i = float(np.isfinite(d[:, :, 0]).all(axis=-1).mean())
                fin_quv = float(np.isfinite(d[:, :, 1:]).all(axis=(2, 3)).mean())
            reason = ("no Q/U/V (Stokes I only)" if fin_i > 0 and fin_quv == 0 else
                      f"only {valid.sum()} valid pixels")
            print(f"scan {k}: {reason}, skipped; data fraction per tuning: "
                  + " ".join(f"{w:.3f}:{f:.2f}" for w, f in zip(wav, frac))
                  + f"; I finite in {fin_i:.2f}, Q/U/V in {fin_quv:.2f} of the field",
                  flush=True)
            skipped[str(k)] = reason
            skip_file.write_text(json.dumps(skipped, indent=1))
            continue
        cont = int(np.flatnonzero(use)[np.argmax(dist_line[use])])   # far-wing point
        d[:, :, :, ~use] = 0.0
        d[~valid] = 0.0
        qs = quiet_sun(d[:, :, 0], d[:, :, 3], cont, valid)
        icq = float(np.median(d[:, :, 0, cont][qs]))
        d /= icq
        noise = [1.4826 * float(np.median(np.abs(d[:, :, s, cont][qs]
                                                  - np.median(d[:, :, s, cont][qs]))))
                 for s in (1, 2, 3)]
        # padding gets a quiet-Sun profile, so the fits there stay tame;
        # its results are written as NaN
        d[~valid] = d[qs].mean(axis=0)
        obs = np.zeros((ny, nx, 4, nfine), dtype=np.float32)
        sig = np.full((4, nfine), 1e32, dtype=np.float32)
        for idx, p in pos:
            obs[:, :, :, p] = d[:, :, :, idx]
            pu = p[use[idx]]                      # only tunings with data
            sig[0, pu] = args.sigma_i
            for s in (1, 2, 3):
                sig[s, pu] = max(noise[s - 1], 1e-4)
        mu = cube.mu(k)
        rows = args.rows or max(1, 150_000 // nx)
        if prev is None or prev.shape[:2] != (ny, nx):
            m0 = me.repeat_model(M0, ny, nx)
            m, chi2 = invert_blocks(me, m0, obs, sig, mu, args.nrandom_first,
                                    args.niter, rows)
            if not args.single_pass:
                m, chi2 = invert_blocks(me, smooth_model(m), obs, sig, mu,
                                        args.nrandom_first, args.niter, rows)
        else:
            m, chi2 = invert_blocks(me, smooth_model(prev), obs, sig, mu,
                                    args.nrandom, args.niter, rows)
        nfix = 0
        if not args.single_pass:
            m, chi2, nfix = repair_outliers(me, m, chi2, obs, sig, mu,
                                            max(4, args.nrandom), args.niter, valid)
        prev = m

        B, inc, azi = m[..., 0], m[..., 1], m[..., 2]
        stack = np.concatenate([m.transpose(2, 0, 1),
                                np.stack([chi2, d[:, :, 0, cont],
                                          B * np.cos(inc), B * np.sin(inc),
                                          np.degrees(inc), np.degrees(azi)])]
                               ).astype(np.float32)
        stack[:, ~valid] = np.nan
        h = fits.Header()
        for i, nm in enumerate(PARAMS + EXTRA):
            h[f"PAR{i:02d}"] = nm
        h["UNITS"] = "G, rad, rad, km/s, A; Blos/Btrans G; *_deg deg"
        h["CUBE"] = os.path.basename(cube.path)
        h["SCAN"] = k
        h["DATE-OBS"] = cube.time(k).isoformat()
        h["MU"] = (fin(mu), "at the field centre")
        h["XCEN"], h["YCEN"] = fin(x_c), fin(y_c)
        h["PIXSCALE"] = (scale * args.bin, "arcsec/px after binning")
        h["BIN"] = args.bin
        if args.crop:
            h["CROP"] = " ".join(map(str, args.crop))
        h["ICQS"] = (fin(icq), "quiet-Sun far-wing I [cube units]")
        h["NVALID"] = (int(valid.sum()), "pixels with data; rest NaN")
        h["WINGWAV"] = (float(wav[cont]), "far-wing wavelength [A]")
        h["NTUNE"] = (int(use.sum()), "tunings used (rest: no data, left out)")
        h["SIGI"], h["SIGQ"], h["SIGU"], h["SIGV"] = args.sigma_i, *map(fin, noise)
        h["PROFILE"] = "none" if args.no_profile else "crisp dual_fpi erh=-0.001"
        h["NITER"] = args.niter
        h["PASSES"] = ("single" if args.single_pass else "two-pass first scan + repair")
        h["NREPAIR"] = (nfix, "outlier pixels re-inverted")
        h["VQS"] = (fin(np.median(m[..., 3][qs])),
                    "QS median vlos [km/s], not subtracted")
        fits.PrimaryHDU(stack, h).writeto(fout, overwrite=True)
        print(f"scan {k} ({cube.time(k):%H:%M:%S}, mu {mu:.2f}): {ny}x{nx} px "
              f"in {time.time() - t0:.0f} s ({(~valid).mean() * 100:.0f} % padding), "
              f"<chi2> {np.median(chi2[valid]):.2f}, "
              f"{nfix} px re-inverted, "
              f"noise QUV {noise[0]:.1e} {noise[1]:.1e} {noise[2]:.1e} -> {fout.name}",
              flush=True)
    cube.close()
    if args.quicklook:
        import argparse as _ap
        import me_quicklook
        try:
            me_quicklook.run_dir(str(out), _ap.Namespace(
                scan=None, all_overviews=False, movie=True, fps=6))
        except (Exception, SystemExit) as e:      # never fail the inversion
            print(f"quicklook failed ({type(e).__name__}: {e}); run me_quicklook.py "
                  "later", flush=True)


if __name__ == "__main__":
    main()
