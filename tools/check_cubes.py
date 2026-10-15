#!/usr/bin/env python
"""Data validity check of SSTRED CRISP cubes (Fe I 6302 and Ca II 8542).

  python tools/check_cubes.py CUBE_OR_FOLDER [...] [--pattern 'nb_*_im.fits']
         [--stride 4] [--workers 4] [--csv DIR] [--plot DIR]

Identical copy in sst-milne/tools and sst-wfa/tools (both repos carry
the same sstred.py). Reads every (scan, Stokes, tuning) map of each cube
(memory-mapped, every --stride-th row and column; --stride 1 reads
everything) and reports, per cube:

  ERROR  EMPTY       map without any data (all NaN / zero; I: no I > 0)
  ERROR  ZERO_POL    Q, U or V exactly 0 in > 50 % of the field where I has data
  ERROR  INF         infinite values
  ERROR  SCAN_SKIP   < 1000 pixels with all Stokes at all used tunings:
                     me_invert.py / wfa_invert.py skip this scan
  ERROR  TIME        scan times not increasing
  WARN   TUNING_DROP tuning with data in < 50 % of the best tuning's field:
                     the inversions leave it out of that scan
  WARN   PARTIAL     map with data in < 90 % of the scan's field (I footprint)
  WARN   NEG_I       I <= 0 (not padding) in > 1 % of the field
  WARN   CONST       map constant over its data
  WARN   I_LEVEL     median I of a scan/tuning > 20 % off that tuning's
                     median over the scans (scaling, clouds, bad scan)
  WARN   NOISE       Q/U/V noise (far-wing tuning) > 3x the cube's median
  WARN   WAVE        tuning wavelengths differ between scans
  INFO   GAP         time gap > 3x the median cadence
  INFO   LINE        which line(s), tunings, field, padding

The SSTRED padding (NaN or zero around the rotated field) is not an error;
the "field" is where Stokes I has data at most tunings. Stokes order is
assumed I, Q, U, V (SSTRED). --csv writes one row per map
(<tag>_maps.csv); --plot a PNG per cube with the data fraction per scan
and tuning for I, Q, U, V. Exit code 1 if any cube has an ERROR.

With --stride N > 1 whole-map and field-wide problems are found reliably,
but isolated bad pixels (a single inf or NaN) can fall between the
sampled rows/columns; --stride 1 for a complete census (reads the whole
cube: ~85 GB for a 2900x2900 px, 50-scan, 12-tuning cube).
"""

import argparse
import csv
import pathlib
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from sstred import SSTRedCube  # noqa: E402

LINES = {"Fe I 6301": 6301.4995, "Fe I 6302": 6302.4931, "Ca II 8542": 8542.091}
STOKES = "IQUV"


def find_cubes(inputs, pattern):
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
    parts = path.name.split("_")
    tag = "_".join(parts[1:3]) if len(parts) > 2 else path.stem
    return tag.replace(":", "-")


def ranges(ks):
    """[0, 1, 2, 5, 7, 8] -> '0-2,5,7-8'"""
    ks, out, i = sorted(ks), [], 0
    while i < len(ks):
        j = i
        while j + 1 < len(ks) and ks[j + 1] == ks[j] + 1:
            j += 1
        out.append(str(ks[i]) if i == j else f"{ks[i]}-{ks[j]}")
        i = j + 1
    return ",".join(out)


def check_scan(data, k, st, nw):
    """Statistics of every map of scan k; returns a dict of (4, nw) arrays
    plus the scan-level counts."""
    shp = (4, nw)
    r = {key: np.zeros(shp) for key in
         ("nvalid", "nnan", "nzero", "ninf", "nneg", "fp_frac", "fp_zero", "med", "mad")}
    r["const"] = np.zeros(shp, bool)
    good = None                     # (nw, ny', nx'): all Stokes usable, as tuning_good
    I_ok = []
    for s in range(4):
        for w in range(nw):
            a = np.asarray(data[k, s, w, ::st, ::st], dtype=np.float32)
            fin = np.isfinite(a)
            r["nnan"][s, w] = np.isnan(a).sum()
            r["ninf"][s, w] = (~fin).sum() - r["nnan"][s, w]
            r["nzero"][s, w] = (a == 0).sum()
            if s == 0:
                ok = fin & (a > 0)
                r["nneg"][s, w] = (fin & (a < 0)).sum()
                I_ok.append(ok)
                if w == nw - 1:     # field: I has data at most tunings
                    cnt = np.sum(I_ok, axis=0)
                    fp = cnt >= max(1, nw // 2)
                    nfp = max(int(fp.sum()), 1)
                    good = np.stack(I_ok)
            else:
                ok = fin & (a != 0)
                good[w] &= fin
            r["nvalid"][s, w] = ok.sum()
            if s > 0:
                r["fp_frac"][s, w] = (ok & fp).sum() / nfp if good is not None else 0
                r["fp_zero"][s, w] = ((a == 0) & fp).sum() / nfp
            v = a[ok]
            if v.size:
                med = float(np.median(v))
                r["med"][s, w] = med
                r["mad"][s, w] = 1.4826 * float(np.median(np.abs(v - med)))
                r["const"][s, w] = v.size > 1 and float(v.max()) == float(v.min())
    for w in range(nw):
        r["fp_frac"][0, w] = (I_ok[w] & fp).sum() / nfp
    r["nfp"] = int(fp.sum())
    # what the inversions do: tunings with < 50 % of the best data fraction
    # are left out, pixels need all Stokes at all remaining tunings
    frac = good.mean(axis=(1, 2))
    use = frac >= 0.5 * frac.max() if frac.max() > 0 else np.zeros(nw, bool)
    r["use"] = use
    r["ninv"] = int(good[use].all(axis=0).sum()) if use.any() else 0
    return r


def check_cube(path, args):
    probs = []                  # (level, code, text)

    def add(level, code, text):
        probs.append((level, code, text))

    try:
        cube = SSTRedCube(path)
    except Exception as e:
        add("ERROR", "READ", f"cannot open as SSTRED cube: {type(e).__name__}: {e}")
        return probs, None
    data, ns, nw = cube.data, cube.nscan, cube.ntune
    st = args.stride
    px = st * st                 # pixels per sampled pixel
    wav = np.array([cube.wavelengths_A(k) for k in range(ns)])
    lines = {n: l for n, l in LINES.items() if np.abs(wav[0] - l).min() < 1.0}
    lam_ref = np.array(list(lines.values())) if lines else np.array([wav[0].mean()])
    far = int(np.argmax(np.min(np.abs(wav[0][:, None] - lam_ref[None]), axis=1)))
    x_c, y_c, scale = cube.pointing(0)
    times = [cube.time(k) for k in range(ns)]
    add("INFO", "LINE", f"{', '.join(lines) or 'no known line'}; {ns} scans "
        f"{times[0]:%Y-%m-%d %H:%M:%S}-{times[-1]:%H:%M:%S}, {nw} tunings "
        f"{wav[0].min():.3f}-{wav[0].max():.3f} A, {cube.nx}x{cube.ny} px at "
        f"{scale:.4f}\"/px, dtype {data.dtype}, BUNIT {cube.header.get('BUNIT', '?')}")
    if not lines:
        add("WARN", "LINE", "no tuning within 1 A of Fe I 6301/6302 or Ca II 8542")
    if ns > 1 and np.abs(wav - wav[0]).max() > 1e-3:
        bad = [k for k in range(ns) if np.abs(wav[k] - wav[0]).max() > 1e-3]
        add("WARN", "WAVE", f"tunings differ from scan 0 in scans {ranges(bad)} "
            f"(up to {np.abs(wav - wav[0]).max() * 1000:.0f} mA)")
    tsec = np.array([(t - times[0]).total_seconds() for t in times])
    if ns > 1:
        dt = np.diff(tsec)
        if (dt <= 0).any():
            add("ERROR", "TIME", f"time not increasing at scan(s) "
                f"{ranges(list(np.flatnonzero(dt <= 0) + 1))}")
        elif ns > 2 and (dt > 3 * np.median(dt)).any():
            gaps = np.flatnonzero(dt > 3 * np.median(dt))
            add("INFO", "GAP", f"cadence {np.median(dt):.0f} s; gaps of "
                + ", ".join(f"{dt[i]:.0f} s before scan {i + 1}" for i in gaps))

    with ThreadPoolExecutor(max(1, args.workers)) as ex:
        res = list(ex.map(lambda k: check_scan(data, k, st, nw), range(ns)))
    cube.close()

    R = {key: np.stack([r[key] for r in res]) for key in res[0]}   # (ns, 4, nw)
    fpfrac = R["nfp"] / (data.shape[3] // st * (data.shape[4] // st) or 1)
    add("INFO", "LINE", f"field (I data at most tunings): "
        f"{np.median(fpfrac) * 100:.0f} % of the frame (rest padding), "
        f"{ranges(list(range(ns)))} scans checked at stride {st}")

    def report(level, code, mask, what, stokes=True):
        """mask (ns, 4, nw) -> grouped lines: Stokes with the same tunings and
        scans together, 'all tunings' when every tuning is affected alike."""
        groups = {}
        for s in range(4):
            per_w = [tuple(np.flatnonzero(mask[:, s, w])) for w in range(nw)]
            by_scans = {}
            for w, ks in enumerate(per_w):
                if ks:
                    by_scans.setdefault(ks, []).append(w)
            for ks, ws in by_scans.items():
                tun = ("all tunings" if len(ws) == nw else
                       f"{len(ws)}/{nw} tunings" if len(ws) > 4 else
                       ", ".join(f"{wav[0][w]:.3f}" for w in ws) + " A")
                groups.setdefault((tun, ks), []).append(STOKES[s])
        for (tun, ks), ss in groups.items():
            who = f"Stokes {''.join(ss)}, {tun}" if stokes else tun
            add(level, code, f"{who}: {what} in scan"
                f"{'s' if len(ks) > 1 else ''} {ranges(ks)}"
                + (f" ({len(ks)}/{ns})" if ns > 1 else ""))

    report("ERROR", "INF", R["ninf"] > 0, "infinite values")
    empty = R["nvalid"] == 0
    tot = R["nnan"] + R["nzero"] + R["ninf"]
    report("ERROR", "EMPTY", empty & (R["nnan"] == tot), "all NaN")
    report("ERROR", "EMPTY", empty & (R["nzero"] == tot), "all 0")
    report("ERROR", "EMPTY", empty & (R["nnan"] < tot) & (R["nzero"] < tot),
           "no data (NaN/0/negative)")
    zero_pol = R["fp_zero"] > 0.5
    zero_pol[:, 0] = False
    report("ERROR", "ZERO_POL", zero_pol & (R["nvalid"] > 0), "exactly 0 in > 50 % of the field")
    part = (R["fp_frac"] < 0.9) & (R["nvalid"] > 0) & ~zero_pol
    report("WARN", "PARTIAL", part, "data in < 90 % of the field")
    neg = np.zeros_like(part)
    neg[:, 0] = R["nneg"][:, 0] > 0.01 * np.maximum(R["nfp"][:, None], 1)
    report("WARN", "NEG_I", neg, "I < 0 in > 1 % of the field")
    report("WARN", "CONST", R["const"], "constant")

    medI = R["med"][:, 0]                                       # (ns, nw)
    ref = np.median(np.where(medI > 0, medI, np.nan), axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        dev = np.abs(medI / ref - 1)
    lvl = np.zeros_like(part)
    lvl[:, 0] = (dev > 0.2) & (medI > 0)
    report("WARN", "I_LEVEL", lvl, "median I > 20 % off the other scans")
    for s in (1, 2, 3):
        n = R["mad"][:, s, far]
        good = n > 0
        if good.sum() > 2:
            ks = np.flatnonzero(good & (n > 3 * np.median(n[good])))
            if ks.size:
                add("WARN", "NOISE", f"Stokes {STOKES[s]} noise at {wav[0][far]:.3f} A > 3x "
                    f"the cube median ({np.median(n[good]):.2e}) in scans {ranges(ks)}")
    skipped = R["ninv"] * px < 1000
    drop = np.zeros_like(part)
    drop[:, 0] = ~R["use"] & ~skipped[:, None]
    report("WARN", "TUNING_DROP", drop, "left out by the inversions (data in "
           "< 50 % of the best tuning's field)", stokes=False)
    skip = np.flatnonzero(skipped)
    if skip.size:
        add("ERROR", "SCAN_SKIP", f"scan{'s' if skip.size > 1 else ''} {ranges(skip)}: "
            "< 1000 usable pixels, the inversions skip "
            f"{'them' if skip.size > 1 else 'it'}")

    if args.csv:
        out = pathlib.Path(args.csv) / f"{cube_tag(path)}_maps.csv"
        with open(out, "w", newline="") as f:
            wr = csv.writer(f)
            wr.writerow(["scan", "stokes", "tuning", "wav_A", "n_valid", "n_nan", "n_zero",
                         "n_inf", "frac_field", "frac_zero_field", "median", "mad", "used"])
            for k in range(ns):
                for s in range(4):
                    for w in range(nw):
                        wr.writerow([k, STOKES[s], w, f"{wav[k][w]:.4f}",
                                     int(R["nvalid"][k, s, w] * px), int(R["nnan"][k, s, w] * px),
                                     int(R["nzero"][k, s, w] * px), int(R["ninf"][k, s, w] * px),
                                     f"{R['fp_frac'][k, s, w]:.4f}",
                                     f"{R['fp_zero'][k, s, w]:.4f}",
                                     f"{R['med'][k, s, w]:.4g}", f"{R['mad'][k, s, w]:.4g}",
                                     int(R["use"][k, w])])
        add("INFO", "CSV", str(out))
    if args.plot:
        png = pathlib.Path(args.plot) / f"{cube_tag(path)}_validity.png"
        plot(R, wav[0], path.name, png)
        add("INFO", "PLOT", str(png))
    return probs, R


def plot(R, wav, name, png):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ns, _, nw = R["fp_frac"].shape
    fig, axs = plt.subplots(1, 4, figsize=(14, 2.5 + 0.08 * ns), sharey=True)
    for s, ax in enumerate(axs):
        f = R["fp_frac"][:, s] * 100
        f = np.where(R["fp_zero"][:, s] > 0.5, -1, f) if s else f
        im = ax.imshow(f, origin="lower", aspect="auto", cmap="RdYlGn", vmin=0, vmax=100,
                       extent=[-0.5, nw - 0.5, -0.5, ns - 0.5], interpolation="nearest")
        ax.set_xticks(range(nw))
        ax.set_xticklabels([f"{w:.2f}" for w in wav], rotation=90, fontsize=7)
        ax.set_title(f"Stokes {STOKES[s]}: % of field with data", fontsize=9)
        ax.set_xlabel("tuning [A]")
    axs[0].set_ylabel("scan")
    fig.colorbar(im, ax=axs, fraction=0.02, pad=0.01)
    fig.suptitle(name, fontsize=9, y=1.04)
    fig.savefig(png, dpi=110, bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="+", help="SSTRED cube(s) and/or folders holding them")
    ap.add_argument("--pattern", default="nb_*_im.fits",
                    help="cube file pattern inside folders (default: 6302 and 8542)")
    ap.add_argument("--stride", type=int, default=4,
                    help="read every N-th row/column (1: every pixel, exact but slow)")
    ap.add_argument("--workers", type=int, default=4, help="scans read in parallel")
    ap.add_argument("--csv", help="folder for per-map statistics (<tag>_maps.csv)")
    ap.add_argument("--plot", help="folder for per-cube validity PNGs")
    ap.add_argument("--quiet", action="store_true", help="only WARN and ERROR lines")
    args = ap.parse_args()
    for d in (args.csv, args.plot):
        if d:
            pathlib.Path(d).mkdir(parents=True, exist_ok=True)

    cubes = find_cubes(args.inputs, args.pattern)
    if not cubes:
        raise SystemExit(f"no cubes matching {args.pattern} in {args.inputs}")
    summary = []
    for n, path in enumerate(cubes, 1):
        print(f"\n=== [{n}/{len(cubes)}] {path}", flush=True)
        try:
            probs, _ = check_cube(path, args)
        except Exception as e:
            probs = [("ERROR", "CHECK", f"check failed: {type(e).__name__}: {e}")]
        for level, code, text in probs:
            if not (args.quiet and level == "INFO"):
                print(f"  {level:5s} {code:11s} {text}", flush=True)
        worst = ("ERROR" if any(p[0] == "ERROR" for p in probs) else
                 "WARN" if any(p[0] == "WARN" for p in probs) else "OK")
        if worst == "OK":
            print("  OK", flush=True)
        summary.append((path, worst, sum(p[0] == "ERROR" for p in probs),
                        sum(p[0] == "WARN" for p in probs)))
    print(f"\n{'status':6s} {'err':>4s} {'warn':>4s}  cube")
    for path, worst, ne, nwn in summary:
        print(f"{worst:6s} {ne:4d} {nwn:4d}  {path}")
    sys.exit(1 if any(s[1] == "ERROR" for s in summary) else 0)


if __name__ == "__main__":
    main()
