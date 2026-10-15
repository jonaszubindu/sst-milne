#!/usr/bin/env python
"""Which cubes of a list are inverted, and which have quicklooks?

  python tools/me_status.py cubes.txt /path/to/me_out [/path/to/me_out_17 ...]
  python tools/me_status.py --array-only cubes.txt /path/to/me_out_17
      # just the comma list of incomplete cubes (empty if none), for sbatch

For every cube in the list (one path per line, as given to the Slurm
array; relative paths are resolved from the current directory) it looks
up the cube's folder (<out>/<cube tag>/) in each output folder and
prints: scans in the cube, scans inverted, overview and movie present and
up to date (newer than the newest scan). At the end it prints the
`--array` list of cubes without a complete inversion (line numbers of
the list, for sbatch) and the cube folders without an up-to-date movie.

Scans that me_invert.py skipped for missing data (e.g. no Q/U/V) are
recorded in <cube folder>/skipped_scans.json; they count as done (column
"skip", reasons listed) and are not offered for a rerun. Cubes inverted
before that record existed show such scans as missing until rerun once.
"""

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from me_invert import cube_tag  # noqa: E402
from sstred import SSTRedCube   # noqa: E402


def main():
    argv = sys.argv[1:]
    array_only = "--array-only" in argv
    argv = [a for a in argv if a != "--array-only"]
    if len(argv) < 2:
        sys.exit(__doc__)
    if array_only:                      # silence the table, keep the list
        import io
        import contextlib
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            inc = report(pathlib.Path(argv[0]), [pathlib.Path(p) for p in argv[1:]])
        print(",".join(map(str, inc)))
        return
    report(pathlib.Path(argv[0]), [pathlib.Path(p) for p in argv[1:]])


def report(lst, outs):
    """Print the table; return the list indices of incomplete cubes."""
    cubes = [l.strip() for l in lst.read_text().splitlines() if l.strip()]
    tags = [cube_tag(pathlib.Path(c)) for c in cubes]
    for t in sorted({t for t in tags if tags.count(t) > 1}):
        print(f"WARNING: cubes {[i for i, x in enumerate(tags) if x == t]} share the "
              f"output folder {t} -- update sst-milne (cube_tag) before running them")
    incomplete, no_movie, unknown = [], [], []
    print(f"{'#':>3s}  {'cube tag':28s} {'scans':>5s} {'done':>5s} {'skip':>4s}  overview  movie   folder")
    import json
    for i, c in enumerate(cubes):
        path = pathlib.Path(c)
        tag = cube_tag(path)
        try:
            cube = SSTRedCube(path)
            nscan = cube.nscan
            cube.close()
        except Exception as e:
            nscan = None
            unknown.append(i)
            print(f"{i:3d}  {tag:28s}   cannot open the cube ({type(e).__name__}) -- "
                  "scan count unknown; run from the folder the list's paths refer to")
        ns = "?" if nscan is None else str(nscan)
        hits = [o / tag for o in outs if (o / tag).is_dir()]
        if not hits:
            print(f"{i:3d}  {tag:28s} {ns:>5s} {0:5d} {0:4d}  -         -       (no results)")
            incomplete.append(i)
            continue
        best = 0
        reasons = {}
        for d in hits:
            files = sorted(d.glob("me_scan*.fits"))
            sk = {}
            if (d / "skipped_scans.json").exists():
                sk = json.loads((d / "skipped_scans.json").read_text())
                done_ids = {int(f.stem[7:]) for f in files}
                sk = {k: v for k, v in sk.items() if int(k) not in done_ids}
            for v in sk.values():
                reasons[v] = reasons.get(v, 0) + 1
            best = max(best, len(files) + len(sk))
            newest = max((f.stat().st_mtime for f in files), default=0)
            ql = d / "quicklook"
            ov = sorted(ql.glob("overview_scan*.png"))
            mv = ql / "movie_me.mp4"
            ov_ok = bool(ov) and all(p.stat().st_mtime > newest for p in ov)
            mv_ok = mv.exists() and mv.stat().st_mtime > newest
            state = lambda ok, exists: "yes" if ok else ("old" if exists else "-")
            print(f"{i:3d}  {tag:28s} {ns:>5s} {len(files):5d} {len(sk):4d}  "
                  f"{state(ov_ok, bool(ov)):8s}  {state(mv_ok, mv.exists()):6s}  {d.parent}")
            if files and not mv_ok:
                no_movie.append(str(d))
        if reasons:
            print(f"{'':5s}skipped scans: " + "; ".join(f"{n}x {r}" for r, n in reasons.items()))
        # complete if any output folder holds all scans of the cube (inverted
        # or recorded as skipped for missing data)
        if nscan is not None and best < nscan:
            incomplete.append(i)
    inc = sorted(set(incomplete))
    done = len(cubes) - len(set(inc) | set(unknown))
    print(f"\n{done}/{len(cubes)} cubes completely inverted"
          + (f"; {len(unknown)} could not be opened (lines {','.join(map(str, unknown))})"
             if unknown else ""))
    if inc:
        print(f"incomplete (resubmit with the same list and output folder):\n"
              f"  --array={','.join(map(str, inc))}")
    if no_movie:
        print(f"{len(no_movie)} cube folder(s) without an up-to-date movie, e.g.:\n"
              f"  python me_quicklook.py {' '.join(sorted({str(pathlib.Path(p).parent) for p in no_movie}))} "
              "--movie --skip-existing")
    return inc


if __name__ == "__main__":
    main()
