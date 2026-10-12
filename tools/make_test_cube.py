"""Small test cube in SSTRED layout from the prepared AR13010 maps (TEST ONLY).

  python tools/make_test_cube.py "<dir with Prepared_map_SIR_AR13010_nb_6302_*_frame_N.fits>" out.fits

Uses frames 0, 11, 19 of 2022-05-16T08:28:21. Those maps hold 12 slots of
the 13-point grid 6302.4931 + (-0.28 ... +0.20) A in 40 mA steps (the grid
in invert_Milne_Daniel.ipynb); slot 1 is a -1 placeholder and is dropped,
leaving 11 tunings. A 200x240 px crop around the darkest area, scaled to
absolute-like units, with WCS-TAB corners at (-519", -227"), 0.059"/px.
"""
import sys, numpy as np
from astropy.io import fits
D = sys.argv[1]; out = sys.argv[2]
frames = [0, 11, 19]
offs = -0.28 + 0.04 * np.arange(12)
keep = [j for j in range(12) if j != 1]                  # slot 1 is the -1 placeholder
wav_nm = (6302.4931 + offs[keep]) / 10.0
maps = [fits.getdata(f"{D}/Prepared_map_SIR_AR13010_nb_6302_2022-05-16T08:28:21_frame_{f}.fits")[:, keep] for f in frames]
ic = maps[0][0, 0]
yy, xx = np.unravel_index(np.argmin(ndi := __import__("scipy.ndimage", fromlist=["x"]).uniform_filter(ic, 15)), ic.shape)
y0, x0 = int(np.clip(yy - 100, 0, ic.shape[0] - 200)), int(np.clip(xx - 120, 0, ic.shape[1] - 240))
cube = np.stack([m[:, :, y0:y0 + 200, x0:x0 + 240] for m in maps]).astype(np.float32) * 3.1e-8   # absolute-like units
ns, nst, nw, ny, nx = cube.shape
scale, xc, yc = 0.059, -519.0, -227.0
t = 8 * 3600 + 28 * 60 + 21 + 20.0 * np.array(frames)
coords = np.zeros((ns, nw, 2, 2, 4))
for k in range(ns):
    for w in range(nw):
        for a in range(2):
            for b in range(2):
                coords[k, w, a, b] = [xc + (b - 0.5) * (nx - 1) * scale, yc + (a - 0.5) * (ny - 1) * scale, wav_nm[w], t[k]]
h = fits.Header()
for k, v in [("CTYPE1", "HPLN-TAB"), ("CTYPE2", "HPLT-TAB"), ("CTYPE3", "WAVE-TAB"), ("CTYPE4", "STOKES"), ("CTYPE5", "UTC--TAB"),
             ("CUNIT3", "nm"), ("CUNIT5", "s"), ("PS3_0", "WCS-TAB"), ("PS3_1", "HPLN+HPLT+WAVE+TIME"), ("PV3_3", 3),
             ("BUNIT", "W m^-2 Hz^-1 sr^-1"), ("DATEREF", "2022-05-16T00:00:00.000000"), ("INSTRUME", "CRISP"),
             ("COMMENT", f"TEST CUBE from prepared maps frames {frames}, crop x{x0} y{y0}")]:
    h[k] = v
prim = fits.PrimaryHDU(cube, h)
col = fits.Column(name="HPLN+HPLT+WAVE+TIME", format=f"{coords.size}D", dim=f"(4,2,2,{nw},{ns})", array=coords[None])
c2 = fits.Column(name="HPLN-INDEX", format="2E", array=np.array([[1, nx]], dtype=np.float32))
c3 = fits.Column(name="HPLT-INDEX", format="2E", array=np.array([[1, ny]], dtype=np.float32))
tab = fits.BinTableHDU.from_columns([col, c2, c3], name="WCS-TAB")
fits.HDUList([prim, tab]).writeto(out, overwrite=True)
print("wrote", out, cube.shape, "crop", x0, y0)
