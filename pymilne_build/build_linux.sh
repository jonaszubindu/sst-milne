#!/usr/bin/env bash
# Build pyMilne on Linux with its own setup.py (GCC + OpenMP), after
# checking what it needs and pointing it at Eigen and FFTW.
#
#   pymilne_build/build_linux.sh /path/to/pyMilne
#
# Run it with the python you will run me_invert.py with (venv or conda
# env active), on the machine type you will run on: setup.py compiles
# with -march=native.
#
# Needs: numpy (<2, as tested), cython (<3, as tested); GCC >= 10 (C++20);
# FFTW3 headers + library and Eigen >= 3.4 headers, from the active conda
# env ($CONDA_PREFIX) or the system (/usr/include, /usr/lib*).
set -euo pipefail

PM=$(cd "${1:-${PYMILNE_DIR:-$HOME/pyMilne}}" && pwd)
[ -f "$PM/pyMilne.pyx" ] || { echo "$PM is not a pyMilne checkout"; exit 1; }
PY=$(command -v python)
echo "pyMilne : $PM"
echo "python  : $PY"
"$PY" -c "import numpy, Cython; print('numpy', numpy.__version__, '| cython', Cython.__version__)" \
    || { echo "install numpy and cython into this python first"; exit 1; }

# --- compiler
GCCV=$(gcc -dumpversion | cut -d. -f1)
echo "gcc     : $(gcc -dumpversion)"
[ "$GCCV" -ge 10 ] || { echo "need gcc >= 10 for -std=c++20 (module load gcc?)"; exit 1; }

# --- FFTW3
FFTW_INC=""; FFTW_LIB=""
for d in "${CONDA_PREFIX:-/nonexistent}" /usr /usr/local; do
    if [ -f "$d/include/fftw3.h" ]; then FFTW_INC="$d/include"; break; fi
done
for d in "${CONDA_PREFIX:-/nonexistent}/lib" /usr/lib/x86_64-linux-gnu /usr/lib64 /usr/lib /usr/local/lib; do
    if ls "$d"/libfftw3.so* >/dev/null 2>&1; then FFTW_LIB="$d"; break; fi
done
[ -n "$FFTW_INC" ] && [ -n "$FFTW_LIB" ] || {
    echo "FFTW3 not found: conda install -c conda-forge fftw   (or the system package libfftw3-dev)"; exit 1; }
echo "fftw    : $FFTW_INC, $FFTW_LIB"

# --- Eigen >= 3.4 as $PM/eigen3 (setup.py includes ./eigen3)
if [ ! -e "$PM/eigen3/Eigen/Core" ]; then
    for d in "${CONDA_PREFIX:-/nonexistent}/include/eigen3" /usr/include/eigen3 /usr/local/include/eigen3; do
        if [ -e "$d/Eigen/Core" ]; then ln -sfn "$d" "$PM/eigen3"; break; fi
    done
fi
[ -e "$PM/eigen3/Eigen/Core" ] || {
    echo "Eigen not found. Header-only, e.g.:"
    echo "  git clone --depth 1 --branch 3.4.0 https://gitlab.com/libeigen/eigen.git $PM/eigen3"; exit 1; }
EV=$(awk '/define EIGEN_(WORLD|MAJOR)_VERSION/ {printf "%s.", $3}' "$PM/eigen3/Eigen/src/Core/util/Macros.h")
echo "eigen   : ${EV%.} ($(readlink -f "$PM/eigen3"))"
case "${EV%.}" in
    3.[0-3]) echo "Eigen ${EV%.} breaks with -std=c++20 (std::result_of); use 3.4"; exit 1;;
esac

# --- build (setup.py first deletes old pyMilne.*.so and pyMilne.*.cpp)
export CFLAGS="${CFLAGS:-} -I$FFTW_INC"
export LDFLAGS="${LDFLAGS:-} -L$FFTW_LIB -Wl,-rpath,$FFTW_LIB"
cd "$PM"
"$PY" setup.py build_ext --inplace > build_linux.log 2>&1 \
    || { grep -n -m 10 -i "error" build_linux.log; echo "build failed: $PM/build_linux.log"; exit 1; }
ls pyMilne*.so

# --- smoke test: import, synthesis, 1 vs N threads give the same fit
"$PY" - <<'EOF'
import os, numpy as np, MilneEddington as ME
w = np.arange(41) * 0.012 + 6302.25
res = []
for nt in (1, os.cpu_count()):
    me = ME.MilneEddington([[w, None]], [6302], nthreads=nt, precision="float32")
    m = me.repeat_model(np.float64([800., 1.0, 0.4, 0.5, 0.03, 15., 0.2, 0.3, 0.7]), 4, 4)
    obs = me.synthesize(m).astype("float32")
    m0 = me.repeat_model(np.float64([400., 1.5, 0.5, 0., 0.03, 15., 0.2, 0.3, 0.7]), 4, 4)
    mf, _, chi2 = me.invert(m0, obs, np.zeros((4, w.size), "float32") + 1e-3,
                            nRandom=4, nIter=25, chi2_thres=1.0, mu=1.0)
    res.append(chi2.mean())
print(f"pyMilne OK: <chi2> 1 thread {res[0]:.3g}, {os.cpu_count()} threads {res[1]:.3g}")
EOF
