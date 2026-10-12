# Build pyMilne on Apple-silicon Macs (arm64) with Apple clang + Homebrew
# libomp, without touching the repository sources.
#
#   conda activate Milne
#   python pymilne_build/setup_macos_arm64.py [/path/to/pyMilne]
#
# pyMilne checkout: the argument, else $PYMILNE_DIR, else
# ~/Documents/GitHub/pyMilne (https://github.com/jaimedelacruz/pyMilne).
# Leaves pyMilne.<tag>.so next to MilneEddington.py in that checkout. Needs Homebrew fftw and
# libomp, eigen3 (here: the conda env's include/eigen3, 3.3.7) and Cython.
#
# Differences from setup.py, and why:
#  * no -march/-mtune=native, -mprefer-vector-width (x86 only; setup.py's own
#    comment says to drop them on M-series Macs);
#  * OpenMP the clang way (-Xpreprocessor -fopenmp), against the conda
#    env's libomp if present (it is the runtime that gets loaded);
#  * -std=c++17 instead of c++20: Eigen 3.3.7 uses std::result_of, which
#    libc++ removes in C++20; the pyMilne sources use no C++20 features;
#  * src/profile.hpp calls std::abs<T>(v) (3x), which libc++ rejects as
#    ambiguous (libstdc++ accepts it). The sources are copied to
#    _build_macos_arm64/ and those calls become std::abs(v) -- identical
#    for T = float/double. The repository files stay unchanged;
#  * -include iostream, for Eigen 3.3.7's std::cerr under libc++;
#  * schedule(dynamic,k) -> schedule(static,k) in the build copy (3 loops):
#    clang >= 19 needs __kmpc_dispatch_deinit for dynamic schedules, absent
#    from the env's llvm-openmp 16 runtime that OpenBLAS loads first;
#    likewise -DEIGEN_DONT_PARALLELIZE (Eigen's own dynamic schedule;
#    pyMilne parallelises over pixels itself).
# (Homebrew GCC 14.2 cannot be used on this machine: it was built for an
# older macOS and its include-fixed headers clash with the installed SDK.)
import os
import pathlib
import shutil
import sys

import numpy
from Cython.Distutils import build_ext
from setuptools import setup
from setuptools.extension import Extension

repo = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else
                    os.environ.get("PYMILNE_DIR",
                                   pathlib.Path.home() / "Documents/GitHub/pyMilne")
                    ).expanduser().resolve()
if not (repo / "pyMilne.pyx").exists():
    sys.exit(f"{repo} is not a pyMilne checkout (no pyMilne.pyx)")
work = repo / "_build_macos_arm64"
env_root = pathlib.Path(sys.executable).parents[1]
brew = pathlib.Path("/opt/homebrew/opt")

if work.exists():
    shutil.rmtree(work)
shutil.copytree(repo / "src", work / "src")
shutil.copy2(repo / "pyMilne.pyx", work / "pyMilne.pyx")
prof = work / "src" / "profile.hpp"
txt = prof.read_text()
n = txt.count("std::abs<T>(")
prof.write_text(txt.replace("std::abs<T>(", "std::abs("))
print(f"patched {n} std::abs<T>( calls in the build copy of profile.hpp")
# Apple clang >= 19 lowers schedule(dynamic) to __kmpc_dispatch_deinit, which
# the conda env's llvm-openmp 16 (loaded first, by OpenBLAS) does not have.
# Round-robin static chunks need only __kmpc_for_static_*; the per-pixel
# work is unchanged, only how pixels are dealt to threads.
nd = 0
for f in (work / "src").glob("*.[hc]pp"):
    t = f.read_text()
    if "schedule(dynamic" in t:
        nd += t.count("schedule(dynamic")
        f.write_text(t.replace("schedule(dynamic", "schedule(static"))
print(f"patched {nd} schedule(dynamic,..) -> schedule(static,..) in the build copy")

os.environ["CC"] = "clang"
os.environ["CXX"] = "clang++"

comp_flags = ["-O3", "-flto", "-g0", "-fstrict-aliasing", "-mcpu=apple-m1",
              "-std=c++17", "-fPIC", "-I./src", "-DNPY_NO_DEPRECATED_API",
              "-DNDEBUG", "-Wall", "-Wno-unused-function",
              "-Wno-unused-but-set-variable", "-Xpreprocessor", "-fopenmp",
              # Eigen 3.3.7 (unsupported/ConstrainedConjGrad.h) uses std::cerr
              # without <iostream>; libstdc++ pulls it in, libc++ does not
              "-include", "iostream",
              # Eigen's own OpenMP (a dynamic schedule in SparseDenseProduct)
              # off: pyMilne parallelises over pixels itself
              "-DEIGEN_DONT_PARALLELIZE",
              # ...after which Eigen no longer includes <omp.h> for pyMilne
              "-include", "omp.h"]
# OpenMP runtime: use the conda env's own libomp (llvm-openmp) when it has
# one -- its OpenBLAS loads that runtime first, and headers from a newer
# Homebrew libomp then reference symbols it lacks (__kmpc_dispatch_deinit)
omp_root = env_root if (env_root / "lib" / "libomp.dylib").exists() \
    else brew / "libomp"
link_flags = ["-flto", f"-L{omp_root / 'lib'}", "-lomp",
              f"-Wl,-rpath,{omp_root / 'lib'}",
              f"-L{brew / 'fftw' / 'lib'}", f"-Wl,-rpath,{brew / 'fftw' / 'lib'}"]

extension = Extension(
    "pyMilne",
    sources=["pyMilne.pyx", "src/wrapper_tools_spatially_coupled.cpp",
             "src/lm_sc.cpp", "src/spatially_coupled_helper.cpp"],
    include_dirs=["./", numpy.get_include(), str(repo / "eigen3"),
                  str(omp_root / "include"),
                  str(env_root / "include"), str(env_root / "include" / "eigen3"),
                  str(brew / "fftw" / "include")],
    language="c++",
    extra_compile_args=comp_flags,
    extra_link_args=link_flags,
    libraries=["fftw3"])
extension.cython_directives = {"language_level": "3"}

os.chdir(work)
setup(name="pyMilne", version="2.0",
      author="J. de la Cruz Rodriguez (ISP-SU 2018 - 2023)",
      ext_modules=[extension], cmdclass={"build_ext": build_ext},
      script_args=["build_ext", "--inplace"])

built = sorted(work.glob("pyMilne*.so"))
if not built:
    sys.exit("build produced no pyMilne*.so")
for so in built:
    shutil.copy2(so, repo / so.name)
    print(f"-> {repo / so.name}")
