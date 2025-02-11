# pysisyphus-addons
Compiled Fortran addons to pysisyphus, that would be too slow in pure Python. The Fortran code
is interfaced by Cython functions that are easily called from Python.

When using procedures from this project that involve handling of bigger matrices please execute
```bash
ulimit -s unlimited
```
in your shell; otherwise the code will segfault.

## Overview
This project provides fast and OMP-parallel Fortran procedures for the following tasks:

- `diabatization/`: Construction of the Coulomb-Tensor in the basis of adiabatic electronic states using density fitting. Required for Edmiston-Ruedenberg diabatization.
- `gdma/`: Evaluation of the density on a given DFT grid, as required for GDMA.
- `grid/`: Evaluation of (multiple) densities on a given cubic grid, as required for the calculation of .cub files.
- `std/`: 2-electron integral tensors for eXact simplified TD-DFT w/ density fitting.

Code found in `common/` includes basic driver procedures to calculate various integrals, to set up integral screeners as well as
simple wrappers for various linear algebra procedures from LAPACK. The `intor/` directory contains procedures for benchmarking the performance of the density fitting integrals and their associated screeners.


## Requirements & Installation
As this project contains mainly Fortran code it is built with `meson` and `meson-python`. The required integral code
for the density fitting (3-center-2-electron- and 2-center-2-electron-integrals), as well as selected 4-center-2-electron-integrals (AB|AB) for Schwarz-screening are generated on-the-fly at build-time by invoking [sympleints](https://github.com/eljost/sympleints), which must be available/installed.
By default, selected 4c2e-integrals are built up to (GG|GG) and 3c2e- and 2c2e integrals for density fitting are built up to (GG|H) and (H|H).
Maximum angular momenta for integral generation can be controlled with the `-Dlmax` and `-Dlauxmax` options.

`pysisyphus-addons` utilizes various `BLAS` and `LAPACK` procedures, so both libraries must be available. Different implementations can be selected  via `-Dblas=[blas|openblas]` and `-Dlapack=[lapack|openblas]`.

A newer GCC version, e.g. 13+,  is required for building this project.

Below you can find an example to check if `pysisyphus-addons` can be built on your system. Depending on your setup you'll have to adapt the compiler version and the path to the `OpenBLAS` pkgconfig-file, as well as the number of cores when building `-j12`.

```bash
PKG_CONFIG_PATH=/home/johannes/programme/openblas/lib/pkgconfig \
FC=gfortran-13 CC=gcc-13 CXX=g++-13	 \
python -m build \
	-Csetup-args=-Dlmax=2 \
  -Csetup-args=-Dlauxmax=2 \
	-Csetup-args=-Dblas=openblas \
	-Csetup-args=-Dlapack=openblas \
	-Ccompile-args="-j12"
```

For an actually useful production build `-Dlmax` and `-Dlauxmax` should be kept at the default values. At their default values, integral generation currently takes about 1 hour and actually compiling them will also take about 1 hour and will require a machine with at least 30 GB free memory.

## Usage
On its own `pysisyphus-addons` is not useful, but it is called by [pysisyphus](https://github.com/eljost/pysisyphus) with the required inputs.
