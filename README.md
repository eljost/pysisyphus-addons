# pysisyphus-addons
Compiled addons for pysisyphus utilizing Fortran & Cython, built by meson/meson-python.

Currently, fast and OMP-parallelized Fortran procedures are provided for the following tasks:

- `common/`: Basic driver procedures to calculate various integrals including screening and to carry out basic linear algebra tasks.
- `diabatization/`: Construction of the Coulomb-Tensor in the basis of adiabatic electronic states using density fitting. Required for Edmiston-Ruedenberg diabatization.
- `gdma/`: Evaluation of the density on a given DFT grid, as required for GDMA.
- `grid/`: Evaluation of the density on cubic grids, as required for the calculation of .cub files.
- `intor/`: Procedure for benchmarking the performance of the density fitting integrals and associated screening procedures.
- `std/`: 2-electron integral tensors for eXact simplified TD-DFT w/ density fitting.

When using procedures from `diabatization/`
```bash
ulimit -s unlimited
```
is probably required.

## Requirements & Installation
As this project contains mainly Fortran code it is built with `meson` and `meson-python`. The required integral code
for the density fitting (3-center-2-electron- and 2-center-2-electron-integrals), as well as selected 4-center-2-electron-integrals (AB|AB) for the Schwarz-screening are generated on-the-fly at compilation time by invoking [sympleints](https://github.com/eljost/sympleints), so [sympleints](https://github.com/eljost/sympleints) must be available/installed at build-time.
By default, selected 4c2e-integrals are built up to (GG|GG) and 3c2e- and 2c2e integrals are built up to (GG|H) and (H|H).
Maximum angular momenta for integral generation can be controlled with the `-Dlmax` and `-Dlauxmax` options.

`pysisyphus-addons` utilizes various `BLAS` and `LAPACK` procedures, so both libraries must be available. Different implementations can be selected  via `-Dblas=[blas|openblas]` and `-Dlapack=[lapack|openblas]`.

A newer GCC version, e.g. 13+,  is required for building this project.

Below you can find an example for building checking if `pysisyphus-addons` can be build on your system. Depending on your setup you'll have to adapt the compiler version and the path to the `OpenBLAS` pkgconfig-file.

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

For a production build that is actually useful `-Dlmax` and `-Dlauxmax` should be kept at the default values. At their default values, integral generation currently takes about 1 hour and actually compiling them will also take about 1 hour and will require a machine with at least 30 GB free memory.

## Usage
On its own `pysisyphus-addons` is not useful, but it is called by [pysisyphus](https://github.com/eljost/pysisyphus).
