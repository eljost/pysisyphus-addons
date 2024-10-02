{ buildPythonPackage
, lib
, gfortran
, cmake
, meson-python
, pkg-config
, sympleints
, which
, blas
, lapack
, cython
, numpy
, test-drive
, lmax ? 4
, lauxmax ? 5
}:

buildPythonPackage rec {
  pname = "pysisyphus-addons";
  version = "0.1.0";

  src = lib.cleanSource ../.;

  pyproject = true;

  build-system = [
    meson-python
    cython
  ];

  dependencies = [
    #sympleints    
    numpy
  ];

  nativeBuildInputs = [
    gfortran
    pkg-config
    sympleints
    which
  ];

  buildInputs = [
    blas
    lapack
  ];

  pypaBuildFlags = [
    "-Csetup-args=-Dlmax=${builtins.toString lmax}"
    "-Csetup-args=-Dlauxmax=${builtins.toString lauxmax}"
  ];

  checkInputs = [ 
    test-drive
  ];

  meta = with lib; {
    description = "Compiled addons for pysisyphus";
    license = licenses.eupl12;
    homepage = "https://github.com/eljost/sympleints";
    maintainers = [ maintainers.sheepforce ];
  };
}
