#!/bin/bash

set -e

while [[ $# -gt 0 ]]
do
  key="$1"
  case $key in
    -h|--help)
      echo "bin/pyodide_build_dependencies.sh [options]"
      echo
      echo "Build local emscripten installs of python-flint's dependencies."
      echo
      echo "Supported options:"
      echo "  --help            - show this help message"
      echo "  --wasm-library-dir <WASM_LIBRARY_DIR> - directory to install libraries"
      echo "  --flint-commit <FLINT_COMMIT> - flint commit to build"
      echo
      exit
    ;;
    --wasm-library-dir)
      # e.g. --wasm-library-dir /path/to/wasm-library-dir
      WASM_LIBRARY_DIR="$2"
      shift
      shift
    ;;
    --flint-commit)
      # e.g. --flint-commit 3.3.1
      FLINT_COMMIT="$2"
      shift
      shift
    ;;
  *)
    2>&1 echo "unrecognised argument:" $key
    exit 1
    ;;
  esac
done


if [ -z "$WASM_LIBRARY_DIR" ]; then
  echo "WASM_LIBRARY_DIR not set"
  exit 1
fi

# Sets versions of GMP, MPFR and FLINT:

source bin/build_variables.sh

# Download mirrored copy of source distributions for GMP and MPFR

git clone https://github.com/oscarbenjamin/gmp_mirror.git
cp gmp_mirror/gmp-$GMPVER.tar.xz .
cp gmp_mirror/mpfr-$MPFRVER.tar.gz .
tar -xf gmp-$GMPVER.tar.xz
tar -xf mpfr-$MPFRVER.tar.gz


# ---------------------------Build GMP ----------------------------------#


cd gmp-$GMPVER

    emconfigure ./configure \
        --disable-dependency-tracking \
        --host none \
        --disable-shared \
        --enable-static \
        --enable-cxx \
        --prefix=$WASM_LIBRARY_DIR

    emmake make -j $(nproc)
    emmake make install

cd ..


# ---------------------------Build MPFR ----------------------------------#


cd mpfr-$MPFRVER

    emconfigure ./configure \
        --disable-dependency-tracking \
        --disable-shared \
        --with-gmp=$WASM_LIBRARY_DIR \
        --prefix=$WASM_LIBRARY_DIR

    emmake make -j $(nproc)
    emmake make install

cd ..


# ---------------------------Build FLINT----------------------------------#


if [ -z "$FLINT_COMMIT" ]; then
    curl -O -L https://github.com/flintlib/flint/releases/download/v$FLINTVER/flint-$FLINTVER.tar.gz
    tar xf flint-$FLINTVER.tar.gz
    cd flint-$FLINTVER
else
    git clone https://github.com/flintlib/flint --branch $FLINT_COMMIT
    cd flint
fi

    ./bootstrap.sh

    emconfigure ./configure \
        --disable-dependency-tracking \
        --disable-shared \
        --prefix=$WASM_LIBRARY_DIR \
        --with-gmp=$WASM_LIBRARY_DIR \
        --with-mpfr=$WASM_LIBRARY_DIR \
        --host=wasm32-unknown-emscripten \
        --disable-assembly \
        --disable-pthread

    emmake make -j $(nproc)
    emmake make install

cd ..
