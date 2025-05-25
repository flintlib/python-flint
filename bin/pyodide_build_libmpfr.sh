#!usr/bin/env bash

set -e

curl -L https://ftp.gnu.org/gnu/mpfr/mpfr-4.2.1.tar.xz -o mpfr-4.2.1.tar.xz
tar -xf mpfr-4.2.1.tar.xz

cd mpfr-4.2.1

emconfigure ./configure \
    --disable-dependency-tracking \
    --disable-shared \
    --with-gmp=$WASM_LIBRARY_DIR \
    --prefix=$WASM_LIBRARY_DIR

emmake make -j $(nproc)
emmake make install
