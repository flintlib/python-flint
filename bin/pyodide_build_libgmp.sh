#!usr/bin/env bash

set -e

curl -L https://ftp.gnu.org/gnu/gmp/gmp-6.3.0.tar.xz -o gmp-6.3.0.tar.xz
tar -xf gmp-6.3.0.tar.xz

cd gmp-6.3.0

emconfigure ./configure \
    --disable-dependency-tracking \
    --host none \
    --disable-shared \
    --enable-static \
    --enable-cxx \
    --prefix=$WASM_LIBRARY_DIR

emmake make -j $(nproc)
emmake make install
