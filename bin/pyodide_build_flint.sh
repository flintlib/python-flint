#!/bin/bash

set -e

curl -L https://github.com/flintlib/flint/releases/download/v3.2.2/flint-3.2.2.tar.xz -o flint-3.2.2.tar.xz
tar -xf flint-3.2.2.tar.xz

cd flint-3.2.2


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
