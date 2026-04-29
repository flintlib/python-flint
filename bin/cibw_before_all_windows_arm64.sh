#!/bin/bash

set -o errexit

pacman -S --noconfirm \
    mingw-w64-clang-aarch64-toolchain \
    mingw-w64-clang-aarch64-llvm-tools \
    mingw-w64-clang-aarch64-tools-git \
    m4 \
    make \
    base-devel \
    autoconf-wrapper \
    automake-wrapper \
    libtool \
    git \
    #

bin/cibw_before_all_windows.sh \
    arm64 \
    --use-gmp-github-mirror \
    --disable-assembly \
    --host aarch64-w64-mingw32 \
    --patch-C23 \
    --patch-ldd \
    --patch-immintrin \
    #
