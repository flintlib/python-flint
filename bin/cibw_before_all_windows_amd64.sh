#!/bin/bash

set -o errexit

pacman -S --noconfirm \
    mingw-w64-ucrt-x86_64-gcc \
    mingw-w64-ucrt-x86_64-llvm-tools \
    mingw-w64-ucrt-x86_64-tools-git \
    m4 \
    make \
    base-devel \
    autoconf-wrapper \
    automake-wrapper \
    libtool \
    git \
    #

bin/cibw_before_all_windows.sh \
    x64 \
    --use-gmp-github-mirror \
    --patch-C23 \
    #
