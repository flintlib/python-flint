#!/bin/bash

set -o errexit

pacman -S --noconfirm \
    mingw-w64-clang-aarch64-toolchain\
    mingw-w64-clang-aarch64-tools-git\
    m4\
    make\
    base-devel\
    autoconf-wrapper\
    automake-wrapper\
    libtool\
    git\
    #

bin/build_dependencies_unix.sh \
    --disable-fat\
    --use-gmp-github-mirror\
    --patch-C23\
    #
