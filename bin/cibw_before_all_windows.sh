#!/bin/bash

set -o errexit

pacman -S --noconfirm \
    mingw-w64-x86_64-gcc\
    mingw-w64-x86_64-tools-git\
    m4\
    make\
    base-devel\
    autoconf-wrapper\
    automake-wrapper\
    libtool\
    #

bin/build_dependencies_unix.sh \
    --use-gmp-github-mirror\
    --patch-C23\
    #
