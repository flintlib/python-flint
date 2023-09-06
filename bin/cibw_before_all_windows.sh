#!/bin/bash

set -o errexit

# Uncomment this to run cibuildwheel locally on Windows:
# export PATH=$PATH:/c/msys64/usr/bin:/c/msys64/mingw64/bin

#
# Make a setup.cfg to specify compiling with mingw64 (even though it says
# mingw32...)
#
echo '[build]' > setup.cfg
echo 'compiler = mingw32' >> setup.cfg
cat setup.cfg

# Install the mingw-w64 toolchain and build tools
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

# This takes ~30mins
bin/build_dependencies_unix.sh --use-gmp-github-mirror
