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

echo "NEED UCRT64 gcc:"
which gcc

echo $PATH

# Install the mingw-w64 toolchain
# pacman -S --noconfirm mingw-w64-x86_64-gcc m4 make mingw-w64-x86_64-tools-git

# Install the mingw-ucrt toolchain
# pacman -S --noconfirm mingw-w64-ucrt-x86_64-gcc m4 make mingw-w64-ucrt-x86_64-tools-git

# This takes ~30mins
bin/build_dependencies_unix.sh
