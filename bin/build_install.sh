#!/bin/bash

set -o errexit

FLINTVER=$1

# Install GMP, MPFR and build tools from Ubuntu repos
sudo apt-get update
sudo apt-get install libgmp-dev libmpfr-dev xz-utils ninja-build

# Build flint and install to /usr/local
curl -O -L https://www.flintlib.org/flint-$FLINTVER.tar.gz
tar xf flint-$FLINTVER.tar.gz
cd flint-$FLINTVER
  ./bootstrap.sh
  ./configure --disable-static
  make -j
  sudo make install
cd ..

# Ensure the the libflint.so is found at runtime
sudo ldconfig /usr/local/lib

echo "Contents of /usr/local/lib:"
ls -l /usr/local/lib
echo

echo "Contents of /usr/local/lib/pkgconfig/flint.pc:"
cat /usr/local/lib/pkgconfig/flint.pc
echo

# Build python-flint using meson-python. This will use pkgconfig to find the
# flint library installed in /usr/local.
pip install .
