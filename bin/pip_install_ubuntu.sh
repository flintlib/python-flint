#!/bin/bash

set -o errexit

# This script should work to install python-flint on Ubuntu from a VCS checkout
#
#  $ git checkout https://github.com/flintlib/python-flint.git
#  $ bin/pip_install_ubuntu.sh .
#
# To install an sdist from PyPI, use
#
#  $ bin/pip_install_ubuntu.sh python-flint
#

# Install runtime and build dependencies

# First install their dependencies and build dependencies
sudo apt-get update
sudo apt-get install libgmp-dev libmpfr-dev xz-utils

# Only Flint 3 or newer will work.
FLINTVER=3.0.0

# This will default to installing in /usr/local. If you want to install in a
# non-standard location then configure flint with
#    ./configure --disable-static --prefix=$PREFIX
# If $PREFIX is not in default search paths, then at build time set
#    export C_INCLUDE_PATH=$PREFIX/include
# and at runtime set
#    export LD_LIBRARY_PATH=$PREFIX/lib

curl -O -L https://www.flintlib.org/flint-$FLINTVER.tar.gz
tar xf flint-$FLINTVER.tar.gz
cd flint-$FLINTVER
  ./bootstrap.sh
  ./configure --disable-static
  make -j
  sudo make install
cd ..

ls -l /usr/local/lib
sudo ldconfig /usr/local/lib

# Python build requirements. Ideally these would be in pyproject.toml, but
# first need to migrate from setup.py to pyproject.toml.
pip install numpy cython setuptools wheel

# Install from checkout (or sdist).
echo -----------------------------------------------------------
echo
echo     Running:
echo        $ pip install --no-binary :all: --no-build-isolation $1
echo
echo -----------------------------------------------------------

pip install --no-binary :all: --no-build-isolation $1
