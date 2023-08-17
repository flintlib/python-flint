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

# The commented commands below would attempt to build python-flint against a
# system installation of Flint and Arb in Ubuntu.
#
# Ubuntu 23.04 has Flint 2.9.0 and Arb 2.23.0, so this script might work there
# (for python-flint 0.4.x). That is untested though (23.04 not available in CI).
#
# With Ubuntu 22.04, this will build but then crashes when running the tests.
# most likely this is because the versions of flint and flint-arb are too old.
# At the time of writing in Ubuntu 22.04 there is Flint 2.8.4 and Arb 2.22. The
# main CI tests and wheels for python-flint 0.4.x are built with Flint 2.9.0
# and Arb 2.23.0.
#
# Link against libflint-arb instead of libarb on Ubuntu
# export PYTHON_FLINT_LIBFLINT_ARB=1
# sudo apt-get update
# sudo apt-get install libflint-dev libflint-arb-dev


# Build Flint and Arb manually
#
# First install their dependencies and build dependencies
sudo apt-get update
sudo apt-get install libgmp-dev libmpfr-dev xz-utils

# Only these *EXACT* versions will work.
FLINTVER=2.9.0
ARBVER=2.23.0

# This will default to installing in /usr/local. If you want to install in a
# non-standard location then configure flint with
#    ./configure --disable-static --prefix=$PREFIX
# and arb with
#    ./configure --disable-static --prefix=$PREFIX --with-flint=$PREFIX
# If $PREFIX is no in default search paths, then at build time set
#    export C_INCLUDE_PATH=$PREFIX/include
# and at runtime set
#    export LD_LIBRARY_PATH=$PREFIX/lib

curl -O -L https://www.flintlib.org/flint-$FLINTVER.tar.gz
tar xf flint-$FLINTVER.tar.gz
cd flint-$FLINTVER
  ./configure --disable-static
  make -j
  sudo make install
cd ..

curl -O -L https://github.com/fredrik-johansson/arb/archive/refs/tags/$ARBVER.tar.gz
mv $ARBVER.tar.gz arb-$ARBVER.tar.gz
tar xf arb-$ARBVER.tar.gz
cd arb-$ARBVER
  ./configure --disable-static
  make -j
  sudo make install
cd ..

ls -l /usr/local/lib
sudo ldconfig /usr/local/lib

# Python build requirements. Ideally these would be in pyprojec.toml, but
# first need to migrate from setup.py to pyproject.toml.
pip install 'cython>=3' numpy wheel

# Install from checkout (or sdist).
echo -----------------------------------------------------------
echo
echo     Running:
echo        $ pip install --no-build-isolation $1
echo
echo -----------------------------------------------------------

pip install --no-build-isolation $1
