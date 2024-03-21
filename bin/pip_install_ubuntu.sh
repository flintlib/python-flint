#!/bin/bash

set -o errexit

# This script should work to install python-flint on Ubuntu from a VCS checkout
#
#  $ git clone https://github.com/flintlib/python-flint.git
#  $ cd python-flint
#  $ bin/pip_install_ubuntu.sh .
#
# To install an sdist from PyPI, use
#
#  $ bin/pip_install_ubuntu.sh python-flint
#
# The Flint version to build can be provided as an argument, e.g.
#
#  $ bin/pip_install_ubuntu.sh python-flint v3.0.1
#
# The version is a tag or branch in the Flint repository. If not provided, the
# script will default to downloading the release tarball hard-coded below. Only
# Flint 3 or newer will work.

if [ -z "$2" ]; then
    echo "Building from release tarball"
    FLINT_GIT=""
    FLINTVER=3.0.1
else
    echo "Building from git: $2"
    FLINT_GIT=$2
fi
# Either . or python-flint. Passed to pip install
PYTHON_FLINT=$1

# Install runtime and build dependencies

# First install their dependencies and build dependencies
sudo apt-get update
sudo apt-get install libgmp-dev libmpfr-dev xz-utils ninja-build

if [ -z "$FLINT_GIT" ]; then
    # Install from release tarball
    echo "Installing Flint $FLINTVER from release tarball"
    curl -O -L https://www.flintlib.org/flint-$FLINTVER.tar.gz
    tar xf flint-$FLINTVER.tar.gz
    cd flint-$FLINTVER
else
    # Install from git
    echo "Installing Flint from git: $FLINT_GIT"
    git clone https://github.com/flintlib/flint.git
    cd flint
    git checkout $FLINT_GIT
fi

# This will default to installing in /usr/local. If you want to install in a
# non-standard location then configure flint with
#    ./configure --disable-static --prefix=$PREFIX
# If $PREFIX is not in default search paths, then at build time set
#    export C_INCLUDE_PATH=$PREFIX/include
# and at runtime set
#    export LD_LIBRARY_PATH=$PREFIX/lib
./bootstrap.sh
./configure --disable-static
make -j
sudo make install
cd ..

ls -l /usr/local/lib
sudo ldconfig /usr/local/lib

# Build python-flint from sdist
echo -----------------------------------------------------------
echo
echo     Running:
echo        $ pip install --no-binary python-flint $PYTHON_FLINT
echo
echo -----------------------------------------------------------

pip install --no-binary python-flint $PYTHON_FLINT
