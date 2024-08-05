#!/bin/bash

set -o errexit

#
# This script should work to build and install Flint from git on Ubuntu
#
#  $ git clone https://github.com/flintlib/python-flint.git
#  $ cd python-flint
#  $ bin/install_flint_ubuntu.sh v3.1.0
#
# The version is a tag or branch in the Flint repository.
#
# Then to install an sdist from PyPI, use
#
#  $ pip install python-flint
#
# To install python-flint from the git checkout
#
#  $ pip install .
#

echo "Building from git: $1"
GIT_REF=$1

# Install runtime and build dependencies

# First install their dependencies and build dependencies
sudo apt-get update
sudo apt-get install libgmp-dev libmpfr-dev xz-utils ninja-build

#
# This will default to installing in /usr/local. If you want to install in a
# non-standard location then configure flint with
#    ./configure --disable-static --prefix=$PREFIX
# If $PREFIX is not in default search paths, then at build time set
#    export C_INCLUDE_PATH=$PREFIX/include
# and at runtime set
#    export LD_LIBRARY_PATH=$PREFIX/lib
#
echo "Installing Flint from git: $GIT_REF"
git clone https://github.com/flintlib/flint.git
cd flint
  git checkout $GIT_REF
  ./bootstrap.sh
  ./configure --disable-static
  make -j
  sudo make install
cd ..

ls -l /usr/local/lib
sudo ldconfig /usr/local/lib
