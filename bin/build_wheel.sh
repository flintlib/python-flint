#!/usr/bin/env bash
#
# Compile a python-flint wheel using the dependencies built by
# build_dependencies_unix.sh (which should be run first).

set -o errexit

source bin/build_variables.sh

python3 -m venv $PREFIX/venv
source $PREFIX/venv/bin/activate
pip install -U pip
pip install numpy cython wheel

C_INCLUDE_PATH=.local/include/ LIBRARY_PATH=.local/lib/ pip wheel .
