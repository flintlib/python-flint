#!/usr/bin/env bash
#
# Compile a python-flint wheel using the dependencies built by
# build_dependencies_unix.sh (which should be run first).

set -o errexit

PYTHONFLINTVER=0.3.0

source bin/build_variables.sh

python3 -m venv $PREFIX/venv
source $PREFIX/venv/bin/activate
pip install -U pip wheel delocate
pip install numpy cython==0.27.3
# N.B. bugs in both older and newer Cython versions...

C_INCLUDE_PATH=.local/include/ LIBRARY_PATH=.local/lib/ pip wheel .

wheelfinal=*.whl

# On OSX bundle the dynamic libraries for the dependencies
mkdir -p wheelhouse
delocate-wheel -w wheelhouse $wheelfinal

echo ------------------------------------------
echo
echo Built wheel: wheelhouse/$wheelfinal
echo
echo Link dependencies:
delocate-listdeps wheelhouse/$wheelfinal
echo
pip install wheelhouse/$wheelfinal
echo
echo Demonstration:
echo
python -c 'import flint; print("(3/2)**2 =", flint.fmpq(3, 2)**2)'
echo
echo Done!
echo
echo ------------------------------------------
