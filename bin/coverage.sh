#!/bin/bash
#
# Arguments to this script are passed to python -m flint.test e.g. to skip
# doctests and run in quiet mode:
#
#    bin/coverage.sh -qt
#
set -o errexit

RC="--rcfile=.coveragerc.meson"

# See https://github.com/cython/cython/issues/6658
# Needed for Python 3.13 only but the plugin does not work with 3.13 anyway...
#pip uninstall -y cython
#pip install git+https://github.com/cython/cython.git@fdbca99

pip uninstall -y cython
pip install --pre cython       # unpinned to pick up new releases in CI
# pip install cython==3.1.0a1  # known working version for Python < 3.13

meson setup build -Dcoverage=true
spin run -- coverage run $RC -m flint.test $@
coverage report $RC -m --sort=cover
coverage html $RC
