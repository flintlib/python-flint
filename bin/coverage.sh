#!/bin/bash
#
# This needs a patched Cython:
#
#   pip install git+https://github.com/oscarbenjamin/cython.git@pr_relative_paths
#
# That patch has been submitted as a pull request:
#
#   https://github.com/cython/cython/pull/6341
#
# Arguments to this script are passed to python -m flint.test e.g. to skip
# doctests and run in quiet mode:
#
#    bin/coverage.sh -qt
#
set -o errexit

RC="--rcfile=.coveragerc.meson"

# See https://github.com/cython/cython/issues/6658
# Needed for Python 3.13 only
#pip uninstall cython
#pip install git+https://github.com/cython/cython.git@fdbca99

meson setup build -Dcoverage=true
spin run -- coverage run $RC -m flint.test $@
coverage report $RC -m --sort=cover
coverage html $RC
