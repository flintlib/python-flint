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

meson setup build -Dcoverage=true
spin run -- coverage run -m flint.test $@
coverage report -m --sort=cover
coverage html
