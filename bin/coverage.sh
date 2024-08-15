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
set -o errexit

meson configure build -Dcoverage=true
spin run -- coverage run -m flint.test
coverage report -m
coverage html
