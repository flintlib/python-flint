#!/usr/bin/env bash
#
# Build the flint._flint module in place.

C_INCLUDE_PATH=.local/include/ LIBRARY_PATH=.local/lib/ python setup.py build_ext --inplace
