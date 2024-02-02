#!/usr/bin/env bash
#
# This script can be used to test cibuildwheel locally on OSX/Linux
#
# It is also worth commenting out the BEFORE_ALL line to build GMP etc after you have
# built those once because that is by far the slowest step.
#

rm -f wheelhouse/*

export CIBW_SKIP='*-win32 *-manylinux_i686 *-musllinux_*'

# export CIBW_ARCHS_MACOS="x86_64"
export CIBW_ARCHS_MACOS="arm64"

# export CIBW_BEFORE_ALL_MACOS=bin/cibw_before_all_macosx_x86_64.sh
export CIBW_BEFORE_ALL_MACOS=bin/cibw_before_all_macosx_arm64.sh
export CIBW_BEFORE_ALL_WINDOWS='C:\\msys64\\usr\\bin\\bash bin/cibw_before_all_windows.sh'

export CIBW_BEFORE_BUILD='pip install numpy cython delvewheel'
export CIBW_BEFORE_BUILD_WINDOWS='C:\\msys64\\usr\\bin\\bash bin/cibw_before_build_windows.sh'

export CIBW_TEST_COMMAND="python -m flint.test"  # override setting in pyproject.toml

# cibuildwheel --platform linux
# cibuildwheel --platform windows
cibuildwheel --platform macos
