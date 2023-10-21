#!/usr/bin/env bash
#
# This script can be used to test cibuildwheel locally on OSX/Linux
#
# It is also worth commenting out the BEFORE_ALL line to build GMP etc after you have
# built those once because that is by far the slowest step.
#

rm -f wheelhouse/*

# bin/build_dependencies_unix.sh places headers and shared libraries under .local
export CIBW_ENVIRONMENT='C_INCLUDE_PATH=$(pwd)/.local/include/ LIBRARY_PATH=$(pwd)/.local/lib/ LD_LIBRARY_PATH=$(pwd)/.local/lib:$LD_LIBRARY_PATH PYTHON_FLINT_MINGW64=true'

export CIBW_BUILD='cp39-* cp310-* cp311-* cp312-*'
# export CIBW_BUILD='cp311-*'
export CIBW_SKIP='*-win32 *-manylinux_i686 *-musllinux_*'

# export CIBW_ARCHS_MACOS="x86_64"
export CIBW_ARCHS_MACOS="arm64"

export CIBW_BEFORE_ALL_LINUX=bin/cibw_before_all_linux.sh
# export CIBW_BEFORE_ALL_MACOS=bin/cibw_before_all_macosx_x86_64.sh
export CIBW_BEFORE_ALL_MACOS=bin/cibw_before_all_macosx_arm64.sh
export CIBW_BEFORE_ALL_WINDOWS='C:\\msys64\\usr\\bin\\bash bin/cibw_before_all_windows.sh'

export CIBW_BEFORE_BUILD='pip install numpy cython delvewheel'
export CIBW_BEFORE_BUILD_WINDOWS='C:\\msys64\\usr\\bin\\bash bin/cibw_before_build_windows.sh'

export CIBW_REPAIR_WHEEL_COMMAND_WINDOWS='bin\cibw_repair_wheel_command_windows.bat {dest_dir} {wheel}'

# export CIBW_TEST_COMMAND="python -c 'import flint; print(str(flint.fmpz(2)))'"
export CIBW_TEST_COMMAND="python -m flint.test"

# cibuildwheel --platform linux
# cibuildwheel --platform windows
cibuildwheel --platform macos
