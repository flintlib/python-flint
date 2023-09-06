#!/bin/bash
#
# Repair Windows wheels. See e.g.:
#
# https://github.com/scipy/scipy/blob/main/tools/wheels/repair_windows.sh

set -o errexit

# Uncomment this to run cibuildwheel locally on Windows:
# export PATH=$PATH:/c/msys64/usr/bin:/c/msys64/mingw64/bin

# We cannot use ordinary command line arguments in CI because msys2 -c mangles
# them. Instead we have a batch file to receive the arguments and convert them
# into environment variables before calling this script. When running locally
# this script could be run directly giving the parameters as command line
# arguments instead.

if [[ -z "${WHEELHOUSE}" ]]; then
  WHEELNAME=$1
fi
if [[ -z "${WHEELNAME}" ]]; then
  WHEELHOUSE=$2
fi

echo WHEELHOUSE=$WHEELHOUSE
echo WHEELNAME=$WHEELNAME

wheeldir=$(dirname $WHEELNAME)
echo $wheeldir

# delvewheel requires DLLs created by mingw64 to be stripped. This strips the
# DLLs for GMP etc that will have been built in the CIBW_BEFORE_ALL step.
#
# Previously the Arb DLLs would have been placed in .local/lib, but as of
# flint 3.0.0 all DLLs aRe in .local/bin.
strip .local/bin/*.dll

# Make sure to leave the wheel in the same directory
wheeldir=$(dirname $WHEELNAME)
pushd $wheeldir
  # Unpack the wheel and strip any .pyd DLLs inside
  wheel unpack $WHEELNAME
  rm $WHEELNAME
  strip python_flint-*/flint/*.pyd
  wheel pack python_flint-*
popd

# Make the wheel relocatable. This will fail with an error message about
# --no-mangle if strip has not been applied to all mingw64-created .dll and
# .pyd files that are needed for the wheel.
delvewheel repair $WHEELNAME 	\
	-w $WHEELHOUSE		\
	--add-path .local/bin
