#!/bin/bash
#
# make_wheels.sh
#
# Build relocatable Windows wheel for python_flint using the mingw-w64
# toolchain in the MSYS2 enironment.
#
# - First install Python
#
#       https://www.python.org/ftp/python/3.10.8/python-3.10.8-amd64.exe
#
# - Then checkout the code:
#
#       $ git clone https://github.com/flintlib/python-flint.git
#
# - Then install msys2
#
#       https://repo.msys2.org/distrib/x86_64/msys2-x86_64-20221028.exe
#
# - Then open msys2, cd into the checked out repo. Make sure setup.py says
#
#       libraries = ["arb", "flint", "mpfr", "gmp"]
#
# - Set the environment variable to the directory containing the installed
#   Python that we want to build a wheel for i.e. the one installed from
#   python.org. If python was on PATH then it would be
#
#        PYTHONDIR=`dirname $(which python)`
#	 PYTHONVER=3.10
#
# - Then run this script.

set -o errexit

#
# In CI this environment variable needs to be set to the directory containing
# the python.org installation of Python. If Python is installed in msys2 then
# it is also necesary to set this environment variable so that it picks up the
# right installation of Python i.e. the one that we want to build a wheel for.
#
if [[ -z "$PYTHONDIR" ]]; then
  PYTHONDIR=`dirname $(which python)`
fi
PYTHON=$PYTHONDIR/python
VER="${PYTHONVER//./}"

WHEELNAME=python_flint-0.6.0-cp$VER-cp$VER-win_amd64.whl

$PYTHON -c 'print("hello world")'

echo $PYTHONDIR
ls $PYTHONDIR
ls $PYTHONDIR/libs

# Install the mingw-w64 toolchain
pacman -S --noconfirm mingw-w64-x86_64-gcc m4 make mingw-w64-x86_64-tools-git

# This takes ~30mins
#bin/build_dependencies_unix.sh

# Add the libpython$VER.a file to Python installation
cd $PYTHONDIR
  gendef python$VER.dll
  dlltool --dllname python$VER.dll --def python$VER.def --output-lib libpython$VER.a
  mv libpython$VER.a libs
cd -

# Make a virtual environment to install the build dependencies
$PYTHON -m venv .local/venv
source .local/venv/Scripts/activate
pip install numpy cython wheel delvewheel psutil

# Pass this flag to setup.py
export PYTHON_FLINT_MINGW64_TMP=true

# Build the wheel
C_INCLUDE_PATH=.local/include/ LIBRARY_PATH=.local/lib/ python setup.py build_ext -cmingw32 -f bdist_wheel

# delvewheel requires DLLs created by mingw64 to be stripped
#
# https://github.com/scipy/scipy/blob/main/tools/wheels/repair_windows.sh
strip .local/bin/*.dll .local/lib/*.dll

# Unpack the wheel and strip the .pyd DLL inside
cd dist
wheel unpack $WHEELNAME
rm $WHEELNAME
strip python_flint-*/flint/*.pyd
wheel pack python_flint-*
cd ..

# Make the wheel relocatable
delvewheel repair dist/python_flint-0.6.0-cp$VER-cp$VER-win_amd64.whl \
        --add-path .local/bin:.local/lib/

# Make a virtual enironment to test the wheel
$PYTHON -m venv test_env
source test_env/Scripts/activate
pip install wheelhouse/$WHEELNAME
python -c 'import flint; print(flint.fmpz(2) + 2)'  # 2 + 2 = 4?
