#!/bin/bash

#
# This script was previously needed to make libpythonXX.a on Windows when using
# MinGW and setuptools. This is no longer needed now that we use meson.
#

set -o errexit

# Uncomment this to run cibuildwheel locally on Windows:
# export PATH=$PATH:/c/msys64/usr/bin:/c/msys64/mingw64/bin

# VER should be set be e.g. 310 for Python 3.10
VER=`python -c 'import sys; print("%s%s" % sys.version_info[:2])'`
echo VER=${VER}

###################################################
#  Find parent Python installation from the venv  #
###################################################

which python
PYTHONBIN=`dirname $(which python)`
PYTHONDIR=`dirname $PYTHONBIN`
cfgfile=$PYTHONDIR/pyvenv.cfg
homeline=`grep home $cfgfile`
homepath=${homeline#*=}

echo ---------------------------------------------------
echo $homepath
echo ---------------------------------------------------

###################################################
#  Find pythonXX.dll and make a .a library        #
###################################################

cd $homepath
gendef python${VER}.dll
dlltool --dllname python${VER}.dll 	\
	--def python${VER}.def 		\
	--output-lib libpython${VER}.a

mv libpython${VER}.a libs

###################################################
#  Install build dependencies                     #
###################################################

pip install cython numpy delvewheel wheel
