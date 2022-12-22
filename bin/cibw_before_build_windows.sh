#!/bin/bash

set -o errexit

# Uncomment this to run cibuildwheel locally on Windows:
# export PATH=$PATH:/c/msys64/usr/bin:/c/msys64/mingw64/bin

# msys2 will not inherit the PATH for the virtual environment
export PATH=$PATH:$VIRTUAL_ENV_BIN
echo PATH=$PATH

# VER should be set be e.g. 310 for Python 3.10
VER=`python -c 'import sys; print("%s%s" % sys.version_info[:2])'`
echo VER=${VER}

###################################################
#  Find parent Python installation from the venv  #
###################################################

which python
PYTHONBIN=$VIRTUAL_ENV_BIN
echo PYTHONBIN=$PYTHONBIN
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

pip install Cython numpy delvewheel wheel
