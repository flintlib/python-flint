#!/usr/bin/env bash
#
# Build local installs of python-flint's dependencies.
# This should be run after bin/download_dependencies.sh
#

set -o errexit
set -o xtrace

# Sets variables like PREFIX and dependency versions
source bin/build_variables.sh

cd $PREFIX
mkdir -p src
cd src

curl -O https://gmplib.org/download/gmp/gmp-$GMPVER.tar.xz
curl -O https://www.mpfr.org/mpfr-current/mpfr-$MPFRVER.tar.gz
curl -O https://www.flintlib.org/flint-$FLINTVER.tar.gz
curl -O -L https://github.com/fredrik-johansson/arb/archive/refs/tags/$ARBVER.tar.gz
mv $ARBVER.tar.gz arb-$ARBVER.tar.gz

tar xf mpfr-$MPFRVER.tar.gz
tar xf flint-$FLINTVER.tar.gz
tar xf arb-$ARBVER.tar.gz
tar xf gmp-$GMPVER.tar.xz

cd gmp-$GMPVER
  ./configure --prefix=$PREFIX --enable-fat --enable-shared=yes --enable-static=no
  make -j3
  make install
cd ..

cd mpfr-$MPFRVER
  ./configure --prefix=$PREFIX --with-gmp=$PREFIX --enable-shared=yes --enable-static=no
  make -j3
  make install
cd ..

cd flint-$FLINTVER
  ./configure --prefix=$PREFIX --with-gmp=$PREFIX --with-mpfr=$PREFIX --disable-static
  make -j3
  make install
cd ..

cd arb-$ARBVER
  ./configure --prefix=$PREFIX --with-flint=$PREFIX --with-gmp=$PREFIX --with-mpfr=$PREFIX --disable-static
  make -j3
  make install
cd ..


echo
echo -----------------------------------------------------------------------
echo
echo Build dependencies for python-flint compiled as shared libraries in:
echo $PREFIX
echo
echo Versions:
echo GMP: $GMPVER
echo MPFR: $MPFRVER
echo Flint: $FLINTVER
echo Arb: $ARBVER
echo
echo -----------------------------------------------------------------------
echo
