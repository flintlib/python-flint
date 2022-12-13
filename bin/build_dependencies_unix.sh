#!/usr/bin/env bash
#
# Build local installs of python-flint's dependencies. This should be run
# before attempting to build python-flint itself.

set -o errexit

# ------------------------------------------------------------------------- #
#                                                                           #
# Supported options:                                                        #
#                                                                           #
# --gmp gmp   - build based on GMP (default)                                #
# --gmp mpir  - build based on MPIR (instead of GMP)                        #
#                                                                           #
# ------------------------------------------------------------------------- #

USE_GMP=gmp

while [[ $# -gt 0 ]]
do
  key="$1"
  case $key in
    -h|--help)
      echo "bin/download_dependencies.sh [--gmp gmp|mpir] [--host HOST]"
      exit
    ;;
    --gmp)
      # e.g. --gmp gmp or --gmp mpir
      USE_GMP="$2"
      if [[ "$USE_GMP" != "gmp" && "$USE_GMP" != "mpir" ]]; then
        echo "--gmp option should be gmp or mpir"
        exit 1
      fi
      shift
      shift
    ;;
    --host)
      # e.g. --host x86_64-unknown-linux-gnu
      # or   --host x86_64-apple-darwin
      HOST_ARG="$2"
      shift
      shift
    ;;
  esac
done

# ------------------------------------------------------------------------- #
#                                                                           #
# The build_variables.sh script sets variables specifying the versions to   #
# use for all dependencies and also the PREFIX variable.                    #
#                                                                           #
# ------------------------------------------------------------------------- #

source bin/build_variables.sh

cd $PREFIX
mkdir -p src
cd src

# ------------------------------------------------------------------------- #
#                                                                           #
# Now build all dependencies.                                               #
#                                                                           #
# ------------------------------------------------------------------------- #

if [ $USE_GMP = "gmp" ]; then

  # ----------------------------------------------------------------------- #
  #                                                                         #
  #                            GMP                                          #
  #                                                                         #
  # ----------------------------------------------------------------------- #

  curl -O https://gmplib.org/download/gmp/gmp-$GMPVER.tar.xz
  tar xf gmp-$GMPVER.tar.xz
  cd gmp-$GMPVER
    # Show the output of configfsf.guess
    ./configfsf.guess
    ./configure --prefix=$PREFIX\
      --enable-fat\
      --enable-shared=yes\
      --enable-static=no\
      --host=$HOSTARG
    make -j3
    make install
  cd ..

  FLINTARB_WITHGMP="--with-gmp=$PREFIX"

else

  # ----------------------------------------------------------------------- #
  #                                                                         #
  #                    YASM (needed to build MPIR)                          #
  #                                                                         #
  # ----------------------------------------------------------------------- #

  curl -O http://www.tortall.net/projects/yasm/releases/yasm-$YASMVER.tar.gz
  tar xf yasm-$YASMVER.tar.gz
  cd yasm-$YASMVER
    ./configure --prefix=$PREFIX
    make -j3
    make install
  cd ..

  # ----------------------------------------------------------------------- #
  #                                                                         #
  #                            MPIR                                         #
  #                                                                         #
  # ----------------------------------------------------------------------- #

  curl -O http://mpir.org/mpir-$MPIRVER.tar.bz2
  tar xf mpir-$MPIRVER.tar.bz2
  cd mpir-$MPIRVER
    ./configure --prefix=$PREFIX\
      --with-yasm=$PREFIX/bin/yasm\
      --enable-fat\
      --enable-shared=yes\
      --enable-static=no\
      --enable-gmpcompat
    make -j3
    make install
  cd ..

  FLINTARB_WITHGMP="--with-mpir=$PREFIX"

fi

# ------------------------------------------------------------------------- #
#                                                                           #
#                              MPFR                                         #
#                                                                           #
# ------------------------------------------------------------------------- #

curl -O https://ftp.gnu.org/gnu/mpfr/mpfr-$MPFRVER.tar.gz
tar xf mpfr-$MPFRVER.tar.gz
cd mpfr-$MPFRVER
  ./configure --prefix=$PREFIX\
    --with-gmp=$PREFIX\
    --enable-shared=yes\
    --enable-static=no
  make -j3
  make install
cd ..

# ------------------------------------------------------------------------- #
#                                                                           #
#                              FLINT                                        #
#                                                                           #
# ------------------------------------------------------------------------- #

curl -O -L https://www.flintlib.org/flint-$FLINTVER.tar.gz
tar xf flint-$FLINTVER.tar.gz
cd flint-$FLINTVER
  ./configure --prefix=$PREFIX\
    $FLINTARB_WITHGMP\
    --with-mpfr=$PREFIX\
    --disable-static
  make -j3
  make install
cd ..

# ------------------------------------------------------------------------- #
#                                                                           #
#                               ARB                                         #
#                                                                           #
# ------------------------------------------------------------------------- #

curl -O -L https://github.com/fredrik-johansson/arb/archive/refs/tags/$ARBVER.tar.gz
mv $ARBVER.tar.gz arb-$ARBVER.tar.gz
tar xf arb-$ARBVER.tar.gz
cd arb-$ARBVER
  ./configure --prefix=$PREFIX\
    --with-flint=$PREFIX\
    $FLINTARB_WITHGMP\
    --with-mpfr=$PREFIX\
    --disable-static
  make -j3
  make install
  #
  # Here make check passes for Linux and OSX but fails for Windows probably
  # because of a linker error or something like that. It would be nice to
  # enable this check when it can work for Windows but for now we disable it
  # because if it fails then we don't get any wheels built.
  #
  # ARB_TEST_MULTIPLIER=0.1 make check
cd ..

# ------------------------------------------------------------------------- #
#                                                                           #
#                              Done!                                        #
#                                                                           #
# ------------------------------------------------------------------------- #

echo
echo -----------------------------------------------------------------------
echo
echo Build dependencies for python-flint compiled as shared libraries in:
echo $PREFIX
echo
echo Versions:
if [[ $USE_GMP = "gmp" ]]; then
  echo GMP: $GMPVER
else
  echo MPIR: $MPIRVER
fi
echo MPFR: $MPFRVER
echo Flint: $FLINTVER
echo Arb: $ARBVER
echo
echo -----------------------------------------------------------------------
echo
