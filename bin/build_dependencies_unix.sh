#!/usr/bin/env bash
#
# Build local installs of python-flint's dependencies. This should be run
# before attempting to build python-flint itself.

set -o errexit

# ------------------------------------------------------------------------- #
#                                                                           #
# Supported options:                                                        #
#                                                                           #
# --gmp gmp         - build based on GMP (default)                          #
# --gmp mpir        - build based on MPIR (no longer works)                 #
# --host <HOST>     - set the host (target) for GMP build                   #
# --patch-gmp-arm64 - apply patch to GMP for OSX arm64                      #
#                                                                           #
# ------------------------------------------------------------------------- #

SKIP_GMP=no
SKIP_MPFR=no

USE_GMP=gmp
PATCH_GMP_ARM64=no
BUILD_ARB=no

while [[ $# -gt 0 ]]
do
  key="$1"
  case $key in
    -h|--help)
      echo "bin/download_dependencies.sh [options]"
      echo
      echo "Build local installs of python-flint's dependencies."
      echo
      echo "Supported options:"
      echo "  --help            - show this help message"
      echo "  --host <HOST>     - set the host (target) for GMP build"
      echo "  --skip-gmp        - skip building GMP"
      echo "  --skip-mpfr       - skip building MPFR"
      echo
      echo "Legacy options:"
      echo "  --gmp gmp         - build based on GMP (default)"
      echo "  --gmp mpir        - build based on MPIR (no longer works)"
      echo "  --patch-gmp-arm64 - apply patch to GMP 6.2.1 for OSX arm64"
      echo "  --arb             - build Arb (only needed for flint < 3.0.0)"
      echo
      exit
    ;;
    --host)
      # e.g. --host x86_64-unknown-linux-gnu
      # or   --host x86_64-apple-darwin
      HOST_ARG="$2"
      shift
      shift
    ;;
    --gmp)
      # e.g. --gmp gmp or --gmp mpir
      # The mpir build no longer works because the download fails.
      USE_GMP="$2"
      if [[ "$USE_GMP" != "gmp" && "$USE_GMP" != "mpir" ]]; then
        echo "--gmp option should be gmp or mpir"
        exit 1
      fi
      shift
      shift
    ;;
    --arb)
      # With flint >= 3.0.0 Arb is included so we do not need to build it
      # separately. Pass --arb if building for older versions of flint.
      BUILD_ARB=yes
      shift
    ;;
    --skip-gmp)
      # If you already have a local install of GMP you can pass --skip-gmp
      # to skip building it.
      SKIP_GMP=yes
      shift
    ;;
    --skip-mpfr)
      # If you already have a local install of MPFR you can pass --skip-mpfr
      # to skip building it.
      SKIP_MPFR=yes
      shift
    ;;
    --patch-gmp-arm64)
      # Needed only for GMP 6.2.1 on OSX arm64 (Apple M1) hardware
      # As of GMP 6.3.0 this patch is no longer needed
      PATCH_GMP_ARM64=yes
      shift
    ;;
    --use-gmp-github-mirror)
      USE_GMP_GITHUB_MIRROR=yes
      shift
    ;;
  *)
    2>&1 echo "unrecognised argument:" $key
    exit 1
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

  if [ $SKIP_GMP = "yes" ]; then
    echo
    echo --------------------------------------------
    echo "           skipping GMP"
    echo --------------------------------------------
    echo
  else
    echo
    echo --------------------------------------------
    echo "           building GMP"
    echo --------------------------------------------
    echo

    if [ $USE_GMP_GITHUB_MIRROR = "yes" ]; then
      # Needed in GitHub Actions because it is blocked from gmplib.org
      git clone https://github.com/oscarbenjamin/gmp_mirror.git
      cp gmp_mirror/gmp-$GMPVER.tar.xz .
    else
      curl -O https://gmplib.org/download/gmp/gmp-$GMPVER.tar.xz
    fi

    tar xf gmp-$GMPVER.tar.xz
    cd gmp-$GMPVER

      #
      # See https://github.com/aleaxit/gmpy/issues/350
      #
      # We need to patch GMP for OSX arm64 (Apple M1) hardware. This patch is
      # from the GMP repo but was applied after the release of GMP 6.2.1.
      # This patch is no longer needed for GMP 6.3.0.
      #
      if [ $PATCH_GMP_ARM64 = "yes" ]; then
        echo
        echo --------------------------------------------
        echo "           patching GMP"
        echo --------------------------------------------
        patch -N -Z -p0 < ../../../bin/patch-arm64.diff
      fi

      # Show the output of configfsf.guess
      chmod +x configfsf.guess
      ./configfsf.guess

      ./configure --prefix=$PREFIX\
        --enable-fat\
        --enable-shared=yes\
        --enable-static=no\
        --host=$HOST_ARG
      make -j6
      make install

    cd ..

  fi

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
    make -j6
    make install
  cd ..

  # ----------------------------------------------------------------------- #
  #                                                                         #
  #                            MPIR                                         #
  #                                                                         #
  # ----------------------------------------------------------------------- #

  #
  # The mpir.org domain has expired and no longer hosts the source code so the
  # call to curl below will fail.
  # We could try to download from https://github.com/wbhart/mpir/releases.
  #
  # Ultimately it seems that MPIR is no longer maintained though so for now
  # this remains unfixed.
  #

  >&2 echo "MPIR build of python_flint is no longer supported"
  exit 1

  curl -O http://mpir.org/mpir-$MPIRVER.tar.bz2
  tar xf mpir-$MPIRVER.tar.bz2
  cd mpir-$MPIRVER
    ./configure --prefix=$PREFIX\
      --with-yasm=$PREFIX/bin/yasm\
      --enable-fat\
      --enable-shared=yes\
      --enable-static=no\
      --enable-gmpcompat
    make -j6
    make install
  cd ..

  FLINTARB_WITHGMP="--with-mpir=$PREFIX"

fi

# ------------------------------------------------------------------------- #
#                                                                           #
#                              MPFR                                         #
#                                                                           #
# ------------------------------------------------------------------------- #

if [ $SKIP_MPFR = "yes" ]; then
  echo
  echo --------------------------------------------
  echo "           skipping MPFR"
  echo --------------------------------------------
  echo
else
  echo
  echo --------------------------------------------
  echo "           building MPFR"
  echo --------------------------------------------
  echo

  curl -O https://ftp.gnu.org/gnu/mpfr/mpfr-$MPFRVER.tar.gz
  tar xf mpfr-$MPFRVER.tar.gz
  cd mpfr-$MPFRVER
    ./configure --prefix=$PREFIX\
      --host=$HOST_ARG\
      --with-gmp=$PREFIX\
      --enable-shared=yes\
      --enable-static=no
    make -j6
    make install
  cd ..
fi

# ------------------------------------------------------------------------- #
#                                                                           #
#                              FLINT                                        #
#                                                                           #
# ------------------------------------------------------------------------- #

echo
echo --------------------------------------------
echo "           building Flint"
echo --------------------------------------------
echo

curl -O -L https://github.com/flintlib/flint/releases/download/v$FLINTVER/flint-$FLINTVER.tar.gz
tar xf flint-$FLINTVER.tar.gz
cd flint-$FLINTVER
  ./bootstrap.sh
  ./configure --prefix=$PREFIX\
    --host=$HOST_ARG\
    $FLINTARB_WITHGMP\
    --with-mpfr=$PREFIX\
    --disable-static\
    --disable-debug
  make -j6
  make install
cd ..

# ------------------------------------------------------------------------- #
#                                                                           #
#                               ARB                                         #
#                                                                           #
# ------------------------------------------------------------------------- #

if [ $BUILD_ARB = "yes" ]; then

  echo
  echo --------------------------------------------
  echo "           building Arb"
  echo --------------------------------------------
  echo

  curl -O -L https://github.com/fredrik-johansson/arb/archive/refs/tags/$ARBVER.tar.gz
  mv $ARBVER.tar.gz arb-$ARBVER.tar.gz
  tar xf arb-$ARBVER.tar.gz
  cd arb-$ARBVER
    ./configure --prefix=$PREFIX\
      --with-flint=$PREFIX\
      $FLINTARB_WITHGMP\
      --with-mpfr=$PREFIX\
      --disable-static
    make -j6
    make install
    #
    # Set PATH so that DLLs are picked up on Windows.
    #
    PATH=$PATH:$PREFIX/lib:$PREFIX/bin \
    ARB_TEST_MULTIPLIER=0.1            \
    # Skip Arb tests now because they are slow.
    # make check
  cd ..
fi

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

if [ $SKIP_GMP = "yes" ]; then
  echo GMP: skipped
else
  if [[ $USE_GMP = "gmp" ]]; then
    echo GMP: $GMPVER
  else
    echo MPIR: $MPIRVER
  fi
fi

if [ $SKIP_MPFR = "yes" ]; then
  echo MPFR: skipped
else
  echo MPFR: $MPFRVER
fi

echo Flint: $FLINTVER

if [ $BUILD_ARB = "yes" ]; then
  echo Arb: $ARBVER
fi
echo
echo -----------------------------------------------------------------------
echo
