#!/usr/bin/env bash
#
# Build local installs of python-flint's dependencies. This should be run
# before attempting to build python-flint itself.

set -o errexit

SKIP_GMP=no
SKIP_MPFR=no
PATCH_GMP_C23=no
PATCH_LDD=no
PATCH_IMMINTRIN=no
GMP_FAT_ARG="--enable-fat"
GMP_ASSEMBLY_ARG=
HOST_ARG=

USE_GMP_GITHUB_MIRROR=no

while [[ $# -gt 0 ]]
do
  key="$1"
  case $key in
    -h|--help)
      echo "bin/build_dependencies_unix.sh [options]"
      echo
      echo "Build local installs of python-flint's dependencies."
      echo
      echo "Supported options:"
      echo "  --help            - show this help message"
      echo "  --host <HOST>     - set the host (target) for GMP build"
      echo "  --skip-gmp        - skip building GMP"
      echo "  --skip-mpfr       - skip building MPFR"
      echo "  --disable-assembly - disable GMP assembly routines"
      echo "  --patch-C23       - apply patch to GMP 6.3.0 for C23 compatibility"
      echo "  --patch-ldd       - patch flint shared linking for mingw on arm64"
      echo "  --patch-immintrin - patch flint arm64 msvc header to avoid immintrin.h"
      echo
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
    --disable-assembly)
      # GMP does not allow --enable-fat together with --disable-assembly.
      GMP_FAT_ARG=
      GMP_ASSEMBLY_ARG="--disable-assembly"
      shift
    ;;
    --patch-C23)
      # Patch GMP 6.3.0 for newer gcc versions
      PATCH_GMP_C23=yes
      shift
    ;;
    --patch-ldd)
      # Needed only for the FLINT shared build on mingw arm64.
      PATCH_LDD=yes
      shift
    ;;
    --patch-immintrin)
      # Needed only for the FLINT headers consumed by MSVC on Windows arm64.
      PATCH_IMMINTRIN=yes
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

mkdir -p "$PREFIX"
cd $PREFIX
mkdir -p src
cd src

# ------------------------------------------------------------------------- #
#                                                                           #
# Now build all dependencies.                                               #
#                                                                           #
# ------------------------------------------------------------------------- #

if [ "$SKIP_GMP" = "yes" ]; then
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

    if [ "$USE_GMP_GITHUB_MIRROR" = "yes" ]; then
      # Needed in GitHub Actions because it is blocked from gmplib.org
      git clone https://github.com/oscarbenjamin/gmp_mirror.git
      cp gmp_mirror/gmp-$GMPVER.tar.xz .
    else
      curl -O https://gmplib.org/download/gmp/gmp-$GMPVER.tar.xz
    fi

    tar xf gmp-$GMPVER.tar.xz
    cd gmp-$GMPVER

      #
      # https://github.com/msys2/MSYS2-packages/issues/5499
      #
      # This patch needed for GMP 6.3.0 building with msys2 or probably just
      # newer gcc versions.
      #
      if [ $PATCH_GMP_C23 = "yes" ]; then
        echo
        echo --------------------------------------------
        echo "           patching GMP"
        echo --------------------------------------------
        patch -N -Z < ../../../bin/patch-C23.diff
        autoreconf -fi
      fi

      # Show the output of configfsf.guess
      chmod +x configfsf.guess
      ./configfsf.guess

      ./configure --prefix=$PREFIX\
        $GMP_FAT_ARG\
        $GMP_ASSEMBLY_ARG\
        --enable-shared=yes\
        --enable-static=no\
        --host=$HOST_ARG
      make -j6
      make install

    cd ..

fi

# ------------------------------------------------------------------------- #
#                                                                           #
#                              MPFR                                         #
#                                                                           #
# ------------------------------------------------------------------------- #

if [ "$SKIP_MPFR" = "yes" ]; then
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

  if [ $USE_GMP_GITHUB_MIRROR = "yes" ]; then
    if [ ! -d "gmp_mirror" ] ; then
      git clone https://github.com/oscarbenjamin/gmp_mirror.git
    fi
    cp gmp_mirror/mpfr-$MPFRVER.tar.gz .
  else
    curl -O https://ftp.gnu.org/gnu/mpfr/mpfr-$MPFRVER.tar.gz
  fi

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
  if [ "$PATCH_LDD" = "yes" ]; then
    echo
    echo --------------------------------------------
    echo "           patching FLINT"
    echo --------------------------------------------
    patch -N -Z -p1 < ../../../bin/patch-flint-windows-arm64-link-$FLINTVER.diff
  fi
  if [ "$PATCH_IMMINTRIN" = "yes" ]; then
    echo
    echo --------------------------------------------
    echo "           patching FLINT"
    echo --------------------------------------------
    patch -N -Z -p1 < ../../../bin/patch-flint-windows-arm64-immintrin.diff
  fi
  ./bootstrap.sh
  ./configure --prefix=$PREFIX\
    --host=$HOST_ARG\
    --with-gmp=$PREFIX\
    --with-mpfr=$PREFIX\
    --disable-static\
    --disable-debug
  make -j6
  make install
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

if [ "$SKIP_GMP" = "yes" ]; then
  echo GMP: skipped
else
  echo GMP: $GMPVER
fi

if [ "$SKIP_MPFR" = "yes" ]; then
  echo MPFR: skipped
else
  echo MPFR: $MPFRVER
fi

echo Flint: $FLINTVER
echo
echo -----------------------------------------------------------------------
echo
