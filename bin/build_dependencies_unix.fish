#!/usr/bin/env fish
#
# Build local installs of python-flint's dependencies. This should be run
# before attempting to build python-flint itself.

# There is no equivalent in fish AFAIK
# set -o errexit

# ------------------------------------------------------------------------- #
#                                                                           #
# Supported options:                                                        #
#                                                                           #
# --jobs <NUM>      - run NUM jobs in parallel [default=3]                  #
# --host <HOST>     - set the host (target) for GMP build                   #
# --patch-gmp-arm64 - apply patch to GMP for OSX arm64                      #
# --gmp gmp         - build based on GMP (default)                          #
# --gmp mpir        - build based on MPIR (no longer works)                 #
# ------------------------------------------------------------------------- #

set NUM_CORES 3
set SKIP_GMP no
set SKIP_MPFR no

set USE_GMP gmp
set PATCH_GMP_ARM64 no
set BUILD_ARB no

while test (count $argv) -gt 0
    set key $argv[1]
    echo "key: $key"
    switch "$key"
        case -h --help
            echo "bin/download_dependencies.sh [options]"
            echo
            echo "Build local installs of python-flint's dependencies."
            echo
            echo "Supported options:"
            echo "  --help            - show this help message"
            echo "  --cores <NUM>     - set number of cores to use (default: 3)"
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
        case -j --jobs --cores
            # e.g. --jobs 4
            # or   --jobs 16
            set NUM_CORES "$argv[2]"
            set argv $argv[3..-1]
        case --host
            # e.g. --host x86_64-unknown-linux-gnu
            # or   --host x86_64-apple-darwin
            set HOST_ARG "$argv[2]"
            set argv $argv[3..-1]
        case --gmp
            # e.g. --gmp gmp or --gmp mpir
            # The mpir build no longer works because the download fails.
            set USE_GMP "$argv[2]"
            if test "$USE_GMP" != "gmp" -a "$USE_GMP" != "mpir"
                echo "--gmp option should be gmp or mpir"
                exit 1
            end
            set argv $argv[3..-1]
        case --arb
            # With flint >= 3.0.0 Arb is included so we do not need to build it
            # separately. Pass --arb if building for older versions of flint.
            set BUILD_ARB yes
            set argv $argv[2..-1]
        case --skip-gmp
            # If you already have a local install of GMP you can pass --skip-gmp
            # to skip building it.
            set SKIP_GMP yes
            set argv $argv[2..-1]
        case --skip-mpfr
            # If you already have a local install of MPFR you can pass --skip-mpfr
            # to skip building it.
            set SKIP_MPFR yes
            set argv $argv[2..-1]
        case --patch-gmp-arm64
            # Needed only for GMP 6.2.1 on OSX arm64 (Apple M1) hardware
            # As of GMP 6.3.0 this patch is no longer needed
            set PATCH_GMP_ARM64 yes
            set argv $argv[2..-1]
        case --use-gmp-github-mirror
            set argv $argv[2..-1]
        case "*"
            echo "unrecognised argument:" $key
            exit 1
    end
end

# ------------------------------------------------------------------------- #
#                                                                           #
# The build_variables.sh script sets variables specifying the versions to   #
# use for all dependencies and also the PREFIX variable.                    #
#                                                                           #
# ------------------------------------------------------------------------- #

source bin/build_variables.fish

echo "=== Parameters ==="
echo "NUM_CORES: $NUM_CORES"
echo "PREFIX: $PREFIX"

cd $PREFIX
mkdir -p src
cd src

# ------------------------------------------------------------------------- #
#                                                                           #
# Now build all dependencies.                                               #
#                                                                           #
# ------------------------------------------------------------------------- #

if test "$USE_GMP" = "gmp"

    # ----------------------------------------------------------------------- #
    #                                                                         #
    #                            GMP                                          #
    #                                                                         #
    # ----------------------------------------------------------------------- #

    if test "$SKIP_GMP" = "yes"
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

        if test "$USE_GMP_GITHUB_MIRROR" = "yes"
            # Needed in GitHub Actions because it is blocked from gmplib.org
            git clone https://github.com/oscarbenjamin/gmp_mirror.git
            cp gmp_mirror/gmp-$GMPVER.tar.xz .
        else
            curl -O https://gmplib.org/download/gmp/gmp-$GMPVER.tar.xz
        end

        tar xf gmp-$GMPVER.tar.xz
        cd gmp-$GMPVER

            #
            # See https://github.com/aleaxit/gmpy/issues/350
            #
            # We need to patch GMP for OSX arm64 (Apple M1) hardware. This patch is
            # from the GMP repo but was applied after the release of GMP 6.2.1.
            # This patch is no longer needed for GMP 6.3.0.
            #
            if test "$PATCH_GMP_ARM64" = "yes"
                echo
                echo --------------------------------------------
                echo "           patching GMP"
                echo --------------------------------------------
                patch -N -Z -p0 < ../../../bin/patch-arm64.diff
            end

            # Show the output of configfsf.guess
            chmod +x configfsf.guess
            ./configfsf.guess

            ./configure --prefix=$PREFIX \
                --enable-fat             \
                --enable-shared=yes      \
                --enable-static=no       \
                --host=$HOST_ARG
            make -j"$NUM_CORES"
            make install

        cd ..

    end

    set FLINTARB_WITHGMP "-DGMP_INCLUDE_DIRS=$PREFIX"

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
        make -j"$NUM_CORES"
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

    echo "MPIR build of python_flint is no longer supported"
    exit 1

    # curl -O http://mpir.org/mpir-$MPIRVER.tar.bz2
    # tar xf mpir-$MPIRVER.tar.bz2
    # cd mpir-$MPIRVER
    #     ./configure --prefix=$PREFIX     \
    #         --with-yasm=$PREFIX/bin/yasm \
    #         --enable-fat                 \
    #         --enable-shared=yes          \
    #         --enable-static=no           \
    #         --enable-gmpcompat
    #     make -j"$NUM_CORES"
    #     make install
    # cd ..
    #
    # set FLINTARB_WITHGMP "--with-mpir=$PREFIX"

end

# ------------------------------------------------------------------------- #
#                                                                           #
#                              MPFR                                         #
#                                                                           #
# ------------------------------------------------------------------------- #

if test "$SKIP_MPFR" = "yes"
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
        ./configure --prefix=$PREFIX \
            --with-gmp=$PREFIX       \
            --enable-shared=yes      \
            --enable-static=no
        make -j"$NUM_CORES"
        make install
    cd ..
end

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

# Use Ninja if it is installed (fish users exclusive feature)
echo -n "Checking if ninja is available..."
if command -q ninja
    set NINJA_FLAG "-GNinja"
    echo "Yes"
else
    echo "No"
end

# Use mold if it is installed (fish users exclusive feature)
# Please make sure the compiler versions are up to date
echo -n "Checking if mold is available..."
if command -q mold
    set MOLD_C_FLAGS "-DCMAKE_C_FLAGS=\"-fuse-ld=mold\""
    set MOLD_CXX_FLAGS "-DCMAKE_CXX_FLAGS=\"-fuse-ld=mold\""
    echo "Yes"
else
    echo "No"
end

curl -O -L https://www.flintlib.org/flint-$FLINTVER.tar.gz
tar xf flint-$FLINTVER.tar.gz
cd flint-$FLINTVER
    # I can't get ./bootstrap to work with fish
    mkdir build
    cd build
    # I don't know how to translate --disable-static
    cmake ..                                \
        -DCMAKE_INSTALL_PREFIX=$PREFIX      \
        -DMPFR_INCLUDE_DIRS=$PREFIX/include \
        $FLINTARB_WITHGMP                   \
        $MOLD_C_FLAGS                       \
        $MOLD_CXX_FLAGS                     \
        $NINJA_FLAG
    cmake --build . -j "$NUM_CORES"
    cmake --build . --target install
cd ..

# ------------------------------------------------------------------------- #
#                                                                           #
#                               ARB                                         #
#                                                                           #
# ------------------------------------------------------------------------- #

if test "$BUILD_ARB" = "yes"

    echo
    echo --------------------------------------------
    echo "           building Arb"
    echo --------------------------------------------
    echo

    curl -O -L https://github.com/fredrik-johansson/arb/archive/refs/tags/$ARBVER.tar.gz
    mv $ARBVER.tar.gz arb-$ARBVER.tar.gz
    tar xf arb-$ARBVER.tar.gz
    cd arb-$ARBVER
        ./configure --prefix=$PREFIX \
            --with-flint=$PREFIX     \
            $FLINTARB_WITHGMP        \
            --with-mpfr=$PREFIX      \
            --disable-static
        make -j"$NUM_CORES"
        make install
        #
        # Set PATH so that DLLs are picked up on Windows.
        #
        fish_add_path $PREFIX/lib
        fish_add_path $PREFIX/bin
        set ARB_TEST_MULTIPLIER 0.1
        # Skip Arb tests now because they are slow.
        # make check
    cd ..
end

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

if test "$SKIP_GMP" = "yes"
    echo "GMP: skipped"
else if test "$USE_GMP" = "gmp"
    echo "GMP: $GMPVER"
else
    echo "MPIR: $MPIRVER"
end

if test "$SKIP_MPFR" = "yes"
    echo "MPFR: skipped"
else
    echo "MPFR: $MPFRVER"
end

echo "Flint: $FLINTVER"

if test "$BUILD_ARB" = "yes"
    echo "Arb: $ARBVER"
end
echo
echo -----------------------------------------------------------------------
echo
