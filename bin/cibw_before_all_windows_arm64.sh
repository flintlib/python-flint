#!/bin/bash

set -o errexit

pacman -S --noconfirm \
    mingw-w64-clang-aarch64-toolchain\
    mingw-w64-clang-aarch64-llvm-tools\
    mingw-w64-clang-aarch64-tools-git\
    m4\
    make\
    base-devel\
    autoconf-wrapper\
    automake-wrapper\
    libtool\
    git\
    #

bin/build_dependencies_unix.sh \
    --use-gmp-github-mirror\
    --host aarch64-w64-mingw32\
    --patch-C23\
    --patch-ldd\
    --patch-immintrin\
    #

mkdir -p .local/lib
cd .local/bin
for dll_file in libgmp-*.dll libmpfr-*.dll libflint*.dll
do
  lib_name=$(basename -s .dll "${dll_file}")
  def_file=${lib_name}.def
  name=$(echo "${lib_name}" | sed 's/^lib//;s/[-.][0-9].*$//')

  gendef "${dll_file}"
  llvm-lib /def:"${def_file}" /out:"../lib/${name}.lib" /machine:arm64 /nologo
  rm "${def_file}"
done
cd ../..
