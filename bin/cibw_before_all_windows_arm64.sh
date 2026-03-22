#!/bin/bash

set -o errexit

pacman -S --noconfirm \
    mingw-w64-clang-aarch64-toolchain\
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
    --disable-fat\
    --use-gmp-github-mirror\
    --host aarch64-pc-windows-gnullvm\
    --patch-C23\
    --patch-ldd\
    #

PATH="$PATH:$(find "/c/Program Files/Microsoft Visual Studio/2022/" -name "HostARM64")/arm64/"

mkdir -p .local/lib
cd .local/bin
for dll_file in libgmp-*.dll libmpfr-*.dll libflint*.dll
do
  lib_name=$(basename -s .dll ${dll_file})
  exports_file=${lib_name}-exports.txt
  def_file=${lib_name}.def
  lib_file=${lib_name}.lib
  name=$(echo ${lib_name}|sed 's/^lib//;s/[-.][0-9].*$//')

  dumpbin //exports ${dll_file} > ${exports_file}

  echo LIBRARY ${lib_name} > ${def_file}
  echo EXPORTS >> ${def_file}
  awk 'NR>19 && $4 != "" {print $4 " @"$1}' ${exports_file} >> ${def_file}
  sed -i 's/$/\r/' ${def_file}

  lib //def:${def_file} //out:${lib_file} //machine:arm64
  rm ${exports_file} ${def_file} ${lib_name}.exp
  mv ${lib_file} ../lib/${name}.lib
done
cd ../..
