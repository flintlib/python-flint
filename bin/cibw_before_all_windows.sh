#!/bin/bash

set -o errexit

pacman -S --noconfirm \
    mingw-w64-ucrt-x86_64-gcc\
    mingw-w64-ucrt-x86_64-tools-git\
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
    --patch-C23\
    #

# Assumes the standard GitHub Actions Windows 2022 runner layout.
PATH="$PATH:$(find "/c/Program Files/Microsoft Visual Studio/2022/" -name "Hostx86")/x64/"

if [ "${RUNNER_ARCH}" = "ARM64" ]
then
  msvc_machine=arm64
else
  msvc_machine=x64
fi

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

  lib //def:${def_file} //out:${lib_file} //machine:${msvc_machine}
  rm ${exports_file} ${def_file} ${lib_name}.exp
  mv ${lib_file} ../lib/${name}.lib
done
cd ../..
