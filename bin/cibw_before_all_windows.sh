#!/bin/bash

set -o errexit

if [ $# -lt 2 ]; then
    echo "usage: $0 <llvm-machine> <build-dependencies-args>..."
    exit 1
fi

llvm_machine="$1"
shift

bin/build_dependencies_unix.sh "$@"

mkdir -p .local/lib
cd .local/bin
for dll_file in libgmp-*.dll libmpfr-*.dll libflint*.dll
do
  lib_name=$(basename -s .dll "${dll_file}")
  def_file=${lib_name}.def
  name=$(echo "${lib_name}" | sed 's/^lib//;s/[-.][0-9].*$//')

  gendef "${dll_file}"
  llvm-lib /def:"${def_file}" /out:"../lib/${name}.lib" /machine:${llvm_machine} /nologo
  rm "${def_file}"
done
cd ../..
