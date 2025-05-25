#!/usr/bin/env bash

set -e

cat wheels/LICENSE_macos_wheels.txt >> LICENSE

export CPPFLAGS=" --target=arm64-apple-macos11"
export LDFLAGS=" -arch arm64"

brew install automake libtool

bin/build_dependencies_unix.sh\
  --gmp gmp\
  --host aarch64-apple-darwin\
  --use-gmp-github-mirror
