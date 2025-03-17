#!/usr/bin/env bash

set -e

cat wheels/LICENSE_macos_wheels.txt >> LICENSE

brew install automake libtool

bin/build_dependencies_unix.sh\
  --gmp gmp\
  --host x86_64-apple-darwin\
  --use-gmp-github-mirror
