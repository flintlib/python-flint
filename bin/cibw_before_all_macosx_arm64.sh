#!/usr/bin/env bash

export CPPFLAGS=" --target=arm64-apple-macos11"
export LDFLAGS=" -arch arm64"

brew install automake

bin/build_dependencies_unix.sh\
  --gmp gmp\
  --host aarch64-apple-darwin\
  --use-gmp-github-mirror
