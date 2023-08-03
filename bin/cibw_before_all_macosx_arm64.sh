#!/usr/bin/env bash

export CPPFLAGS=" --target=arm64-apple-macos11"
export LDFLAGS=" -arch arm64"

bin/build_dependencies_unix.sh\
  --gmp gmp\
  --host aarch64-apple-darwin\
  --patch-gmp-arm64\
  --use-gmp-github-mirror
