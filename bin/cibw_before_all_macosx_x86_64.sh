#!/usr/bin/env bash

brew install automake libtool

bin/build_dependencies_unix.sh\
  --gmp gmp\
  --host x86_64-apple-darwin\
  --use-gmp-github-mirror
