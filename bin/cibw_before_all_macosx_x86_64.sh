#!/usr/bin/env bash

set -e

brew install automake libtool

bin/build_dependencies_unix.sh\
  --host x86_64-apple-darwin\
  --use-gmp-github-mirror
