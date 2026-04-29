#!/usr/bin/env bash

apt-get install xz-utils

bin/build_dependencies_unix.sh\
  --gmp gmp\
  --host x86_64-pc-linux-gnu\
  --use-gmp-github-mirror
