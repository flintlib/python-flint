#!/usr/bin/env bash

yum install -y xz
bin/build_dependencies_unix.sh\
  --gmp gmp\
  --host aarch64-pc-linux-gnu\
  --use-gmp-github-mirror
