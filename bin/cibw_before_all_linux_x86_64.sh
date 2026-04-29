#!/usr/bin/env bash

set -e

yum install -y xz
bin/build_dependencies_unix.sh\
  --gmp gmp\
  --host x86_64-pc-linux-gnu\
  --use-gmp-github-mirror
