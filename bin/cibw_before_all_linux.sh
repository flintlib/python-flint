#!/usr/bin/env bash

#yum install -y xz
apt-get install -y xz-utils

bin/build_dependencies_unix.sh\
  --gmp gmp\
  --host x86_64-pc-linux-gnu\
  --use-gmp-github-mirror
