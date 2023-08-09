#!/usr/bin/env bash

set -xe

#yum install -y xz
sudo apt-get update
sudo apt-get install -y xz-utils

bin/build_dependencies_unix.sh\
  --gmp gmp\
  --host x86_64-pc-linux-gnu\
  --use-gmp-github-mirror
