#!/usr/bin/env bash

yum install -y xz
bin/build_dependencies_unix.sh\
  --gmp gmp\
  --host x86_64-unknown-linux-gnu
