#!/usr/bin/env bash

yum install -y xz
bin/build_dependencies_unix.sh\
  --gmp mpir\
  --host x86_64-pc-linux-gnu
