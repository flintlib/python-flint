#!/usr/bin/env bash

set -e

cat wheels/LICENSE_linux_wheels.txt >> LICENSE

yum install -y xz
bin/build_dependencies_unix.sh\
  --gmp gmp\
  --host aarch64-pc-linux-gnu\
  --use-gmp-github-mirror
