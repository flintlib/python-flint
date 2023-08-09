#!/usr/bin/env bash

set -xe

#yum install -y xz

# The echo lines are needed for manylinux_2_24
# https://github.com/pypa/manylinux/issues/1369#issuecomment-1546594841
echo "deb http://archive.debian.org/debian stretch main" > /etc/apt/sources.list
echo "deb http://archive.debian.org/debian-security stretch/updates main" >> /etc/apt/sources.list
apt-get update
apt-get install -y xz-utils

bin/build_dependencies_unix.sh\
  --gmp gmp\
  --host x86_64-pc-linux-gnu\
  --use-gmp-github-mirror
