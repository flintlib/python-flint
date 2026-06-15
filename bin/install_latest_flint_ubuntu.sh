#!/usr/bin/env bash

set -e

source bin/build_variables.sh

sudo apt-get update
sudo apt-get install libgmp-dev libmpfr-dev xz-utils ninja-build

curl -O -L https://github.com/flintlib/flint/releases/download/v$FLINTVER/flint-$FLINTVER.tar.gz
tar -xzf flint-$FLINTVER.tar.gz
cd flint-$FLINTVER && ./configure --disable-static && make -j$(expr $(nproc) + 1) && sudo make install

ls -l /usr/local/lib
sudo ldconfig /usr/local/lib
