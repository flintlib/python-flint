#!/usr/bin/env bash

set -e

source bin/build_variables.sh

sudo apt-get update
sudo apt-get install libgmp-dev libmpfr-dev xz-utils ninja-build

git clone https://github.com/flintlib/flint
cd flint && ./bootstrap.sh && ./configure --disable-static && make -j$(expr $(nproc) + 1) && sudo make install

ls -l /usr/local/lib
sudo ldconfig /usr/local/lib
