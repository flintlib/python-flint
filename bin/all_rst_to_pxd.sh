#!/usr/bin/env bash

FLINT_DOC_DIR=$1

set -e

modules="fmpz"

for module in $modules; do
    echo "Processing $module"
    bin/rst_to_pxd.py flint/$module --flint-doc-dir=$FLINT_DOC_DIR > src/flint/flintlib/$module.pxd
done
