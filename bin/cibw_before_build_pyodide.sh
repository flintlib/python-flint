#!/bin/bash

set -euo pipefail

python -m pip install wheel auditwheel_emscripten

EMSCRIPTEN_DIR="$HOME/.cache/cibuildwheel/emsdk-$PYODIDE_EMSCRIPTEN_VERSION/emsdk-$PYODIDE_EMSCRIPTEN_VERSION/upstream/emscripten"

for patch_file in pyodide-patches/emsdk/patches/*.patch; do
  echo "Applying Pyodide Emscripten patch $(basename "$patch_file")"
  patch -p1 --verbose -d "$EMSCRIPTEN_DIR" < "$patch_file"
done
