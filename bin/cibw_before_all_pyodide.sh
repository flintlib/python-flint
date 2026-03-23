#!/bin/bash

set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
PROJECT_DIR=$(cd -- "$SCRIPT_DIR/.." && pwd)
WASM_LIBRARY_DIR=${WASM_LIBRARY_DIR:-"$PROJECT_DIR/wasm-library-dir"}

if [ -f "$WASM_LIBRARY_DIR/lib/pkgconfig/flint.pc" ]; then
  echo "Using cached Pyodide dependencies from $WASM_LIBRARY_DIR"
  exit 0
fi

mkdir -p "$WASM_LIBRARY_DIR"

export CFLAGS="-fPIC ${CFLAGS:-}"

"$SCRIPT_DIR/pyodide_build_dependencies.sh" --wasm-library-dir "$WASM_LIBRARY_DIR"
