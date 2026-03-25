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

# The Pyodide cibuildwheel environment for building python-flint itself sets
# include/library flags that reference libgmp/libmpfr/libflint. Those libraries
# do not exist yet while bootstrapping the static dependencies here, so letting
# them leak into configure would break link tests such as "does the C compiler
# work?" for GMP. Use a minimal environment for the bootstrap step.
unset PKG_CONFIG_PATH
unset LDFLAGS
export CFLAGS="-fPIC"

"$SCRIPT_DIR/pyodide_build_dependencies.sh" --wasm-library-dir "$WASM_LIBRARY_DIR"
