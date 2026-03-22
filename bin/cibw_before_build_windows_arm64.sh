#!/bin/bash

set -o errexit

env_file=.local/cibw_before_build_windows_arm64.env
if [ ! -f "$env_file" ]; then
    echo "Could not find $env_file" >&2
    exit 1
fi

# shellcheck disable=SC1090
. "$env_file"

lib_dir_msys="$LIB_DIR"
pkgconfig_dir_msys="$PKGCONFIG_DIR"
dll_path_msys="$DLL_PATH"
if command -v cygpath >/dev/null 2>&1; then
    lib_dir_msys=$(cygpath -u "$LIB_DIR")
    pkgconfig_dir_msys=$(cygpath -u "$PKGCONFIG_DIR")
    dll_path_msys=$(cygpath -u "$DLL_PATH")
fi
mkdir -p "$pkgconfig_dir_msys"

if ! command -v gendef >/dev/null 2>&1; then
    echo "Could not find gendef on PATH" >&2
    exit 1
fi

dlltool=
if command -v cc >/dev/null 2>&1; then
    dlltool=$(cc -print-prog-name=dlltool)
fi
if [ -n "${dlltool:-}" ] && [ -x "$dlltool" ]; then
    :
elif [ -n "${dlltool:-}" ] && command -v "$dlltool" >/dev/null 2>&1; then
    dlltool=$(command -v "$dlltool")
elif command -v dlltool >/dev/null 2>&1; then
    dlltool=$(command -v dlltool)
elif command -v llvm-dlltool >/dev/null 2>&1; then
    dlltool=$(command -v llvm-dlltool)
else
    echo "Could not find dlltool or llvm-dlltool on PATH" >&2
    exit 1
fi

dll_stem=${DLL_NAME%.dll}
def_path_msys="$lib_dir_msys/$dll_stem.def"
import_lib_msys="$lib_dir_msys/lib$dll_stem.dll.a"
pc_path_msys="$pkgconfig_dir_msys/python-$PKG_VERSION.pc"

rm -f "$def_path_msys" "$import_lib_msys"
(
    cd "$lib_dir_msys"
    gendef "$dll_path_msys"
)
echo "Generating $import_lib_msys with $dlltool"
"$dlltool" -d "$def_path_msys" -D "$DLL_NAME" -l "$import_lib_msys"
if [ ! -f "$import_lib_msys" ]; then
    echo "Failed to create $import_lib_msys" >&2
    exit 1
fi

cat > "$pc_path_msys" <<EOF
prefix=$REPO_ROOT
exec_prefix=\${prefix}
libdir=$LIB_DIR
includedir=$INCLUDE_DIR

Name: Python
Description: CPython import library for MinGW wheel builds
Version: $PKG_VERSION
Libs: $LIB_DIR/lib$dll_stem.dll.a
Cflags: -DMS_WIN64 -I\${includedir}
EOF

echo "Generated $import_lib_msys"
echo "Generated $pc_path_msys"
