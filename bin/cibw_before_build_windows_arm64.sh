#!/bin/bash

set -o errexit

repo_root=$(python -c 'from pathlib import Path; print(Path.cwd().resolve().as_posix())')
lib_dir="$repo_root/.local/lib"
pkgconfig_dir="$lib_dir/pkgconfig"
mkdir -p "$pkgconfig_dir"

dll_name=$(python -c 'import sysconfig; print(sysconfig.get_config_var("LDLIBRARY") or "")')
pkg_version=$(python -c 'import sysconfig; print(sysconfig.get_config_var("LDVERSION") or sysconfig.get_python_version())')
include_dir=$(python -c 'import sysconfig; print(sysconfig.get_config_var("INCLUDEPY") or "")')
base_prefix=$(python -c 'import sys; print(getattr(sys, "base_prefix", sys.prefix))')
python_bin=$(python -c 'import os, sys; print(os.path.dirname(sys.executable))')

include_dir=${include_dir//\\//}
base_prefix=${base_prefix//\\//}
python_bin=${python_bin//\\//}

if [ -z "$dll_name" ] || [ -z "$include_dir" ]; then
    echo "Could not determine Python DLL or include dir" >&2
    exit 1
fi

dll_path=""
for candidate in \
    "$python_bin/$dll_name" \
    "$base_prefix/$dll_name" \
    "$base_prefix/DLLs/$dll_name" \
    "$base_prefix/libs/$dll_name"
do
    if [ -f "$candidate" ]; then
        dll_path="$candidate"
        break
    fi
done

if [ -z "$dll_path" ]; then
    echo "Could not find $dll_name" >&2
    exit 1
fi

if ! command -v gendef >/dev/null 2>&1; then
    echo "Could not find gendef on PATH" >&2
    exit 1
fi

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

dll_stem=${dll_name%.dll}
def_path="$lib_dir/$dll_stem.def"
import_lib="$lib_dir/lib$dll_stem.dll.a"
pc_path="$pkgconfig_dir/python-$pkg_version.pc"

rm -f "$def_path" "$import_lib"
(
    cd "$lib_dir"
    gendef "$dll_path"
)
"$dlltool" -d "$def_path" -D "$dll_name" -l "$import_lib"

printf 'prefix=%s\n' "$repo_root" > "$pc_path"
printf 'exec_prefix=${prefix}\n' >> "$pc_path"
printf 'libdir=%s\n' "$lib_dir" >> "$pc_path"
printf 'includedir=%s\n\n' "$include_dir" >> "$pc_path"
printf 'Name: Python\n' >> "$pc_path"
printf 'Description: CPython import library for MinGW wheel builds\n' >> "$pc_path"
printf 'Version: %s\n' "$pkg_version" >> "$pc_path"
printf 'Libs: -L${libdir} -l%s\n' "$dll_stem" >> "$pc_path"
printf 'Cflags: -I${includedir}\n' >> "$pc_path"

if command -v pkg-config >/dev/null 2>&1; then
    pkg-config --exists "python-$pkg_version"
fi

echo "Generated $import_lib"
echo "Generated $pc_path"
