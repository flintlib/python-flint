#!/usr/bin/env python3

from __future__ import annotations

import shlex
import sys
import sysconfig
from pathlib import Path


def _normalize_path(path: Path | str) -> str:
    return str(path).replace('\\', '/')


def _default_dll_name() -> str:
    vernum = sysconfig.get_config_var('py_version_nodot')
    if not vernum:
        vernum = f"{sys.version_info.major}{sys.version_info.minor}"
    is_freethreaded = bool(sysconfig.get_config_var('Py_GIL_DISABLED'))
    if not is_freethreaded:
        abiflags = sysconfig.get_config_var('ABIFLAGS') or getattr(sys, 'abiflags', '') or ''
        is_freethreaded = 't' in abiflags
    suffix = 't' if is_freethreaded else ''
    return f'python{vernum}{suffix}.dll'


def _normalize_dll_name(name: str) -> str:
    if name.endswith('.dll'):
        return name
    if name.startswith('lib') and name.endswith('.dll.a'):
        return _default_dll_name()
    return _default_dll_name()


def _find_dll(dll_name: str) -> Path:
    candidates: list[Path] = []
    libdir = sysconfig.get_config_var('LIBDIR')
    if libdir:
        candidates.append(Path(libdir) / dll_name)
    exe_dir = Path(sys.executable).resolve().parent
    candidates.append(exe_dir / dll_name)
    base_prefix = Path(getattr(sys, 'base_prefix', sys.prefix))
    candidates.append(base_prefix / dll_name)
    candidates.append(base_prefix / 'DLLs' / dll_name)
    candidates.append(base_prefix / 'libs' / dll_name)

    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()

    paths = '\n'.join(_normalize_path(path) for path in candidates)
    raise SystemExit(f'Could not find {dll_name} in:\n{paths}')


def main() -> None:
    repo_root = Path.cwd().resolve()
    lib_dir = repo_root / '.local' / 'lib'
    pkgconfig_dir = lib_dir / 'pkgconfig'
    lib_dir.mkdir(parents=True, exist_ok=True)
    pkgconfig_dir.mkdir(parents=True, exist_ok=True)

    raw_name = sysconfig.get_config_var('DLLLIBRARY') or sysconfig.get_config_var('LDLIBRARY') or ''
    dll_name = _normalize_dll_name(raw_name)
    include_dir = sysconfig.get_config_var('INCLUDEPY')
    if not include_dir:
        raise SystemExit('Could not determine Python include dir')
    pkg_version = sysconfig.get_config_var('LDVERSION') or sysconfig.get_python_version()
    dll_path = _find_dll(dll_name)

    values = {
        'REPO_ROOT': _normalize_path(repo_root),
        'LIB_DIR': _normalize_path(lib_dir),
        'PKGCONFIG_DIR': _normalize_path(pkgconfig_dir),
        'DLL_NAME': dll_name,
        'DLL_PATH': _normalize_path(dll_path),
        'PKG_VERSION': pkg_version,
        'INCLUDE_DIR': _normalize_path(include_dir),
    }

    env_path = repo_root / '.local' / 'cibw_before_build_windows_arm64.env'
    env_text = ''.join(f'{key}={shlex.quote(value)}\n' for key, value in values.items())
    env_path.write_text(env_text, encoding='utf-8')

    python_site = repo_root / '.local' / 'python-site'
    python_site.mkdir(parents=True, exist_ok=True)
    sitecustomize = python_site / 'sitecustomize.py'
    sitecustomize.write_text(
        'import sysconfig\n'
        f"sysconfig.get_config_vars()['LIBPC'] = {values['PKGCONFIG_DIR']!r}\n",
        encoding='utf-8',
    )

    print(f'Generated {env_path}')
    print(f'Generated {sitecustomize}')


if __name__ == '__main__':
    main()
