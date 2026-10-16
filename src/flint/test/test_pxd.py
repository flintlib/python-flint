"""Installed wheels must ship .pxd so third-party Cython can cimport acb_t."""
from pathlib import Path

import flint


def test_cython_pxd_files_are_installed() -> None:
    root = Path(flint.__file__).resolve().parent
    needed = [
        root / "pyflint.pxd",
        root / "types" / "acb.pxd",
        root / "flint_base" / "flint_base.pxd",
        root / "flintlib" / "functions" / "acb.pxd",
        root / "flintlib" / "types" / "acb.pxd",
        root / "flintlib" / "functions" / "acb_hypgeom.pxd",
    ]
    missing = [str(p.relative_to(root)) for p in needed if not p.is_file()]
    assert missing == [], missing
