from __future__ import annotations

from typing import Callable

from flint import acb, arb


def raises(f: Callable[[], object], exception: type[Exception]) -> bool:
    try:
        f()
    except exception:
        return True
    return False


def is_close_arb(
    x: arb,
    y: int | float | str | arb,
    *,
    tol: int | float | str | arb = 1e-10,
    rel_tol: int | float | str | arb = 1e-10,
    max_width: int | float | str | arb = 1e-10,
) -> bool:
    y = arb(y)
    tol_arb = arb(tol)
    rel_tol_arb = arb(rel_tol)
    max_width_arb = arb(max_width)
    return (
        isinstance(x, arb)
        and x.rad() < max_width_arb
        and y.rad() < max_width_arb
        and abs(x - y) <= tol_arb + rel_tol_arb * max(abs(x), abs(y))
    )


def is_close_acb(
    x: acb,
    y: int | float | complex | str | acb,
    *,
    tol: int | float | str | arb = 1e-10,
    rel_tol: int | float | str | arb = 1e-10,
    max_width: int | float | str | arb = 1e-10,
) -> bool:
    y = acb(y)
    return (
        isinstance(x, acb)
        and is_close_arb(x.real, y.real, tol=tol, rel_tol=rel_tol, max_width=max_width)
        and is_close_arb(x.imag, y.imag, tol=tol, rel_tol=rel_tol, max_width=max_width)
    )
