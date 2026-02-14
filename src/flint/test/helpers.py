from __future__ import annotations

from typing import Callable

from flint import acb, arb, arb_poly


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


def is_close_arb_poly(
    x: arb_poly,
    y: arb_poly | int | float | arb,
    *,
    tol: int | float | str | arb = 1e-10,
    rel_tol: int | float | str | arb = 1e-10,
    max_width: int | float | str | arb = 1e-10,
) -> bool:
    if not isinstance(x, arb_poly):
        return False
    if not isinstance(y, arb_poly):
        y = arb_poly(y)
    d = x - y
    for i in range(d.length()):
        if not is_close_arb(d[i], 0, tol=tol, rel_tol=rel_tol, max_width=max_width):
            return False
    return True
