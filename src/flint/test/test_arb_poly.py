"""Tests for python-flint's `arb_poly` type."""

from __future__ import annotations

from flint import acb, acb_poly, arb, arb_poly, fmpq, fmpq_poly, fmpz, fmpz_poly
from flint.test.helpers import is_close_acb, is_close_arb, is_close_arb_poly as is_close, raises


def test_arb_poly_constructor_and_basic() -> None:
    x = arb_poly([0, 1])
    p0 = arb_poly()
    assert len(p0) == 0
    assert p0.length() == 0
    assert p0.degree() == -1

    p = arb_poly([1, 2, 3])
    assert len(p) == 3
    assert p.length() == 3
    assert p.degree() == 2
    assert is_close(p, 3*x**2 + 2*x + 1)
    assert is_close_arb(p[-1], 0)
    assert p.repr() == "arb_poly([1.00000000000000, 2.00000000000000, 3.00000000000000])"

    q = arb_poly(p)
    assert q is not p
    assert is_close(q, p)

    assert is_close(arb_poly(5), arb_poly([5]))
    assert is_close(arb_poly(fmpz_poly([1, 2])), 2*x + 1)
    assert is_close(arb_poly(fmpq_poly([1, fmpq(1, 2)])), 0.5*x + 1)

    mixed = arb_poly([1, 2.5, arb(3), fmpz(4), fmpq(5, 2)])
    assert is_close(mixed, 2.5*x**4 + 4*x**3 + 3*x**2 + 2.5*x + 1)

    assert raises(lambda: arb_poly([1, object()]), TypeError)  # type: ignore[list-item]


def test_arb_poly_setitem() -> None:
    p = arb_poly([1])
    p[3] = 5
    assert p.degree() == 3
    assert is_close_arb(p[3], 5)
    p[2] = arb(1.25)
    assert is_close_arb(p[2], 1.25)
    assert raises(lambda: p.__setitem__(-1, 1), ValueError)


def test_arb_poly_from_roots_and_complex_roots() -> None:
    x = arb_poly([0, 1])
    p = arb_poly.from_roots([0, 1, 2])
    assert is_close(p, x*(x - 1)*(x - 2))
    assert is_close_arb(p(0), 0)
    assert is_close_arb(p(1), 0)
    assert is_close_arb(p(2), 0)

    one = arb_poly.from_roots([])
    assert one.degree() == 0
    assert is_close_arb(one[0], 1)

    roots = arb_poly([1, 0, 1]).complex_roots()
    assert len(roots) == 2
    assert any(is_close_acb(r, acb(0, 1), tol=1e-5) for r in roots)
    assert any(is_close_acb(r, acb(0, -1), tol=1e-5) for r in roots)


def test_arb_poly_evaluate() -> None:
    p = arb_poly([1, 2, 3])
    ys_fast = p.evaluate([0, 1, 2], algorithm="fast")
    ys_iter = p.evaluate([0, 1, 2], algorithm="iter")

    assert len(ys_fast) == 3
    assert is_close_arb(ys_fast[0], 1)
    assert is_close_arb(ys_fast[1], 6)
    assert is_close_arb(ys_fast[2], 17)

    assert len(ys_iter) == 3
    assert is_close_arb(ys_iter[0], 1)
    assert is_close_arb(ys_iter[1], 6)
    assert is_close_arb(ys_iter[2], 17)

    assert p.evaluate([], algorithm="fast") == []
    assert raises(lambda: p.evaluate([1], algorithm="bad"), AssertionError)


def test_arb_poly_interpolate() -> None:
    x = arb_poly([0, 1])
    xs = [0, 1, 2]
    ys = [1, 6, 17]

    p_fast = arb_poly.interpolate(xs, ys, algorithm="fast")
    p_newton = arb_poly.interpolate(xs, ys, algorithm="newton")
    p_bary = arb_poly.interpolate(xs, ys, algorithm="barycentric")

    assert is_close(p_fast, 3*x**2 + 2*x + 1)
    assert is_close(p_newton, 3*x**2 + 2*x + 1)
    assert is_close(p_bary, 3*x**2 + 2*x + 1)

    assert raises(lambda: arb_poly.interpolate([0], [1, 2]), ValueError)
    assert raises(lambda: arb_poly.interpolate([0], [1], algorithm="bad"), AssertionError)


def test_arb_poly_calculus_and_unary() -> None:
    x = arb_poly([0, 1])
    p = 3*x**2 + 2*x + 1
    assert is_close(p.derivative(), 6*x + 2)
    assert is_close(p.integral(), x**3 + x**2 + x)

    assert (+p) is p
    assert is_close(-p, -(3*x**2 + 2*x + 1))


def test_arb_poly_arithmetic_real_coercions() -> None:
    x = arb_poly([0, 1])
    p = 3*x**2 + 2*x + 1

    assert is_close(p + 2, 3*x**2 + 2*x + 3)
    assert is_close(2 + p, 3*x**2 + 2*x + 3)
    assert is_close(p - 2, 3*x**2 + 2*x - 1)
    assert is_close(2 - p, -3*x**2 - 2*x + 1)
    assert is_close(p*2, 6*x**2 + 4*x + 2)
    assert is_close(2*p, 6*x**2 + 4*x + 2)

    assert is_close(p + fmpz(2), 3*x**2 + 2*x + 3)
    assert is_close(p + fmpq(1, 2), 3*x**2 + 2*x + 1.5)
    assert is_close(p + arb(1.25), 3*x**2 + 2*x + 2.25)
    assert is_close(p + fmpz_poly([1, 1]), 3*x**2 + 3*x + 2)
    assert is_close(p + fmpq_poly([1, fmpq(1, 2)]), 3*x**2 + 2.5*x + 2)
    assert is_close(p // 2, 1.5*x**2 + x + 0.5)
    assert is_close(p % 2, 0)
    q, r = divmod(p, 2)
    assert is_close(q, 1.5*x**2 + x + 0.5)
    assert is_close(r, 0)


def test_arb_poly_arithmetic_divmod_and_errors() -> None:
    x = arb_poly([0, 1])
    p = 3*x**2 + 2*x + 1
    d = x + 1

    q = p // d
    r = p % d
    q2, r2 = divmod(p, d)
    assert is_close(q, 3*x - 1)
    assert is_close(r, 2)
    assert is_close(q2, 3*x - 1)
    assert is_close(r2, 2)

    assert is_close(3 // d, 0)
    assert is_close(3 % d, 3)
    rq, rr = divmod(3, d)
    assert is_close(rq, 0)
    assert is_close(rr, 3)

    z = arb_poly([0])
    assert raises(lambda: p // z, ZeroDivisionError)
    assert raises(lambda: p % z, ZeroDivisionError)
    assert raises(lambda: divmod(p, z), ZeroDivisionError)


def test_arb_poly_arithmetic_notimplemented_paths() -> None:
    p = arb_poly([1, 2, 3])

    assert raises(lambda: p + object(), TypeError)  # type: ignore[operator]
    assert raises(lambda: p - object(), TypeError)  # type: ignore[operator]
    assert raises(lambda: p * object(), TypeError)  # type: ignore[operator]
    assert raises(lambda: p // object(), TypeError)  # type: ignore[operator]
    assert raises(lambda: p % object(), TypeError)  # type: ignore[operator]
    assert raises(lambda: divmod(p, object()), TypeError)  # type: ignore[operator]

    assert raises(lambda: object() + p, TypeError)  # type: ignore[operator]
    assert raises(lambda: object() - p, TypeError)  # type: ignore[operator]
    assert raises(lambda: object() * p, TypeError)  # type: ignore[operator]
    assert raises(lambda: object() // p, TypeError)  # type: ignore[operator]
    assert raises(lambda: object() % p, TypeError)  # type: ignore[operator]
    assert raises(lambda: divmod(object(), p), TypeError)  # type: ignore[operator]


def test_arb_poly_complex_coercion_branch() -> None:
    p = arb_poly([1, 2, 3])
    z = 2 + 3j
    a = acb(2, 3)

    r1 = p + z
    r2 = z + p
    r3 = p - a
    r4 = a - p
    r5 = p * z
    r6 = z * p

    assert isinstance(r1, acb_poly)
    assert isinstance(r2, acb_poly)
    assert isinstance(r3, acb_poly)
    assert isinstance(r4, acb_poly)
    assert isinstance(r5, acb_poly)
    assert isinstance(r6, acb_poly)


def test_arb_poly_truncate_shifts_pow() -> None:
    x = arb_poly([0, 1])
    p = 3*x**2 + 2*x + 1

    assert p.truncate(0).length() == 0
    assert is_close(p.truncate(2), 2*x + 1)
    assert is_close(p.truncate(5), p)

    assert is_close(p.left_shift(0), p)
    assert is_close(p.left_shift(2), x**2*p)

    assert is_close(p.right_shift(0), p)
    assert is_close(p.right_shift(2), 3)
    assert p.right_shift(6).length() == 0

    assert raises(lambda: p.left_shift(-1), ValueError)
    assert raises(lambda: p.right_shift(-1), ValueError)

    p3 = p**3
    assert is_close(p3, p*p*p)
    assert raises(lambda: pow(p, 2, 3), NotImplementedError)  # type: ignore[misc]


def test_arb_poly_call_paths() -> None:
    p = arb_poly([1, 2, 3])

    cp = p(arb_poly([1, 1]))
    x = arb_poly([0, 1])
    assert is_close(cp, 3*x**2 + 8*x + 6)

    assert is_close_arb(p(arb(2)), 17)
    assert is_close_arb(p(2), 17)
    assert is_close_arb(p(2.0), 17)
    assert is_close_arb(p(fmpz(2)), 17)
    assert is_close_arb(p(fmpq(1, 2)), 2.75)

    fp = p(fmpz_poly([1, 1]))
    fqp = p(fmpq_poly([1, 1]))
    assert is_close(fp, 3*x**2 + 8*x + 6)
    assert is_close(fqp, 3*x**2 + 8*x + 6)

    ac = p(acb(2, 3))
    cz = p(2 + 3j)
    assert is_close_acb(ac, acb(-10, 42))
    assert is_close_acb(cz, acb(-10, 42))

    assert raises(lambda: p(object()), TypeError)  # type: ignore[call-overload]


def test_arb_poly_unique_fmpz_poly() -> None:
    p = arb_poly([1, 2, 3])
    q = arb_poly([1.1, 2, 3])

    up = p.unique_fmpz_poly()
    assert up is not None
    assert up == fmpz_poly([1, 2, 3])
    assert q.unique_fmpz_poly() is None
