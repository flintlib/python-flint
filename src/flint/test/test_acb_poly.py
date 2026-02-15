"""Tests for python-flint's `acb_poly` type."""

from __future__ import annotations

from flint import acb, acb_poly, arb, arb_poly, fmpq, fmpq_poly, fmpz, fmpz_poly
from flint.test.helpers import is_close_acb, is_close_acb_poly as is_close, is_close_arb, raises


def test_acb_poly_constructor_and_basic() -> None:
    x = acb_poly([0, 1])
    p0 = acb_poly()
    assert len(p0) == 0
    assert p0.length() == 0
    assert p0.degree() == -1

    p = acb_poly([1, 2, 3])
    assert len(p) == 3
    assert p.length() == 3
    assert p.degree() == 2
    assert is_close(p, 3*x**2 + 2*x + 1)
    assert is_close_acb(p[-1], 0)
    assert p.repr() == "acb_poly([1.00000000000000, 2.00000000000000, 3.00000000000000])"

    q = acb_poly(p)
    assert q is not p
    assert is_close(q, p)

    assert is_close(acb_poly(5), acb_poly([5]))
    assert is_close(acb_poly(2 + 3j), acb_poly([acb(2, 3)]))
    assert is_close(acb_poly(arb_poly([1, 2])), 2*x + 1)
    assert is_close(acb_poly(fmpz_poly([1, 2])), 2*x + 1)
    assert is_close(acb_poly(fmpq_poly([1, fmpq(1, 2)])), 0.5*x + 1)

    mixed = acb_poly([1, 2.5, 3 + 4j, arb(4), fmpz(5), fmpq(7, 2), acb(1, -2)])
    assert is_close(mixed, (1 - 2j)*x**6 + 3.5*x**5 + 5*x**4 + 4*x**3 + (3 + 4j)*x**2 + 2.5*x + 1)

    assert raises(lambda: acb_poly([1, object()]), TypeError)  # type: ignore[list-item]


def test_acb_poly_setitem() -> None:
    p = acb_poly([1])
    p[3] = 5
    assert p.degree() == 3
    assert is_close_acb(p[3], 5)
    p[2] = acb(1.25, -0.5)
    assert is_close_acb(p[2], acb(1.25, -0.5))
    assert raises(lambda: p.__setitem__(-1, 1), ValueError)


def test_acb_poly_from_roots() -> None:
    x = acb_poly([0, 1])
    p = acb_poly.from_roots([0, 1, 2])
    assert is_close(p, x*(x - 1)*(x - 2))
    assert is_close_acb(p(0), 0)
    assert is_close_acb(p(1), 0)
    assert is_close_acb(p(2), 0)

    q = acb_poly.from_roots([1j, -1j])
    assert is_close(q, x**2 + 1)
    assert is_close_acb(q(1j), 0)
    assert is_close_acb(q(-1j), 0)

    one = acb_poly.from_roots([])
    assert one.degree() == 0
    assert is_close_acb(one[0], 1)


def test_acb_poly_evaluate() -> None:
    p = acb_poly([1, 2, 3])
    ys_fast = p.evaluate([0, 1, 2], algorithm="fast")
    ys_iter = p.evaluate([0, 1, 2], algorithm="iter")

    assert len(ys_fast) == 3
    assert is_close_acb(ys_fast[0], 1)
    assert is_close_acb(ys_fast[1], 6)
    assert is_close_acb(ys_fast[2], 17)

    assert len(ys_iter) == 3
    assert is_close_acb(ys_iter[0], 1)
    assert is_close_acb(ys_iter[1], 6)
    assert is_close_acb(ys_iter[2], 17)

    assert p.evaluate([], algorithm="fast") == []
    assert raises(lambda: p.evaluate([1], algorithm="bad"), AssertionError)


def test_acb_poly_interpolate() -> None:
    x = acb_poly([0, 1])
    xs = [0, 1, 2]
    ys = [1, 6, 17]

    p_fast = acb_poly.interpolate(xs, ys, algorithm="fast")
    p_newton = acb_poly.interpolate(xs, ys, algorithm="newton")
    p_bary = acb_poly.interpolate(xs, ys, algorithm="barycentric")

    assert is_close(p_fast, 3*x**2 + 2*x + 1)
    assert is_close(p_newton, 3*x**2 + 2*x + 1)
    assert is_close(p_bary, 3*x**2 + 2*x + 1)

    assert raises(lambda: acb_poly.interpolate([0], [1, 2]), ValueError)
    assert raises(lambda: acb_poly.interpolate([0], [1], algorithm="bad"), AssertionError)


def test_acb_poly_calculus_and_unary() -> None:
    x = acb_poly([0, 1])
    p = 3*x**2 + 2*x + 1
    assert is_close(p.derivative(), 6*x + 2)
    assert is_close(p.integral(), x**3 + x**2 + x)
    assert (+p) is p
    assert is_close(-p, -(3*x**2 + 2*x + 1))


def test_acb_poly_arithmetic_real_and_complex_coercions() -> None:
    x = acb_poly([0, 1])
    p = 3*x**2 + 2*x + 1

    assert is_close(p + 2, 3*x**2 + 2*x + 3)
    assert is_close(2 + p, 3*x**2 + 2*x + 3)
    assert is_close(p - 2, 3*x**2 + 2*x - 1)
    assert is_close(2 - p, -3*x**2 - 2*x + 1)
    assert is_close(p*2, 6*x**2 + 4*x + 2)
    assert is_close(2*p, 6*x**2 + 4*x + 2)

    assert is_close(p + 1j, 3*x**2 + 2*x + (1 + 1j))
    assert is_close((1 + 1j) + p, 3*x**2 + 2*x + (2 + 1j))
    assert is_close(p - 1j, 3*x**2 + 2*x + (1 - 1j))
    assert is_close((1 + 1j) - p, -3*x**2 - 2*x + 1j)
    assert is_close(p*(1 + 1j), (3 + 3j)*x**2 + (2 + 2j)*x + (1 + 1j))
    assert is_close((1 + 1j)*p, (3 + 3j)*x**2 + (2 + 2j)*x + (1 + 1j))

    assert is_close(p + fmpz(2), 3*x**2 + 2*x + 3)
    assert is_close(p + fmpq(1, 2), 3*x**2 + 2*x + 1.5)
    assert is_close(p + arb(1.25), 3*x**2 + 2*x + 2.25)
    assert is_close(p + acb(1, -2), 3*x**2 + 2*x + acb(2, -2))
    assert is_close(p + fmpz_poly([1, 1]), 3*x**2 + 3*x + 2)
    assert is_close(p + fmpq_poly([1, fmpq(1, 2)]), 3*x**2 + 2.5*x + 2)
    assert is_close(p + arb_poly([1, 1]), 3*x**2 + 3*x + 2)


def test_acb_poly_arithmetic_divmod_and_errors() -> None:
    x = acb_poly([0, 1])
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

    assert is_close(p // acb(2), 1.5*x**2 + x + 0.5)
    assert is_close(p % acb(2), 0)
    q3, r3 = divmod(p, acb(2))
    assert is_close(q3, 1.5*x**2 + x + 0.5)
    assert is_close(r3, 0)

    z = acb_poly([0])
    assert raises(lambda: p // z, ZeroDivisionError)
    assert raises(lambda: p % z, ZeroDivisionError)
    assert raises(lambda: divmod(p, z), ZeroDivisionError)


def test_acb_poly_arithmetic_unsupported_operands() -> None:
    p = acb_poly([1, 2, 3])
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


def test_acb_poly_truncate_shifts_pow() -> None:
    x = acb_poly([0, 1])
    p = 3*x**2 + 2*x + 1

    assert p.truncate(0).length() == 0
    assert is_close(p.truncate(-1), 0)
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


def test_acb_poly_call_paths() -> None:
    p = acb_poly([1, 2, 3])
    x = acb_poly([0, 1])

    assert is_close(p(acb_poly([1, 1])), 3*x**2 + 8*x + 6)
    assert is_close(p(fmpz_poly([1, 1])), 3*x**2 + 8*x + 6)
    assert is_close(p(fmpq_poly([1, 1])), 3*x**2 + 8*x + 6)
    assert is_close(p(arb_poly([1, 1])), 3*x**2 + 8*x + 6)

    assert is_close_acb(p(acb(2, 3)), acb(-10, 42))
    assert is_close_acb(p(2 + 3j), acb(-10, 42))
    assert is_close_acb(p(arb(2)), 17)
    assert is_close_acb(p(2), 17)
    assert is_close_acb(p(2.0), 17)
    assert is_close_acb(p(fmpz(2)), 17)
    assert is_close_acb(p(fmpq(1, 2)), 2.75)

    assert raises(lambda: p(object()), TypeError)  # type: ignore[call-overload]


def test_acb_poly_unique_fmpz_poly() -> None:
    p = acb_poly([1, 2, 3])
    q = acb_poly([1.1, 2, 3])
    assert p.unique_fmpz_poly() == fmpz_poly([1, 2, 3])
    assert q.unique_fmpz_poly() is None


def test_acb_poly_roots_and_complex_roots() -> None:
    x = acb_poly([0, 1])
    p = x**2 + 1
    rs0 = p.roots()
    assert len(rs0) == 2
    assert any(is_close_acb(r, 1j, tol=1e-5) for r in rs0)
    assert any(is_close_acb(r, -1j, tol=1e-5) for r in rs0)

    rs = p.roots(tol=1e-30, maxprec=128)
    assert len(rs) == 2
    assert any(is_close_acb(r, 1j, tol=arb("1e-20"), rel_tol=arb("1e-20")) for r in rs)
    assert any(is_close_acb(r, -1j, tol=arb("1e-20"), rel_tol=arb("1e-20")) for r in rs)
    assert is_close_acb(p(rs[0]), 0, tol=arb("1e-20"), rel_tol=arb("1e-20"))
    assert is_close_acb(p(rs[1]), 0, tol=arb("1e-20"), rel_tol=arb("1e-20"))

    rs2 = p.complex_roots(tol=1e-30, maxprec=128)
    assert len(rs2) == 2
    assert any(is_close_acb(r, 1j, tol=arb("1e-20"), rel_tol=arb("1e-20")) for r in rs2)
    assert any(is_close_acb(r, -1j, tol=arb("1e-20"), rel_tol=arb("1e-20")) for r in rs2)

    z = acb_poly()
    assert z.roots() == []
    assert z.complex_roots() == []

    squareful = x**2
    assert raises(lambda: squareful.roots(tol=1e-20, maxprec=128), ValueError)
    assert raises(lambda: (x**3 - 1).roots(tol=1e-30, maxprec=40), ValueError)


def test_acb_poly_root_bound() -> None:
    x = acb_poly([0, 1])
    p = x**2 + 1
    b = p.root_bound()
    assert is_close_arb(b, 1.4142135623730951, tol=1e-7, rel_tol=1e-7)
    assert (b >= 1) is True
