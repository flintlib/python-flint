from __future__ import annotations

import sys
from unittest.mock import patch

from flint import _has_acb_theta, acb, acb_mat, arb_mat, ctx, fmpq_mat, fmpz_mat
from flint.test.helpers import is_close_acb, is_close_acb_mat as is_close, is_close_arb_mat, raises


class _DummyMatrix:
    rows = 2
    cols = 2

    def __getitem__(self, index: tuple[int, int]) -> complex:
        i, j = index
        return complex(i + 2 * j + 1, i - j)


def test_acb_mat_constructor() -> None:
    z = fmpz_mat([[1, 2], [3, 4]])
    q = fmpq_mat([[1, 2], [3, 4]])
    a = arb_mat([[1, 2], [3, 4]])
    b = acb_mat([[1, 2], [3, 4]])

    assert is_close(acb_mat(b), [[1, 2], [3, 4]])
    assert is_close(acb_mat(a), [[1, 2], [3, 4]])
    assert is_close(acb_mat(z), [[1, 2], [3, 4]])
    assert is_close(acb_mat(q), [[1, 2], [3, 4]])
    assert is_close(acb_mat(_DummyMatrix()), [[1, 3 - 1j], [2 + 1j, 4]])

    assert isinstance(acb_mat.convert(z), acb_mat)
    assert isinstance(acb_mat.convert(q), acb_mat)
    assert isinstance(acb_mat.convert(a), acb_mat)
    assert raises(lambda: acb_mat.convert(object()), TypeError)  # type: ignore[arg-type]

    assert is_close(acb_mat(2, 2), [[0, 0], [0, 0]])
    assert is_close(acb_mat(2, 2, [1, 2, 3, 4]), [[1, 2], [3, 4]])
    assert is_close(acb_mat(3, 3, 5), [[5, 0, 0], [0, 5, 0], [0, 0, 5]])

    assert raises(lambda: acb_mat([1, 2]), TypeError)  # type: ignore[arg-type,list-item]
    assert raises(lambda: acb_mat([[1], [2, 3]]), ValueError)
    assert raises(lambda: acb_mat(object()), TypeError)  # type: ignore[call-overload]
    assert raises(lambda: acb_mat(2, 2, [1, 2, 3]), ValueError)
    assert raises(lambda: acb_mat(1, 2, 3, 4), ValueError)  # type: ignore[call-overload]


def test_acb_mat_basics_and_indexing() -> None:
    a = acb_mat([[1, 2], [3, 4]])
    assert a.nrows() == 2
    assert a.ncols() == 2
    assert is_close_acb(a[0, 1], 2)
    a[0, 1] = 1 + 2j
    assert is_close_acb(a[0, 1], 1 + 2j)

    assert raises(lambda: a[2, 0], ValueError)  # type: ignore[index]
    assert raises(lambda: a[0, 2], ValueError)  # type: ignore[index]

    def set_oob_1() -> None:
        a[2, 0] = 1

    def set_oob_2() -> None:
        a[0, 2] = 1

    def set_bad() -> None:
        a[0, 0] = object()  # type: ignore[assignment]

    assert raises(set_oob_1, ValueError)
    assert raises(set_oob_2, ValueError)
    assert raises(set_bad, TypeError)

    assert is_close(a.transpose(), [[1, 3], [1 + 2j, 4]])
    assert is_close(a.conjugate(), [[1, 1 - 2j], [3, 4]])
    assert is_close(+a, a)
    assert is_close(-a, [[-1, -1 - 2j], [-3, -4]])
    assert raises(lambda: bool(a), NotImplementedError)

    oldpretty = ctx.pretty
    try:
        ctx.pretty = False
        assert "acb_mat(" in repr(a)
    finally:
        ctx.pretty = oldpretty
    assert "[" in str(a)


def test_acb_mat_add_sub() -> None:
    a = acb_mat([[1, 2], [3, 4]])
    b = acb_mat([[4, 5], [6, 7]])
    z = fmpz_mat([[4, 5], [6, 7]])
    q = fmpq_mat([[4, 5], [6, 7]])
    r = arb_mat([[4, 5], [6, 7]])

    assert is_close(a + b, [[5, 7], [9, 11]])
    assert is_close(a + z, [[5, 7], [9, 11]])
    assert is_close(a + q, [[5, 7], [9, 11]])
    assert is_close(a + r, [[5, 7], [9, 11]])
    assert is_close(z + a, [[5, 7], [9, 11]])
    assert is_close(q + a, [[5, 7], [9, 11]])
    assert is_close(r + a, [[5, 7], [9, 11]])

    assert is_close(a - b, [[-3, -3], [-3, -3]])
    assert is_close(a - z, [[-3, -3], [-3, -3]])
    assert is_close(a - q, [[-3, -3], [-3, -3]])
    assert is_close(a - r, [[-3, -3], [-3, -3]])
    assert is_close(z - a, [[3, 3], [3, 3]])
    assert is_close(q - a, [[3, 3], [3, 3]])
    assert is_close(r - a, [[3, 3], [3, 3]])

    assert is_close(a + 2, [[3, 2], [3, 6]])
    assert is_close(2 + a, [[3, 2], [3, 6]])
    assert is_close(a - 2, [[-1, 2], [3, 2]])
    assert is_close(2 - a, [[1, -2], [-3, -2]])

    assert is_close(a + (1 + 2j), [[2 + 2j, 2], [3, 5 + 2j]])
    assert is_close((1 + 2j) + a, [[2 + 2j, 2], [3, 5 + 2j]])
    assert is_close(a - (1 + 2j), [[-2j, 2], [3, 3 - 2j]])
    assert is_close((1 + 2j) - a, [[2j, -2], [-3, -3 + 2j]])

    assert raises(lambda: a + acb_mat([[1, 2, 3]]), ValueError)
    assert raises(lambda: a - acb_mat([[1, 2, 3]]), ValueError)
    assert raises(lambda: a + object(), TypeError)  # type: ignore[operator]
    assert raises(lambda: object() + a, TypeError)  # type: ignore[operator]
    assert raises(lambda: a - object(), TypeError)  # type: ignore[operator]
    assert raises(lambda: object() - a, TypeError)  # type: ignore[operator]


def test_acb_mat_mul_div() -> None:
    a = acb_mat([[1, 2], [3, 4]])
    b = acb_mat([[4, 5], [6, 7]])
    z = fmpz_mat([[4, 5], [6, 7]])
    q = fmpq_mat([[4, 5], [6, 7]])
    r = arb_mat([[4, 5], [6, 7]])

    assert is_close(a * b, [[16, 19], [36, 43]])
    assert is_close(a * z, [[16, 19], [36, 43]])
    assert is_close(a * q, [[16, 19], [36, 43]])
    assert is_close(a * r, [[16, 19], [36, 43]])
    assert is_close(z * a, [[19, 28], [27, 40]])
    assert is_close(q * a, [[19, 28], [27, 40]])
    assert is_close(r * a, [[19, 28], [27, 40]])

    assert is_close(a * 2, [[2, 4], [6, 8]])
    assert is_close(2 * a, [[2, 4], [6, 8]])
    assert is_close(a * 0.5, [[0.5, 1], [1.5, 2]])
    assert is_close(a / 2, [[0.5, 1], [1.5, 2]])
    assert is_close(a * (1 + 2j), [[1 + 2j, 2 + 4j], [3 + 6j, 4 + 8j]])
    assert is_close((1 + 2j) * a, [[1 + 2j, 2 + 4j], [3 + 6j, 4 + 8j]])

    assert raises(lambda: a * acb_mat([[1, 2, 3]]), ValueError)
    assert raises(lambda: a * object(), TypeError)  # type: ignore[operator]
    assert raises(lambda: object() * a, TypeError)  # type: ignore[operator]
    assert raises(lambda: a / object(), TypeError)  # type: ignore[operator]


def test_acb_mat_pow_inv_solve() -> None:
    a = acb_mat([[1, 2], [3, 4]])
    assert is_close(a**2, [[7, 10], [15, 22]])
    assert raises(lambda: pow(a, 2, 3), TypeError)  # type: ignore[misc]
    assert raises(lambda: acb_mat([[1, 2, 3]]) ** 2, ValueError)

    ai = a.inv()
    assert is_close(a * ai, [[1, 0], [0, 1]], tol=1e-8, rel_tol=1e-8, max_width=1e-8)
    assert raises(lambda: acb_mat([[1, 2], [2, 4]]).inv(), ZeroDivisionError)
    assert raises(lambda: acb_mat([[1, 2, 3]]).inv(), ValueError)

    inv_ns = acb_mat([[1, 2], [2, 4]]).inv(nonstop=True)
    assert inv_ns[0, 0].is_finite() is False
    assert inv_ns[1, 1].is_finite() is False

    x = acb_mat([[1], [2]])
    b = a * x
    assert is_close(a.solve(b), x, tol=1e-8, rel_tol=1e-8, max_width=1e-8)
    assert is_close(a.solve(b, algorithm="lu"), x, tol=1e-8, rel_tol=1e-8, max_width=1e-8)
    assert is_close(a.solve(b, algorithm="precond"), x, tol=1e-8, rel_tol=1e-8, max_width=1e-8)
    assert is_close(a.solve(b, algorithm="approx"), x, tol=1e-8, rel_tol=1e-8, max_width=1e-8)

    assert raises(lambda: a.solve(b, algorithm="bad"), ValueError)  # type: ignore[arg-type]
    assert raises(lambda: acb_mat([[1, 2], [2, 4]]).solve(b), ZeroDivisionError)

    solve_ns = acb_mat([[1, 2], [2, 4]]).solve(b, nonstop=True)
    assert solve_ns[0, 0].is_finite() is False
    assert solve_ns[1, 0].is_finite() is False

    assert raises(lambda: acb_mat([[1, 2, 3]]).solve(acb_mat([[1], [2], [3]])), ValueError)
    assert raises(lambda: a.solve([[1], [2]]), TypeError)  # type: ignore[arg-type]


def test_acb_mat_special_methods() -> None:
    a = acb_mat([[1, 2], [3, 4]])
    assert is_close_acb(a.det(), -2)
    assert is_close_acb(a.trace(), 5)
    assert is_close(a.mid(), a)

    e = acb_mat(2, 2, [1, 4, -2, 1]).exp()
    assert is_close(
        e,
        [
            [-2.58607310345045, 1.18429895089106],
            [-0.592149475445530, -2.58607310345045],
        ],
        tol=1e-12,
        rel_tol=1e-12,
        max_width=1e-12,
    )

    assert raises(lambda: acb_mat([[1, 2, 3]]).det(), ValueError)
    assert raises(lambda: acb_mat([[1, 2, 3]]).trace(), ValueError)
    assert raises(lambda: acb_mat([[1, 2, 3]]).exp(), ValueError)

    p = acb_mat([[1, 1], [1, 0]]).charpoly()
    assert is_close_acb(p[0], -1)
    assert is_close_acb(p[1], -1)
    assert is_close_acb(p[2], 1)
    assert raises(lambda: acb_mat([[1, 2, 3]]).charpoly(), ValueError)

    d = acb_mat.dft(3)
    assert d.nrows() == 3
    assert d.ncols() == 3
    assert is_close_acb(d[0, 0], 0.5773502691896257)
    assert is_close_acb(d[1, 1], acb(-0.28867513459481287, -0.5), tol=1e-12, rel_tol=1e-12)

    d2 = acb_mat.dft(2, 3)
    assert d2.nrows() == 2
    assert d2.ncols() == 3


def test_acb_mat_contains_overlap_chop_cmp_real_imag() -> None:
    a = acb_mat([[1, 2], [3, 4]])
    b = (a / 3) * 3

    assert b.contains(a) is True
    assert a.contains(b) is False
    assert b.contains(fmpz_mat([[1, 2], [3, 4]])) is True
    assert (a / 3).contains(fmpq_mat([[1, 2], [3, 4]]) / 3) is True
    assert ((a / 3) * 3).contains(arb_mat([[1, 2], [3, 4]])) is True
    assert raises(lambda: a.contains(object()), TypeError)  # type: ignore[arg-type]

    assert b.overlaps(a) is True
    assert (a + 100).overlaps(a) is False

    c = acb_mat([[1e-20 + 1e-20j, 2], [3j, -1e-20 + 1e-20j]])
    chopped = c.chop(1e-10)
    assert is_close(chopped, [[0, 2], [3j, 0]])

    assert (a == acb_mat([[1, 2], [3, 4]])) is True
    assert (a != acb_mat([[1, 2], [3, 4]])) is False
    assert (a == fmpz_mat([[1, 2], [3, 4]])) is True
    assert (a != fmpz_mat([[1, 2], [3, 5]])) is True
    assert raises(lambda: a < acb_mat([[1, 2], [3, 4]]), ValueError)  # type: ignore[operator]
    assert raises(lambda: a <= acb_mat([[1, 2], [3, 4]]), ValueError)  # type: ignore[operator]
    assert raises(lambda: a > acb_mat([[1, 2], [3, 4]]), ValueError)  # type: ignore[operator]
    assert raises(lambda: a >= acb_mat([[1, 2], [3, 4]]), ValueError)  # type: ignore[operator]
    assert (a == object()) is False
    assert (a != object()) is True

    d = acb_mat.dft(3)
    assert is_close_arb_mat(d.real, [[0.5773502691896257, 0.5773502691896257, 0.5773502691896257], [0.5773502691896257, -0.28867513459481287, -0.28867513459481287], [0.5773502691896257, -0.28867513459481287, -0.28867513459481287]])
    assert is_close_arb_mat(d.imag, [[0, 0, 0], [0, -0.5, 0.5], [0, 0.5, -0.5]])


def test_acb_mat_eig_theta_and_helper() -> None:
    a = acb_mat([[1, 0], [0, 2]])

    vals = a.eig()
    assert len(vals) == 2
    assert any(v.real.contains(1) for v in vals)
    assert any(v.real.contains(2) for v in vals)

    vals_l, left = a.eig(left=True)
    assert len(vals_l) == 2
    assert isinstance(left, acb_mat)

    vals_r, right = a.eig(right=True)
    assert len(vals_r) == 2
    assert isinstance(right, acb_mat)

    vals_lr, left2, right2 = a.eig(left=True, right=True)
    assert len(vals_lr) == 2
    assert isinstance(left2, acb_mat)
    assert isinstance(right2, acb_mat)

    vals_approx = acb_mat.dft(4).eig(algorithm="approx")
    assert len(vals_approx) == 4

    vals_rump = a.eig(algorithm="rump")
    assert len(vals_rump) == 2

    vals_vm = a.eig(algorithm="vdhoeven_mourrain")
    assert len(vals_vm) == 2

    vals_tol = a.eig(tol=1e-12)
    assert len(vals_tol) == 2

    vals_multi = acb_mat.dft(4).eig(multiple=True)
    assert len(vals_multi) == 4

    vals_multi_rump = acb_mat.dft(4).eig(multiple=True, algorithm="rump")
    assert len(vals_multi_rump) == 4

    assert raises(lambda: acb_mat.dft(4).eig(multiple=True, right=True), NotImplementedError)
    assert raises(lambda: acb_mat.dft(4).eig(), ValueError)
    assert raises(lambda: acb_mat([[1, 2, 3]]).eig(), ValueError)

    assert acb_mat(0, 0).eig() == []

    tau = acb_mat([[1j]])
    if _has_acb_theta():
        theta_vals = acb_mat.theta(tau, acb_mat([[0]]))
        assert isinstance(theta_vals, acb_mat)
        assert theta_vals.nrows() == 1
        assert theta_vals.ncols() == 4
    else:
        assert raises(lambda: acb_mat.theta(tau, acb_mat([[0]])), NotImplementedError)

    with patch.dict(sys.modules, {"flint.types.acb_theta": None}):
        assert raises(lambda: acb_mat.theta(tau, acb_mat([[0]])), NotImplementedError)

    assert is_close(a, [[1, 0], [0, 2]]) is True
    assert is_close(a, acb_mat([[1, 0], [0, 2]])) is True
    assert is_close(a, [[1, 0, 0], [0, 2, 0]]) is False
    assert is_close(a, [[1, 0], [0, 3]]) is False
    assert is_close(object(), [[1]]) is False  # type: ignore[arg-type]
