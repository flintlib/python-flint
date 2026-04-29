from __future__ import annotations

from flint import acb_mat, arb, arb_mat, ctx, fmpq_mat, fmpz_mat
from flint.test.helpers import is_close_arb, is_close_arb_mat as is_close, raises


class _DummyMatrix:
    rows = 2
    cols = 2

    def __getitem__(self, index: tuple[int, int]) -> float:
        i, j = index
        return float(i + 2 * j + 1)


def test_arb_mat_constructor() -> None:
    z = fmpz_mat([[1, 2], [3, 4]])
    q = fmpq_mat([[1, 2], [3, 4]])
    a = arb_mat([[1, 2], [3, 4]])
    b = arb_mat(a)
    assert is_close(b, [[1, 2], [3, 4]])
    assert is_close(arb_mat(z), [[1, 2], [3, 4]])
    assert is_close(arb_mat(q), [[1, 2], [3, 4]])
    assert is_close(arb_mat(_DummyMatrix()), [[1, 3], [2, 4]])
    assert isinstance(arb_mat.convert(z), arb_mat)
    assert isinstance(arb_mat.convert(q), arb_mat)
    assert raises(lambda: arb_mat.convert(object()), TypeError)  # type: ignore[arg-type]

    assert is_close(arb_mat(2, 2), [[0, 0], [0, 0]])
    assert is_close(arb_mat(2, 2, [1, 2, 3, 4]), [[1, 2], [3, 4]])
    assert is_close(arb_mat(3, 3, 5), [[5, 0, 0], [0, 5, 0], [0, 0, 5]])

    assert raises(lambda: arb_mat([1, 2]), TypeError)  # type: ignore[arg-type,list-item]
    assert raises(lambda: arb_mat([[1], [2, 3]]), ValueError)
    assert raises(lambda: arb_mat(object()), TypeError)  # type: ignore[call-overload]
    assert raises(lambda: arb_mat(2, 2, [1, 2, 3]), ValueError)
    assert raises(lambda: arb_mat(1, 2, 3, 4), ValueError)  # type: ignore[call-overload]


def test_arb_mat_basics_and_indexing() -> None:
    a = arb_mat([[1, 2], [3, 4]])
    assert a.nrows() == 2
    assert a.ncols() == 2
    assert is_close_arb(a[0, 1], 2)
    a[0, 1] = 1.5
    assert is_close_arb(a[0, 1], 1.5)
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

    assert is_close(a.transpose(), [[1, 3], [1.5, 4]])
    assert is_close(+a, a)
    assert is_close(-a, [[-1, -1.5], [-3, -4]])
    assert raises(lambda: bool(a), NotImplementedError)

    oldpretty = ctx.pretty
    try:
        ctx.pretty = False
        assert "arb_mat(" in repr(a)
    finally:
        ctx.pretty = oldpretty
    assert "[" in str(a)


def test_arb_mat_add_sub() -> None:
    a = arb_mat([[1, 2], [3, 4]])
    b = arb_mat([[4, 5], [6, 7]])
    z = fmpz_mat([[4, 5], [6, 7]])
    q = fmpq_mat([[4, 5], [6, 7]])

    assert is_close(a + b, [[5, 7], [9, 11]])
    assert is_close(a + z, [[5, 7], [9, 11]])
    assert is_close(a + q, [[5, 7], [9, 11]])
    assert is_close(z + a, [[5, 7], [9, 11]])
    assert is_close(q + a, [[5, 7], [9, 11]])

    assert is_close(a - b, [[-3, -3], [-3, -3]])
    assert is_close(a - z, [[-3, -3], [-3, -3]])
    assert is_close(z - a, [[3, 3], [3, 3]])

    assert is_close(a + 2, [[3, 2], [3, 6]])
    assert is_close(2 + a, [[3, 2], [3, 6]])
    assert is_close(a - 2, [[-1, 2], [3, 2]])
    assert is_close(2 - a, [[1, -2], [-3, -2]])

    c = a + 1j
    d = 1j + a
    e = a - 1j
    f = 1j - a
    assert isinstance(c, acb_mat)
    assert isinstance(d, acb_mat)
    assert isinstance(e, acb_mat)
    assert isinstance(f, acb_mat)

    assert raises(lambda: a + arb_mat([[1, 2, 3]]), ValueError)
    assert raises(lambda: a - arb_mat([[1, 2, 3]]), ValueError)
    assert raises(lambda: a + object(), TypeError)  # type: ignore[operator]
    assert raises(lambda: object() + a, TypeError)  # type: ignore[operator]
    assert raises(lambda: a - object(), TypeError)  # type: ignore[operator]
    assert raises(lambda: object() - a, TypeError)  # type: ignore[operator]


def test_arb_mat_mul_div() -> None:
    a = arb_mat([[1, 2], [3, 4]])
    b = arb_mat([[4, 5], [6, 7]])
    z = fmpz_mat([[4, 5], [6, 7]])
    q = fmpq_mat([[4, 5], [6, 7]])

    assert is_close(a * b, [[16, 19], [36, 43]])
    assert is_close(a * z, [[16, 19], [36, 43]])
    assert is_close(a * q, [[16, 19], [36, 43]])
    assert is_close(z * a, [[19, 28], [27, 40]])
    assert is_close(q * a, [[19, 28], [27, 40]])

    assert is_close(a * 2, [[2, 4], [6, 8]])
    assert is_close(2 * a, [[2, 4], [6, 8]])
    assert is_close(a * 0.5, [[0.5, 1], [1.5, 2]])
    assert is_close(a / 2, [[0.5, 1], [1.5, 2]])

    c = a * (1 + 2j)
    d = (1 + 2j) * a
    assert isinstance(c, acb_mat)
    assert isinstance(d, acb_mat)

    assert raises(lambda: a * arb_mat([[1, 2, 3]]), ValueError)
    assert raises(lambda: a * object(), TypeError)  # type: ignore[operator]
    assert raises(lambda: object() * a, TypeError)  # type: ignore[operator]
    nan_mat = a / 0
    assert nan_mat[0, 0].is_nan() is True
    assert nan_mat[1, 1].is_nan() is True
    assert raises(lambda: a / object(), TypeError)  # type: ignore[operator]


def test_arb_mat_pow_inv_solve() -> None:
    a = arb_mat([[1, 2], [3, 4]])
    assert is_close(a**2, [[7, 10], [15, 22]])
    assert raises(lambda: pow(a, 2, 3), NotImplementedError)  # type: ignore[misc]
    assert raises(lambda: arb_mat([[1, 2, 3]])**2, ValueError)

    ai = a.inv()
    eye = a * ai
    assert is_close(eye, [[1, 0], [0, 1]], tol=1e-8, rel_tol=1e-8, max_width=1e-8)
    assert raises(lambda: arb_mat([[1, 2], [2, 4]]).inv(), ZeroDivisionError)
    assert raises(lambda: arb_mat([[1, 2, 3]]).inv(), ValueError)
    inv_ns = arb_mat([[1, 2], [2, 4]]).inv(nonstop=True)
    assert inv_ns[0, 0].is_nan() is True
    assert inv_ns[1, 1].is_nan() is True

    x = arb_mat([[1], [2]])
    b = a * x
    assert is_close(a.solve(b), x, tol=1e-8, rel_tol=1e-8, max_width=1e-8)
    assert is_close(a.solve(b, algorithm="lu"), x, tol=1e-8, rel_tol=1e-8, max_width=1e-8)
    assert is_close(a.solve(b, algorithm="precond"), x, tol=1e-8, rel_tol=1e-8, max_width=1e-8)
    assert is_close(a.solve(b, algorithm="approx"), x, tol=1e-8, rel_tol=1e-8, max_width=1e-8)
    assert raises(lambda: a.solve(b, algorithm="bad"), ValueError)  # type: ignore[arg-type]
    assert raises(lambda: arb_mat([[1, 2], [2, 4]]).solve(b), ZeroDivisionError)
    solve_ns = arb_mat([[1, 2], [2, 4]]).solve(b, nonstop=True)
    assert solve_ns[0, 0].is_nan() is True
    assert solve_ns[1, 0].is_nan() is True
    assert raises(lambda: arb_mat([[1, 2, 3]]).solve(arb_mat([[1], [2], [3]])), ValueError)
    assert raises(lambda: a.solve([[1], [2]]), TypeError)  # type: ignore[arg-type]


def test_arb_mat_special_methods() -> None:
    a = arb_mat([[1, 2], [3, 4]])
    assert is_close_arb(a.det(), -2)
    assert is_close_arb(a.trace(), 5)
    assert is_close(a.mid(), a)
    assert is_close(
        a.exp(),
        [[51.9689561987050, 74.7365645670032], [112.104846850505, 164.073803049210]],
        tol=1e-7,
        rel_tol=1e-7,
        max_width=1e-7,
    )
    assert raises(lambda: arb_mat([[1, 2, 3]]).det(), ValueError)
    assert raises(lambda: arb_mat([[1, 2, 3]]).trace(), ValueError)
    assert raises(lambda: arb_mat([[1, 2, 3]]).exp(), ValueError)

    p = a.charpoly()
    assert is_close_arb(p[0], -2)
    assert is_close_arb(p[1], -5)
    assert is_close_arb(p[2], 1)
    assert raises(lambda: arb_mat([[1, 2, 3]]).charpoly(), ValueError)

    h = arb_mat.hilbert(2, 2)
    assert is_close(h, [[1, 1/2], [1/2, 1/3]], tol=1e-12, rel_tol=1e-12, max_width=1e-12)

    ps = arb_mat.pascal(3, 4)
    assert is_close(ps, [[1, 1, 1, 1], [1, 2, 3, 4], [1, 3, 6, 10]])
    pu = arb_mat.pascal(3, 4, 1)
    assert is_close(pu, [[1, 1, 1, 1], [0, 1, 2, 3], [0, 0, 1, 3]])
    pl = arb_mat.pascal(3, 4, -1)
    assert is_close(pl, [[1, 0, 0, 0], [1, 1, 0, 0], [1, 2, 1, 0]])

    st0 = arb_mat.stirling(4, 3, 0)
    assert is_close(st0, [[1, 0, 0], [0, 1, 0], [0, 1, 1], [0, 2, 3]])
    st1 = arb_mat.stirling(4, 3, 1)
    assert is_close(st1, [[1, 0, 0], [0, 1, 0], [0, -1, 1], [0, 2, -3]])
    st2 = arb_mat.stirling(4, 3, 2)
    assert is_close(st2, [[1, 0, 0], [0, 1, 0], [0, 1, 1], [0, 1, 3]])
    assert raises(lambda: arb_mat.stirling(2, 2, 5), ValueError)

    d = arb_mat.dct(2)
    assert d.nrows() == 2
    assert d.ncols() == 2
    d2 = arb_mat.dct(2, 3)
    assert d2.nrows() == 2
    assert d2.ncols() == 3


def test_arb_mat_contains_overlap_chop_cmp_eig() -> None:
    a = arb_mat([[1, 2], [3, 4]])
    b = (a / 3) * 3
    assert b.contains(a) is True
    assert a.contains(b) is False
    assert b.contains(fmpz_mat([[1, 2], [3, 4]])) is True
    assert (a / 3).contains(fmpq_mat([[1, 2], [3, 4]]) / 3) is True
    assert raises(lambda: a.contains(object()), TypeError)  # type: ignore[arg-type]

    assert b.overlaps(a) is True
    assert (a + 100).overlaps(a) is False

    c = arb_mat([[1e-20, 2], [3, -1e-20]])
    chopped = c.chop(1e-10)
    assert is_close(chopped, [[0, 2], [3, 0]])

    assert (a == arb_mat([[1, 2], [3, 4]])) is True
    assert (a != arb_mat([[1, 2], [3, 4]])) is False
    assert (a == fmpz_mat([[1, 2], [3, 4]])) is True
    assert (a != fmpz_mat([[1, 2], [3, 5]])) is True
    assert raises(lambda: a < arb_mat([[1, 2], [3, 4]]), ValueError)  # type: ignore[operator]
    assert (a == object()) is False
    assert (a != object()) is True

    eigvals = arb_mat([[1, 0], [0, 2]]).eig()
    assert len(eigvals) == 2
    assert any(v.real.contains(1) for v in eigvals)
    assert any(v.real.contains(2) for v in eigvals)


def test_is_close_arb_mat() -> None:
    x = arb_mat([[1, 2], [3, 4]])
    assert is_close(x, [[1, 2], [3, 4]]) is True
    assert is_close(x, arb_mat([[1, 2], [3, 4]])) is True
    assert is_close(x, [[1, 2, 3]]) is False
    assert is_close(x, [[1, 2], [3, 5]]) is False
    assert is_close(object(), [[1]]) is False  # type: ignore[arg-type]
