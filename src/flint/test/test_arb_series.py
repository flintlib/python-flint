"""Tests for python-flint's `arb_series` type."""

from __future__ import annotations

from flint import acb_series, arb, arb_poly, arb_series, ctx, fmpq, fmpq_poly, fmpq_series, fmpz, fmpz_poly, fmpz_series
from flint.test.helpers import is_close_arb, is_close_arb_series as is_close, raises


def test_arb_series_constructor_and_basic() -> None:
    old_cap = ctx.cap
    try:
        ctx.cap = 8
        p0 = arb_series()
        assert len(p0) == 0
        assert p0.length() == 0
        assert p0.prec == 8

        p = arb_series([1, 2], prec=3)
        assert len(p) == 2
        assert p.length() == 2
        assert p.prec == 3
        assert p.repr() == "arb_series([1.00000000000000, 2.00000000000000], prec=3)"
        assert str(p) == "1.00000000000000 + 2.00000000000000*x + O(x^3)"

        q = arb_series(p)
        assert q is not p
        assert is_close(q, p)
        assert q.prec == p.prec

        assert is_close(arb_series(fmpz_series([1, 2], prec=4)), arb_series([1, 2], prec=4))
        assert is_close(arb_series(fmpz_poly([1, 2])), arb_series([1, 2], prec=8))
        assert is_close(arb_series(arb_poly([1, 2])), arb_series([1, 2], prec=8))
        assert is_close(arb_series(5), arb_series([5], prec=8))
        assert is_close(arb_series(fmpq(1, 2)), arb_series([0.5], prec=8))

        assert str(arb_series([1], prec=0)) == "O(x^0)"
        assert str(arb_series([1], prec=-1)) == "(invalid power series)"

        x = arb_series([1], prec=4)
        x[3] = 5
        assert is_close_arb(x[3], 5)
        assert is_close_arb(x[-1], 0)
        assert raises(lambda: x.__setitem__(-1, 1), ValueError)
    finally:
        ctx.cap = old_cap


def test_arb_series_arithmetic_and_valuation() -> None:
    old_cap = ctx.cap
    try:
        ctx.cap = 8
        x = arb_series([0, 1], prec=8)
        one = arb_series([1], prec=8)

        assert (+x) is x
        assert is_close(-x, arb_series([0, -1], prec=8))

        assert is_close(x + 1, one + x)
        assert is_close(1 + x, one + x)
        assert is_close(x - 1, x + (-1))
        assert is_close(1 - x, one - x)
        assert is_close((one + x)*(one - x), one - x*x)
        assert is_close(2*x, x + x)
        assert is_close(x*2, x + x)
        assert is_close(x + fmpz(2), x + 2)
        assert is_close(x + fmpq(1, 2), x + 0.5)
        assert is_close(x + arb(1.25), x + 1.25)
        assert is_close(x + fmpz_poly([1, 1]), arb_series([1, 2], prec=8))
        assert is_close(x + fmpq_poly([1, fmpq(1, 2)]), arb_series([1, 1.5], prec=8))
        assert is_close(x + fmpz_series([1, 1], prec=8), arb_series([1, 2], prec=8))
        assert is_close(x + fmpq_series([1, fmpq(1, 2)], prec=8), arb_series([1, 1.5], prec=8))

        assert x.valuation() == 1
        assert arb_series([0], prec=8).valuation() == -1
        assert arb_series([0, 0, 2], prec=8).valuation() == 2

        assert raises(lambda: x + object(), TypeError)  # type: ignore[operator]
        assert raises(lambda: x - object(), TypeError)  # type: ignore[operator]
        assert raises(lambda: x * object(), TypeError)  # type: ignore[operator]
        assert raises(lambda: object() + x, TypeError)  # type: ignore[operator]
        assert raises(lambda: object() - x, TypeError)  # type: ignore[operator]
        assert raises(lambda: object() * x, TypeError)  # type: ignore[operator]
    finally:
        ctx.cap = old_cap


def test_arb_series_mixed_complex_paths() -> None:
    old_cap = ctx.cap
    try:
        ctx.cap = 8
        x = arb_series([0, 1], prec=8)
        r1 = x + 1j  # type: ignore[operator]
        r2 = 1j + x  # type: ignore[operator]
        r3 = x - 1j  # type: ignore[operator]
        r4 = 1j - x  # type: ignore[operator]
        r5 = x * 1j  # type: ignore[operator]
        r6 = 1j * x  # type: ignore[operator]
        assert isinstance(r1, acb_series)
        assert isinstance(r2, acb_series)
        assert isinstance(r3, acb_series)
        assert isinstance(r4, acb_series)
        assert isinstance(r5, acb_series)
        assert isinstance(r6, acb_series)
    finally:
        ctx.cap = old_cap


def test_arb_series_division_and_pow() -> None:
    old_cap = ctx.cap
    try:
        ctx.cap = 8
        x = arb_series([0, 1], prec=8)
        one = arb_series([1], prec=8)

        assert is_close(x/(one + x), x - x**2 + x**3 - x**4 + x**5 - x**6 + x**7)
        assert ((x*x)/(x*(one + x))).prec == 7
        assert is_close((x*x)/(x*(one + x)), arb_series([0, 1, -1, 1, -1, 1, -1], prec=7))
        assert (arb_series([0], prec=8)/(one + x)).prec == 8
        assert is_close(x/2, arb_series([0, 0.5], prec=8))
        assert is_close(2/(one + x), arb_series([2, -2, 2, -2, 2, -2, 2, -2], prec=8))
        assert raises(lambda: x / arb_series([0], prec=8), ZeroDivisionError)
        assert raises(lambda: one / x, ValueError)
        assert raises(lambda: 2 / x, ValueError)
        bad = arb_series([arb(0, 1), 1], prec=8)
        assert raises(lambda: x / bad, ValueError)
        assert raises(lambda: x / object(), TypeError)  # type: ignore[operator]
        assert raises(lambda: object() / x, TypeError)  # type: ignore[operator]

        assert is_close((one + x)**2, one + 2*x + x*x)
        y = 2**x
        assert is_close_arb(y[0], 1)
        assert is_close_arb(y[1], 0.6931471805599453, tol=1e-12, rel_tol=1e-12)
        assert raises(lambda: pow(x, 2, 3), NotImplementedError)  # type: ignore[misc]
        assert raises(lambda: x ** object(), TypeError)  # type: ignore[operator]
        assert raises(lambda: object() ** x, TypeError)  # type: ignore[operator]
    finally:
        ctx.cap = old_cap


def test_arb_series_call_reversion_inv() -> None:
    old_cap = ctx.cap
    try:
        ctx.cap = 8
        x = arb_series([0, 1], prec=8)
        one = arb_series([1], prec=8)
        assert is_close((one + x)(x), one + x)
        assert raises(lambda: (one + x)(one + x), ValueError)
        assert raises(lambda: (one + x)(1), TypeError)  # type: ignore[arg-type]

        r = (x + x*x).reversion()
        assert is_close((x + x*x)(r), x)
        assert raises(lambda: (one + x).reversion(), ValueError)

        inv = (one + x).inv()
        assert is_close((one + x)*inv, one)
        assert raises(lambda: arb_series([0], prec=8).inv(), ZeroDivisionError)
    finally:
        ctx.cap = old_cap


def test_arb_series_elementary_functions() -> None:
    old_cap = ctx.cap
    try:
        ctx.cap = 8
        x = arb_series([0, 1], prec=8)
        one = arb_series([1], prec=8)

        assert is_close((one + x + x*x).derivative(), arb_series([1, 2], prec=7))
        assert is_close((one + x).integral(), x + 0.5*x*x)
        assert is_close((one + x).sqrt()*(one + x).sqrt(), one + x)
        assert is_close((one + x).sqrt()*(one + x).rsqrt(), one)
        assert is_close(x.exp().log(), x)
        assert is_close((one + x).log().exp(), one + x)
        assert is_close(x.atan().tan(), x)
        assert is_close(x.asin().sin(), x)
        assert is_close(x.acos().cos(), x)

        s, c = x.sin_cos()
        assert is_close(s, x.sin())
        assert is_close(c, x.cos())
        assert is_close(s*s + c*c, one)

        sp, cp = x.sin_cos_pi()
        assert is_close(sp, x.sin_pi())
        assert is_close(cp, x.cos_pi())
        assert is_close(x.sin()/x.cos(), x.tan())
        assert (1 + x).cot_pi().prec == 8
    finally:
        ctx.cap = old_cap


def test_arb_series_special_functions() -> None:
    old_cap = ctx.cap
    try:
        ctx.cap = 8
        x = arb_series([0, 1], prec=8)
        one = arb_series([1], prec=8)

        assert is_close((one + x).gamma()*(one + x).rgamma(), one)
        assert is_close((one + x).lgamma().exp(), (one + x).gamma())
        assert is_close(x.rising(3), 2*x + 3*x*x + x**3)
        assert is_close_arb((2 + x).zeta()[0], arb.zeta(arb(2)), tol=1e-8, rel_tol=1e-8)
        assert (1 + x).riemann_siegel_theta().prec == 8
        assert (14 + x).riemann_siegel_z().prec == 8

        assert is_close(x.erf() + x.erfc(), one)
        assert x.erfi().prec == 8
        fs, fc = x.fresnel()
        assert is_close(fs, x.fresnel_s())
        assert is_close(fc, x.fresnel_c())

        assert (one + x).ei().prec == 8
        assert is_close_arb(x.si()[0], 0)
        assert is_close_arb((one + x).ci()[0], arb(0.337403922900968), tol=1e-8, rel_tol=1e-8)
        assert is_close_arb(x.shi()[0], 0)
        assert is_close_arb((one + x).chi()[0], arb(0.837866940980208), tol=1e-8, rel_tol=1e-8)
        assert (2 + x).li().prec == 8
        assert (2 + x).li(offset=True).prec == 8
    finally:
        ctx.cap = old_cap


def test_arb_series_airy_coulomb_gamma_beta_lambertw() -> None:
    old_cap = ctx.cap
    try:
        ctx.cap = 8
        x = arb_series([0, 1], prec=8)

        ai, aip, bi, bip = x.airy()
        assert is_close(ai, x.airy_ai())
        assert is_close(aip, x.airy_ai_prime())
        assert is_close(bi, x.airy_bi())
        assert is_close(bip, x.airy_bi_prime())

        f, g = x.coulomb(0, 0)
        assert is_close(f, x.coulomb_f(0, 0))
        assert str(g) == str(x.coulomb_g(0, 0))

        up = arb_series.gamma_upper(2, x)
        lo = arb_series.gamma_lower(2, x)
        assert is_close(up + lo, arb_series([1], prec=8))
        assert arb_series.gamma_upper(2, x, regularized=1).prec == 8
        assert arb_series.gamma_lower(2, x, regularized=1).prec == 8
        assert arb_series.beta_lower(2, 3, x).prec == 8
        assert arb_series.beta_lower(2, 3, x, regularized=1).prec == 8

        w0 = x.lambertw()
        wm1 = x.lambertw(-1)
        assert is_close(w0.exp()*w0, x)
        assert wm1.prec == 8
        assert raises(lambda: x.lambertw(2), ValueError)
    finally:
        ctx.cap = old_cap


def test_arb_series_find_roots() -> None:
    roots = arb_series.find_roots(lambda t: t.sin(), -1, 1, maxn=200)
    assert roots == [(arb(-1), arb(1))]
    assert arb_series.find_roots(lambda t: t + 2, -1, 1, maxn=200) == []
    assert len(arb_series.find_roots(lambda t: t*t - 0.1, -1, 1, maxn=200)) == 2
    assert raises(lambda: arb_series.find_roots(lambda t: t, 0, 1, maxn=20), ValueError)
    assert raises(lambda: arb_series.find_roots(lambda t: t.sin(), -1, 1, maxn=1), ValueError)
    assert raises(lambda: arb_series.find_roots(lambda t: (t - 1)*(t - 33 / 32), 0, 2, maxn=200), ValueError)
