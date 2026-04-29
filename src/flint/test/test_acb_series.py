"""Tests for python-flint's `acb_series` type."""

from __future__ import annotations

from flint import (
    acb,
    acb_poly,
    acb_series,
    arb,
    arb_series,
    ctx,
    dirichlet_char,
    fmpq,
    fmpq_poly,
    fmpq_series,
    fmpz,
    fmpz_poly,
    fmpz_series,
)
from flint.test.helpers import is_close_acb, is_close_acb_series as is_close, raises


def test_acb_series_constructor_and_basic() -> None:
    old_cap = ctx.cap
    try:
        ctx.cap = 8
        p0 = acb_series()
        assert len(p0) == 0
        assert p0.length() == 0
        assert p0.prec == 8

        p = acb_series([1, 2], prec=3)
        assert len(p) == 2
        assert p.length() == 2
        assert p.prec == 3
        assert p.repr() == "acb_series([1.00000000000000, 2.00000000000000], prec=3)"
        assert str(p) == "1.00000000000000 + 2.00000000000000*x + O(x^3)"

        q = acb_series(p)
        assert q is not p
        assert is_close(q, p)
        assert q.prec == p.prec

        assert is_close(acb_series(arb_series([1, 2], prec=4)), acb_series([1, 2], prec=4))
        assert is_close(acb_series(fmpz_series([1, 2], prec=4)), acb_series([1, 2], prec=4))
        assert is_close(acb_series(fmpq_series([1, fmpq(1, 2)], prec=4)), acb_series([1, 1/2], prec=4))
        assert is_close(acb_series(fmpz_poly([1, 2])), acb_series([1, 2], prec=8))
        assert is_close(acb_series(fmpq_poly([1, fmpq(1, 2)])), acb_series([1, 1/2], prec=8))
        assert is_close(acb_series(acb_poly([1, 2])), acb_series([1, 2], prec=8))
        assert is_close(acb_series(5, prec=1), [5])

        assert str(acb_series([1], prec=0)) == "O(x^0)"
        assert str(acb_series([1], prec=-1)) == "(invalid power series)"

        x = acb_series([1], prec=4)
        x[3] = 5
        assert is_close_acb(x[3], 5)
        assert is_close_acb(x[-1], 0)
        assert raises(lambda: x.__setitem__(-1, 1), ValueError)
    finally:
        ctx.cap = old_cap


def test_acb_series_arithmetic_and_valuation() -> None:
    old_cap = ctx.cap
    try:
        ctx.cap = 8
        x = acb_series([0, 1], prec=8)
        one = acb_series([1], prec=8)
        x2 = acb_series([0, 1], prec=2)
        one2 = acb_series([1], prec=2)

        assert (+x) is x
        assert is_close(-x2, [0, -1])

        assert is_close(x + 1, one + x)
        assert is_close(1 + x, one + x)
        assert is_close(x - 1, x + (-1))
        assert is_close(1 - x, one - x)
        assert is_close((one + x) * (one - x), one - x * x)
        assert is_close(2 * x, x + x)
        assert is_close(x * 2, x + x)
        assert is_close(x2 + 2j, [2j, 1])
        assert is_close(2j + x2, [2j, 1])
        assert is_close(x + fmpz(2), x + 2)
        assert is_close(x + fmpq(1, 2), x + 0.5)
        assert is_close(x + arb(1.25), x + 1.25)
        assert is_close(x + acb(1.25, -0.5), x + (1.25 - 0.5j))
        assert is_close(x2 + fmpz_poly([1, 1]), [1, 2])
        assert is_close(x2 + fmpq_poly([1, fmpq(1, 2)]), [1, 3/2])
        assert is_close(x2 + fmpz_series([1, 1], prec=2), [1, 2])
        assert is_close(x2 + fmpq_series([1, fmpq(1, 2)], prec=2), [1, 3/2])

        x8 = acb_series([0, 1], prec=8)
        assert x8.valuation() == 1
        assert acb_series([0], prec=8).valuation() == -1
        assert acb_series([0, 0, 2], prec=8).valuation() == 2

        assert raises(lambda: x8 + object(), TypeError)  # type: ignore[operator]
        assert raises(lambda: x8 - object(), TypeError)  # type: ignore[operator]
        assert raises(lambda: x8 * object(), TypeError)  # type: ignore[operator]
        assert raises(lambda: object() + x8, TypeError)  # type: ignore[operator]
        assert raises(lambda: object() - x8, TypeError)  # type: ignore[operator]
        assert raises(lambda: object() * x8, TypeError)  # type: ignore[operator]
    finally:
        ctx.cap = old_cap


def test_acb_series_division_and_pow() -> None:
    old_cap = ctx.cap
    try:
        ctx.cap = 8
        x = acb_series([0, 1], prec=8)
        one = acb_series([1], prec=8)

        assert is_close(x / (one + x), x - x**2 + x**3 - x**4 + x**5 - x**6 + x**7)
        assert ((x * x) / (x * (one + x))).prec == 7
        assert is_close((x * x) / (x * (one + x)), [0, 1, -1, 1, -1, 1, -1])
        assert (acb_series([0], prec=8) / (one + x)).prec == 8
        x2 = acb_series([0, 1], prec=2)
        assert is_close(x2 / 2, [0, 1/2])
        assert is_close(2 / (one + x), [2, -2, 2, -2, 2, -2, 2, -2])
        assert raises(lambda: x / acb_series([0], prec=8), ZeroDivisionError)
        assert raises(lambda: one / x, ValueError)
        assert raises(lambda: 2 / x, ValueError)
        bad = acb_series([acb(arb(0, 1), 0), 1], prec=8)
        assert raises(lambda: x / bad, ValueError)
        assert raises(lambda: x / object(), TypeError)  # type: ignore[operator]
        assert raises(lambda: object() / x, TypeError)  # type: ignore[operator]

        assert is_close((one + x) ** 2, one + 2 * x + x * x)
        y = 2**x
        assert is_close_acb(y[0], 1)
        assert is_close_acb(y[1], 0.6931471805599453, tol=1e-12, rel_tol=1e-12)
        assert raises(lambda: pow(x, 2, 3), NotImplementedError)  # type: ignore[misc]
        assert raises(lambda: x**object(), TypeError)  # type: ignore[operator]
        assert raises(lambda: object() ** x, TypeError)  # type: ignore[operator]
    finally:
        ctx.cap = old_cap


def test_acb_series_call_reversion_inv() -> None:
    old_cap = ctx.cap
    try:
        ctx.cap = 8
        x = acb_series([0, 1], prec=8)
        one = acb_series([1], prec=8)

        assert is_close((one + x)(x), one + x)
        assert raises(lambda: (one + x)(one + x), ValueError)
        assert raises(lambda: (one + x)(1), TypeError)  # type: ignore[arg-type]

        r = (x + x * x).reversion()
        assert is_close((x + x * x)(r), x)
        assert raises(lambda: (one + x).reversion(), ValueError)

        inv = (one + x).inv()
        assert is_close((one + x) * inv, one)
        assert raises(lambda: acb_series([0], prec=8).inv(), ZeroDivisionError)
    finally:
        ctx.cap = old_cap


def test_acb_series_elementary_functions() -> None:
    old_cap = ctx.cap
    try:
        ctx.cap = 8
        x = acb_series([0, 1], prec=8)
        one = acb_series([1], prec=8)
        x2 = acb_series([0, 1], prec=2)
        one2 = acb_series([1], prec=2)
        x3 = acb_series([0, 1], prec=3)
        one3 = acb_series([1], prec=3)

        x7 = acb_series([0, 1], prec=7)
        assert is_close((one3 + x3 + x3 * x3).derivative(), [1, 2])
        assert is_close((one2 + x2).integral(), [0, 1, 1/2])
        assert is_close((one + x).sqrt(), [1, 1/2, -1/8, 1/16, -5/128, 7/256, -21/1024, 33/2048])
        assert is_close((one + x).rsqrt(), [1, -1/2, 3/8, -5/16, 35/128, -63/256, 231/1024, -429/2048])
        assert is_close(x.exp(), [1, 1, 1/2, 1/6, 1/24, 1/120, 1/720, 1/5040])
        assert is_close((one + x).log(), [0, 1, -1/2, 1/3, -1/4, 1/5, -1/6, 1/7])
        assert is_close(x.atan(), [0, 1, 0, -1/3, 0, 1/5, 0, -1/7])
        assert is_close(x.sin(), [0, 1, 0, -1/6, 0, 1/120, 0, -1/5040])
        assert is_close(x7.cos(), [1, 0, -1/2, 0, 1/24, 0, -1/720])
        assert is_close(x.tan(), [0, 1, 0, 1/3, 0, 2/15, 0, 17/315])

        s, c = x.sin_cos()
        assert is_close(s, [0, 1, 0, -1/6, 0, 1/120, 0, -1/5040])
        assert is_close(c, x.cos())
        assert is_close(s * s + c * c, one)

        pi = acb.pi()
        assert is_close(x.sin_pi(), [0, pi, 0, -(pi**3)/6, 0, (pi**5)/120, 0, -(pi**7)/5040])
        assert is_close(x7.cos_pi(), [1, 0, -(pi**2)/2, 0, (pi**4)/24, 0, -(pi**6)/720])
        sp, cp = x.sin_cos_pi()
        assert is_close(sp, [0, pi, 0, -(pi**3)/6, 0, (pi**5)/120, 0, -(pi**7)/5040])
        assert is_close(cp, x.cos_pi())
        assert is_close(x.sin() / x.cos(), x.tan())
        assert (1 + x).cot_pi().prec == 8
    finally:
        ctx.cap = old_cap


def test_acb_series_special_and_classmethods() -> None:
    old_cap = ctx.cap
    try:
        ctx.cap = 8
        x = acb_series([0, 1], prec=8)
        one = acb_series([1], prec=8)

        assert is_close((one + x).gamma() * (one + x).rgamma(), one)
        assert is_close((one + x).lgamma().exp(), (one + x).gamma())
        assert is_close(x.rising(3), 2 * x + 3 * x * x + x**3)
        assert is_close_acb((2 + x).zeta()[0], acb(2).zeta(), tol=1e-8, rel_tol=1e-8)
        chi = dirichlet_char(3, 1)
        assert is_close((2 + x).dirichlet_l(chi), (2 + x).dirichlet_l((3, 1)))
        assert (2 + x).dirichlet_l((3, 1), deflate=True).prec == 8

        p = acb_series.polylog(1 + x, 0.5)
        assert p.prec == 8

        assert (2 + x).agm().prec == 8
        assert is_close((2 + x).agm(3 + x), ((2 + x) / (3 + x)).agm() * (3 + x))
        assert (x * x).elliptic_k().prec == 8
        assert x.elliptic_p(1j).prec == 8

        assert is_close(x.erf() + x.erfc(), one)
        assert x.erfi().prec == 8

        up = acb_series.gamma_upper(2, x)
        lo = acb_series.gamma_lower(2, x)
        assert is_close(acb_series([up[0] + lo[0]], prec=1), [1])
        assert acb_series.gamma_upper(2, x, regularized=1).prec == 8
        assert acb_series.gamma_lower(2, x, regularized=1).prec == 8
        assert acb_series.beta_lower(2, 3, x).prec == 8
        assert acb_series.beta_lower(2, 3, x, regularized=1).prec == 8

        assert acb_series.hypgeom([1 + x], [2 + x], x).prec == 8
        assert acb_series.hypgeom([1 + x], [2 + x], x, n=4, regularized=True).prec == 8

        w0 = x.lambertw()
        wm1 = x.lambertw(-1)
        assert is_close(w0.exp() * w0, x)
        assert wm1.prec == 8
    finally:
        ctx.cap = old_cap


def test_acb_series_airy_modular_coulomb_fresnel_and_misc() -> None:
    old_cap = ctx.cap
    try:
        ctx.cap = 8
        x = acb_series([0, 1], prec=8)

        ai, aip, bi, bip = x.airy()
        assert is_close(ai, x.airy_ai())
        assert is_close(aip, x.airy_ai_prime())
        assert is_close(bi, x.airy_bi())
        assert is_close(bip, x.airy_bi_prime())

        t1, t2, t3, t4 = x.modular_theta(1j)
        assert t1.prec == 8
        assert t2.prec == 8
        assert t3.prec == 8
        assert t4.prec == 8

        f, g, hpos, hneg = x.coulomb(0, 0)
        assert is_close(f, x.coulomb_f(0, 0))
        assert str(g) == str(x.coulomb_g(0, 0))
        assert hpos.prec == 8
        assert hneg.prec == 8

        fs, fc = x.fresnel()
        assert is_close(fs, x.fresnel_s())
        assert is_close(fc, x.fresnel_c())
        assert x.fresnel(normalized=False)[0].prec == 8
        assert x.fresnel_s(normalized=False).prec == 8
        assert x.fresnel_c(normalized=False).prec == 8

        assert x.ei().prec == 8
        assert is_close_acb(x.si()[0], 0)
        assert x.ci().prec == 8
        assert is_close_acb(x.shi()[0], 0)
        assert x.chi().prec == 8
        assert (2 + x).li().prec == 8
        assert (2 + x).li(offset=True).prec == 8
    finally:
        ctx.cap = old_cap
