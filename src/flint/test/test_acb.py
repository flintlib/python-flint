"""Tests for python-flint's `acb` constructor."""

import math

from flint import acb, arb, dirichlet_char, fmpz_poly
from flint.test.helpers import is_close_acb, is_close_arb, raises


def test_acb_constructor() -> None:
    """Cover all branches of acb.__init__."""

    z0 = acb()
    assert z0.real == 0
    assert z0.imag == 0

    z1 = acb(acb(1, 2))
    assert z1.real == 1
    assert z1.imag == 2

    z2 = acb(3)
    assert z2.real == 3
    assert z2.imag == 0

    z3 = acb(4 + 5j)
    assert z3.real == 4
    assert z3.imag == 5

    class DummyMpc:
        def __init__(self) -> None:
            self._mpc_ = acb(6, 7)._mpc_

    z4 = acb(DummyMpc())
    assert z4.real == 6
    assert z4.imag == 7

    z5 = acb(8, 9)
    assert z5.real == 8
    assert z5.imag == 9

    assert raises(lambda: acb(1 + 2j, 3), ValueError)
    assert raises(lambda: acb(object()), TypeError)  # type: ignore[call-overload]
    assert raises(lambda: acb(1, object()), TypeError)  # type: ignore[call-overload]


def test_fmpz_poly_call_with_acb() -> None:
    p = fmpz_poly([1, 2])
    assert is_close_acb(p(acb(1, 2)), acb(3, 4))
    assert is_close_acb(p(acb(1)), acb(3))


def test_acb_is_zero() -> None:
    assert acb().is_zero() is True
    assert acb(0).is_zero() is True
    assert acb(0, 0).is_zero() is True
    assert acb(1).is_zero() is False
    assert acb(0, 1).is_zero() is False


def test_acb_is_finite_and_exact() -> None:
    assert acb(1, 2).is_finite() is True
    assert acb("1 +/- 0.1").is_finite() is True
    assert acb(float("inf")).is_finite() is False
    assert acb(float("nan")).is_finite() is False

    assert acb(1, 2).is_exact() is True
    assert acb("1 +/- 0.1").is_exact() is False
    assert acb(1, "0 +/- 0.1").is_exact() is False


def test_acb_pos() -> None:
    x = acb(2, 3)
    y = +x
    assert y == x
    assert y is not x


def test_acb_neg_dunder() -> None:
    x = acb(2, 3)
    assert -x == acb(-2, -3)


def test_acb_neg_method() -> None:
    x = acb(2, 3)
    assert x.neg(exact=True) == acb(2, 3)
    assert x.neg(exact=False) == acb(2, 3)


def test_acb_conjugate() -> None:
    x = acb(2, 3)
    assert x.conjugate(exact=True) == acb(2, -3)
    assert x.conjugate(exact=False) == acb(2, -3)


def test_acb_abs() -> None:
    x = acb(3, 4)
    y = abs(x)
    assert isinstance(y, arb)
    assert y == arb(5)


def test_acb_add() -> None:
    x = acb(2, 3)
    assert x + acb(5, 7) == acb(7, 10)
    assert x + 1 == acb(3, 3)
    assert raises(lambda: x + object(), TypeError)  # type: ignore[operator]


def test_acb_radd() -> None:
    x = acb(2, 3)
    assert 1 + x == acb(3, 3)
    assert acb(5, 7) + x == acb(7, 10)
    assert raises(lambda: object() + x, TypeError)  # type: ignore[operator]


def test_acb_sub() -> None:
    x = acb(2, 3)
    assert x - acb(5, 7) == acb(-3, -4)
    assert x - 1 == acb(1, 3)
    assert raises(lambda: x - object(), TypeError)  # type: ignore[operator]


def test_acb_rsub() -> None:
    x = acb(2, 3)
    assert 1 - x == acb(-1, -3)
    assert acb(5, 7) - x == acb(3, 4)
    assert raises(lambda: object() - x, TypeError)  # type: ignore[operator]


def test_acb_mul() -> None:
    x = acb(2, 3)
    assert x * acb(5, 7) == acb(-11, 29)
    assert x * 2 == acb(4, 6)
    assert raises(lambda: x * object(), TypeError)  # type: ignore[operator]


def test_acb_rmul() -> None:
    x = acb(2, 3)
    assert 2 * x == acb(4, 6)
    assert acb(5, 7) * x == acb(-11, 29)
    assert raises(lambda: object() * x, TypeError)  # type: ignore[operator]


def test_acb_truediv() -> None:
    x = acb(2, 3)
    assert x / 2 == acb(1, 1.5)
    assert is_close_acb(x / acb(2, -1), acb(0.2, 1.6))
    assert raises(lambda: x / object(), TypeError)  # type: ignore[operator]


def test_acb_rtruediv() -> None:
    x = acb(2, 3)
    assert is_close_acb(2 / x, acb(4 / 13, -6 / 13))
    assert is_close_acb(acb(2, -1) / x, acb(0.07692307692307693, -0.6153846153846154))
    assert raises(lambda: object() / x, TypeError)  # type: ignore[operator]


def test_acb_pow() -> None:
    x = acb(2, 3)
    assert x ** 2 == acb(-5, 12)
    assert x ** acb(2) == acb(-5, 12)
    assert raises(lambda: x ** object(), TypeError)  # type: ignore[operator]
    assert raises(lambda: pow(x, 2, 3), ValueError)  # type: ignore[misc]


def test_acb_rpow() -> None:
    x = acb(2, 3)
    y = 2 ** x
    assert is_close_acb(y, acb(2) ** x)
    assert raises(lambda: object() ** x, TypeError)  # type: ignore[operator]
    assert raises(lambda: pow(2, x, 3), ValueError)  # type: ignore[misc]


def test_acb_eq_ne() -> None:
    x = acb(2, 3)
    assert (x == acb(2, 3)) is True
    assert (x != acb(2, 3)) is False
    assert (x == acb(2, -3)) is False
    assert (x != acb(2, -3)) is True
    assert (x == (2 + 3j)) is True
    assert (x != (2 + 3j)) is False
    assert (x == 2) is False
    assert (x != 2) is True


def test_acb_eq_ne_intervals() -> None:
    a = acb("1 +/- 0.1", "2 +/- 0.2")
    b = acb("1 +/- 0.1", "2 +/- 0.2")
    c = acb(1, 2)
    d = acb("1 +/- 0.2", "2 +/- 0.2")
    e = acb("10 +/- 0.1", "2 +/- 0.2")

    assert (a == b) is False
    assert (a != b) is False
    assert (a == c) is False
    assert (a != c) is False
    assert (a == d) is False
    assert (a != d) is False
    assert (a == e) is False
    assert (a != e) is True


def test_acb_eq_ne_notimplemented_path() -> None:
    class Unsupported:
        pass

    x = acb(2, 3)
    assert (x == Unsupported()) is False
    assert (x != Unsupported()) is True


def test_acb_ordering_raises() -> None:
    x = acb(2, 3)
    y = acb(3, 4)
    assert raises(lambda: x < y, TypeError)  # type: ignore[operator]
    assert raises(lambda: x <= y, TypeError)  # type: ignore[operator]
    assert raises(lambda: x > y, TypeError)  # type: ignore[operator]
    assert raises(lambda: x >= y, TypeError)  # type: ignore[operator]


def test_acb_contains_methods() -> None:
    outer = acb("1 +/- 0.5", "2 +/- 0.5")
    inner = acb("1 +/- 0.25", "2 +/- 0.25")
    touching = acb("1.5 +/- 0.5", "2 +/- 0.5")
    disjoint = acb(4, 5)

    assert (inner in outer) is True
    assert (outer in inner) is False
    assert outer.contains(inner) is True
    assert outer.contains(1 + 2j) is True
    assert outer.contains_interior(inner) is True
    assert outer.contains_interior(outer) is False
    assert outer.overlaps(touching) is True
    assert outer.overlaps(disjoint) is False
    assert acb(2).contains_integer() is True
    assert acb(1.1).contains_integer() is False
    assert raises(lambda: outer.contains(object()), TypeError)  # type: ignore[arg-type]


def test_acb_union() -> None:
    a = acb(1, 2)
    b = acb("2 +/- 0.5", "3 +/- 0.5")
    u = a.union(b)
    assert u.contains(a)
    assert u.contains(b)


def test_acb_mid_rad_complex_rad() -> None:
    x = acb(arb(1, (1, -2)), arb(2, (1, -1)))
    assert x.mid() == acb(1, 2)
    assert is_close_arb(x.rad(), math.hypot(0.25, 0.5), tol=1e-8, rel_tol=1e-8)
    assert is_close_acb(x.complex_rad(), acb(0.25, 0.5), tol=1e-8, rel_tol=1e-8)


def test_acb_repr_str() -> None:
    x = acb(2)
    y = acb(0, 3)
    z = acb(2, -3)
    w = acb(2, 3)

    assert x.repr().startswith("acb(")
    assert ", " in z.repr()

    assert x.str() == x.real.str()
    assert y.str().endswith("j")
    assert " + " in w.str()
    assert " - " in z.str()


def test_acb_complex_dunder() -> None:
    x = acb(2, 3)
    assert complex(x) == complex(2, 3)


def test_acb_abs_methods() -> None:
    x = acb(3, -4)
    y1 = abs(x)
    y2 = x.__abs__()
    assert is_close_arb(y1, 5.0)
    assert is_close_arb(y2, 5.0)
    assert is_close_arb(x.abs_lower(), 5.0)
    assert is_close_arb(x.abs_upper(), 5.0)


def test_acb_csgn_sgn() -> None:
    assert is_close_arb(acb(2, 3).csgn(), 1.0)
    assert is_close_arb(acb(-1).csgn(), -1.0)

    assert is_close_acb(acb(-1).sgn(), acb(-1))
    assert is_close_acb(acb(0).sgn(), acb(0))
    assert is_close_acb(acb(5, 5).sgn(), acb(math.sqrt(0.5), math.sqrt(0.5)))


def test_acb_arg() -> None:
    assert is_close_arb(acb(3.3).arg(), 0.0)
    assert is_close_arb(acb(-1).arg(), math.pi)


def test_acb_elementary_transcendentals() -> None:
    x = acb(1, 2)
    y = acb(2, 3)
    q = acb(0.25)

    assert is_close_acb(x.pow(3), x ** 3)
    assert is_close_acb(x.log().exp(), x, tol=1e-8, rel_tol=1e-8)
    assert is_close_acb(x.log1p(), acb(1.039720770839918, 0.7853981633974483))
    assert is_close_acb(acb(2).asin(), acb(1.5707963267948966, -1.3169578969248166))
    assert is_close_acb(acb(2).acos(), acb(0, 1.3169578969248166))
    assert is_close_acb(x.atan(), acb(1.3389725222944935, 0.4023594781085251))
    assert is_close_acb(y.asinh(), acb(1.9686379257930963, 0.9646585044076028))
    assert is_close_acb(y.acosh(), acb(1.9833870299165354, 1.0001435424737972))
    assert is_close_acb(y.atanh(), acb(0.14694666622552975, 1.3389725222944935))
    assert is_close_acb(acb.pi(), acb(math.pi))
    assert is_close_acb(x.sqrt(), acb(1.272019649514069, 0.7861513777574233))
    assert is_close_acb(x.rsqrt(), acb(0.5688644810057831, -0.35157758425414293))
    assert is_close_acb(x.exp(), acb(-1.1312043837568137, 2.4717266720048188))
    assert is_close_acb(x.exp_pi_i(), acb("-0.001867442731707988814430213"))
    assert is_close_acb(acb("1e-8").expm1(), acb("1e-8"), tol=1e-12, rel_tol=1e-12)

    s1, c1 = x.sin_cos()
    assert is_close_acb(s1, x.sin())
    assert is_close_acb(c1, x.cos())
    sp, cp = x.sin_cos_pi()
    assert is_close_acb(sp, x.sin_pi())
    assert is_close_acb(cp, x.cos_pi())
    sh, ch = x.sinh_cosh()
    assert is_close_acb(sh, x.sinh())
    assert is_close_acb(ch, x.cosh())

    assert is_close_acb(x.cot(), 1 / x.tan(), tol=1e-8, rel_tol=1e-8)
    assert is_close_acb(x.tan_pi(), x.sin_pi() / x.cos_pi(), tol=1e-8, rel_tol=1e-8)
    assert is_close_acb(x.cot_pi(), x.cos_pi() / x.sin_pi(), tol=1e-8, rel_tol=1e-8)
    assert is_close_acb(y.sec(), 1 / y.cos(), tol=1e-8, rel_tol=1e-8)
    assert is_close_acb(y.csc(), 1 / y.sin(), tol=1e-8, rel_tol=1e-8)
    assert is_close_acb(y.tanh(), y.sinh() / y.cosh(), tol=1e-8, rel_tol=1e-8)
    assert is_close_acb(y.coth(), y.cosh() / y.sinh(), tol=1e-8, rel_tol=1e-8)
    assert is_close_acb(y.sech(), 1 / y.cosh(), tol=1e-8, rel_tol=1e-8)
    assert is_close_acb(y.csch(), 1 / y.sinh(), tol=1e-8, rel_tol=1e-8)

    assert is_close_acb(y.sinc(), y.sin() / y, tol=1e-8, rel_tol=1e-8)
    assert is_close_acb(q.sinc_pi(), q.sin_pi() / (acb.pi() * q), tol=1e-8, rel_tol=1e-8)


def test_acb_special_functions() -> None:
    x = acb(1, 2)
    y = acb(2, 3)

    assert is_close_acb(x.agm(), x.agm(acb(1)))
    assert is_close_acb(x.agm(acb(1.5)), acb(1.5).agm(x))
    assert is_close_acb(x.gamma(), acb(0.15190400267003614, 0.01980488016185498))
    assert is_close_acb(x.rgamma(), acb(6.4730736260191345, -0.8439438407732021))
    assert is_close_acb(x.lgamma(), acb(-1.8760787864309293, 0.12964631630978832))
    assert is_close_acb(x.digamma(), acb(0.7145915153739775, 1.3208072826422302))
    assert is_close_acb(acb(0.5, 1000).zeta(), acb(0.35633436719439606, 0.9319978312329937))
    assert is_close_acb(x.zeta(acb(2, 3)), acb(-2.953059572088557, 3.4109625245120507))
    assert is_close_acb(x.lerch_phi(3, 4), acb(0.006872751459699249, 0.011125353146863518))
    assert is_close_acb(y.erf(), acb(-20.829461427614568, 8.687318271470163))
    assert is_close_acb(acb("77.7").erfc(), acb("7.929310690520378873143053e-2625"), tol=arb("1e-2630"), rel_tol=1e-8)
    assert is_close_acb(acb(10).erfi(), acb("1.524307422708669699360547e+42"), tol=1e20, rel_tol=1e-8, max_width=1e30)
    assert is_close_acb(y.gamma_upper(1 + 2j), acb(0.02614303924198793, -0.0007536537278463329))
    assert is_close_acb(y.gamma_lower(2.5), acb(1.6460770101348767, 1.140585862703101))
    assert is_close_acb(y.beta_lower(1, 2.5), acb(0.2650137734913867, -7.111836702381955))
    assert is_close_acb(y.expint(1 + 2j), acb(-0.014426614955270803, 0.019423483729866873))
    assert is_close_acb(acb(10).ei(), acb("2492.228976241877759138440"), tol=1e-8, rel_tol=1e-8)
    assert is_close_acb(acb(10).si(), acb("1.658347594218874049330972"))
    assert is_close_acb(acb(10).ci(), acb("-0.04545643300445537263453283"))
    assert is_close_acb(acb(10).shi(), acb("1246.114490199423344411882"), tol=1e-8, rel_tol=1e-8)
    assert is_close_acb(acb(10).chi(), acb("1246.114486042454414726558"), tol=1e-8, rel_tol=1e-8)
    assert is_close_acb(acb(10).li(), acb("6.165599504787297937522982"))
    assert is_close_acb(acb(10).li(offset=True), acb("5.120435724669805152678393"))
    assert is_close_acb(x.rising(5), acb(-540, -100))
    r1, r2 = x.rising2(5)
    assert is_close_acb(r1, acb(-540, -100))
    assert is_close_acb(r2, acb(-666, 420))
    assert is_close_acb(acb(3).polylog(2), acb(2.3201804233130984, -3.4513922952232027))
    z = acb(-1, 1)
    assert is_close_acb(z.airy_ai(), acb(0.8221174265552726, -0.11996634266442434))
    assert is_close_acb(z.airy_bi(), acb(0.21429040153487357, 0.6739169237227052))
    ai, aip, bi, bip = y.airy()
    assert is_close_acb(ai, y.airy_ai())
    assert is_close_acb(aip, y.airy_ai(derivative=1))
    assert is_close_acb(bi, y.airy_bi())
    assert is_close_acb(bip, y.airy_bi(derivative=1))
    assert is_close_acb(acb(1).lambertw(), acb("0.5671432904097838729999687"))
    assert is_close_acb(acb(3).fresnel_s(), acb("0.4963129989673750360976123"))
    assert is_close_acb(acb(3).fresnel_c(), acb("0.6057207892976856295561611"))
    assert is_close_acb(acb(2).bessel_j(1), acb("0.5767248077568733872024482"))
    assert is_close_acb(acb(2).bessel_y(1), acb("-0.1070324315409375468883708"))
    assert is_close_acb(acb(2).bessel_k(1), acb("0.1398658818165224272845988"))
    assert is_close_acb(acb(2).bessel_i(1), acb("1.590636854637329063382254"))


def test_acb_accuracy_bits() -> None:
    x = acb(1, 2)
    assert isinstance(x.rel_accuracy_bits(), int)
    assert isinstance(x.rel_one_accuracy_bits(), int)
    assert x.bits() >= 1
    assert acb("2047/2048").bits() == 11
    assert acb(float("nan")).bits() == 0


def test_acb_real_methods() -> None:
    assert is_close_acb(acb(-5, 2).real_abs(), acb(5, -2))
    assert is_close_acb(acb(-5, 2).real_sgn(), acb(-1))
    assert is_close_acb(acb(0).real_sgn(), acb(0))
    assert acb(0).real_sgn(analytic=True).is_finite() is False

    assert is_close_acb(acb(-5, 2).real_heaviside(), acb(0))
    assert is_close_acb(acb(5, 2).real_heaviside(), acb(1))
    assert is_close_acb(acb(0).real_heaviside(), acb(0.5))
    assert acb(0).real_heaviside(analytic=True).is_finite() is False

    assert is_close_acb(acb(2.7).real_floor(), acb(2))
    assert is_close_acb(acb(2.2).real_ceil(), acb(3))
    assert is_close_acb(acb(2, 10).real_max(acb(3, -1)), acb(3, -1))
    assert is_close_acb(acb(2, 10).real_min(acb(3, -1)), acb(2, 10))
    assert is_close_acb(acb(9).real_sqrt(), acb(3))


def test_acb_elliptic_methods() -> None:
    assert is_close_acb(2 * acb(0).elliptic_k(), acb.pi())
    assert is_close_acb(2 * acb(0).elliptic_e(), acb.pi())

    assert is_close_acb(acb.elliptic_rf(1, 2 + 3j, 3 + 4j), acb(0.5577655465453922, -0.2202042457195556))
    assert is_close_acb(acb.elliptic_rc(1, 2 + 3j), acb(0.5952169239306156, -0.23879819090905094))
    assert is_close_acb(acb.elliptic_rj(1, 2, 1 + 2j, 2 + 3j), acb(0.1604659632144333, -0.2502751672723324))
    assert is_close_acb(acb.elliptic_rd(1, 2, 1 + 2j), acb(0.20437225103026298, -0.3559745898273716))
    assert is_close_acb(acb.elliptic_rg(1, 2, 1 + 2j), acb(1.2065571680567230, 0.27521766887077397))

    assert is_close_acb(acb.elliptic_f(2, 0.75), acb(2.9525696736557795))
    assert is_close_acb(acb.elliptic_e_inc(2, 0.75), acb(1.4434330690994616))
    assert is_close_acb(acb.elliptic_pi(0.25, 0.125), acb(1.8793494518796038))
    assert is_close_acb(acb.elliptic_pi_inc(0.25, 0.5, 0.125), acb(0.5128718023282086))
    assert is_close_acb(acb.elliptic_pi_inc(0.25, 0.5, 0.125, pi=True), acb.elliptic_pi(0.25, 0.125))

    z = acb("1/3", "1/5")
    wp = z.elliptic_p(1j)
    assert is_close_acb(wp, acb(3.6863806460788798, -4.5914983714972594))
    assert is_close_acb(z.elliptic_zeta(1j), acb(2.2196803395084187, -1.5049479257552417))
    assert is_close_acb(z.elliptic_sigma(1j), acb(0.33965494971368865, 0.19706907623509313))
    assert is_close_acb(wp.elliptic_inv_p(1j), z, tol=1e-9, rel_tol=1e-9, max_width=1e-9)

    e1, e2, e3 = acb(0.5 + 1j).elliptic_roots()
    assert is_close_acb(e1 + e2 + e3, acb(0), tol=1e-9, rel_tol=1e-9, max_width=1e-9)
    assert is_close_acb(e1, acb(6.2853881186694))

    g2, g3 = acb(0.5 + 1j).elliptic_invariants()
    assert is_close_acb(g2, acb(72.64157667926128))
    assert is_close_acb(g3, acb(536.6642788346023))


def test_acb_hypgeom_functions() -> None:
    z = acb(5)
    assert is_close_acb(z.hypgeom_2f1(1, 2, 3), acb(-0.5109035488895912, -0.25132741228718347))
    assert is_close_acb(z.hypgeom_2f1(1, 2, 3, regularized=True), acb(-0.2554517744447956, -0.12566370614359174))

    w = acb(0.2)
    assert is_close_acb(w.hypgeom_0f1(2.5), acb(1.0823198864627344))
    assert is_close_acb(w.hypgeom_0f1(2.5, regularized=True), acb(0.8141781413451533))
    assert is_close_acb(w.hypgeom([1, 2], [3]), acb(1.1571775657104878))
    assert is_close_acb(w.hypgeom([1, 2], [3], regularized=True), acb(0.5785887828552439))
    assert is_close_acb(w.hypgeom_u(1, 2.5), acb(11.378859082050016))
    assert is_close_acb(w.hypgeom_1f1(1, 2.5), acb(1.0847822248780032))
    assert is_close_acb(w.hypgeom_1f1(1, 2.5, regularized=True), acb(0.8160304422585721))


def test_acb_modular_functions() -> None:
    z = acb(1 + 1j)
    tau = acb(1.25 + 3j)

    t1, t2, t3, t4 = acb.modular_theta(z, tau)
    assert is_close_acb(t1, acb(1.8202359101249896, -1.216251950154478))
    assert is_close_acb(t2, acb(-1.2207902675769677, -1.8270555167911547))
    assert is_close_acb(t3, acb(0.9694430387796704, -0.030556961208168033))
    assert is_close_acb(t4, acb(1.0305569611960065, 0.030556961208168033))

    assert is_close_acb(acb(1 + 1j).modular_eta(), acb(0.7420487758365647, 0.19883137022991072))
    assert is_close_acb((1 + acb(-163).sqrt() / 2).modular_j(), acb(262537412640769488.0), max_width=1e5)
    assert is_close_acb(acb(0.25 + 5j).modular_lambda(), acb(1.7049954156680393e-6, 1.7049925086620794e-6))
    assert is_close_acb(acb(0.25 + 5j).modular_delta(), acb(1.2378960150102817e-26, 2.271101068324094e-14))

    u1, u2, u3, u4 = acb.modular_theta(z, tau, 1)
    assert is_close_acb(u1, acb(-3.8353056542516, -5.7398107897127), tol=1e-9, rel_tol=1e-9)
    assert is_close_acb(u2, acb(-5.7184931625874, 3.8208882734627), tol=1e-9, rel_tol=1e-9)
    assert is_close_acb(u3, acb(-0.19199371059495, 0.19199371074778), tol=1e-9, rel_tol=1e-9)
    assert is_close_acb(u4, acb(0.19199371059495, -0.19199371044212), tol=1e-9, rel_tol=1e-9)


def test_acb_orthogonal_and_related_functions() -> None:
    x = acb(1) / 3
    assert is_close_acb(x.chebyshev_t(3), acb(-0.8518518518518519))
    assert is_close_acb(x.chebyshev_u(3), acb(-1.037037037037037))
    assert is_close_acb(x.jacobi_p(5, 0.25, 0.5), acb(0.4131944444444444))
    assert is_close_acb(x.gegenbauer_c(5, 0.25), acb(0.13218557098765432))
    assert is_close_acb(x.laguerre_l(5, 0.25), acb(0.03871323490012003))
    assert is_close_acb(x.hermite_h(5), acb(34.20576131687243))
    assert is_close_acb(x.legendre_p(5), acb(1 / 3))
    assert is_close_acb(x.legendre_q(5), acb(0.1655245300933242))
    assert is_close_acb(
        acb(3).legendre_q(5, 1.5, type=3), acb(0, -0.0003010942389043591), tol=1e-8, rel_tol=1e-8, max_width=1e-8
    )
    assert is_close_acb(acb.spherical_y(5, 3, 0.25, 0.75), acb(0.02451377199072374, -0.03036343496553117))
    assert raises(lambda: x.legendre_p(5, type=4), ValueError)
    assert raises(lambda: x.legendre_q(5, type=4), ValueError)


def test_acb_misc_special_functions() -> None:
    assert is_close_acb(acb(5).dirichlet_eta(), acb(0.9721197704469093))
    assert is_close_acb(acb(1).dirichlet_eta(), acb(math.log(2.0)))
    assert is_close_acb(acb(2).polygamma(1), acb(0.6449340668482264))
    assert is_close_acb(acb(3).log_barnes_g(), acb(0))
    assert is_close_acb(acb(3).barnes_g(), acb(1))
    assert is_close_acb(acb.stieltjes(1), acb(-0.07281584548367672))
    assert raises(lambda: acb.stieltjes(-1), ValueError)
    assert is_close_acb(acb(0.25 + 0.25j).bernoulli_poly(5), acb(-0.05859375, 0.006510416666666667))
    assert is_close_acb(acb(5 + 2j).log_sin_pi(), acb(5.590034639271204, -14.13716694115407))

    chi = dirichlet_char(3, 1)
    assert is_close_acb(acb(2).dirichlet_l(chi), acb(1.46216361497620))
    assert is_close_acb(acb(2).dirichlet_l((3, 1)), acb(1.46216361497620))

    vals = acb.dft(range(1, 12))
    back = acb.dft(vals, inverse=True)
    for i, v in enumerate(back, start=1):
        assert is_close_acb(v, acb(i), tol=1e-9, rel_tol=1e-9, max_width=1e-9)
    assert acb.dft([]) == []

    assert acb("5 +/- 0.1").unique_fmpz() is not None
    assert acb(5.5, 0.1).unique_fmpz() is None


def test_acb_more_special_function_branches() -> None:
    z = acb(-1, 1)
    assert is_close_acb(z.airy_ai(derivative=1), acb(-0.3790604792268335, -0.6045001308622461))
    assert is_close_acb(z.airy_bi(derivative=1), acb(0.8344734885227826, -0.3465260632668285))
    assert raises(lambda: z.airy_ai(derivative=2), ValueError)
    assert raises(lambda: z.airy_bi(derivative=2), ValueError)

    assert is_close_acb(acb(-5, "+/- 1e-20").lambertw(left=True), acb(0.844844605432170, 1.97500875488903))
    assert is_close_acb(acb(-0.25, "+/- 1e-20").lambertw(middle=True), acb(-2.15329236411035))

    assert is_close_acb(acb(5).bessel_k(1, scaled=True), acb(0.6002738587883126))
    assert is_close_acb(acb(5).bessel_i(1, scaled=True), acb(0.16397226694454236))
    assert is_close_acb(acb(-1).root(3), acb(0.5, 0.8660254037844386))


def test_acb_zeta_zero_methods() -> None:
    z1 = acb.zeta_zero(1)
    assert is_close_acb(z1, acb(0.5, 14.134725141734694))
    zs = acb.zeta_zeros(1000, 3)
    assert len(zs) == 3
    assert is_close_acb(zs[0], acb(0.5, 1419.42248094600))
    assert is_close_acb(zs[1], acb(0.5, 1420.41652632375))
    assert is_close_acb(zs[2], acb(0.5, 1421.85056718705))
    assert raises(lambda: acb.zeta_zero(0), ValueError)
    assert raises(lambda: acb.zeta_zeros(0, 1), ValueError)
    assert raises(lambda: acb.zeta_zeros(1, -1), ValueError)


def test_acb_integral() -> None:
    res = acb.integral(lambda z, d: z, 0, 1)
    assert is_close_acb(res, acb(0.5), tol=1e-9, rel_tol=1e-9, max_width=1e-9)
    assert is_close_acb(acb.integral(lambda z, d: z, 0, 1, rel_tol=arb("1e-20"), abs_tol=arb("1e-20")), acb(0.5))
    assert is_close_acb(
        acb.integral(lambda z, d: z, 0, 1, deg_limit=8, eval_limit=1000, depth_limit=30, use_heap=True, verbose=False),
        acb(0.5),
    )
    assert raises(lambda: acb.integral(lambda z, d: 1, 0, 1), TypeError)  # type: ignore[arg-type]


def test_acb_coulomb_functions() -> None:
    x = acb(1)
    l = acb(0.5)
    eta = acb(0.25)
    f, g, hp, hn = x.coulomb(l, eta)
    assert is_close_acb(f, x.coulomb_f(l, eta), tol=1e-9, rel_tol=1e-9, max_width=1e-9)
    assert is_close_acb(g, x.coulomb_g(l, eta), tol=1e-9, rel_tol=1e-9, max_width=1e-9)
    assert is_close_acb(hp, g + 1j * f, tol=1e-9, rel_tol=1e-9, max_width=1e-9)
    assert is_close_acb(hn, g - 1j * f, tol=1e-9, rel_tol=1e-9, max_width=1e-9)
    assert is_close_acb(f, acb(0.4283180781043542))
    assert is_close_acb(g, acb(1.218454487206368))


def test_acb_hypgeom_remaining_branches() -> None:
    x = acb("11/10")
    a = acb(2).sqrt()
    b = acb("1/2")
    c = a + acb("3/2")

    assert x.hypgeom_2f1(a, b, c, ab=True).is_finite() is False
    assert x.hypgeom_2f1(a, b, c, ac=True).is_finite() is False
    assert x.hypgeom_2f1(a, b, c, bc=True).is_finite() is False
    assert is_close_acb(x.hypgeom_2f1(a, b, c, abc=True), acb(1.8017826594800542, -0.3114019850045849))

    w = acb(0.2)
    assert is_close_acb(w.hypgeom([1, 2], [3], n=5), acb(1.1571775657104878), tol=1e-3, rel_tol=1e-3, max_width=1e-3)
    assert raises(lambda: w.hypgeom([1, 2], [3], regularized=True, n=5), NotImplementedError)
    assert is_close_acb(acb(-30).hypgeom_u(1 + 1j, 2 + 3j, n=30, asymp=True), acb(0.7808944974, -0.2674783065), tol=1e-8, rel_tol=1e-8)
