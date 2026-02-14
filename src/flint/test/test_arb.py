"""Test for python-flint's `arb` type."""

from __future__ import annotations

from typing import Callable

import math

from flint import arb, ctx

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
    tol: float = 1e-10,
    rel_tol: float = 1e-10,
    max_width: float = 1e-10,
) -> bool:
    y = arb(y)
    return (
        isinstance(x, arb)
        and x.rad() < max_width
        and y.rad() < max_width
        and abs(x - y) <= max(rel_tol * max(abs(x), abs(y)), tol)
    )

def assert_almost_equal(x: float | int, y: float | int, places: int = 7) -> None:
    """Helper method for approximate comparisons."""
    assert round(x-y, ndigits=places) == 0

def test_from_int() -> None:
    """Tests instantiating `arb`s from ints."""
    for val in [
        -42 * 10**9,
        -42 * 10**7,
        -42,
        0,
        42,
        42 * 10**7,
        42 * 10**9,
        42 * 10**11,
    ]:
        x = arb(val)
        man, exp = x.man_exp()
        assert (man * 2**exp) == val

def test_from_float() -> None:
    """Tests instantiating `arb`s from floats."""
    for val in [0.0, 1.1, -1.1, 9.9 * 0.123, 99.12]:
        x = arb(val)
        man, exp = x.man_exp()
        assert (int(man) * 2 ** int(exp)) == val

def test_from_float_inf() -> None:
    """Tests `arb` works with +/- inf."""
    posinf = arb(float("inf"))
    neginf = arb(float("-inf"))

    assert posinf.is_finite() is False
    assert neginf.is_finite() is False
    assert float(posinf) == float("inf")
    assert float(neginf) == float("-inf")

def test_from_man_exp() -> None:
    """Tests instantiating `arb`s with mantissa and exponent."""
    for man, exp in [(2, 30), (4, 300), (5 * 10**2, 7**8)]:
        x = arb(mid=(man, exp))
        m, e = x.man_exp()
        assert (m * 2**e) == (man * 2**exp)

def test_from_midpoint_radius() -> None:
    """Tests instantiating `arb`s with midpoint and radius."""
    for mid, rad in [(10, 1), (10000, 5), (10, 1), (10, 1)]:
        mid_arb = arb(mid)
        rad_arb = arb(rad)
        x = arb(mid_arb, rad_arb)
        assert x.mid() == mid_arb
        actual_radius = float(x.rad())
        assert_almost_equal(actual_radius, rad)

def test_is_exact() -> None:
    """Tests `arb.is_exact`."""
    for arb_val, exact in [
        (arb(10), True),
        (arb(0.01), True),
        (arb(-float("inf")), True),
        (arb(1, 0), True),
        (arb(1, 1), False),
    ]:
        assert arb_val.is_exact() == exact

def test_is_finite() -> None:
    """Tests `arb.is_finite`."""
    assert not (arb(-float("inf")).is_finite())
    assert not (arb(float("inf")).is_finite())
    assert (arb(10).is_finite())

def test_is_nan() -> None:
    """Tests `arb.is_nan`."""
    assert (arb(float("nan")).is_nan())
    assert not (arb(0.0).is_nan())

def test_lower() -> None:
    """Tests `arb.lower`."""
    with ctx.workprec(100):
        arb_val = arb(1, 0.5)
        assert_almost_equal(float(arb_val.lower()), 0.5)

def test_upper() -> None:
    """Tests `arb.upper`."""
    with ctx.workprec(100):
        arb_val = arb(1, 0.5)
        assert_almost_equal(float(arb_val.upper()), 1.5)

def test_contains() -> None:
    """`y.__contains__(x)` returns True iff every number in `x` is also in `y`."""
    for x, y, expected in [
        (
            arb(mid=9, rad=1),
            arb(mid=10, rad=2),
            True,
        ),
        (
            arb(mid=10, rad=2),
            arb(mid=9, rad=1),
            False,
        ),
        (arb(10), arb(mid=9, rad=1), True),
        (arb(10.1), arb(mid=9, rad=1), False),
    ]:
        assert (x in y) == expected

def test_hash() -> None:
    """`x` and `y` hash to the same value if they have the same midpoint and radius.

    Args:
        x: An arb.
        y: An arb.
        expected: Whether `x` and `y` should hash to the same value.
    """
    def arb_pi(prec: int) -> arb:
        """Helper to calculate arb to a given precision."""
        with ctx.workprec(prec):
            return arb.pi()
    for x, y, expected in [
        (arb(10), arb(10), True),
        (arb(10), arb(11), False),
        (arb(10.0), arb(10), True),
    ]:
        assert (hash(x) == hash(y)) == expected

    for x in [
        arb(mid=10, rad=2),
        arb_pi(100),
    ]:
        try:
            hash(x)
        except ValueError:
            pass
        else:
            assert False, f"Expected {x} to raise an error if hashed, but succeeded."

def test_arb_constructor() -> None:
    """Cover constructor conversion failures and string parsing branches."""
    assert raises(lambda: arb(object()), TypeError)  # type: ignore[arg-type]
    assert raises(lambda: arb((1, "bad")), TypeError)  # type: ignore[arg-type]
    assert raises(lambda: arb("not-a-number"), ValueError)

    x = arb("1.5 ± 0.25")
    assert x.contains(arb(1.5))
    assert x.contains(arb(1.25))
    assert x.contains(arb(1.75))

def test_arb_contains() -> None:
    """Cover contains/overlaps helper methods and real/imag properties."""
    outer = arb(10, 2)      # [8, 12]
    inner = arb(10, 1)      # [9, 11]
    touching = arb(12, 1)   # [11, 13]
    disjoint = arb(20, 0.5) # [19.5, 20.5]

    assert outer.contains(inner)
    assert outer.contains_interior(arb(10, 0.5))
    assert not outer.contains_interior(outer)
    assert outer.overlaps(touching)
    assert not outer.overlaps(disjoint)

    assert arb(1.1, 0.2).contains_integer()
    assert not arb(1.1, 0.05).contains_integer()

    x = arb(3, 0.1)
    assert x.real is x
    assert x.imag == arb(0)

def test_arb_comparison() -> None:
    x = arb(2)
    y = arb(3)
    assert (x < y) is True
    assert (x > y) is False
    assert (x <= y) is True
    assert (x >= y) is False
    assert (x < x) is False
    assert (x <= x) is True
    assert (y < x) is False
    assert (y <= x) is False
    assert (y > x) is True
    assert (y >= x) is True
    assert (x != y) is True
    assert (x == arb(2)) is True

    xi = arb('2 +/- 1')
    assert (xi < x) is False
    assert (xi > x) is False

def test_arb_comparison_scalars() -> None:
    x = arb(2)
    assert (x == 2) is True
    assert (x == 2.0) is True
    assert (x < 3) is True
    assert (x <= 2) is True
    assert (x > 1) is True
    assert (x >= 2.0) is True
    assert (3 > x) is True
    assert (2 == x) is True

def test_arb_comparison_invalid_type() -> None:
    x = arb(2)
    assert raises(lambda: x < "bad", TypeError)  # type: ignore[operator]
    assert raises(lambda: x <= "bad", TypeError)  # type: ignore[operator]
    assert raises(lambda: x > "bad", TypeError)  # type: ignore[operator]
    assert raises(lambda: x >= "bad", TypeError)  # type: ignore[operator]
    assert raises(lambda: "bad" < x, TypeError)  # type: ignore[operator]
    assert raises(lambda: "bad" <= x, TypeError)  # type: ignore[operator]
    assert raises(lambda: "bad" > x, TypeError)  # type: ignore[operator]
    assert raises(lambda: "bad" >= x, TypeError)  # type: ignore[operator]

def test_arb_pos_floor_ceil() -> None:
    x = arb(2.25, 0.25)  # [2.0, 2.5]
    assert x in (+x)
    assert (+x) in x
    assert (arb(2.9).floor() == arb(2)) is True
    assert (arb(2.1).ceil() == arb(3)) is True

def test_arb_abs_bounds() -> None:
    x = arb(-5, 2)  # [-7, -3]
    lower = x.abs_lower()
    upper = x.abs_upper()
    assert float(lower) <= 3.0
    assert float(upper) >= 7.0
    assert (lower <= upper) is True

def test_arb_exact_conversions_and_string_options() -> None:
    x = arb(3)
    assert x.is_integer() is True
    assert int(x.fmpz()) == 3
    assert x.fmpq() == 3

    y = arb(1, 1)
    assert raises(lambda: y.man_exp(), ValueError)
    assert raises(lambda: y.fmpq(), ValueError)
    assert raises(lambda: y.fmpz(), ValueError)

    assert raises(lambda: arb(float("inf")).man_exp(), ValueError)
    assert raises(lambda: arb(float("nan")).fmpq(), ValueError)

    z = arb("1.5 +/- 0.25")
    assert z.repr().startswith("arb(")
    assert "+/-" in z.str(10)
    assert "+/-" not in z.str(10, radius=False)
    assert isinstance(z.str(10, more=True), str)
    assert isinstance(z.str(10, condense=2), str)

    oldunicode = ctx.unicode
    try:
        ctx.unicode = True
        assert "±" in z.str(10)
    finally:
        ctx.unicode = oldunicode

def test_arb_mpf_mid_rad_and_interval_ops() -> None:
    t_pos = arb(3)._mpf_
    t_neg = arb(-3)._mpf_
    t_inf = arb(float("inf"))._mpf_
    assert t_pos[0] == 0
    assert t_neg[0] == 1
    assert t_inf[3] == -1

    mid, rad, exp = (arb(1) / 3).mid_rad_10exp(12)
    assert int(rad) >= 0
    assert isinstance(int(mid), int)
    assert isinstance(int(exp), int)

    a = arb(1, 1)  # [0,2]
    b = arb(2, 1)  # [1,3]
    c = a.intersection(b)
    assert c.contains(arb(1.5))
    assert raises(lambda: arb(2).intersection(3), ValueError)

    nn = arb(-1, 2).nonnegative_part()
    assert nn.contains(0)
    assert nn.contains(1)

def test_arb_arithmetic() -> None:
    x = arb(2)
    assert (x + 1) == arb(3)
    assert (1 + x) == arb(3)
    assert (x - 1) == arb(1)
    assert (5 - x) == arb(3)
    assert (x * 3) == arb(6)
    assert (3 * x) == arb(6)
    assert (x / 2) == arb(1)
    assert (6 / x) == arb(3)
    assert (x ** 2) == arb(4)
    assert (2 ** x) == arb(4)

    assert raises(lambda: x + "bad", TypeError)  # type: ignore[operator]
    assert raises(lambda: "bad" + x, TypeError)  # type: ignore[operator]
    assert raises(lambda: x ** "bad", TypeError)  # type: ignore[operator]
    assert raises(lambda: "bad" ** x, TypeError)  # type: ignore[operator]

def test_arb_functions() -> None:
    x = arb(0.5)
    xpi = arb(0.25)
    y = arb(2)

    with ctx.workprec(53):
        assert is_close_arb(y.sqrt(), 1.4142135623730951)
        assert is_close_arb(y.rsqrt(), 1.0 / math.sqrt(2.0))
        assert is_close_arb(x.expm1(), math.expm1(0.5))
        assert is_close_arb(y.log(), math.log(2.0))
        assert is_close_arb(x.log1p(), math.log1p(0.5))

        assert is_close_arb(x.sin(), math.sin(0.5))
        assert is_close_arb(x.cos(), math.cos(0.5))
        s, c = x.sin_cos()
        assert is_close_arb(s, math.sin(0.5))
        assert is_close_arb(c, math.cos(0.5))
        assert is_close_arb(xpi.sin_pi(), math.sin(math.pi * 0.25))
        assert is_close_arb(xpi.cos_pi(), math.cos(math.pi * 0.25))
        s, c = xpi.sin_cos_pi()
        assert is_close_arb(s, math.sin(math.pi * 0.25))
        assert is_close_arb(c, math.cos(math.pi * 0.25))
        assert is_close_arb(x.tan(), math.tan(0.5))
        assert is_close_arb(x.cot(), 1.0 / math.tan(0.5))
        assert is_close_arb(xpi.tan_pi(), math.tan(math.pi * 0.25))
        assert is_close_arb(xpi.cot_pi(), 1.0 / math.tan(math.pi * 0.25))
        assert is_close_arb(x.sec(), 1.0 / math.cos(0.5))
        assert is_close_arb(x.csc(), 1.0 / math.sin(0.5))
        assert is_close_arb(x.sinc(), math.sin(0.5) / 0.5)
        assert is_close_arb(x.sinc_pi(), math.sin(math.pi * 0.5) / (math.pi * 0.5))

        assert is_close_arb(x.atan(), math.atan(0.5))
        assert is_close_arb(x.acos(), math.acos(0.5))
        assert is_close_arb(x.asin(), math.asin(0.5))
        assert is_close_arb(x.atanh(), math.atanh(0.5))
        assert is_close_arb(x.asinh(), math.asinh(0.5))
        assert is_close_arb(y.acosh(), math.acosh(2.0))
        assert is_close_arb(x.sinh(), math.sinh(0.5))
        assert is_close_arb(x.cosh(), math.cosh(0.5))
        s, c = x.sinh_cosh()
        assert is_close_arb(s, math.sinh(0.5))
        assert is_close_arb(c, math.cosh(0.5))
        assert is_close_arb(x.tanh(), math.tanh(0.5))
        assert is_close_arb(x.coth(), 1.0 / math.tanh(0.5))
        assert is_close_arb(x.sech(), 1.0 / math.cosh(0.5))
        assert is_close_arb(x.csch(), 1.0 / math.sinh(0.5))

def test_arb_constants_and_special_values() -> None:
    with ctx.workprec(53):
        assert is_close_arb(arb.const_e(), math.e)
        assert is_close_arb(arb.const_log2(), math.log(2.0))
        assert is_close_arb(arb.const_log10(), math.log(10.0))

    posinf = arb.pos_inf()
    neginf = arb.neg_inf()
    nan = arb.nan()

    assert (posinf.is_finite()) is False
    assert (neginf.is_finite()) is False
    assert (nan.is_nan()) is True
    assert float(posinf) == float("inf")
    assert float(neginf) == float("-inf")

def test_arb_unique_fmpz_and_bits() -> None:
    u1 = arb("5 +/- 0.1").unique_fmpz()
    u2 = arb("5 +/- 0.9").unique_fmpz()
    assert u1 is not None
    assert u2 is not None
    assert int(u1) == 5
    assert int(u2) == 5
    assert arb("5.1 +/- 0.9").unique_fmpz() is None
    assert arb("2047/2048").bits() == 11
    assert arb.nan().bits() == 0

def test_arb_atan2_log_base_and_lambertw() -> None:
    with ctx.workprec(53):
        assert is_close_arb(arb.atan2(arb(3), arb(2)), math.atan2(3.0, 2.0))
        assert is_close_arb(arb(8).log_base(2), 3.0)

        w = arb(1).lambertw()
        assert is_close_arb(w * w.exp(), 1.0, tol=1e-9, rel_tol=1e-9, max_width=1e-9)

    wneg = arb("-0.1").lambertw(-1)
    assert (wneg < -1) is True
    assert raises(lambda: arb(1).lambertw(2), ValueError)

def test_arb_special_functions() -> None:
    with ctx.workprec(53):
        assert is_close_arb(arb(5).gamma(), 24.0)
        assert is_close_arb(arb.gamma_fmpq(arb(1).fmpq() / 2), math.sqrt(math.pi))
        assert is_close_arb(arb(5).rgamma(), 1.0 / 24.0)
        assert is_close_arb(arb(5).lgamma(), math.lgamma(5.0))
        assert is_close_arb(arb(1).digamma(), -0.5772156649015329, tol=1e-9, rel_tol=1e-9, max_width=1e-9)

        assert is_close_arb(arb(3).rising(3), 60.0)
        assert is_close_arb(arb.rising_fmpq_ui(arb(1).fmpq() / 2, 2), 0.75)
        u, v = arb(3).rising2(5)
        assert is_close_arb(u, 2520.0)
        assert is_close_arb(v, 2754.0)

        assert is_close_arb(arb(2).zeta(), math.pi**2 / 6)
        assert is_close_arb(arb(2).zeta(3), "0.3949340668482264364724152")
        assert is_close_arb(arb(2).agm(3), "2.474680436236304462606658")

        assert is_close_arb(arb.bernoulli(2), 1.0 / 6.0)
        assert is_close_arb(arb.bell_number(6), 203.0)
        assert is_close_arb(arb.partitions_p(5), 7.0)
        assert is_close_arb(arb(arb(1) / 3).bernoulli_poly(3), "0.03703703703703703703703716")

        assert is_close_arb(arb(5).fac(), 120.0)
        assert is_close_arb(arb.fac_ui(6), 720.0)
        assert is_close_arb(arb(10).bin(3), 120.0)
        assert is_close_arb(arb.bin_uiui(10, 3), 120.0)
        assert is_close_arb(arb.fib(10), 55.0)
        assert is_close_arb(arb(0.5).polylog(2), "0.5822405264650125059026562")

        assert is_close_arb(arb(1).airy_ai(), "0.1352924163128814155241473")
        assert is_close_arb(arb(1).airy_ai(derivative=1), "-0.1591474412967932127875002")
        assert raises(lambda: arb(1).airy_ai(derivative=2), ValueError)

        assert is_close_arb(arb(1).airy_bi(), "1.207423594952871259436378")
        assert is_close_arb(arb(1).airy_bi(derivative=1), "0.9324359333927756329594508")
        assert raises(lambda: arb(1).airy_bi(derivative=2), ValueError)

        ai, aip, bi, bip = arb(1).airy()
        assert is_close_arb(ai, "0.1352924163128814155241473")
        assert is_close_arb(aip, "-0.1591474412967932127875002")
        assert is_close_arb(bi, "1.207423594952871259436378")
        assert is_close_arb(bip, "0.9324359333927756329594508")

        assert is_close_arb(arb.airy_ai_zero(1), "-2.338107410459767038489195")
        assert is_close_arb(arb.airy_ai_zero(1, derivative=1), "-1.018792971647471089017324")
        assert raises(lambda: arb.airy_ai_zero(0), ValueError)
        assert raises(lambda: arb.airy_ai_zero(1, derivative=2), ValueError)

        assert is_close_arb(arb.airy_bi_zero(1), "-1.173713222709127924919979")
        assert is_close_arb(arb.airy_bi_zero(1, derivative=1), "-2.294439682614123246622457")
        assert raises(lambda: arb.airy_bi_zero(0), ValueError)
        assert raises(lambda: arb.airy_bi_zero(1, derivative=2), ValueError)

        x = arb(1) / 3
        assert is_close_arb(x.chebyshev_t(3), "-0.8518518518518518518518510")
        assert is_close_arb(x.chebyshev_u(3), "-1.037037037037037037037036")
        assert is_close_arb(x.jacobi_p(3, 0.25, 0.5), "-0.4722222222222222222222269")
        assert is_close_arb(x.gegenbauer_c(3, 0.25), "-0.1736111111111111111111118")
        assert is_close_arb(x.laguerre_l(3, 0.25), "0.4790702160493827160493819")
        assert is_close_arb(x.hermite_h(3), "-3.703703703703703703703700")
        assert is_close_arb(x.legendre_p(3), "-0.4074074074074074074074072")
        assert is_close_arb(x.legendre_q(3), "0.2476922409970481777113051")
        assert raises(lambda: x.legendre_p(3, type=1), ValueError)
        assert raises(lambda: x.legendre_q(3, type=4), ValueError)

        root = arb.legendre_p_root(5, 2)
        assert is_close_arb(root, 0.0)
        root_w, w = arb.legendre_p_root(5, 2, weight=True)
        assert is_close_arb(root_w, 0.0)
        assert is_close_arb(w, "0.5688888888888888888888887")
        assert raises(lambda: arb.legendre_p_root(5, 5), ValueError)

        x = arb(2)
        assert is_close_arb(arb(0.5).erfcinv(), "0.4769362762044698733814183")
        assert is_close_arb(x.erfi(), "18.56480241457555259870427")
        assert is_close_arb(x.fresnel_s(), "0.3434156783636982421953005")
        assert is_close_arb(x.fresnel_s(normalized=False), "0.8047764893437561102962750")
        assert is_close_arb(x.fresnel_c(), "0.4882534060753407545002233")
        assert is_close_arb(x.fresnel_c(normalized=False), "0.4614614624332163728664738")
        assert is_close_arb(x.ei(), "4.954234356001890163379494")
        assert is_close_arb(x.si(), "1.605412976802694848576720")
        assert is_close_arb(x.ci(), "0.4229808287748649956985648")
        assert is_close_arb(x.shi(), "2.501567433354975641473370")
        assert is_close_arb(x.chi(), "2.452666922646914521906127")
        assert is_close_arb(arb(10).li(), "6.165599504787297937522937")
        assert is_close_arb(arb(10).li(offset=True), "5.120435724669805152678346")

        assert is_close_arb(x.bessel_j(1), "0.5767248077568733872024483")
        assert is_close_arb(x.bessel_y(1), "-0.1070324315409375468883694")
        assert is_close_arb(x.bessel_k(1), "0.1398658818165224272845938")
        assert is_close_arb(x.bessel_k(1, scaled=True), "1.033476847068688573175318")
        assert is_close_arb(x.bessel_i(1), "1.590636854637329063382253")
        assert is_close_arb(x.bessel_i(1, scaled=True), "0.2152692892489376591585047")

        assert is_close_arb(x.gamma_upper(2.5), "0.7303608140431147358169869")
        assert is_close_arb(x.gamma_upper(2.5, regularized=1), "0.5494159513527802326058330")
        assert is_close_arb(x.gamma_lower(2.5), "0.5989795741360222846566369")
        assert is_close_arb(x.gamma_lower(2.5, regularized=1), "0.4505840486472197673941662")
        assert is_close_arb(x.gamma_lower(2.5, regularized=2), "0.07965275907323457828282656")
        assert is_close_arb(arb("0.9").beta_lower(1, 2.5), "0.3987350889359326482671957")
        assert is_close_arb(arb("0.9").beta_lower(1, 2.5, regularized=True), "0.9968377223398316206679875")
        assert is_close_arb(x.expint(2), "0.03753426182049045275952322")

        z = arb(0.2)
        assert is_close_arb(z.hypgeom([1, 2], [3]), "1.157177565710487798620111")
        assert is_close_arb(z.hypgeom([1, 2], [3], regularized=True), "0.5785887828552438993100557")
        assert is_close_arb(z.hypgeom_u(1, 2.5), "11.37885908205001566617169")
        assert is_close_arb(z.hypgeom_1f1(1, 2.5), "1.084782224878003156149143")
        assert is_close_arb(z.hypgeom_1f1(1, 2.5, regularized=True), "0.8160304422585721469158928")
        assert is_close_arb(z.hypgeom_0f1(2.5), "1.082319886462734398651672")
        assert is_close_arb(z.hypgeom_0f1(2.5, regularized=True), "0.8141781413451533169519156")
        assert is_close_arb(z.hypgeom_2f1(1, 2, 3), "1.157177565710487798620111")
        assert is_close_arb(z.hypgeom_2f1(1, 2, 3, regularized=True), "0.5785887828552438993100557")
        assert is_close_arb(z.hypgeom_2f1(1, 2, 3, ab=True, ac=True, bc=True, abc=True), "1.157177565710487798620111")

        assert is_close_arb(arb.const_sqrt_pi(), math.sqrt(math.pi))
        assert is_close_arb(arb.const_euler(), 0.5772156649015329, tol=1e-9, rel_tol=1e-9, max_width=1e-9)
        assert is_close_arb(arb.const_catalan(), 0.915965594177219, tol=1e-9, rel_tol=1e-9, max_width=1e-9)
        assert is_close_arb(arb.const_khinchin(), 2.6854520010653062, tol=1e-9, rel_tol=1e-9, max_width=1e-9)
        assert is_close_arb(arb.const_glaisher(), 1.2824271291006226, tol=1e-9, rel_tol=1e-9, max_width=1e-9)

        assert is_close_arb(arb(27).root(3), 3.0)
        assert is_close_arb(arb.gram_point(0), "17.84559954041086081682633")
        assert is_close_arb(arb(100).zeta_nzeros(), 29.0)
        assert is_close_arb(arb(123).backlund_s(), "0.4757920863536796196115762")
        F, G = arb(1).coulomb(0.5, 0.25)
        assert is_close_arb(F, "0.4283180781043541845555929")
        assert is_close_arb(G, "1.218454487206367973745601")
        assert is_close_arb(arb(1).coulomb_f(0.5, 0.25), "0.4283180781043541845555929")
        assert is_close_arb(arb(1).coulomb_g(0.5, 0.25), "1.218454487206367973745601")

def test_arb_internal_branches() -> None:
    import flint
    assert is_close_arb(arb(flint.arf(1.25)), 1.25)
    assert is_close_arb(arb(flint.fmpz(7)), 7.0)
    assert is_close_arb(arb(flint.fmpq(3, 4)), 0.75)
    assert is_close_arb(arb(1, (1, -2)).rad(), 0.25)
    assert (arb(0).is_zero()) is True

    class DummyMpf:
        def __init__(self, mpf: tuple[int, int, int, int]) -> None:
            self._mpf_ = mpf

    assert is_close_arb(arb(DummyMpf((0, 0, 0, 0))), 0.0)
    assert (arb(DummyMpf((0, 0, 1, 0))).is_nan()) is True
    assert is_close_arb(arb(DummyMpf((0, 3, -1, 2))), 1.5)
    assert is_close_arb(arb(DummyMpf((1, 3, -1, 2))), -1.5)

    assert raises(lambda: arb(1).contains(object()), TypeError)  # type: ignore[arg-type]
    p = flint.fmpz_poly([1, 2])  # 1 + 2*x
    assert is_close_arb(p(arb(3)), 7.0)
    assert is_close_arb(p(1.5), 4.0)
    assert raises(lambda: p(object()), TypeError)  # type: ignore

    assert arb(2).repr() == "arb((0x1, 0x1))"

    x = arb(2)
    bad = object()
    assert raises(lambda: x - bad, TypeError)  # type: ignore[operator]
    assert raises(lambda: bad - x, TypeError)  # type: ignore[operator]
    assert raises(lambda: x * bad, TypeError)  # type: ignore[operator]
    assert raises(lambda: bad * x, TypeError)  # type: ignore[operator]
    assert raises(lambda: x / bad, TypeError)  # type: ignore[operator]
    assert raises(lambda: bad / x, TypeError)  # type: ignore[operator]
    assert raises(lambda: x ** bad, TypeError)  # type: ignore[operator]
    assert raises(lambda: bad ** x, TypeError)  # type: ignore[operator]
    assert raises(lambda: pow(x, 2, 3), TypeError)  # type: ignore[misc]
    assert raises(lambda: pow(2, x, 3), TypeError)  # type: ignore[misc]
    assert x.neg().is_finite() is True

    q = flint.fmpq(3, 4)
    assert is_close_arb(arb.sin_pi_fmpq(q), math.sin(math.pi * 0.75))
    assert is_close_arb(arb.cos_pi_fmpq(q), math.cos(math.pi * 0.75))
    s, c = arb.sin_cos_pi_fmpq(q)
    assert is_close_arb(s, math.sin(math.pi * 0.75))
    assert is_close_arb(c, math.cos(math.pi * 0.75))

    assert isinstance(int(arb(1).rel_accuracy_bits()), int)
    assert isinstance(int(arb(1).rel_one_accuracy_bits()), int)
    assert isinstance(arb(1).mid_rad_10exp()[0], flint.fmpz)

# Tests for arithmetic functions in `flint.arb`.

# NOTE: Correctness of an arb function `F` is specified as follows:

#     If `f` is the corresponding real-valued arithmetic function, `F` is correct
#     only if, for any Arb X and any real number x in the interval X,
#     `f(x)` is in `F(X)`.

# These tests assume arb.__contains__ is correct.

def test_arb_sub() -> None:
    """`arb.__sub__` works as expected."""
    arb1 = arb(2, 0.5)
    arb2 = arb(1, 1)
    with ctx.workprec(100):
        actual = arb1 - arb2
    # Smallest value in diff => 1.5 - 2 = -0.5
    # Largest value in diff => 2.5 - 0 = 2.5
    true_interval = arb(1, 1.5)  # [-0.5, 2.5]
    assert true_interval in actual

def test_arb_add() -> None:
    """`arb.__add__` works as expected."""
    arb1 = arb(2, 1)
    arb2 = arb(1, 1)
    with ctx.workprec(100):
        actual = arb1 + arb2
    true_interval = arb(3, 2)  # [1, 5]
    assert true_interval in actual

def test_arb_mul() -> None:
    """`arb.__mul__` works as expected."""
    arb1 = arb(2, 1)
    arb2 = arb(1, 1)
    with ctx.workprec(100):
        actual = arb1 * arb2
    true_interval = arb(3, 3)  # [0, 6]
    assert true_interval in actual

def test_arb_div() -> None:
    """`arb.__div__` works as expected."""
    arb1 = arb(4, 1)
    arb2 = arb(2, 1)
    with ctx.workprec(100):
        actual = arb1 / arb2
    true_interval = arb(4, 1)  # [3, 5]
    assert true_interval in actual

def test_arb_log() -> None:
    """`arb.log` works as expected."""
    midpoint = (1 + math.exp(10)) / 2
    arb_val = arb(midpoint, midpoint - 1)  # [1, exp(10)]
    with ctx.workprec(100):
        actual = arb_val.log()
    true_interval = arb(5, 5)  # [0,10]
    assert true_interval in actual

def test_arb_exp() -> None:
    """`arb.exp` works as expected."""
    midpoint = math.log(9) / 2
    arb_val = arb(midpoint, midpoint)  # [0, log(9)]
    with ctx.workprec(100):
        actual = arb_val.exp()
    true_interval = arb(5, 4)  # [1,9]
    assert true_interval in actual

def test_arb_max() -> None:
    """`arb.max` works as expected."""
    arb1 = arb(1.5, 0.5)  # [1, 2]
    arb2 = arb(1, 2)  # [-1, 3]
    with ctx.workprec(100):
        actual = arb1.max(arb2)
    true_interval = arb(2, 1)  # [1, 3]
    assert true_interval in actual

def test_arb_min() -> None:
    """`arb.min` works as expected."""
    arb1 = arb(1.5, 0.5)  # [1, 2]
    arb2 = arb(1, 2)  # [-1, 3]
    with ctx.workprec(100):
        actual = arb1.min(arb2)
    true_interval = arb(0.5, 1.5)  # [-1, 2]
    assert true_interval in actual

def test_arb_abs() -> None:
    """`arb.__abs__` works as expected."""
    arb_val = arb(1, 2)  # [-1,3]
    actual = abs(arb_val)
    true_interval = arb(1.5, 1.5)
    assert true_interval in actual

def test_arb_neg() -> None:
    """`arb.neg` works as expected."""
    arb_val = arb(1, 2)  # [-1,3]
    actual = arb_val.neg(exact=True)
    true_interval = arb(-2, 1)  # [-3,1]
    assert true_interval in actual

def test_arb_neg_dunder() -> None:
    """`arb.__neg__` works as expected."""
    arb_val = arb(1, 2)  # [-1,3]
    actual = -arb_val
    true_interval = arb(-2, 1)  # [-3,1]
    assert true_interval in actual

def test_arb_sgn() -> None:
    """`arb.sgn` works as expected."""
    arb1 = arb(1, 0.5)  # [0.5,1.5]
    arb2 = arb(-1, 0.5)  # [-1.5,-0.5]
    arb3 = arb(1, 2)  # [-1,3]
    assert_almost_equal(float(arb1.sgn()), 1)
    assert_almost_equal(float(arb2.sgn()), -1)
    # arb3 contains both positive and negative numbers
    # So, arb_sgn returns [0, 1]
    assert_almost_equal(float(arb3.sgn().mid()), 0)
    assert_almost_equal(float(arb3.sgn().rad()), 1)

def test_arb_erfinv() -> None:
    """`arb.erfinv` works as expected."""
    midpoint = (math.erf(1 / 8) + math.erf(1 / 16)) / 2
    radius = midpoint - math.erf(1 / 16)
    arb_val = arb(midpoint, radius)
    with ctx.workprec(100):
        actual = arb_val.erfinv()
    true_interval = arb(3 / 32, 1 / 32)  # [1/16, 1/8]
    assert true_interval in actual

def test_arb_erf() -> None:
    """`arb.erf` works as expected."""
    arb_val = arb(2, 1)
    with ctx.workprec(100):
        actual = arb_val.erf()
    true_interval = arb(
        (math.erf(1) + math.erf(3)) / 2,
        (math.erf(1) + math.erf(3)) / 2 - math.erf(1)
    )
    assert true_interval in actual

def test_arb_erfc() -> None:
    """`arb.erfc` works as expected."""
    arb_val = arb(2, 1)
    with ctx.workprec(100):
        actual = arb_val.erfc()
    true_interval = arb(
        (math.erfc(1) + math.erfc(3)) / 2,
        (math.erfc(1) + math.erfc(3)) / 2 - math.erfc(3)
    )
    assert true_interval in actual

def test_arb_const_pi() -> None:
    """`arb.pi` works as expected."""
    with ctx.workprec(100):
        actual = arb.pi()
    interval_around_pi = arb(math.pi, 1e-10)
    assert actual in interval_around_pi

def test_arb_union() -> None:
    """`arb.union` works as expected."""
    arb1 = arb(1, 0.5)  # [0.5,1.5]
    arb2 = arb(3, 0.5)  # [2.5,3.5]
    with ctx.workprec(100):
        actual = arb1.union(arb2)
    true_interval = arb(2, 1.5)  # [0.5, 3.5]
    assert true_interval in actual

def test_arb_sum() -> None:
    """`arb.__sum__` works as expected."""
    arb1 = arb(1, 0.5)  # [0.5,1.5]
    arb2 = arb(2, 0.5)  # [1.5,2.5]
    arb3 = arb(3, 0.5)  # [2.5,3.5]
    with ctx.workprec(100):
        actual = arb1 + arb2 + arb3
    true_interval = arb(6, 1.5)  # [4.5, 7.5]
    assert true_interval in actual

def test_no_tests_missing() -> None:
    """Make sure all arb tests are included in all_tests."""
    test_funcs = {f for name, f in globals().items() if name.startswith("test_")}
    untested = test_funcs - set(all_tests)
    assert not untested, f"Untested functions: {untested}"

all_tests = [
    test_no_tests_missing,
    test_from_int,
    test_from_float,
    test_from_float_inf,
    test_from_man_exp,
    test_from_midpoint_radius,
    test_is_exact,
    test_is_finite,
    test_is_nan,
    test_lower,
    test_upper,
    test_contains,
    test_hash,
    test_arb_constructor,
    test_arb_contains,
    test_arb_comparison,
    test_arb_comparison_scalars,
    test_arb_comparison_invalid_type,
    test_arb_pos_floor_ceil,
    test_arb_abs_bounds,
    test_arb_exact_conversions_and_string_options,
    test_arb_mpf_mid_rad_and_interval_ops,
    test_arb_arithmetic,
    test_arb_functions,
    test_arb_constants_and_special_values,
    test_arb_unique_fmpz_and_bits,
    test_arb_atan2_log_base_and_lambertw,
    test_arb_special_functions,
    test_arb_internal_branches,
    test_arb_sub,
    test_arb_add,
    test_arb_mul,
    test_arb_div,
    test_arb_log,
    test_arb_exp,
    test_arb_max,
    test_arb_min,
    test_arb_abs,
    test_arb_neg,
    test_arb_neg_dunder,
    test_arb_sgn,
    test_arb_erfinv,
    test_arb_erf,
    test_arb_erfc,
    test_arb_const_pi,
    test_arb_union,
    test_arb_sum,
]
