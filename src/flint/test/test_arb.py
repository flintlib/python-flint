"""Test for python-flint's `arb` type."""

import math

from flint import arb, ctx

def assert_almost_equal(x, y, places=7):
    """Helper method for approximate comparisons."""
    assert round(x-y, ndigits=places) == 0

def test_from_int():
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

def test_from_float():
    """Tests instantiating `arb`s from floats."""
    for val in [0.0, 1.1, -1.1, 9.9 * 0.123, 99.12]:
        x = arb(val)
        man, exp = x.man_exp()
        assert (int(man) * 2 ** int(exp)) == val

def test_from_float_inf():
    """Tests `arb` works with +/- inf."""
    posinf = arb(float("inf"))
    neginf = arb(float("-inf"))

    assert not posinf.is_finite()
    assert not neginf.is_finite()
    assert float(posinf) == float("inf")
    assert float(neginf) == float("-inf")

def test_from_man_exp():
    """Tests instantiating `arb`s with mantissa and exponent."""
    for man, exp in [(2, 30), (4, 300), (5 * 10**2, 7**8)]:
        x = arb(mid=(man, exp))
        m, e = x.man_exp()
        assert (m * 2**e) == (man * 2**exp)

def test_from_midpoint_radius():
    """Tests instantiating `arb`s with midpoint and radius."""
    for mid, rad in [(10, 1), (10000, 5), (10, 1), (10, 1)]:
        mid_arb = arb(mid)
        rad_arb = arb(rad)
        x = arb(mid_arb, rad_arb)
        assert x.mid() == mid_arb
        actual_radius = float(x.rad())
        assert_almost_equal(actual_radius, rad)

def test_is_exact():
    """Tests `arb.is_exact`."""
    for arb_val, exact in [
        (arb(10), True),
        (arb(0.01), True),
        (arb(-float("inf")), True),
        (arb(1, 0), True),
        (arb(1, 1), False),
    ]:
        assert arb_val.is_exact() == exact

def test_is_finite():
    """Tests `arb.is_finite`."""
    assert not (arb(-float("inf")).is_finite())
    assert not (arb(float("inf")).is_finite())
    assert (arb(10).is_finite())

def test_is_nan():
    """Tests `arb.is_nan`."""
    assert (arb(float("nan")).is_nan())
    assert not (arb(0.0).is_nan())

def test_lower():
    """Tests `arb.lower`."""
    with ctx.workprec(100):
        arb_val = arb(1, 0.5)
        assert_almost_equal(float(arb_val.lower()), 0.5)

def test_upper():
    """Tests `arb.upper`."""
    with ctx.workprec(100):
        arb_val = arb(1, 0.5)
        assert_almost_equal(float(arb_val.upper()), 1.5)

def test_contains():
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

def test_hash():
    """`x` and `y` hash to the same value if they have the same midpoint and radius.

    Args:
        x: An arb.
        y: An arb.
        expected: Whether `x` and `y` should hash to the same value.
    """
    def arb_pi(prec):
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



# Tests for arithmetic functions in `flint.arb`.

# NOTE: Correctness of an arb function `F` is specified as follows:

#     If `f` is the corresponding real-valued arithmetic function, `F` is correct
#     only if, for any Arb X and any real number x in the interval X,
#     `f(x)` is in `F(X)`.

# These tests assume arb.__contains__ is correct.

def test_arb_sub():
    """`arb.__sub__` works as expected."""
    arb1 = arb(2, 0.5)
    arb2 = arb(1, 1)
    with ctx.workprec(100):
        actual = arb1 - arb2
    # Smallest value in diff => 1.5 - 2 = -0.5
    # Largest value in diff => 2.5 - 0 = 2.5
    true_interval = arb(1, 1.5)  # [-0.5, 2.5]
    assert true_interval in actual

def test_arb_add():
    """`arb.__add__` works as expected."""
    arb1 = arb(2, 1)
    arb2 = arb(1, 1)
    with ctx.workprec(100):
        actual = arb1 + arb2
    true_interval = arb(3, 2)  # [1, 5]
    assert true_interval in actual

def test_arb_mul():
    """`arb.__mul__` works as expected."""
    arb1 = arb(2, 1)
    arb2 = arb(1, 1)
    with ctx.workprec(100):
        actual = arb1 * arb2
    true_interval = arb(3, 3)  # [0, 6]
    assert true_interval in actual

def test_arb_div():
    """`arb.__div__` works as expected."""
    arb1 = arb(4, 1)
    arb2 = arb(2, 1)
    with ctx.workprec(100):
        actual = arb1 / arb2
    true_interval = arb(4, 1)  # [3, 5]
    assert true_interval in actual

def test_arb_log():
    """`arb.log` works as expected."""
    midpoint = (1 + math.exp(10)) / 2
    arb_val = arb(midpoint, midpoint - 1)  # [1, exp(10)]
    with ctx.workprec(100):
        actual = arb_val.log()
    true_interval = arb(5, 5)  # [0,10]
    assert true_interval in actual

def test_arb_exp():
    """`arb.exp` works as expected."""
    midpoint = math.log(9) / 2
    arb_val = arb(midpoint, midpoint)  # [0, log(9)]
    with ctx.workprec(100):
        actual = arb_val.exp()
    true_interval = arb(5, 4)  # [1,9]
    assert true_interval in actual

def test_arb_max():
    """`arb.max` works as expected."""
    arb1 = arb(1.5, 0.5)  # [1, 2]
    arb2 = arb(1, 2)  # [-1, 3]
    with ctx.workprec(100):
        actual = arb1.max(arb2)
    true_interval = arb(2, 1)  # [1, 3]
    assert true_interval in actual

def test_arb_min():
    """`arb.min` works as expected."""
    arb1 = arb(1.5, 0.5)  # [1, 2]
    arb2 = arb(1, 2)  # [-1, 3]
    with ctx.workprec(100):
        actual = arb1.min(arb2)
    true_interval = arb(0.5, 1.5)  # [-1, 2]
    assert true_interval in actual

def test_arb_abs():
    """`arb.__abs__` works as expected."""
    arb_val = arb(1, 2)  # [-1,3]
    actual = abs(arb_val)
    true_interval = arb(1.5, 1.5)
    assert true_interval in actual

def test_arb_neg():
    """`arb.neg` works as expected."""
    arb_val = arb(1, 2)  # [-1,3]
    actual = arb_val.neg(exact=True)
    true_interval = arb(-2, 1)  # [-3,1]
    assert true_interval in actual

def test_arb_neg_dunder():
    """`arb.__neg__` works as expected."""
    arb_val = arb(1, 2)  # [-1,3]
    actual = -arb_val
    true_interval = arb(-2, 1)  # [-3,1]
    assert true_interval in actual

def test_arb_sgn():
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

def test_arb_erfinv():
    """`arb.erfinv` works as expected."""
    midpoint = (math.erf(1 / 8) + math.erf(1 / 16)) / 2
    radius = midpoint - math.erf(1 / 16)
    arb_val = arb(midpoint, radius)
    with ctx.workprec(100):
        actual = arb_val.erfinv()
    true_interval = arb(3 / 32, 1 / 32)  # [1/16, 1/8]
    assert true_interval in actual

def test_arb_erf():
    """`arb.erf` works as expected."""
    arb_val = arb(2, 1)
    with ctx.workprec(100):
        actual = arb_val.erf()
    true_interval = arb(
        (math.erf(1) + math.erf(3)) / 2,
        (math.erf(1) + math.erf(3)) / 2 - math.erf(1)
    )
    assert true_interval in actual

def test_arb_erfc():
    """`arb.erfc` works as expected."""
    arb_val = arb(2, 1)
    with ctx.workprec(100):
        actual = arb_val.erfc()
    true_interval = arb(
        (math.erfc(1) + math.erfc(3)) / 2,
        (math.erfc(1) + math.erfc(3)) / 2 - math.erfc(3)
    )
    assert true_interval in actual

def test_arb_const_pi():
    """`arb.pi` works as expected."""
    with ctx.workprec(100):
        actual = arb.pi()
    interval_around_pi = arb(math.pi, 1e-10)
    assert actual in interval_around_pi

def test_arb_union():
    """`arb.union` works as expected."""
    arb1 = arb(1, 0.5)  # [0.5,1.5]
    arb2 = arb(3, 0.5)  # [2.5,3.5]
    with ctx.workprec(100):
        actual = arb1.union(arb2)
    true_interval = arb(2, 1.5)  # [0.5, 3.5]
    assert true_interval in actual

def test_arb_sum():
    """`arb.__sum__` works as expected."""
    arb1 = arb(1, 0.5)  # [0.5,1.5]
    arb2 = arb(2, 0.5)  # [1.5,2.5]
    arb3 = arb(3, 0.5)  # [2.5,3.5]
    with ctx.workprec(100):
        actual = arb1 + arb2 + arb3
    true_interval = arb(6, 1.5)  # [4.5, 7.5]
    assert true_interval in actual

def test_no_tests_missing():
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