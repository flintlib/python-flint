import sys
import math
import operator
import pickle
import doctest
import platform

from flint.utils.flint_exceptions import DomainError

import flint

if sys.version_info[0] >= 3:
    long = int

PYPY = platform.python_implementation() == "PyPy"

ctx = flint.ctx

def raises(f, exception):
    try:
        f()
    except exception:
        return True
    return False


_default_ctx_string = """\
pretty = True      # pretty-print repr() output
unicode = False    # use unicode characters in output
prec = 53          # real/complex precision (in bits)
dps = 15           # real/complex precision (in digits)
cap = 10           # power series precision
threads = 1        # max number of threads used internally
"""

def test_pyflint():

    assert flint.__version__ == "0.5.0"

    ctx = flint.ctx
    assert str(ctx) == repr(ctx) == _default_ctx_string
    assert ctx.prec == 53
    oldprec = ctx.prec
    try:
        f1 = flint.arb(2).sqrt()
        ctx.prec = 10
        f2 = flint.arb(2).sqrt()
        assert f1.rad() < f2.rad()
        assert 1e-17 < f1.rad() < 1e-15
        assert 1e-3 < f2.rad() < 1e-2
    finally:
        ctx.prec = oldprec

    assert ctx.dps == 15
    olddps = ctx.dps
    try:
        f1 = flint.arb(2).sqrt()
        ctx.dps = 3
        f2 = flint.arb(2).sqrt()
        assert f1.rad() < f2.rad()
        assert 1e-17 < f1.rad() < 1e-15
        assert 1e-4 < f2.rad() < 1e-3
    finally:
        ctx.prec = oldprec

    assert ctx.cap == 10
    oldcap = ctx.cap
    try:
        f1 = flint.fmpz_series([1,1])
        ctx.cap = 5
        f2 = flint.fmpz_series([1,1])
        assert f1._equal_repr(flint.fmpz_series([1,1],10))
        assert f2._equal_repr(flint.fmpz_series([1,1],5))
    finally:
        ctx.cap = oldcap

    assert raises(lambda: setattr(ctx, "cap", -1), ValueError)
    assert raises(lambda: setattr(ctx, "prec", -1), ValueError)
    assert raises(lambda: setattr(ctx, "dps", -1), ValueError)

def test_fmpz():
    assert flint.fmpz() == flint.fmpz(0)
    L = [0, 1, 2, 3, 2**31-1, 2**31, 2**63-1, 2**63, 2**64-1, 2**64]
    L += [-x for x in L]
    for i in L:
        f = flint.fmpz(i)
        assert int(f) == i
        assert flint.fmpz(f) == f
        assert flint.fmpz(str(i)) == f
    assert raises(lambda: flint.fmpz("qwe"), ValueError)
    assert raises(lambda: flint.fmpz([]), TypeError)
    for s in L:
        for t in L:
            for ltype in (flint.fmpz, int, long):
                for rtype in (flint.fmpz, int, long):

                    assert (ltype(s) == rtype(t)) == (s == t)
                    assert (ltype(s) != rtype(t)) == (s != t)
                    assert (ltype(s) < rtype(t)) == (s < t)
                    assert (ltype(s) <= rtype(t)) == (s <= t)
                    assert (ltype(s) > rtype(t)) == (s > t)
                    assert (ltype(s) >= rtype(t)) == (s >= t)

                    assert ltype(s) + rtype(t) == s + t
                    assert ltype(s) - rtype(t) == s - t
                    assert ltype(s) * rtype(t) == s * t
                    assert ltype(s) & rtype(t) == s & t
                    assert ltype(s) | rtype(t) == s | t
                    assert ltype(s) ^ rtype(t) == s ^ t
                    assert ~ltype(s) == ~s

                    if t == 0:
                        assert raises(lambda: ltype(s) // rtype(t), ZeroDivisionError)
                        assert raises(lambda: ltype(s) % rtype(t), ZeroDivisionError)
                        assert raises(lambda: divmod(ltype(s), rtype(t)), ZeroDivisionError)
                    else:
                        assert ltype(s) // rtype(t) == s // t
                        assert ltype(s) % rtype(t) == s % t
                        assert divmod(ltype(s), rtype(t)) == divmod(s, t)

                    if 0 <= t < 10:
                        assert (ltype(s) ** rtype(t)) == (s ** t)
                        assert ltype(s) << rtype(t) == s << t
                        assert ltype(s) >> rtype(t) == s >> t
                    elif -10 <= t < 0:
                        assert raises(lambda: ltype(s) << rtype(t), ValueError)
                        assert raises(lambda: ltype(s) >> rtype(t), ValueError)

    assert 2 ** flint.fmpz(2) == 4
    assert type(2 ** flint.fmpz(2)) == flint.fmpz
    assert raises(lambda: () ** flint.fmpz(1), TypeError)
    assert raises(lambda: flint.fmpz(1) ** (), TypeError)
    assert raises(lambda: flint.fmpz(1) ** -1, ValueError)

    mega = flint.fmpz(2) ** 8000000
    assert raises(lambda: mega ** mega, OverflowError)

    pow_mod_examples = [
        (2, 2, 3, 1),
        (2, -1, 5, 3),
        (2, 0, 5, 1),
        (2, 5, 1000, 32),
    ]
    for a, b, c, ab_mod_c in pow_mod_examples:
        assert pow(a, b, c) == ab_mod_c
        assert pow(flint.fmpz(a), b, c) == ab_mod_c
        assert pow(flint.fmpz(a), flint.fmpz(b), c) == ab_mod_c
        assert pow(flint.fmpz(a), b, flint.fmpz(c)) == ab_mod_c
        assert pow(flint.fmpz(a), flint.fmpz(b), flint.fmpz(c)) == ab_mod_c

        # 3-arg pow cannot be made to work with fmpz on PyPy
        # https://github.com/flintlib/python-flint/issues/74
        if not PYPY:
            assert pow(a, flint.fmpz(b), c) == ab_mod_c
            assert pow(a, b, flint.fmpz(c)) == ab_mod_c
            assert pow(a, flint.fmpz(b), flint.fmpz(c)) == ab_mod_c

    assert raises(lambda: pow(flint.fmpz(2), 2, 0), ValueError)
    # XXX: Handle negative modulus like int?
    assert raises(lambda: pow(flint.fmpz(2), 2, -1), ValueError)

    f = flint.fmpz(2)
    assert f.numerator == f
    assert type(f.numerator) is flint.fmpz
    assert f.denominator == 1
    assert type(f.denominator) is flint.fmpz

    assert int(f) == 2
    assert type(int(f)) is int
    assert operator.index(f) == 2
    assert type(operator.index(f)) is int
    assert float(f) == 2.0
    assert type(float(f)) is float
    assert round(f) == 2
    assert type(round(f)) is flint.fmpz
    assert round(f, 1) == 2
    assert type(round(f, 1)) is flint.fmpz
    assert round(f, -1) == 0
    assert type(round(f, -1)) is flint.fmpz
    assert math.trunc(f) == 2
    assert type(math.trunc(f)) is flint.fmpz
    assert math.floor(f) == 2
    assert type(math.floor(f)) is flint.fmpz
    assert math.ceil(f) == 2
    assert type(math.ceil(f)) is flint.fmpz

    assert flint.fmpz(2) != []
    assert +flint.fmpz(0) == 0
    assert +flint.fmpz(1) == 1
    assert +flint.fmpz(-1) == -1
    assert -flint.fmpz(0) == 0
    assert -flint.fmpz(1) == -1
    assert -flint.fmpz(-1) == 1
    assert abs(flint.fmpz(0)) == 0
    assert abs(flint.fmpz(1)) == 1
    assert abs(flint.fmpz(-1)) == 1

    assert bool(flint.fmpz(0)) == False
    assert bool(flint.fmpz(1)) == True

    assert flint.fmpz(2).bit_length() == 2
    assert flint.fmpz(-2).bit_length() == 2
    assert flint.fmpz(2).height_bits() == 2
    assert flint.fmpz(-2).height_bits() == 2
    assert flint.fmpz(2).height_bits(signed=True) == 2
    assert flint.fmpz(-2).height_bits(signed=True) == -2

    f1 = flint.fmpz(1)
    f2 = flint.fmpz(2)
    f3 = flint.fmpz(3)
    f8 = flint.fmpz(8)

    assert f2 << 2 == 8
    assert f2 << f2 == 8
    assert 2 << f2 == 8
    assert raises(lambda: f2 << -1, ValueError)
    assert raises(lambda: 2 << -f1, ValueError)

    assert f8 >> 2 == f2
    assert f8 >> f2 == f2
    assert 8 >> f2 == f2
    assert raises(lambda: f2 >> -1, ValueError)
    assert raises(lambda: 2 >> -f1, ValueError)

    assert f2 & 3 == 2
    assert f2 & f3 == 2
    assert 2 & f3 == 2
    assert f2 | 3 == 3
    assert f2 | f3 == 3
    assert 2 | f3 == 3
    assert f2 ^ 3 == 1
    assert f2 ^ f3 == 1
    assert 2 ^ f3 == 1

    assert raises(lambda: f2 << (), TypeError)
    assert raises(lambda: () << f2, TypeError)
    assert raises(lambda: f2 >> (), TypeError)
    assert raises(lambda: () >> f2, TypeError)
    assert raises(lambda: f2 & (), TypeError)
    assert raises(lambda: () & f2, TypeError)
    assert raises(lambda: f2 | (), TypeError)
    assert raises(lambda: () | f2, TypeError)
    assert raises(lambda: f2 ^ (), TypeError)
    assert raises(lambda: () ^ f2, TypeError)

    l = [1, 2, 3]
    l[flint.fmpz(1)] = -2
    assert l == [1, -2, 3]
    d = {flint.fmpz(2): 3}
    d[flint.fmpz(2)] = -1

    assert d == {flint.fmpz(2): -1}
    ctx.pretty = False
    assert repr(flint.fmpz(0)) == "fmpz(0)"
    assert repr(flint.fmpz(-27)) == "fmpz(-27)"
    ctx.pretty = True
    assert repr(flint.fmpz(0)) == "0"
    assert repr(flint.fmpz(-27)) == "-27"
    bigstr = '1' * 100
    big = flint.fmpz(bigstr)
    assert big.str() == bigstr
    assert big.str(condense=10) == '1111111111{...80 digits...}1111111111'

def test_fmpz_factor():
    assert flint.fmpz(6).gcd(flint.fmpz(9)) == 3
    assert flint.fmpz(6).gcd(9) == 3
    assert raises(lambda: flint.fmpz(2).gcd('asd'), TypeError)
    assert flint.fmpz(6).lcm(flint.fmpz(9)) == 18
    assert flint.fmpz(6).lcm(9) == 18
    assert raises(lambda: flint.fmpz(2).lcm('asd'), TypeError)
    assert flint.fmpz(25).factor() == [(5, 2)]
    n = flint.fmpz(10**100 + 1)
    assert n.factor() == [
        (73, 1), (137, 1), (401, 1), (1201, 1), (1601, 1), (1676321, 1), (5964848081, 1),
        (129694419029057750551385771184564274499075700947656757821537291527196801, 1)]
    assert n.factor(trial_limit=100) == [
        (73, 1), (137, 1), (401, 1),
        (2493516234411471571047384039650897753117456334167082044912715710972543643391271845384040149601, 1)]
    assert n.factor_smooth() == [
        (73, 1), (137, 1), (401, 1), (1201, 1), (1601, 1),
        (1296814508839693536173209832765271992846610925502473758289451540212712414540699659186801, 1)]

def test_fmpz_functions():
    T, F, VE, OE = True, False, ValueError, OverflowError
    cases = [
        # (f, [f(-1), f(0), f(1), f(2), ... f(10)]),
        (lambda n: flint.fmpz(n).is_prime(),
            [0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0]),
        (lambda n: flint.fmpz(n).is_probable_prime(),
            [0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0]),
        (lambda n: flint.fmpz(n).is_perfect_power(),
            [T, T, T, F, F, T, F, F, F, T, T, F]),
        (lambda n: flint.fmpz(n).partitions_p(),
            [0, 1, 1, 2, 3, 5, 7, 11, 15, 22, 30, 42]),
        (lambda n: flint.fmpz(n).moebius_mu(),
            [1, 0, 1, -1, -1, 0, -1, 1, -1, 0, 0, 1]),
        (lambda n: flint.fmpz.fac_ui(n),
            [OE, 1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800]),
        (lambda n: flint.fmpz.primorial_ui(n),
            [OE, 1, 1, 2, 6, 6, 30, 30, 210, 210, 210, 210, 2310, 2310]),
        (lambda n: flint.fmpz.fib_ui(n),
            [OE, 0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55]),
        (lambda n: flint.fmpz(n).rising(0),
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]),
        (lambda n: flint.fmpz(n).rising(1),
            [-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]),
        (lambda n: flint.fmpz(n).rising(2),
            [0, 0, 2, 6, 12, 20, 30, 42, 56, 72, 90, 110]),
        (lambda n: flint.fmpz.bin_uiui(n, 0),
            [OE, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]),
        (lambda n: flint.fmpz.bin_uiui(n, 1),
            [OE, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]),
        (lambda n: flint.fmpz.bin_uiui(n, 2),
            [OE, 0, 0, 1, 3, 6, 10, 15, 21, 28, 36, 45]),
        (lambda n: flint.fmpz.bell_number(n),
            [OE, 1, 1, 2, 5, 15, 52, 203, 877, 4140, 21147, 115975]),
        (lambda n: flint.fmpz.euler_number(n),
            [OE, 1, 0, -1, 0, 5, 0, -61, 0, 1385, 0, -50521]),
        (lambda n: flint.fmpz.stirling_s1(n, 1),
            [OE, 0, 1, -1, 2, -6, 24, -120, 720, -5040, 40320, -362880]),
        (lambda n: flint.fmpz.stirling_s2(n, 2),
            [OE, 0, 0, 1, 3, 7, 15, 31, 63, 127, 255, 511]),
        (lambda n: flint.fmpz(n).divisor_sigma(2),
            [1, 0, 1, 5, 10, 21, 26, 50, 50, 85, 91, 130]),
        (lambda n: flint.fmpz(n).euler_phi(),
            [0, 0, 1, 1, 2, 2, 4, 2, 6, 4, 6, 4]),
        (lambda n: flint.fmpz(n).isqrt(),
            [VE, 0, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3]),
        (lambda n: flint.fmpz(n).sqrtrem(),
            [VE, (0, 0), (1, 0), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2), (2, 3), (2, 4), (3, 0), (3, 1)]),
        (lambda n: flint.fmpz(n).sqrtmod(11),
            [VE, 0, 1, VE, 5, 2, 4, VE, VE, VE, 3, VE]),
        (lambda n: flint.fmpz(n).root(3),
            [VE, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2]),
        (lambda n: flint.fmpz(n).jacobi(3),
            [-1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1]),
        (lambda n: flint.fmpz(2).jacobi(n),
            [1, 1, 1, -1, -1, -1, -1, 1, 1, 1, 1, -1]),
    ]
    is_exception = lambda v: isinstance(v, type) and issubclass(v, Exception)

    for func, values in cases:
        for n, val in enumerate(values, -1):
            if is_exception(val):
                assert raises(lambda: func(n), val)
            else:
                assert func(n) == val

    assert raises(lambda: flint.fmpz(1).root(-1), ValueError)
    assert raises(lambda: flint.fmpz(1).jacobi('bad'), TypeError)

def test_fmpz_poly():
    Z = flint.fmpz_poly
    assert Z() == Z([])
    assert Z() == Z([0])
    assert Z() == Z([0,flint.fmpz(0),0])
    assert Z() == Z([0,0,0])
    assert Z() != Z([1])
    assert Z([1]) == Z([1])
    assert Z([1]) == Z([flint.fmpz(1)])
    assert Z(Z([1,2])) == Z([1,2])
    assert raises(lambda: Z([1,2,[]]), TypeError)
    assert raises(lambda: Z({}), TypeError)
    # XXX: This should probably be made to work:
    assert raises(lambda: Z((1,2,3)), TypeError)
    for ztype in [int, long, flint.fmpz]:
        assert Z([1,2,3]) + ztype(5) == Z([6,2,3])
        assert ztype(5) + Z([1,2,3]) == Z([6,2,3])
        assert Z([1,2,3]) - ztype(5) == Z([-4,2,3])
        assert ztype(5) - Z([1,2,3]) == Z([4,-2,-3])
        assert Z([1,2,3]) * ztype(5) == Z([5,10,15])
        assert ztype(5) * Z([1,2,3]) == Z([5,10,15])
        assert Z([11,6,2]) // ztype(5) == Z([2,1])
        assert ztype(5) // Z([-2]) == Z([-3])
        assert ztype(5) // Z([1,2]) == 0
        assert Z([11,6,2]) % ztype(5) == Z([1,1,2])
        assert ztype(5) % Z([-2]) == Z([-1])
        assert ztype(5) % Z([1,2]) == 5
        assert Z([1,2,3]) ** ztype(0) == 1
        assert Z([1,2,3]) ** ztype(1) == Z([1,2,3])
        assert Z([1,2,3]) ** ztype(2) == Z([1,4,10,12,9])
    assert divmod(Z([11,6,2]), Z([1,2])) == (Z([2,1]), Z([9,1]))
    assert raises(lambda: pow(Z([1,2]), 2, 3), NotImplementedError)
    assert +Z([1,2]) == Z([1,2])
    assert -Z([1,2]) == Z([-1,-2])
    assert raises(lambda: Z([1,2]) + [], TypeError)
    assert raises(lambda: Z([1,2]) - [], TypeError)
    assert raises(lambda: Z([1,2]) * [], TypeError)
    assert raises(lambda: Z([1,2]) / [], TypeError)
    assert raises(lambda: Z([1,2]) // [], TypeError)
    assert raises(lambda: Z([1,2]) % [], TypeError)
    assert raises(lambda: divmod(Z([1,2]), []), TypeError)
    assert raises(lambda: [] + Z([1,2]), TypeError)
    assert raises(lambda: [] - Z([1,2]), TypeError)
    assert raises(lambda: [] * Z([1,2]), TypeError)
    assert raises(lambda: [] / Z([1,2]), TypeError)
    assert raises(lambda: [] // Z([1,2]), TypeError)
    assert raises(lambda: [] % Z([1,2]), TypeError)
    assert raises(lambda: divmod([], Z([1,2])), TypeError)
    assert raises(lambda: Z([1,2,3]) ** -1, (OverflowError, ValueError))
    assert raises(lambda: Z([1,2,3]) ** Z([1,2]), TypeError)
    assert raises(lambda: Z([1,2]) // Z([]), ZeroDivisionError)
    assert raises(lambda: Z([]) // Z([]), ZeroDivisionError)
    assert raises(lambda: Z([1,2]) % Z([]), ZeroDivisionError)
    assert raises(lambda: divmod(Z([1,2]), Z([])), ZeroDivisionError)
    assert raises(lambda: Z([1,2]) < Z([1,2]), TypeError)
    assert raises(lambda: Z([1,2]) < [], TypeError)
    assert raises(lambda: [] < Z([1,2]), TypeError)
    assert Z([]).degree() == -1
    assert Z([]).length() == 0
    p = Z([1,2])
    assert p.length() == 2
    assert p.degree() == 1
    assert p[0] == 1
    assert p[1] == 2
    assert p[2] == 0
    assert p[-1] == 0
    assert raises(lambda: p.__setitem__(-1, 1), ValueError)
    p[0] = 3
    assert p[0] == 3
    p[4] = 7
    assert p.degree() == 4
    assert p[4] == 7
    assert p[3] == 0
    p[4] = 0
    assert p.degree() == 1
    assert p.coeffs() == [3,2]
    assert Z([]).coeffs() == []
    assert bool(Z([])) == False
    assert bool(Z([1])) == True
    ctx.pretty = False
    assert repr(Z([1,2])) == "fmpz_poly([1, 2])"
    ctx.pretty = True
    assert str(Z([1,2])) == "2*x + 1"
    assert str(Z([])) == "0"
    assert str(Z([1,0,2])) == "2*x^2 + 1"
    assert str(Z([-1,0,2])) == "2*x^2 + (-1)"
    assert str(Z([-1,0,1])) == "x^2 + (-1)"
    p = Z([3,4,5])
    assert p(2) == 31
    assert p(flint.fmpq(2,3)) == flint.fmpq(71,9)
    assert p(Z([1,-1])) == Z([12,-14,5])
    assert p(flint.fmpq_poly([2,3],5)) == flint.fmpq_poly([27,24,9],5)
    assert p(flint.arb("1.1")).overlaps(flint.arb("13.45"))
    assert p(flint.acb("1.1", "1.2")).overlaps(flint.acb("6.25", "18.00"))
    assert raises(lambda: p(None), TypeError)
    assert Z([1,2,3]).derivative() == Z([2,6])
    assert Z([1,2,-4]).height_bits() == 3
    assert Z([1,2,-4]).height_bits(signed=True) == -3
    assert Z([1,2,1]).sqrt() == Z([1,1])
    assert raises(lambda: Z([1,2,2]).sqrt(), ValueError)
    assert Z([1,0,2,0,3]).deflation() == (Z([1,2,3]), 2)
    assert Z([]).deflation() == (Z([]), 1)
    assert Z([1,1]).deflation() == (Z([1,1]), 1)
    [(r,m)] = Z([1,1]).complex_roots()
    assert m == 1
    assert r.overlaps(-1)
    assert Z([]).complex_roots() == []
    assert Z([1]).complex_roots() == []
    [(r,m)] = Z([1,1]).roots()
    assert m == 1
    assert r == -1
    assert Z([]).roots() == []
    assert Z([1]).roots() == []
    assert Z([1, 2]).roots() == []

def test_fmpz_poly_factor():
    Z = flint.fmpz_poly
    assert Z([1,2]).gcd(Z([3,4])) == 1
    assert Z([1,2,1]).gcd(Z([1,1])) == Z([1,1])
    assert raises(lambda: Z([1,2,1]).gcd([]), TypeError)
    assert Z([1,2,1]).factor() == (1, [(Z([1,1]), 2)])

def test_fmpz_poly_functions():
    Z = flint.fmpz_poly
    assert Z.cyclotomic(4) == Z([1,0,1])
    assert Z.cos_minpoly(10) == Z([-1,-1,1])
    assert Z.chebyshev_u(4) == Z([1,0,-12,0,16])
    assert Z.chebyshev_t(4) == Z([1,0,-8,0,8])
    assert Z.swinnerton_dyer(2) == Z([1,0,-10,0,1])
    assert Z.swinnerton_dyer(2, use_arb=False) == Z([1,0,-10,0,1])
    assert raises(lambda: Z.swinnerton_dyer(21), OverflowError)
    assert Z.hilbert_class_poly(-7) == Z([3375,1])
    assert raises(lambda: Z.hilbert_class_poly(-2), ValueError)
    assert Z([1]).is_cyclotomic() == 0
    assert Z([-1,1]).is_cyclotomic() == 1
    assert Z([1,1]).is_cyclotomic() == 2
    assert Z([1,2]).is_cyclotomic() == 0
    assert Z([1,2,1]).is_cyclotomic() == 0
    assert Z([2,2,1,1]).is_cyclotomic() == 0
    assert Z([2,2,1]).is_cyclotomic() == 0
    assert Z([1,2,2]).is_cyclotomic() == 0
    assert Z([1,1,1]).is_cyclotomic() == 3

def test_fmpz_mat():
    M = flint.fmpz_mat
    a = M(2,3,[1,2,3,4,5,6])
    b = M(2,3,[4,5,6,7,8,9])
    assert a == a
    assert a == M(a)
    assert a != b
    assert a.nrows() == 2
    assert a.ncols() == 3
    assert a.entries() == [1,2,3,4,5,6]
    assert a.table() == [[1,2,3],[4,5,6]]
    assert (a + b).entries() == [5,7,9,11,13,15]
    assert (a - b).entries() == [-3,-3,-3,-3,-3,-3]
    assert a.transpose() == M(3,2,[1,4,2,5,3,6])
    assert raises(lambda: a + 1, TypeError)
    assert raises(lambda: 1 + a, TypeError)
    assert raises(lambda: a - 1, TypeError)
    assert raises(lambda: 1 - a, TypeError)
    # XXX: Maybe there should be a ShapeError or something?
    assert raises(a.det, ValueError)
    assert +a == a
    assert -a == M(2,3,[-1,-2,-3,-4,-5,-6])
    c = M(2,2,[1,2,3,4])
    assert c.det() == -2
    assert raises(lambda: a + c, ValueError)
    assert (a * 3).entries() == [3,6,9,12,15,18]
    assert (3 * a).entries() == [3,6,9,12,15,18]
    assert (a * long(3)).entries() == [3,6,9,12,15,18]
    assert (long(3) * a).entries() == [3,6,9,12,15,18]
    assert (a * flint.fmpz(3)).entries() == [3,6,9,12,15,18]
    assert (flint.fmpz(3) * a).entries() == [3,6,9,12,15,18]
    assert M.randrank(5,7,3,10).rank() == 3
    A = M.randbits(5,3,2)
    B = M.randtest(3,7,3)
    C = M.randtest(7,2,4)
    assert (A.nrows(),A.ncols()) == (5,3)
    assert (B.nrows(),B.ncols()) == (3,7)
    assert (C.nrows(),C.ncols()) == (7,2)
    assert A*(B*C) == (A*B)*C
    assert raises(lambda: A*C, ValueError)
    assert bool(M(2,2,[0,0,0,0])) == False
    assert bool(M(2,2,[0,0,0,1])) == True
    ctx.pretty = False
    assert repr(M(2,2,[1,2,3,4])) == 'fmpz_mat(2, 2, [1, 2, 3, 4])'
    ctx.pretty = True
    assert str(M(2,2,[1,2,3,4])) == '[1, 2]\n[3, 4]'
    assert M(1,2,[3,4]) * flint.fmpq(1,3) == flint.fmpq_mat(1, 2, [1, flint.fmpq(4,3)])
    assert flint.fmpq(1,3) * M(1,2,[3,4]) == flint.fmpq_mat(1, 2, [1, flint.fmpq(4,3)])
    assert M(1,2,[3,4]) / 3 == flint.fmpq_mat(1, 2, [1, flint.fmpq(4,3)])
    assert M(2,2,[1,2,3,4]).inv().det() == flint.fmpq(1) / M(2,2,[1,2,3,4]).det()
    assert M(2,2,[1,2,3,4]).inv().inv() == M(2,2,[1,2,3,4])
    assert raises(lambda: M.randrank(4,3,4,1), ValueError)
    assert raises(lambda: M.randrank(3,4,4,1), ValueError)
    assert M(1,1,[3]) ** 5 == M(1,1,[3**5])
    assert raises(lambda: pow(M([[1]]), 2, 3), NotImplementedError)
    assert raises(lambda: M(1,2) ** 3, ValueError)
    assert raises(lambda: M(1,1) ** M(1,1), TypeError)
    assert raises(lambda: 1 ** M(1,1), TypeError)
    assert raises(lambda: M([1]), TypeError)
    assert raises(lambda: M([[1],[2,3]]), ValueError)
    assert raises(lambda: M(None), TypeError)
    assert raises(lambda: M(2,2,[1,2,3]), ValueError)
    assert raises(lambda: M(2,2,2,2), ValueError)
    assert M([[1,2,3],[4,5,6]]) == M(2,3,[1,2,3,4,5,6])
    assert raises(lambda: M([[1]]) < M([[2]]), TypeError)
    assert (M([[1]]) == 1) is False
    assert (1 == M([[1]])) is False
    assert (M([[1]]) != 1) is True
    assert (1 != M([[1]])) is True
    D = M([[1,2],[3,4]])
    assert (D[0,0],D[0,1],D[1,0],D[1,1]) == (1,2,3,4)
    D[0,0] = 3
    assert D == M([[3,2],[3,4]])
    def set_bad(i,j):
        D[i,j] = -1
    # XXX: Should be IndexError
    raises(lambda: set_bad(2,0), ValueError)
    raises(lambda: set_bad(0,2), ValueError)
    raises(lambda: D[0,2], ValueError)
    raises(lambda: D[0,2], ValueError)
    # XXX: Negative indices?
    raises(lambda: set_bad(-1,0), ValueError)
    raises(lambda: set_bad(0,-1), ValueError)
    raises(lambda: D[-1,0], ValueError)
    raises(lambda: D[0,-1], ValueError)
    assert M.hadamard(2) == M([[1,1],[1,-1]])
    assert raises(lambda: M.hadamard(3), ValueError)
    assert M.hadamard(2).is_hadamard() is True
    assert M([[1,2],[3,4]]).is_hadamard() is False
    M1 = M([[1,0],[1,1]])
    M2 = M([[1,0],[-1,1]])
    x = M([[3],[4]])
    b = M([[3],[7]])
    assert M1.inv() == M2
    assert M2.inv() == M1
    assert M1.inv(integer=True) == M2
    assert M2.inv(integer=True) == M1
    assert M1*x == b
    assert M1.solve(b) == x
    assert M1.solve(b, integer=True) == x
    assert raises(lambda: M1.solve([]), TypeError)
    assert raises(lambda: M1.solve(M([[1]])), ValueError)
    assert raises(lambda: M([[1,1],[1,1]]).solve(b), ZeroDivisionError)
    assert raises(lambda: M([[1,2],[3,4],[5,6]]).solve(b), ValueError)
    assert M([[1,0],[1,2]]).solve(b) == flint.fmpq_mat([[3],[2]])
    assert raises(lambda: M([[1,0],[1,2]]).solve(b, integer=True), ValueError)
    assert raises(lambda: M([[1,2,3],[4,5,6]]).inv(), ValueError)
    assert raises(lambda: M([[1,1],[1,1]]).inv(), ZeroDivisionError)
    assert raises(lambda: M([[1,0],[1,2]]).inv(integer=True), ValueError)
    half = flint.fmpq(1,2)
    assert M([[1,0],[1,2]]).inv() == flint.fmpq_mat([[1, 0], [-half, half]])
    M3 = M([[1,2,3],[4,5,6],[7,8,9]])
    M3_copy = M(M3)
    M3r = M([[-3,0,3],[0,-3,-6],[0,0,0]])
    assert M3.rref() == (M3r, -3, 2)
    assert M3 != M3r
    assert M3.rref(inplace=True) == (M3r, -3, 2)
    assert M3 == M3r
    M3 = M3_copy
    M3n = M([[3,0,0],[-6,0,0],[3,0,0]])
    assert M3.nullspace() == (M3n, 1)
    assert M3 * M3.nullspace()[0] == M(3,3,[0]*9)
    # XXX: lll core dumps on a singular matrix
    M4 = M([[1,2,3],[4,5,6],[7,8,10]])
    L4 = M([[0,0,1],[-1,1,0],[2,1,0]])
    T4 = M([[1,-2,1],[0,5,-3],[-2,1,0]])
    assert L4 == T4 * M4
    assert M4.lll() == L4
    assert M4.lll(transform=True) == (L4, T4)
    # XXX: rep="gram" consumes all memory in the system and core dumps
    #for rep in "zbasis", "gram":
    rep = "zbasis"
    for gram in "approx", "exact":
        assert M4.lll(rep=rep, gram=gram) == L4
        assert M4.lll(rep=rep, gram=gram, transform=True) == (L4, T4)
    assert raises(lambda: M4.lll(rep="bad"), ValueError)
    assert raises(lambda: M4.lll(gram="bad"), ValueError)
    M5 = M([[1,2,3],[4,5,6]])
    H5 = M([[1,2,3],[0,3,6]])
    T5 = M([[1,0],[4,-1]])
    assert H5 == T5 * M5
    assert M5.hnf() == H5
    assert M5.hnf(transform=True) == (H5, T5)
    assert M5.is_hnf() is False
    assert H5.is_hnf() is True
    S5 = M([[1,0,0],[0,3,0]])
    assert M5.snf() == S5
    assert M5.is_snf() is False
    assert S5.is_snf() is True
    M6 = M([[2,0,0],[0,2,1],[0,0,2]])
    assert M6.charpoly() == flint.fmpz_poly([-8,12,-6,1])
    assert M6.minpoly() == flint.fmpz_poly([4,-4,1])
    assert list(M6) == [2,0,0,0,2,1,0,0,2]

def test_fmpz_series():
    Zp = flint.fmpz_poly
    Z = flint.fmpz_series
    ctx = flint.ctx
    assert ctx.cap == 10
    s1 = Z([1,2])
    s2 = Z([1,2])
    s3 = Z([1,1])
    s4 = Z([1,2],11)
    p1 = Zp([1,2])
    assert s1._equal_repr(s1) is True
    assert s1._equal_repr(s2) is True
    assert s1._equal_repr(s3) is False
    assert s1._equal_repr(s4) is False
    assert s1._equal_repr(p1) is False
    assert s1._equal_repr(flint.fmpz(1)) is False
    assert s1._equal_repr(1) is False
    assert Z([])._equal_repr(Z([0])) is True
    assert Z([1])._equal_repr(Z([0])) is False
    # XXX: this gives a core dump:
    # s = Z([1,2])
    # s[10**10] = 1
    assert Z(Z([1]))._equal_repr(Z([1]))
    assert Z(Zp([1]))._equal_repr(Z([1]))
    assert Z(1)._equal_repr(Z([1]))
    assert raises(lambda: Z([1]) < Z([1]), TypeError)
    assert len(Z([1,2])) == 2
    assert Z([1,2]).length() == 2
    s5 = Z([1,2])
    assert s5[1] == 2
    assert s5[2] == 0
    assert s5[-1] == 0
    # XXX: This goes beyond cap. Should it give an error?
    assert s5[100] == 0
    s5[2] = -1
    assert s5[2] == -1
    assert s5._equal_repr(Z([1,2,-1]))
    def set_bad():
        s5[-1] = 3
    assert raises(set_bad, ValueError)
    assert Z([1,2,0,4]).str() == "1 + 2*x + 4*x^3 + O(x^10)"
    assert Z([1,2,0,4]).repr() == "fmpz_series([1, 2, 0, 4], prec=10)"
    assert Z([],0).str() == "O(x^0)"
    assert Z([],-1).str() == "(invalid power series)"
    assert (+Z([1,2]))._equal_repr(Z([1,2]))
    assert (-Z([1,2]))._equal_repr(Z([-1,-2]))
    assert (Z([1,2]) + Z([3,4,5]))._equal_repr(Z([4,6,5]))
    assert (Z([1,2]) + 1)._equal_repr(Z([2,2]))
    assert (Z([1,2]) + Zp([3,4,5]))._equal_repr(Z([4,6,5]))
    assert (1 + Z([1,2]))._equal_repr(Z([2,2]))
    assert (Zp([1,2]) + Z([3,4,5]))._equal_repr(Z([4,6,5]))
    assert raises(lambda: Z([1,2]) + [], TypeError)
    assert raises(lambda: [] + Z([1,2]), TypeError)
    assert (Z([1,2]) - Z([3,5]))._equal_repr(Z([-2,-3]))
    assert (Z([1,2]) - 1)._equal_repr(Z([0,2]))
    assert (1 - Z([1,2]))._equal_repr(Z([0,-2]))
    assert raises(lambda: [] - Z([1,2]), TypeError)
    assert raises(lambda: Z([1,2]) - [], TypeError)
    assert (Z([1,2]) * Z([1,2]))._equal_repr(Z([1,4,4]))
    assert (2 * Z([1,2]))._equal_repr(Z([2,4]))
    assert (Z([1,2]) * 2)._equal_repr(Z([2,4]))
    assert raises(lambda: [] * Z([1,2]), TypeError)
    assert raises(lambda: Z([1,2]) * [], TypeError)
    assert Z([1,2]).valuation() == 0
    assert Z([0,2]).valuation() == 1
    assert Z([0,0]).valuation() == -1
    assert (Z([1,1]) / Z([1,-1]))._equal_repr(Z([1,2,2,2,2,2,2,2,2,2]))
    assert raises(lambda: Z([1,1]) / Z([]), ZeroDivisionError)
    assert (Z([]) / Z([1,-1]))._equal_repr(Z([]))
    # quotient would not be a power series
    assert raises(lambda: Z([1,1]) / Z([0,1]), ValueError)
    # leading term in denominator is not a unit
    assert raises(lambda: Z([1,1]) / Z([2,1]), ValueError)
    assert (Z([0,1,1]) / Z([0,1,-1]))._equal_repr(Z([1,2,2,2,2,2,2,2,2],9))
    # XXX: This gives leading term not a unit:
    # assert Z([2,4]) / 2 == Z([1,2])
    assert (Z([1,1]) / -1)._equal_repr(Z([-1,-1]))
    assert raises(lambda: Z([1,1]) / [], TypeError)
    assert raises(lambda: [] / Z([1,1]), TypeError)
    assert (Z([1,2]) ** 2)._equal_repr(Z([1,4,4]))
    assert raises(lambda: pow(Z([1,2]), 3, 5), NotImplementedError)
    assert (Z([1,2])(Z([0,1,2])))._equal_repr(Z([1,2,4]))
    assert raises(lambda: Z([1,2])(Z([1,2])), ValueError)
    assert raises(lambda: Z([1,2])([]), TypeError)
    coeffs = [0, 1, -2, 8, -40, 224, -1344, 8448, -54912, 366080]
    assert Z([0,1,2]).reversion()._equal_repr(Z(coeffs))
    # power series reversion must have valuation 1
    assert raises(lambda: Z([1,1]).reversion(), ValueError)
    # leading term is not a unit
    assert raises(lambda: Z([0,2,1]).reversion(), ValueError)


def test_fmpq():
    Q = flint.fmpq
    assert Q() == Q(0)
    assert Q(0) != Q(1)
    assert Q(0) == 0
    assert 0 == Q(0)
    assert Q(2) != 1
    assert 1 != Q(2)
    assert Q(1) != ()
    assert Q(1,2) != 1
    assert Q(2,3) == Q(flint.fmpz(2),long(3))
    assert Q(-2,-4) == Q(1,2)
    assert Q("1") == Q(1)
    assert Q("1/2") == Q(1,2)
    assert raises(lambda: Q("1.0"), ValueError)
    assert raises(lambda: Q("1.5"), ValueError)
    assert raises(lambda: Q("1/2/3"), ValueError)
    assert raises(lambda: Q([]), TypeError)
    assert raises(lambda: Q(1, []), TypeError)
    assert raises(lambda: Q([], 1), TypeError)
    assert bool(Q(0)) == False
    assert bool(Q(1)) == True
    assert Q(1,3) + Q(2,3) == 1
    assert Q(1,3) - Q(2,3) == Q(-1,3)
    assert Q(1,3) * Q(2,3) == Q(2,9)
    assert Q(1,3) + 2 == Q(7,3)
    assert 2 + Q(1,3) == Q(7,3)
    assert Q(1,3) - 2 == Q(-5,3)
    assert 2 - Q(1,3) == Q(5,3)
    assert Q(1,3) * 2 == Q(2,3)
    assert 2 * Q(1,3) == Q(2,3)
    assert Q(2,3) / Q(4,5) == Q(5,6)
    assert Q(2,3) / 5 == Q(2,15)
    assert Q(2,3) / flint.fmpz(5) == Q(2,15)
    assert 5 / Q(2,3) == Q(15,2)
    assert flint.fmpz(5) / Q(2,3) == Q(15,2)
    assert Q(2,3)/5 == Q(2,15)
    assert Q(1,2) ** 2 == Q(1,4)
    assert Q(1,2) ** -2 == Q(4)
    assert raises(lambda: Q(0) ** -1, ZeroDivisionError)
    assert raises(lambda: Q(1,2) ** Q(1,2), TypeError)
    assert raises(lambda: Q(1,2) ** [], TypeError)
    assert raises(lambda: [] ** Q(1,2), TypeError)
    # XXX: This should NotImplementedError or something.
    assert raises(lambda: pow(Q(1,2),2,3), AssertionError)

    megaz = flint.fmpz(2) ** 8000000
    megaq = Q(megaz)
    assert raises(lambda: megaq ** megaz, OverflowError)

    assert raises(lambda: Q(1,2) + [], TypeError)
    assert raises(lambda: Q(1,2) - [], TypeError)
    assert raises(lambda: Q(1,2) * [], TypeError)
    assert raises(lambda: Q(1,2) / [], TypeError)
    assert raises(lambda: [] + Q(1,2), TypeError)
    assert raises(lambda: [] - Q(1,2), TypeError)
    assert raises(lambda: [] * Q(1,2), TypeError)
    assert raises(lambda: [] / Q(1,2), TypeError)
    assert (Q(1,2) == 1) is False
    assert (Q(1,2) != 1) is True
    assert (Q(1,2) <  1) is True
    assert (Q(1,2) <= 1) is True
    assert (Q(1,2) >  1) is False
    assert (Q(1,2) >= 1) is False
    assert (Q(1,2) == Q(3,4)) is False
    assert (Q(1,2) != Q(3,4)) is True
    assert (Q(1,2) <  Q(3,4)) is True
    assert (Q(1,2) <= Q(3,4)) is True
    assert (Q(1,2) >  Q(3,4)) is False
    assert (Q(1,2) >= Q(3,4)) is False
    assert (Q(1,2) == Q(1,2)) is True
    assert (Q(1,2) != Q(1,2)) is False
    assert (Q(1,2) <  Q(1,2)) is False
    assert (Q(1,2) <= Q(1,2)) is True
    assert (Q(1,2) >  Q(1,2)) is False
    assert (Q(1,2) >= Q(1,2)) is True
    assert raises(lambda: Q(1,2) > [], TypeError)
    assert raises(lambda: [] < Q(1,2), TypeError)

    ctx.pretty = False
    assert repr(Q(-2,3)) == "fmpq(-2,3)"
    assert repr(Q(3)) == "fmpq(3)"
    ctx.pretty = True
    assert str(Q(-2,3)) == "-2/3"
    assert str(Q(3)) == "3"

    assert Q(2,3).p == Q(2,3).numer() == Q(2,3).numerator == 2
    assert Q(2,3).q == Q(2,3).denom() == Q(2,3).denominator == 3
    assert +Q(5,7) == Q(5,7)
    assert -Q(5,7) == Q(-5,7)
    assert -Q(-5,7) == Q(5,7)
    assert abs(Q(5,7)) == Q(5,7)
    assert abs(-Q(5,7)) == Q(5,7)
    assert raises(lambda: Q(1,0), ZeroDivisionError)
    assert raises(lambda: Q(1,2) / Q(0), ZeroDivisionError)
    assert raises(lambda: Q(1,2) / 0, ZeroDivisionError)

    assert Q(5,3).floor() == flint.fmpz(1)
    assert Q(-5,3).floor() == flint.fmpz(-2)
    assert Q(5,3).ceil() == flint.fmpz(2)
    assert Q(-5,3).ceil() == flint.fmpz(-1)

    assert type(int(Q(5,3))) is int
    assert type(math.floor(Q(5,3))) is flint.fmpz
    assert type(math.ceil(Q(5,3))) is flint.fmpz
    assert type(math.trunc(Q(5,3))) is flint.fmpz
    assert type(round(Q(5,3))) is flint.fmpz
    assert type(round(Q(5,3))) is flint.fmpz
    assert type(round(Q(5,3), 0)) is flint.fmpq
    assert type(round(Q(5,3), 1)) is flint.fmpq

    assert int(Q(5,3)) == 1
    assert math.floor(Q(5,3)) == flint.fmpz(1)
    assert math.ceil(Q(5,3)) == flint.fmpz(2)
    assert math.trunc(Q(5,3)) == flint.fmpz(1)
    assert round(Q(5,3)) == flint.fmpz(2)

    assert int(Q(-5,3)) == flint.fmpz(-1)
    assert math.floor(Q(-5,3)) == flint.fmpz(-2)
    assert math.ceil(Q(-5,3)) == flint.fmpz(-1)
    assert math.trunc(Q(-5,3)) == flint.fmpz(-1)
    assert round(Q(-5,3)) == -2

    assert round(Q(100,3), 2) == Q(3333,100)
    assert round(Q(100,3), 0) == Q(33,1)
    assert round(Q(100,3), -1) == Q(30,1)
    assert round(Q(100,3), -2) == Q(0)

    d = {}
    d[Q(1,2)] = 3
    d[Q(1,2)] = 4
    assert d == {Q(1,2):4}

    assert Q(-5,3).height_bits() == 3
    assert Q(-5,3).height_bits(signed=True) == -3

    cases = [
        (lambda q: q.next(),
            [Q(0), Q(1), Q(-1), Q(1,2), Q(-1,2), Q(2), Q(-2), Q(1,3), Q(-1,3), Q(3)]),
        (lambda q: q.next(signed=False),
            [Q(0), Q(1), Q(1,2), Q(2), Q(1,3), Q(3), Q(2,3), Q(3,2), Q(1,4), Q(4)]),
        (lambda q: q.next(minimal=False),
            [Q(0), Q(1), Q(-1), Q(1,2), Q(-1,2), Q(2), Q(-2), Q(1,3), Q(-1,3), Q(3,2)]),
        (lambda q: q.next(signed=False, minimal=False),
            [Q(0), Q(1), Q(1,2), Q(2), Q(1,3), Q(3,2), Q(2,3), Q(3), Q(1,4), Q(4,3)]),
    ]
    for func, values in cases:
        for val1, val2 in zip(values[:-1], values[1:]):
            assert func(val1) == val2
    raises(lambda: Q(-1).next(signed=False), ValueError)

    OE = OverflowError
    cases = [
        (flint.fmpq.bernoulli,
            [OE, Q(1), Q(-1,2), Q(1,6), Q(0), Q(-1,30)]),
        (lambda n: flint.fmpq.bernoulli(n, cache=True),
            [OE, Q(1), Q(-1,2), Q(1,6), Q(0), Q(-1,30)]),
        (flint.fmpq.harmonic,
            [OE, Q(0), Q(1), Q(3,2), Q(11, 6), Q(25, 12)]),
        (lambda n: flint.fmpq.dedekind_sum(n, 3),
            [-Q(1,18), 0, Q(1,18), -Q(1,18), 0, Q(1,18), -Q(1,18)]),
    ]
    is_exception = lambda v: isinstance(v, type) and issubclass(v, Exception)

    for func, values in cases:
        for n, val in enumerate(values, -1):
            if is_exception(val):
                assert raises(lambda: func(n), val)
            else:
                assert func(n) == val

def test_fmpq_poly():
    Q = flint.fmpq_poly
    Z = flint.fmpz_poly
    assert Q() == Q([]) == Q([0]) == Q([0,0])
    assert Q() != Q([1])
    assert Q([1]) == Q([1])
    assert Q([1], 2) == Q([flint.fmpq(1,2)])
    assert raises(lambda: Q([1,[]]), TypeError)
    assert raises(lambda: Q({}), TypeError)
    assert raises(lambda: Q([1], []), TypeError)
    assert raises(lambda: Q([1], 0), ZeroDivisionError)
    assert bool(Q()) == False
    assert bool(Q([1])) == True
    assert Q(Q([1,2])) == Q([1,2])
    assert Q(Z([1,2])) == Q([1,2])
    assert Q([1,2]) + 3 == Q([4,2])
    assert 3 + Q([1,2]) == Q([4,2])
    assert Q([1,2]) - 3 == Q([-2,2])
    assert 3 - Q([1,2]) == Q([2,-2])
    assert raises(lambda: Q([1]) + [], TypeError)
    assert raises(lambda: [] + Q([1]), TypeError)
    assert raises(lambda: Q([1]) - [], TypeError)
    assert raises(lambda: [] - Q([1]), TypeError)
    assert raises(lambda: Q([1]) * [], TypeError)
    assert raises(lambda: [] * Q([1]), TypeError)
    assert Q([1,2,1]) // Q([1,1]) == Q([1,1])
    assert Q([1,2,1]) % Q([1,1]) == 0
    assert divmod(Q([1,2,1]), Q([1,1])) == (Q([1,1]), 0)
    assert raises(lambda: Q([1,2,1]) // [], TypeError)
    assert raises(lambda: [] // Q([1,2,1]), TypeError)
    assert raises(lambda: Q([1,2,1]) % [], TypeError)
    assert raises(lambda: [] % Q([1,2,1]), TypeError)
    assert raises(lambda: divmod(Q([1,2,1]), []), TypeError)
    assert raises(lambda: divmod([], Q([1,2,1])), TypeError)
    assert raises(lambda: Q([1,2,1]) / 0, ZeroDivisionError)
    assert raises(lambda: Q([1,2,1]) // 0, ZeroDivisionError)
    assert raises(lambda: Q([1,2,1]) % 0, ZeroDivisionError)
    assert raises(lambda: divmod(Q([1,2,1]), 0), ZeroDivisionError)
    assert +Q([1,2]) == Q([1,2])
    assert -Q([1,2]) == Q([-1,-2])
    assert Q([flint.fmpq(1,2),1]) * 2 == Q([1,2])
    assert Q([1,2]) == Z([1,2])
    assert Z([1,2]) == Q([1,2])
    assert Q([1,2]) != Z([3,2])
    assert Z([1,2]) != Q([3,2])
    assert Q([1,2]) != []
    assert raises(lambda: Q([1,2]) < Q([1,2]), TypeError)
    assert Q([1,2,3])*Q([1,2]) == Q([1,4,7,6])
    assert Q([1,2,3])*Z([1,2]) == Q([1,4,7,6])
    assert Q([1,2,3]) * 3 == Q([3,6,9])
    assert 3 * Q([1,2,3]) == Q([3,6,9])
    assert Q([1,2,3]) * flint.fmpq(2,3) == (Q([1,2,3]) * 2) / 3
    assert flint.fmpq(2,3) * Q([1,2,3]) == (Q([1,2,3]) * 2) / 3
    assert raises(lambda: Q([1,2]) / Q([1,2]), TypeError)
    assert Q([1,2,3]) / flint.fmpq(2,3) == Q([1,2,3]) * flint.fmpq(3,2)
    assert Q([1,2,3]) ** 2 == Q([1,2,3]) * Q([1,2,3])
    assert raises(lambda: pow(Q([1,2]), 3, 5), NotImplementedError)
    assert Q([1,2,flint.fmpq(1,2)]).coeffs() == [1,2,flint.fmpq(1,2)]
    assert Q().coeffs() == []
    assert Q().degree() == -1
    assert Q([1]).degree() == 0
    assert Q([1,2]).degree() == 1
    assert Q().length() == 0
    assert Q([1]).length() == 1
    assert Q([1,2]).length() == 2
    assert (Q([1,2,3]) / 5).numer() == (Q([1,2,3]) / 5).p == Z([1,2,3])
    assert (Q([1,2,3]) / 5).denom() == (Q([1,2,3]) / 5).q == 5
    ctx.pretty = False
    assert repr(Q([15,20,10])) == "fmpq_poly([15, 20, 10])"
    assert repr(Q([15,20,10]) / 25) == "fmpq_poly([3, 4, 2], 5)"
    ctx.pretty = True
    assert str(Q([3,4,2],5)) == "2/5*x^2 + 4/5*x + 3/5"
    a = Q([2,2,3],4)
    assert a[2] == flint.fmpq(3,4)
    assert a[-1] == flint.fmpq(0)
    a[2] = 4
    def set_bad():
        a[-1] = 2
    assert raises(set_bad, ValueError)
    assert a == Q([1,1,8],2)
    p = Q([3,4,5],7)
    assert p(2) == flint.fmpq(31,7)
    assert p(flint.fmpq(2,3)) == flint.fmpq(71,63)
    assert p(Z([1,-1])) == Q([12,-14,5],7)
    assert p(flint.fmpq_poly([2,3],5)) == flint.fmpq_poly([27,24,9],7*5)
    assert raises(lambda: p([]), TypeError)
    assert Q([1,2,3]).derivative() == Q([2,6])
    assert Q([1,2,3]).integral() == Q([0,1,1,1])
    assert Q([1,2,1]).gcd(Q([1,1])) == Q([1,1])
    assert Q([1,2,1]).xgcd(Q([1,1])) == (Q([1,1]), 0, 1)
    assert raises(lambda: Q([1,2,1]).gcd([]), TypeError)
    assert raises(lambda: Q([1,2,1]).xgcd([]), TypeError)
    assert Q([1,2,1]).factor() == (1, [(Q([1,1]), 2)])
    assert Q.bernoulli_poly(3) == Q([0,1,-3,2],2)
    assert Q.euler_poly(3) == Q([1,0,-6,4],4)
    assert Q.legendre_p(3) == Q([0,-3,0,5],2)
    assert Q([]).complex_roots() == []
    assert Q([1]).complex_roots() == []
    [(r,m)] = Q([1,1]).complex_roots()
    assert m == 1
    assert r.overlaps(-1)
    assert str(Q([1,2]).roots()) == "[(-1/2, 1)]"
    assert Q([2,1]).roots() == [(-2, 1)]

def test_fmpq_mat():
    Q = flint.fmpq_mat
    Z = flint.fmpz_mat
    assert Q(1,2,[3,4]) == Z(1,2,[3,4])
    assert Q(1,2,[3,4]) != Z(1,2,[5,4])
    assert Q(1,2,[3,4]) != Q(1,2,[5,4])
    assert Q(Q(1,2,[3,4])) == Q(1,2,[3,4])
    assert Q(Z(1,2,[3,4])) == Q(1,2,[3,4])
    assert Q(2,3,[1,2,3,4,5,6]) + Q(2,3,[4,5,6,7,8,9]) == Q(2,3,[5,7,9,11,13,15])
    assert Q(2,3,[1,2,3,4,5,6]) - Q(2,3,[4,5,6,7,8,9]) == Q(2,3,[-3,-3,-3,-3,-3,-3])
    assert Q(2,3,[1,2,3,4,5,6]) * Q(3,2,[4,5,6,7,8,9]) == Q(2,2,[40,46,94,109])
    assert Q(2,3,[1,2,3,4,5,6]) * Z(3,2,[4,5,6,7,8,9]) == Q(2,2,[40,46,94,109])
    assert Z(2,3,[1,2,3,4,5,6]) * Q(3,2,[4,5,6,7,8,9]) == Q(2,2,[40,46,94,109])
    assert -Q(2,1,[2,5]) == Q(2,1,[-2,-5])
    assert +Q(1,1,[3]) == Z(1,1,[3])
    assert (Q(1,1,[1]) == 1) is False
    assert (Q(1,1,[1]) != 1) is True
    assert (1 == Q(1,1,[1])) is False
    assert (1 != Q(1,1,[1])) is True
    assert Q(1,2,[3,4]) * 2 == Q(1,2,[6,8])
    assert Q(1,2,[3,4]) * flint.fmpq(1,3) == Q(1,2,[1,flint.fmpq(4,3)])
    assert Q(1,2,[3,4]) * flint.fmpq(5,3) == Q(1,2,[5,flint.fmpq(20,3)])
    assert 2 * Q(1,2,[3,4]) == Q(1,2,[6,8])
    assert flint.fmpq(1,3) * Q(1,2,[3,4]) == Q(1,2,[1,flint.fmpq(4,3)])
    M = Q([[1,2],[3,4]])
    assert M ** -1 == Q([[-4,2],[3,-1]]) / 2
    assert M ** 0 == Q([[1,0],[0,1]])
    assert M ** 1 == M
    assert M ** 2 == Q([[7,10],[15,22]])
    assert M ** 12 == Q([[138067399, 201223170],[301834755, 439902154]])
    M = Q([[1,2],[2,4]])
    assert raises(lambda: M ** -1, ZeroDivisionError)
    assert Q(1,2,[3,4]) / 2 == Q(1,2,[flint.fmpq(3,2),2])
    assert Q(1,2,[3,4]) / flint.fmpq(2,3) == Q(1,2,[flint.fmpq(9,2),6])
    assert Q(3,2,range(6)).table() == Z(3,2,range(6)).table()
    assert Q(3,2,range(6)).entries() == Z(3,2,range(6)).entries()
    assert Q(3,2,range(6)).nrows() == 3
    assert Q(3,2,range(6)).ncols() == 2
    assert Q(2,2,[3,7,4,5]).det() == -13
    assert (Q(2,2,[3,7,4,5]) / 5).det() == flint.fmpq(-13,25)
    assert raises(lambda: Q(1,2,[1,2]).det(), ValueError)
    assert Q(2,2,[1,2,3,4]).inv().inv() == Q(2,2,[1,2,3,4])
    assert raises(lambda: Q(2,2,[1,1,1,1]).inv(), ZeroDivisionError)
    assert raises(lambda: Q(2,1,[1,1]).inv(), ValueError)
    assert raises(lambda: Q([1]), TypeError)
    assert raises(lambda: Q([[1],[2,3]]), ValueError)
    assert raises(lambda: Q(None), TypeError)
    assert Q([[1,2,3],[4,5,6]]) == Q(2,3,[1,2,3,4,5,6])
    assert raises(lambda: Q(2,3,[1,2,3,4,5]), ValueError)
    assert raises(lambda: Q([[1,2,3],[4,[],6]]), TypeError)
    assert raises(lambda: Q(2,3,[1,2,3,4,[],6]), TypeError)
    assert raises(lambda: Q(2,3,[1,2],[3,4]), ValueError)
    assert bool(Q([[1]])) is True
    assert bool(Q([[0]])) is False
    assert raises(lambda: Q([[1]]) < Q([[0]]), TypeError)
    M = Q([[1,2],[3,4]])
    assert M[0,1] == 2
    M[0,1] = -1
    assert M[0,1] == -1
    # XXX: Negative indices should probably be allowed
    def set_bad(i):
        M[i,0] = -1
    raises(lambda: M[-1,0], ValueError)
    raises(lambda: M[0,-1], ValueError)
    raises(lambda: set_bad(-1), ValueError)
    # XXX: Should be IndexError
    raises(lambda: M[2,0], ValueError)
    raises(lambda: M[0,2], ValueError)
    raises(lambda: set_bad(2), ValueError)
    assert Q([[1,2,3],[4,5,6]]).transpose() == Q([[1,4],[2,5],[3,6]])
    raises(lambda: M + [], TypeError)
    raises(lambda: M - [], TypeError)
    raises(lambda: M * [], TypeError)
    raises(lambda: M / [], TypeError)
    raises(lambda: [] / M, TypeError)
    raises(lambda: [] + M, TypeError)
    raises(lambda: [] - M, TypeError)
    raises(lambda: [] * M, TypeError)
    # XXX: Maybe a ShapeError?
    raises(lambda: Q(1,2,[3,4]) + Q(1,3,[5,6,7]), ValueError)
    raises(lambda: Q(1,2,[3,4]) - Q(1,3,[5,6,7]), ValueError)
    raises(lambda: Q(1,2,[3,4]) * Q(1,3,[5,6,7]), ValueError)
    raises(lambda: Q(1,2,[3,4]) * Z(1,3,[5,6,7]), ValueError)
    raises(lambda: Z(1,2,[3,4]) * Q(1,3,[5,6,7]), ValueError)
    A = Q([[3,4],[5,7]]) / 11
    X = Q([[1,2],[3,4]])
    B = A*X
    assert A.solve(B) == X
    for algorithm in None, "fflu", "dixon":
        assert A.solve(B, algorithm=algorithm) == X
    assert raises(lambda: A.solve(B, algorithm="invalid"), ValueError)
    assert raises(lambda: A.solve(None), TypeError)
    assert raises(lambda: A.solve([1,2]), TypeError)
    assert raises(lambda: A.solve(Q([[1,2]])), ValueError)
    assert raises(lambda: Q([[1,2],[2,4]]).solve(Q([[1],[2]])), ZeroDivisionError)
    M = Q([[1,2,3],[flint.fmpq(1,2),5,6]])
    Mcopy = Q(M)
    Mrref = Q([[1,0,flint.fmpq(3,4)],[0,1,flint.fmpq(9,8)]])
    assert M.rref() == (Mrref, 2)
    assert M != Mrref
    assert M == Mcopy
    assert M.rref(inplace=True) == (Mrref, 2)
    assert M == Mrref
    assert M != Mcopy
    assert Q.hilbert(1, 2) == Q([[1,flint.fmpq(1,2)]])
    M2 = Q([[2, 4, 6], [8, 10, 12]]) / 4
    assert M2.numer_denom() == (Q([[1,2,3],[4,5,6]]), 2)
    half = flint.fmpq(1,2)
    M3 = Q([[half,0,0],[0,half,1],[0,0,half]])
    assert M3.charpoly() == flint.fmpq_poly([-1,6,-12,8]) / 8
    assert M3.minpoly() == flint.fmpq_poly([1,-4,4]) / 4

def test_fmpq_series():
    Qp = flint.fmpq_poly
    Q = flint.fmpq_series
    Zp = flint.fmpz_poly
    Z = flint.fmpz_series
    ctx = flint.ctx
    assert ctx.cap == 10
    s1 = Q([1,2])
    s2 = Q([1,2])
    s3 = Q([1,1])
    s4 = Q([1,2],1,11)
    p1 = Qp([1,2])
    sz1 = Z([1,2])
    sz2 = Z([1,1])
    sz3 = Z([1,1],11)
    assert s1._equal_repr(s1) is True
    assert s1._equal_repr(s2) is True
    assert s1._equal_repr(s3) is False
    assert s1._equal_repr(s4) is False
    assert s1._equal_repr(p1) is False
    assert s1._equal_repr(sz1) is False
    assert s1._equal_repr(sz2) is False
    assert s1._equal_repr(sz3) is False
    assert Q([1])._equal_repr(flint.fmpq(1)) is False
    assert Q([1])._equal_repr(1) is False
    # XXX: this gives a core dump:
    # s = Q([1,2])
    # s[10**10] = 1
    assert Q([])._equal_repr(Q([0]))
    assert not Q([1])._equal_repr(Q([0]))
    assert Q(Q([1]))._equal_repr(Q([1]))
    assert Q(Qp([1]))._equal_repr(Q([1]))
    assert Q(Zp([1]))._equal_repr(Q([1]))
    assert Q(Z([1]))._equal_repr(Q([1]))
    assert Q(1)._equal_repr(Q([1]))
    assert Q([1],1)._equal_repr(Q([1]))
    assert Q([1],2)._equal_repr(Q([flint.fmpq(1,2)]))
    assert not Q([1],1,10)._equal_repr(Q([1],1,11))
    assert Q([1,2],3)._equal_repr(Q([flint.fmpq(1,3), flint.fmpq(2,3)]))
    assert raises(lambda: Q([1],[]), TypeError)
    assert raises(lambda: Q([1],0), ZeroDivisionError)
    assert raises(lambda: Q([1]) < Q([1]), TypeError)
    assert len(Q([1,2])) == 2
    assert Q([1,2]).length() == 2
    s5 = Q([1,2])
    assert s5[1] == 2
    assert s5[2] == 0
    assert s5[-1] == 0
    # XXX: This goes beyond cap. Should it give an error?
    assert s5[100] == 0
    s5[2] = -1
    assert s5[2] == -1
    assert s5._equal_repr(Q([1,2,-1]))
    def set_bad():
        s5[-1] = 3
    assert raises(set_bad, ValueError)
    assert Q([1,2,1]).coeffs() == list(Q([1,2,1])) == [1,2,1]
    assert Q([1,2,1],2).coeffs() == [flint.fmpq(1,2),1,flint.fmpq(1,2)]
    assert Q([1,2,0,4]).str() == "1 + 2*x + 4*x^3 + O(x^10)"
    assert Q([1,2,0,4]).repr() == "fmpq_series([1, 2, 0, 4], 1, prec=10)"
    assert Q([],1,0).str() == "O(x^0)"
    assert Q([],1,-1).str() == "(invalid power series)"
    assert (+Q([1,2]))._equal_repr(Q([1,2]))
    assert (-Q([1,2]))._equal_repr(Q([-1,-2]))
    assert (Q([1,2]) + Q([3,4,5]))._equal_repr(Q([4,6,5]))
    assert (Q([1,2]) + 1)._equal_repr(Q([2,2]))
    assert (Q([1,2]) + Qp([3,4,5]))._equal_repr(Q([4,6,5]))
    assert (1 + Q([1,2]))._equal_repr(Q([2,2]))
    assert (Qp([1,2]) + Q([3,4,5]))._equal_repr(Q([4,6,5]))
    assert raises(lambda: Q([1,2]) + [], TypeError)
    assert raises(lambda: [] + Q([1,2]), TypeError)
    assert (Q([1,2]) - Q([3,5]))._equal_repr(Q([-2,-3]))
    assert (Q([1,2]) - 1)._equal_repr(Q([0,2]))
    assert (1 - Q([1,2]))._equal_repr(Q([0,-2]))
    assert raises(lambda: [] - Q([1,2]), TypeError)
    assert raises(lambda: Q([1,2]) - [], TypeError)
    assert (Q([1,2]) * Q([1,2]))._equal_repr(Q([1,4,4]))
    assert (2 * Q([1,2]))._equal_repr(Q([2,4]))
    assert (Q([1,2]) * 2)._equal_repr(Q([2,4]))
    assert raises(lambda: [] * Q([1,2]), TypeError)
    assert raises(lambda: Q([1,2]) * [], TypeError)
    assert Q([1,2]).valuation() == 0
    assert Q([0,2]).valuation() == 1
    assert Q([0,0]).valuation() == -1
    assert (Q([1,1]) / Q([1,-1]))._equal_repr(Q([1,2,2,2,2,2,2,2,2,2]))
    assert raises(lambda: Q([1,1]) / Q([]), ZeroDivisionError)
    assert (Q([],1) / Q([1,-1]))._equal_repr(Q([],1))
    # quotient would not be a power series
    assert raises(lambda: Q([1,1]) / Q([0,1]), ValueError)
    assert (Q([1,1]) / Q([2,1]))._equal_repr(Q([512, 256, -128, 64, -32, 16, -8, 4, -2, 1], 1024))
    assert (Q([0,1,1]) / Q([0,1,-1]))._equal_repr(Q([1,2,2,2,2,2,2,2,2],1,9))
    assert (Q([2,4]) / 2)._equal_repr(Q([1,2]))
    assert (Q([1,4]) / 2)._equal_repr(Q([1,4],2))
    assert (Q([1,1]) / -1)._equal_repr(Q([-1,-1]))
    assert raises(lambda: Q([1,1]) / [], TypeError)
    assert raises(lambda: [] / Q([1,1]), TypeError)
    q = Q([1,2],3)
    assert (q ** 0)._equal_repr(Q([1]))
    assert (q ** 1)._equal_repr(q)
    assert (q ** 2)._equal_repr(Q([1,4,4],9))
    assert (q ** 3)._equal_repr(Q([1,6,12,8],27))
    assert (q ** 4)._equal_repr(Q([1, 8, 24, 32, 16], 81))
    assert (Q([1,2]) ** 2)._equal_repr(Q([1,4,4]))
    assert (Q([1,2]) ** 2)._equal_repr(Q([1,4,4]))
    assert raises(lambda: pow(Q([1,2]), 3, 5), NotImplementedError)
    assert (Q([1,2])(Q([0,1,2])))._equal_repr(Q([1,2,4]))
    assert raises(lambda: Q([1,2])(Q([1,2])), ValueError)
    assert raises(lambda: Q([1,2])([]), TypeError)
    coeffs = [0, 1, -2, 8, -40, 224, -1344, 8448, -54912, 366080]
    assert Q([0,1,2]).reversion()._equal_repr(Q(coeffs))
    # power series reversion must have valuation 1
    assert raises(lambda: Q([1,1]).reversion(), ValueError)
    assert Q([0,2,1]).reversion()._equal_repr( \
        Q([0,32768,-8192,4096,-2560,1792,-1344,1056,-858,715],65536,10))
    x = Q([0,1])
    expx = x.exp()
    assert (expx)._equal_repr(Q([362880,362880,181440,60480,15120,3024,504,72,9,1],362880))
    assert (expx.inv())._equal_repr(Q([362880,-362880,181440,-60480,15120,-3024,504,-72,9,-1], 362880))
    assert (expx.derivative())._equal_repr(Q([40320,40320,20160,6720,1680,336,56,8,1],40320,prec=9))
    assert (Q([1,1]).integral())._equal_repr(Q([0,2,1],2))
    assert (expx.sqrt())._equal_repr(\
        Q([185794560,92897280,23224320,3870720,483840,48384,4032,288,18,1],185794560))
    assert (expx.rsqrt())._equal_repr(\
        Q([185794560,-92897280,23224320,-3870720,483840,-48384,4032,-288,18,-1],185794560))
    assert (expx.log())._equal_repr(x)
    zero = Q()
    one = Q([1])
    assert (zero.exp())._equal_repr(one)
    assert (one.log())._equal_repr(zero)
    assert (x.atan())._equal_repr(Q([0,315,0,-105,0,63,0,-45,0,35],315))
    assert (x.atanh())._equal_repr(Q([0,315,0,105,0,63,0,45,0,35],315))
    assert (x.asin())._equal_repr(Q([0,40320,0,6720,0,3024,0,1800,0,1225],40320))
    assert (x.asinh())._equal_repr(Q([0,40320,0,-6720,0,3024,0,-1800,0,1225],40320))
    assert (x.sin())._equal_repr(Q([0,362880,0,-60480,0,3024,0,-72,0,1],362880))
    assert (x.cos())._equal_repr(Q([40320,0,-20160,0,1680,0,-56,0,1],40320))
    assert (x.tan())._equal_repr(Q([0,2835,0,945,0,378,0,153,0,62],2835))
    assert (x.sinh())._equal_repr(Q([0,362880,0,60480,0,3024,0,72,0,1],362880))
    assert (x.cosh())._equal_repr(Q([40320,0,20160,0,1680,0,56,0,1],40320))
    assert (x.tanh())._equal_repr(Q([0,2835,0,-945,0,378,0,-153,0,62],2835))
    # Constant term must be nonzero
    assert raises(lambda: Q([0,1]).inv(), ValueError)
    # Constant term must be 1
    assert raises(lambda: Q([2,1]).sqrt(), ValueError)
    assert raises(lambda: Q([2,1]).rsqrt(), ValueError)
    assert raises(lambda: Q([2,1]).log(), ValueError)
    assert raises(lambda: Q([]).sqrt(), ValueError)
    assert raises(lambda: Q([]).rsqrt(), ValueError)
    assert raises(lambda: Q([]).log(), ValueError)
    # Constant term must be 0
    assert raises(lambda: Q([1,1]).exp(), ValueError)
    assert raises(lambda: Q([1,1]).atan(), ValueError)
    assert raises(lambda: Q([1,1]).atanh(), ValueError)
    assert raises(lambda: Q([1,1]).asin(), ValueError)
    assert raises(lambda: Q([1,1]).asinh(), ValueError)
    assert raises(lambda: Q([1,1]).sin(), ValueError)
    assert raises(lambda: Q([1,1]).cos(), ValueError)
    assert raises(lambda: Q([1,1]).tan(), ValueError)
    assert raises(lambda: Q([1,1]).sinh(), ValueError)
    assert raises(lambda: Q([1,1]).cosh(), ValueError)
    assert raises(lambda: Q([1,1]).tanh(), ValueError)

def test_nmod():
    G = flint.nmod
    assert G(0,2) == G(2,2) == G(-2,2)
    assert G(1,2) != G(0,2)
    assert G(0,2) != G(0,3)
    assert G(3,5) == G(8,5)
    assert G(1,2) != (1,2)
    assert isinstance(hash(G(3, 5)), int)
    assert raises(lambda: G([], 3), TypeError)
    #assert G(3,5) == 8        # do we want this?
    #assert 8 == G(3,5)
    assert G(3,5) != 7
    assert 7 != G(3,5)
    assert raises(lambda: G(3,5) < G(2,5), TypeError)
    assert bool(G(0,5)) is False
    assert bool(G(2,5)) is True
    assert G(-3,5) == -G(3,5) == G(2,5) == +G(2,5)
    assert G(2,5) + G(1,5) == G(3,5)
    assert G(2,5) + 1 == G(3,5)
    assert 1 + G(2,5) == G(3,5)
    assert G(2,5) - G(3,5) == G(4,5)
    assert G(2,5) - 3 == G(4,5)
    assert 3 - G(2,5) == G(1,5)
    assert G(2,5) * G(3,5) == G(1,5)
    assert G(2,5) * 3 == G(1,5)
    assert 3 * G(2,5) == G(1,5)
    assert G(3,17) / G(2,17) == G(10,17)
    assert G(3,17) / 2 == G(10,17)
    assert 3 / G(2,17) == G(10,17)
    assert G(0,3) / G(1,3) == G(0,3)
    assert G(3,17) * flint.fmpq(11,5) == G(10,17)
    assert G(3,17) / flint.fmpq(11,5) == G(6,17)
    assert raises(lambda: G(flint.fmpq(2, 3), 3), ZeroDivisionError)
    assert raises(lambda: G(2,5) / G(0,5), ZeroDivisionError)
    assert raises(lambda: G(2,5) / 0, ZeroDivisionError)
    assert G(1,6) / G(5,6) == G(5,6)
    assert raises(lambda: G(1,6) / G(3,6), ZeroDivisionError)
    assert G(1,3) ** 2 == G(1,3)
    assert G(2,3) ** flint.fmpz(2) == G(1,3)
    assert ~G(2,7) == G(2,7) ** -1 == G(4,7)
    assert raises(lambda: G(3,6) ** -1, ZeroDivisionError)
    assert raises(lambda: ~G(3,6), ZeroDivisionError)
    assert raises(lambda: pow(G(1,3), 2, 7), TypeError)
    assert G(flint.fmpq(2, 3), 5) == G(4,5)
    assert raises(lambda: G(2,5) ** G(2,5), TypeError)
    assert raises(lambda: flint.fmpz(2) ** G(2,5), TypeError)
    assert raises(lambda: G(2,5) + G(2,7), ValueError)
    assert raises(lambda: G(2,5) - G(2,7), ValueError)
    assert raises(lambda: G(2,5) * G(2,7), ValueError)
    assert raises(lambda: G(2,5) / G(2,7), ValueError)
    assert raises(lambda: G(2,5) + [], TypeError)
    assert raises(lambda: G(2,5) - [], TypeError)
    assert raises(lambda: G(2,5) * [], TypeError)
    assert raises(lambda: G(2,5) / [], TypeError)
    assert raises(lambda: G(2,5) ** [], TypeError)
    assert raises(lambda: [] + G(2,5), TypeError)
    assert raises(lambda: [] - G(2,5), TypeError)
    assert raises(lambda: [] * G(2,5), TypeError)
    assert raises(lambda: [] / G(2,5), TypeError)
    assert raises(lambda: [] ** G(2,5), TypeError)
    assert G(3,17).modulus() == 17
    assert str(G(3,5)) == "3"
    assert G(3,5).repr() == "nmod(3, 5)"

def test_nmod_poly():
    N = flint.nmod
    P = flint.nmod_poly
    Z = flint.fmpz_poly
    assert P([],17) == P([0],17)
    assert P([1,2,3],17) == P([1,2,3],17)
    assert P([1,2,3],17) != P([1,2,3],15)
    assert P([1,2,3],17) != P([1,2,4],15)
    assert P([1,2,3],17) != 1
    assert P([1,2,3],17) != Z([1,2,3])
    assert raises(lambda: P([1,2],3) < P([1,2],3), TypeError)
    assert P(Z([1,2,3]),17) == P([1,2,3],17)
    assert P([1,2,N(3,17)],17) == P([1,2,3],17)
    assert P(P([1,2],17),17) == P([1,2],17)
    assert raises(lambda: P(P([1,2],17),13), ValueError)
    assert raises(lambda: P([1,2,[]],17), TypeError)
    assert raises(lambda: P([1,2,flint.nmod(3,15)],17), ValueError)
    assert raises(lambda: P([1,2],0), ValueError)
    assert raises(lambda: P({},3), TypeError)
    assert P([1,2,3],17).degree() == 2
    assert P([1,2,3],17).length() == 3
    assert len(P([1,2,3],17)) == 3
    assert P([1,2,3],5).coeffs() == [N(1,5),N(2,5),N(3,5)]
    assert P([1,2,3],17) + 2 == P([3,2,3],17)
    assert 2 + P([1,2,3],17) == P([3,2,3],17)
    assert P([1,2,3],17) + P([3,4,5],17) == P([4,6,8],17)
    assert P([1,2,3],17) + P([3,4,5],17) == P([4,6,8],17)
    assert P([1,2,3],17) - 2 == P([16,2,3],17)
    assert 2 - P([1,2,3],17) == -P([16,2,3],17)
    assert +P([1,2,3],17) == P([1,2,3],17)
    assert P([1,2,3],17) - P([3,4,6],17) == P([15,15,14],17)
    assert P([1,2,3],17) * 2 == P([2,4,6],17)
    assert 2 * P([1,2,3],17) == P([2,4,6],17)
    assert P([1,2,3],17) * P([1,2,3],17) == P([1,4,10,12,9], 17)
    assert P([1,2,3],17) * Z([1,2,3]) == P([1,4,10,12,9], 17)
    assert Z([1,2,3]) * P([1,2,3],17) == P([1,4,10,12,9], 17)
    assert P([1,2,3,4,5],17) % P([2,3,4],17) == P([12,12],17)
    assert P([1,2,3,4,5],17) // P([2,3,4],17) == P([3,16,14],17)
    assert raises(lambda: P([1,2],5) // P([],5), ZeroDivisionError)
    assert raises(lambda: P([1,2],5) % P([],5), ZeroDivisionError)
    assert raises(lambda: divmod(P([1,2],5), P([],5)), ZeroDivisionError)
    assert P([1,2,3,4,5],17) ** 2 == P([1,2,3,4,5],17) * P([1,2,3,4,5],17)
    assert P([1,2,3],17) * N(3,17) == P([3,6,9],17)
    s = P([1,2,3],17)
    s2 = P([1,2,3],5)
    assert raises(lambda: s + s2, ValueError)
    assert raises(lambda: s - s2, ValueError)
    assert raises(lambda: s * s2, ValueError)
    assert raises(lambda: s // s2, ValueError)
    assert raises(lambda: s % s2, ValueError)
    assert raises(lambda: s2 + s, ValueError)
    assert raises(lambda: s2 - s, ValueError)
    assert raises(lambda: s2 * s, ValueError)
    assert raises(lambda: s2 // s, ValueError)
    assert raises(lambda: s2 % s, ValueError)
    assert raises(lambda: s + [], TypeError)
    assert raises(lambda: s - [], TypeError)
    assert raises(lambda: s * [], TypeError)
    assert raises(lambda: s // [], TypeError)
    assert raises(lambda: s % [], TypeError)
    assert raises(lambda: [] + s, TypeError)
    assert raises(lambda: [] - s, TypeError)
    assert raises(lambda: [] * s, TypeError)
    assert raises(lambda: [] // s, TypeError)
    assert raises(lambda: [] % s, TypeError)
    assert raises(lambda: pow(P([1,2],3), 3, 4), NotImplementedError)
    assert str(P([1,2,3],17)) == "3*x^2 + 2*x + 1"
    assert P([1,2,3],17).repr() == "nmod_poly([1, 2, 3], 17)"
    p = P([3,4,5],17)
    assert p(14) == N(2,17)
    assert p(P([1,2,3],17)) == P([12,11,11,9,11],17)
    assert raises(lambda: p({}), TypeError)
    p2 = P([3,4,5],17)
    assert p2[1] == N(4,17)
    assert p2[-1] == N(0,17)
    p2[1] = N(2,17)
    assert p2 == P([3,2,5],17)
    p2[1] = 6
    assert p2 == P([3,6,5],17)
    def set_bad1():
        p2[-1] = 3
    def set_bad2():
        p2[2] = []
    assert raises(set_bad1, ValueError)
    assert raises(set_bad2, TypeError)
    assert bool(P([], 5)) is False
    assert bool(P([1], 5)) is True
    assert P([1,2,1],3).gcd(P([1,1],3)) == P([1,1],3)
    raises(lambda: P([1,2],3).gcd([]), TypeError)
    raises(lambda: P([1,2],3).gcd(P([1,2],5)), ValueError)
    p3 = P([1,2,3,4,5,6],7)
    f3 = (N(6,7), [(P([6, 1],7), 5)])
    assert p3.factor() == f3
    # XXX: factor ignores an invalid algorithm string
    for alg in [None, 'berlekamp', 'cantor-zassenhaus']:
        assert p3.factor(alg) == f3
        assert p3.factor(algorithm=alg) == f3
    assert P([1], 11).roots() == []
    assert P([1, 2, 3], 11).roots() == [(8, 1), (6, 1)]
    assert P([1, 6, 1, 8], 11).roots() == [(5, 3)]

def test_nmod_mat():
    M = flint.nmod_mat
    G = flint.nmod
    Z = flint.fmpz_mat
    a = M(2,3,[1,2,3,4,5,6],17)
    b = M(2,3,[4,5,6,7,8,9],17)
    assert a == a
    assert a == M(a)
    assert a != b
    assert a.nrows() == 2
    assert a.ncols() == 3
    assert a.entries() == [G(x,17) for x in [1,2,3,4,5,6]]
    assert a.table() == [[G(x,17) for x in [1,2,3]], [G(x,17) for x in [4,5,6]]]
    assert (a + b).entries() == [G(x,17) for x in [5,7,9,11,13,15]]
    assert raises(a.det, ValueError)
    assert +a == a
    assert -a == M(2,3,[-1,-2,-3,-4,-5,-6],17)
    c = M(2,2,[1,2,3,4],17)
    assert c.det() == G(-2,17)
    assert raises(lambda: a + c, ValueError)
    assert (a * 3).entries() == [G(x,17) for x in [3,6,9,12,15,18]]
    assert (3 * a).entries() == [G(x,17) for x in [3,6,9,12,15,18]]
    assert (a * long(3)).entries() == [G(x,17) for x in [3,6,9,12,15,18]]
    assert (long(3) * a).entries() == [G(x,17) for x in [3,6,9,12,15,18]]
    assert (a * flint.fmpz(3)).entries() == [G(x,17) for x in [3,6,9,12,15,18]]
    assert (flint.fmpz(3) * a).entries() == [G(x,17) for x in [3,6,9,12,15,18]]
    assert M(2,2,[1,1,2,2],17).rank() == 1
    assert M(Z(2,2,[1,2,3,4]),17) == M(2,2,[1,2,3,4],17)
    A = M(5,3,Z.randbits(5,3,5).entries(),17)
    B = M(3,7,Z.randtest(3,7,5).entries(),17)
    C = M(7,2,Z.randtest(7,2,5).entries(),17)
    assert A*(B*C) == (A*B)*C
    assert bool(M(2,2,[0,0,0,0],17)) == False
    assert bool(M(2,2,[0,0,0,1],17)) == True
    ctx.pretty = False
    assert repr(M(2,2,[1,2,3,4],17)) == 'nmod_mat(2, 2, [1, 2, 3, 4], 17)'
    ctx.pretty = True
    assert str(M(2,2,[1,2,3,4],17)) == '[1, 2]\n[3, 4]'
    assert repr(M(2,2,[1,2,3,4],17)) == '[1, 2]\n[3, 4]'
    assert M(1,2,[3,4],17) / 3 == M(1,2,[3,4],17) * (~G(3,17))
    assert M(2,2,[1,2,3,4], 17).inv().det() == ~(M(2,2,[1,2,3,4], 17).det())
    assert M(2,2,[1,2,3,4], 17).inv().inv() == M(2,2,[1,2,3,4], 17)
    assert M(2,2,[0,1,2,3],17) * M(2, 2, [2,3,4,5], 17) == M(2,2,[4,5,16,4],17)
    assert raises(lambda: M([1], 5), TypeError)
    assert raises(lambda: M([[1],[2,3]], 5), ValueError)
    assert raises(lambda: M([[1],[2]], 0), ValueError)
    assert raises(lambda: M(None), TypeError)
    assert raises(lambda: M(None,17), TypeError)
    assert M(2,3,17) == M(2,3,[0,0,0,0,0,0],17)
    assert raises(lambda: M(2,3,[0,0,0,0,0],17), ValueError)
    assert raises(lambda: M(2,3,[0,1],[1,2],17), ValueError)
    assert M([[1,2,3],[4,5,6]], 5) == M(2,3,[1,2,3,4,5,6], 5)
    assert raises(lambda: M([[0]],13) < M([[1]],13), TypeError)
    assert (M([[1]],17) == M([[1]],13)) is False
    assert (M([[1]],17) != M([[1]],13)) is True
    assert (M([[1]],17) == None) is False
    assert (M([[1]],17) != None) is True
    M2 = M.randtest(3,4,5)
    assert all(0 <= int(x) < 5 for x in M2.entries())
    assert (M2.nrows(), M2.ncols()) == (3, 4)
    M3 = M(2,2,[1,2,3,4],17)
    assert M3[0,1] == G(2,17)
    M3_copy = M(M3)
    M3[0,1] = -1
    assert M3[0,1] == G(-1,17)
    def set_bad(i,j):
        M3[i,j] = 2
    # XXX: negative indices should be allowed
    assert raises(lambda: M3[-1,0], ValueError)
    assert raises(lambda: M3[0,-1], ValueError)
    assert raises(lambda: set_bad(-1,0), ValueError)
    assert raises(lambda: set_bad(0,-1), ValueError)
    # XXX: Should be IndexError
    assert raises(lambda: M3[2,0], ValueError)
    assert raises(lambda: M3[0,2], ValueError)
    assert raises(lambda: set_bad(2,0), ValueError)
    assert raises(lambda: set_bad(0,2), ValueError)
    def set_bad2():
        M3[0,0] = 1.5
    assert raises(set_bad2, ValueError)
    assert raises(lambda: M3 + [], TypeError)
    assert raises(lambda: M3 - [], TypeError)
    assert raises(lambda: M3 * [], TypeError)
    assert raises(lambda: M3 / [], TypeError)
    assert raises(lambda: [] + M3, TypeError)
    assert raises(lambda: [] - M3, TypeError)
    assert raises(lambda: [] * M3, TypeError)
    assert raises(lambda: [] / M3, TypeError)
    assert raises(lambda: M([[1]],3) + M([[1]],5), ValueError)
    assert raises(lambda: M([[1]],3) - M([[1]],5), ValueError)
    assert raises(lambda: M([[1]],3) * M([[1]],5), ValueError)
    assert Z([[1]]) + M([[1]],3) == M([[2]],3)
    assert M([[1]],3) + Z([[1]]) == M([[2]],3)
    assert Z([[1]]) - M([[1]],3) == M([[0]],3)
    assert M([[1]],3) - Z([[1]]) == M([[0]],3)
    assert Z([[1]]) * M([[1]],3) == M([[1]],3)
    assert M([[1]],3) * Z([[1]]) == M([[1]],3)
    assert raises(lambda: M([[1]],3) - M([[1,2]],3), ValueError)
    assert raises(lambda: M([[1,2]],3) * M([[1,2]],3), ValueError)
    M4 = M([[1,2],[3,4]],17)
    assert M4.inv() == M([[15,1],[10,8]],17)
    assert raises(lambda: M([[1,2]],17).inv(), ValueError)
    assert raises(lambda: M([[1,2],[2,4]],17).inv(), ZeroDivisionError)
    assert M([[1,2,3],[4,5,6]],17).transpose() == M([[1,4],[2,5],[3,6]],17)
    M5 = M([[1,2],[3,4]],17)
    X = M([[1],[2]],17)
    b = M5*X
    assert M5.solve(b) == X
    assert raises(lambda: M5.solve([]), TypeError)
    assert raises(lambda: b.solve(M5), ValueError)
    assert raises(lambda: M([[1,2],[2,4]],17).solve(b), ZeroDivisionError)
    M6 = M([[1,2,3],[4,5,6]],17)
    M6_rref = M([[1,0,16],[0,1,2]],17)
    M6_copy = M(M6)
    assert M6.rref() == (M6_rref, 2)
    assert M6 == M6_copy
    assert M6.rref(inplace=True) == (M6_rref, 2)
    assert M6 == M6_rref
    M6 = M6_copy
    assert M6.nullspace() == (M([[1,15,1],[0,0,0],[0,0,0]],17).transpose(), 1)

def test_nmod_series():
    # XXX: currently no code in nmod_series.pyx
    pass

def test_arb():
    A = flint.arb
    assert A(3) > A(2.5)
    assert A(3) >= A("2.5")
    assert A(3) < A((3,1))
    assert A(3) <= A("inf")
    assert A(3) == A(3)
    assert A(3) != A(2)
    assert not (A("1.1") == A("1.1"))

def test_pickling():
    objects = [
        flint.fmpz(1),
        flint.fmpq(1,2),
        # XXX: Add pickling for everything else
    ]
    for obj in objects:
        s = pickle.dumps(obj)
        obj2 = pickle.loads(s)
        assert obj == obj2

def test_fmpz_mod():
    from flint import fmpz_mod_ctx, fmpz, fmpz_mod
    
    p_sml = 163
    p_med = 2**127 - 1
    p_big = 2**255 - 19

    F_cmp = fmpz_mod_ctx(10)
    F_sml = fmpz_mod_ctx(p_sml)
    F_med = fmpz_mod_ctx(p_med)
    F_big = fmpz_mod_ctx(p_big)

    assert F_sml.is_prime() is True
    assert F_med.is_prime() is True
    assert F_big.is_prime() is True
    assert F_cmp.is_prime() is False


    # Context tests
    assert raises(lambda: fmpz_mod_ctx("AAA"), TypeError)
    assert raises(lambda: fmpz_mod_ctx(-1), ValueError)
    assert F_sml.modulus() == p_sml
    assert F_med.modulus() == p_med
    assert F_big.modulus() == p_big

    F_big_copy = fmpz_mod_ctx(p_big)
    assert F_big_copy == F_big
    assert F_big != F_sml
    assert hash(F_big_copy) == hash(F_big)
    assert hash(F_big) != hash(F_sml)
    assert F_big_copy != F_sml
    assert F_big_copy != "A"

    assert repr(F_sml) == "fmpz_mod_ctx(163)"
    assert str(F_sml) == "Context for fmpz_mod with modulus: 163"

    # Type tests
    assert raises(lambda: fmpz_mod(1, "AAA"), TypeError)

    # Test for small, medium and large char.
    for F_test in [F_sml, F_med, F_big]:
        test_mod = int(F_test.modulus())
        test_x = (-123) % test_mod # canonical value
        test_y = ((-456) % test_mod)**2 # non-canoncial value

        test_z = F_test.random_element()
        assert int(test_z) < F_test.modulus()
        assert int(test_z) >= 0

        assert raises(lambda: F_test(F_cmp(1)), ValueError)
        assert raises(lambda: F_test("abc"), NotImplementedError)

        F_test_copy = fmpz_mod_ctx(test_mod)
        F_other = fmpz_mod_ctx(11)

        assert raises(lambda: F_test(test_x) > 0, TypeError)
        assert raises(lambda: F_test(test_x) >= 0, TypeError)
        assert raises(lambda: F_test(test_x) < 0, TypeError)
        assert raises(lambda: F_test(test_x) <= 0, TypeError)

        assert (test_x == F_test(test_x)) is True, f"{test_x}, {F_test(test_x)}"
        assert (124 != F_test(test_x)) is True
        assert (F_test(test_x) == test_x) is True
        assert (F_test(test_x) == test_x + test_mod) is True
        assert (F_test(test_x) == 1) is False
        assert (F_test(test_x) != 1) is True
        assert (F_test(test_x) == F_test(test_x)) is True
        assert (F_test(test_x) == F_test(test_x + test_mod)) is True
        assert (F_test(test_x) == F_test(1)) is False
        assert (F_test(test_x) != F_test(1)) is True
        assert (F_test(test_x) != "abc") is True

        assert (hash(F_test(test_x)) == hash(test_x)) is True
        assert (hash(F_test(F_test(test_x))) == hash(test_x)) is True
        assert (hash(F_test(test_x)) == hash(1)) is False
        assert (hash(F_test(test_x)) != hash(1)) is True
        assert (hash(F_test(test_x)) == hash(F_test(test_x))) is True
        assert (hash(F_test(test_x)) == hash(F_test(test_x + test_mod))) is True
        assert (hash(F_test(test_x)) == hash(F_test(1))) is False
        assert (hash(F_test(test_x)) != hash(F_test(1))) is True

        # Is one, zero
        assert (F_test(0) == 0) is True
        assert F_test(0).is_zero() is True
        assert not F_test(0)
        assert not F_test(test_mod)
        assert F_test(1).is_one() is True
        assert F_test(test_mod + 1).is_one() is True
        assert F_test(1).is_one() is True
        assert F_test(2).is_one() is False

        # int, str, repr
        assert str(F_test(11)) == "11"
        assert str(F_test(-1)) == str(test_mod - 1)
        assert repr(F_test(11)) == f"fmpz_mod(11, {test_mod})"
        assert repr(F_test(-1)) == f"fmpz_mod({test_mod - 1}, {test_mod})"

        assert +F_test(5) == F_test(5)

        # Arithmetic tests

        # Negation
        assert -F_test(test_x) == F_test(-test_x) == (-test_x % test_mod)
        assert -F_test(1) == F_test(-1) == F_test(test_mod - 1)

        # Addition
        assert F_test(test_x) + F_test(test_y) == F_test(test_x + test_y)
        assert F_test(test_x) + F_test_copy(test_y) == F_test(test_x + test_y)
        assert F_test(test_x) + F_test(test_y) == F_test_copy(test_x + test_y)
        assert raises(lambda: F_test(test_x) + "AAA", TypeError)

        assert F_test(test_x) + F_test(test_y) == F_test(test_y) + F_test(test_x)
        assert F_test(test_x) + test_y == F_test(test_x + test_y)
        assert test_y + F_test(test_x) == F_test(test_x + test_y)
        assert F_test(test_x) + fmpz(test_y) == F_test(test_y) + F_test(test_x)
        assert raises(lambda: F_test(test_x) + F_other(test_y), ValueError)

        # Subtraction

        assert F_test(test_x) - F_test(test_y) == F_test(test_x - test_y)
        assert F_test(test_x) - test_y == F_test(test_x - test_y)
        assert F_test(test_x) - test_y == F_test(test_x) - F_test(test_y)
        assert F_test(test_y) - test_x == F_test(test_y) - F_test(test_x)
        assert test_x - F_test(test_y) == F_test(test_x) - F_test(test_y)
        assert test_y - F_test(test_x) == F_test(test_y) - F_test(test_x)
        assert F_test(test_x) - fmpz(test_y) == F_test(test_x) - F_test(test_y)
        assert raises(lambda: F_test(test_x) - F_other(test_y), ValueError)
        assert raises(lambda: F_test(test_x) - "AAA", TypeError)
        assert raises(lambda: "AAA" - F_test(test_x), TypeError)


        # Multiplication

        assert F_test(test_x) * F_test(test_y) == (test_x * test_y) % test_mod
        assert F_test(test_x) * test_y == (test_x * test_y) % test_mod
        assert test_y * F_test(test_x) == (test_x * test_y) % test_mod

        assert F_test(1) * F_test(test_x) == F_test(1 * test_x)
        assert F_test(2) * F_test(test_x) == F_test(2 * test_x)
        assert F_test(3) * F_test(test_x) == F_test(3 * test_x)
        assert 1 * F_test(test_x) == F_test(1 * test_x)
        assert 2 * F_test(test_x) == F_test(2 * test_x)
        assert 3 * F_test(test_x) == F_test(3 * test_x)
        assert F_test(test_x) * 1 == F_test(1 * test_x)
        assert F_test(test_x) * 2 == F_test(2 * test_x)
        assert F_test(test_x) * 3 == F_test(3 * test_x)
        assert fmpz(1) * F_test(test_x) == F_test(1 * test_x)
        assert fmpz(2) * F_test(test_x) == F_test(2 * test_x)
        assert fmpz(3) * F_test(test_x) == F_test(3 * test_x)
        assert raises(lambda: F_test(test_x) * "AAA", TypeError)
        assert raises(lambda: F_test(test_x) * F_other(test_x), ValueError)

        # Exponentiation 

        assert F_test(0)**0 == pow(0, 0, test_mod)
        assert F_test(0)**1 == pow(0, 1, test_mod)
        assert F_test(0)**2 == pow(0, 2, test_mod)
        assert raises(lambda: F_test(0)**(-1), ZeroDivisionError)
        assert raises(lambda: F_test(0)**("AA"), NotImplementedError)

        assert F_test(test_x)**fmpz(0) == pow(test_x, 0, test_mod)
        assert F_test(test_x)**fmpz(1) == pow(test_x, 1, test_mod)
        assert F_test(test_x)**fmpz(2) == pow(test_x, 2, test_mod)
        assert F_test(test_x)**fmpz(3) == pow(test_x, 3, test_mod)

        assert F_test(test_x)**0 == pow(test_x, 0, test_mod)
        assert F_test(test_x)**1 == pow(test_x, 1, test_mod)
        assert F_test(test_x)**2 == pow(test_x, 2, test_mod)
        assert F_test(test_x)**3 == pow(test_x, 3, test_mod)
        assert F_test(test_x)**100 == pow(test_x, 100, test_mod)

        assert F_test(test_x)**(-1) == pow(test_x, -1, test_mod)
        assert F_test(test_x)**(-2) == pow(test_x, -2, test_mod)
        assert F_test(test_x)**(-3) == pow(test_x, -3, test_mod)
        assert F_test(test_x)**(-4) == pow(test_x, -4, test_mod)

        # Inversion

        assert raises(lambda: ~F_test(0), ZeroDivisionError)
        assert ~F_test(test_x) == pow(test_x, -1, test_mod)
        assert ~F_test(1) == pow(1, -1, test_mod)
        assert ~F_test(2) == pow(2, -1, test_mod), f"Broken!! {~F_test(2)}, {pow(2, -1, test_mod)}"

        assert F_test(1).inverse(check=False) == pow(1, -1, test_mod)
        assert F_test(2).inverse(check=False) == pow(2, -1, test_mod)
        assert F_test(test_x).inverse(check=False) == pow(test_x, -1, test_mod)

        # Division
        assert raises(lambda: F_test(1) / F_test(0), ZeroDivisionError)
        assert F_test(test_x) / F_test(test_y) == (test_x * pow(test_y, -1, test_mod)) % test_mod
        assert F_test(test_x) / fmpz(test_y) == (test_x * pow(test_y, -1, test_mod)) % test_mod
        assert F_test(test_x) / test_y == (test_x * pow(test_y, -1, test_mod)) % test_mod
        assert raises(lambda: F_test(test_x) / "AAA", TypeError)
        assert raises(lambda: "AAA" / F_test(test_x), TypeError)
        assert raises(lambda: F_other(test_x) / F_test(test_x), ValueError)
        assert raises(lambda: F_test(test_x) // F_test(test_x), TypeError)
        assert 1 / F_test(2) ==  pow(2, -1, test_mod)
        assert 1 / F_test(test_x) ==  pow(test_x, -1, test_mod)
        assert 1 / F_test(test_y) ==  pow(test_y, -1, test_mod)

        assert fmpz(test_y) / F_test(test_x) == (test_y * pow(test_x, -1, test_mod)) % test_mod
        assert test_y / F_test(test_x) == (test_y * pow(test_x, -1, test_mod)) % test_mod
            
def test_fmpz_mod_dlog():
    from flint import fmpz, fmpz_mod_ctx

    # Input modulus must be prime
    F = fmpz_mod_ctx(4)
    g, a = F(1), F(2)
    assert raises(lambda: g.discrete_log(a), NotImplementedError)

    # Moduli must match
    F1, F2 = fmpz_mod_ctx(2), fmpz_mod_ctx(3)
    g = F1(2)
    a = F2(4)
    assert raises(lambda: g.discrete_log(a), ValueError)

    # Need to use either fmpz_mod or something which can be case to
    # fmpz
    assert raises(lambda: g.discrete_log("A"), TypeError)

    F = fmpz_mod_ctx(163)
    g = F(2)
    a = g**123
    
    assert 123 == g.discrete_log(a)

    a_int = pow(2, 123, 163)
    a_fmpz = fmpz(a_int)
    assert 123 == g.discrete_log(a_int)
    assert 123 == g.discrete_log(a_fmpz)

    # Randomised testing with smooth large modulus
    e2, e3 = 92, 79
    p = 2**e2 * 3**e3 + 1
    F = fmpz_mod_ctx(p)

    import random
    for _ in range(10):
        g = F(random.randint(0,p))
        for _ in range(10):
            i = random.randint(0,p)
            a = g**i
            x = g.discrete_log(a)
            assert g**x == a

def test_fmpz_mod_poly():
    from flint import fmpz_poly, fmpz_mod_poly, fmpz_mod_poly_ctx, fmpz_mod_ctx, fmpz

    # fmpz_mod_poly_ctx tests
    F = fmpz_mod_ctx(11)
    R1 = fmpz_mod_poly_ctx(F)
    R2 = fmpz_mod_poly_ctx(11)
    R3 = fmpz_mod_poly_ctx(13)

    assert raises(lambda: fmpz_mod_ctx("AAA"), TypeError)
    assert raises(lambda: fmpz_mod_ctx(-1), ValueError)
    assert (R1 == R1) is True
    assert (R1 != R1) is False
    assert (R1 == R2) is True
    assert (R1 != R2) is False
    assert (R1 != R3) is True
    assert (R1 == R3) is False
    assert (R1 != "AAA") is True
    assert (R1 == "AAA") is False

    assert (hash(R1) == hash(R1)) is True
    assert (hash(R1) == hash(R2)) is True
    assert (hash(R1) != hash(R3)) is True

    assert str(R1) == "Context for fmpz_mod_poly with modulus: 11"
    assert str(R1) == str(R2)
    assert repr(R3) == "fmpz_mod_poly_ctx(13)"

    assert R1.modulus() == 11
    
    assert R1.is_prime()
    assert R1.zero() == 0
    assert R1.one() == 1
    assert R1.gen() == R1([0,1])

    # Random testing
    f = R1.random_element()
    assert f.degree() == 3
    f = R1.random_element(degree=5, monic=True)
    assert f.degree() == 5
    assert f.is_monic()
    f = R1.random_element(degree=100, irreducible=True)
    assert f.degree() == 100
    assert f.is_irreducible()
    f = R1.random_element(degree=1, monic=True, irreducible=True)
    assert f.degree() == 1
    assert f.is_irreducible()
    assert f.is_monic()
    assert raises(lambda: R1.random_element(degree=-123), ValueError)
    assert raises(lambda: R1.random_element(monic="A"), ValueError)
    assert raises(lambda: R1.random_element(irreducible="A"), ValueError)


    # Conversion tests
    F = fmpz_mod_ctx(11)
    F_other = fmpz_mod_ctx(10)
    R = fmpz_mod_poly_ctx(F)
    R_other = fmpz_mod_poly_ctx(F_other)

    assert raises(lambda: fmpz_mod_poly(1, "A"), TypeError) # Need a valid context
    assert raises(lambda: R(R_other([1,2,3])), ValueError), f"{R(R_other([1,2,3])) = }" # moduli must match
    assert raises(lambda: R(F_other(2)), ValueError) # moduli must match
    assert raises(lambda: R([F(1), F_other(2)]), ValueError) # moduli must match
    assert raises(lambda: R([F(1), "A"]), TypeError) # need to be able to cast to fmpz_mod

    f1 = R([int(-1),int(-2),int(-3)])
    f2 = R([fmpz(-1),fmpz(-2),fmpz(-3)])
    f3 = R([F(-1),F(-2),F(-3)])
    f4 = R(fmpz_poly([-1, -2, -3]))
    f5 = R(f4)

    assert str(f1) == "8*x^2 + 9*x + 10"
    assert str(f2) == "8*x^2 + 9*x + 10"
    assert str(f3) == "8*x^2 + 9*x + 10"
    assert str(f4) == "8*x^2 + 9*x + 10"
    assert str(f5) == "8*x^2 + 9*x + 10"

    f1 = R(5)
    f2 = R(fmpz(6))
    f3 = R(F(7))
    assert str(f1) == "5"
    assert str(f2) == "6"
    assert str(f3) == "7"

    # Printing
    f = R([5, 6, 7, 8])
    assert str(f) == "8*x^3 + 7*x^2 + 6*x + 5"
    # assert repr(f) == "fmpz_mod_poly([5, 6, 7, 8], fmpz_mod_poly_ctx(11))"

    # Get and Set tests
    f = R([5, 6, 7, 8])
    assert f[0] == 5
    assert repr(f[0]) == "fmpz_mod(5, 11)"
    f[0] = 7
    assert repr(f[0]) == "fmpz_mod(7, 11)"
    assert str(f) == "8*x^3 + 7*x^2 + 6*x + 7"

    # TODO: currently repr does pretty printing
    # just like str, we should address this. Mainly, 
    # the issue is we want nice `repr` behaviour in
    # interactive shells, which currently is why this
    # choice has been made
    assert str(f) == repr(f)

    assert f[-1] == 0
    assert raises(lambda: f.__setitem__(-1, 1), ValueError)
    assert raises(lambda: f.__setitem__(1, "A"), TypeError)


    # Comparisons
    f1 = R([1,2,3])
    f2 = R([12,13,14])
    f3 = R([4,5,6])
    f4 = R([3])

    assert (f1 == f2) is True
    assert (f1 != f3) is True
    assert (f1 != "1") is True
    assert (f4 == 3) is True
    assert (hash(f1) == hash(f2)) is True
    assert raises(lambda: f1 > f2, TypeError)
    assert raises(lambda: f1 >= f2, TypeError)
    assert raises(lambda: f1 < f2, TypeError)
    assert raises(lambda: f1 <= f2, TypeError)

    assert len(f1) == f1.length() == 3
    assert f1.degree() == 2

    f1 = R([0])
    f2 = R([1])
    f3 = R([0, 1])

    assert f1.is_zero() is True
    assert f2.is_one() is True
    assert f3.is_gen() is True

    # Arithmetic
    p_sml = 163
    p_med = 2**127 - 1
    p_big = 2**255 - 19

    F_sml = fmpz_mod_ctx(p_sml)
    F_med = fmpz_mod_ctx(p_med)
    F_big = fmpz_mod_ctx(p_big)
    
    R_sml = fmpz_mod_poly_ctx(F_sml)
    R_med = fmpz_mod_poly_ctx(F_med)
    R_big = fmpz_mod_poly_ctx(F_big)

    F_cmp = fmpz_mod_ctx(10)
    R_cmp = fmpz_mod_poly_ctx(F_cmp)
    f_cmp = R_cmp([1,2,3,4,5])
    f_bad = R_cmp([2,2,2,2,2])

    for (F_test, R_test) in [(F_sml, R_sml), (F_med, R_med), (F_big, R_big)]:
        
        f = R_test([-1,-2])
        g = R_test([-3,-4])

        # pos, neg
        assert f is +f
        assert -f == R_test([1,2])
        
        # add
        assert raises(lambda: f + f_cmp, ValueError)
        assert raises(lambda: f + "AAA", TypeError)
        assert raises(lambda: "AAA" + f, TypeError)
        assert f + g == R_test([-4,-6])
        assert f + 1 ==  R_test([0,-2])
        assert f + fmpz(1) ==  R_test([0,-2])
        assert f + F_test(1) ==  R_test([0,-2])
        assert 1 + f ==  R_test([0,-2])
        assert fmpz(1) + f ==  R_test([0,-2])
        assert F_test(1) + f ==  R_test([0,-2])

        # sub
        assert raises(lambda: f - f_cmp, ValueError)
        assert raises(lambda: f - "AAA", TypeError)
        assert raises(lambda: "AAA" - f, TypeError)
        assert f - g == R_test([2, 2])
        assert f - 1 ==  R_test([-2,-2])
        assert f - fmpz(1) ==  R_test([-2,-2])
        assert f - F_test(1) ==  R_test([-2,-2])
        assert 1 - f ==  R_test([2, 2])
        assert fmpz(1) - f ==  R_test([2, 2])
        assert F_test(1) - f ==  R_test([2, 2])

        # mul
        assert raises(lambda: f * f_cmp, ValueError)
        assert raises(lambda: f * "AAA", TypeError)
        assert raises(lambda: "AAA" * f, TypeError)
        assert f * g == R_test([3, 4 + 6, 8])
        assert f * 2 ==  R_test([-2,-4])
        assert f * fmpz(2) ==  R_test([-2,-4])
        assert f * F_test(2) ==  R_test([-2,-4])
        assert 2 * f ==  R_test([-2,-4])
        assert fmpz(2) * f ==  R_test([-2,-4])
        assert F_test(2) * f ==  R_test([-2,-4])

        # scalar_mul
        assert 2 * f == f.scalar_mul(2)
        assert raises(lambda: f.scalar_mul("AAA"), TypeError)

        # Exact division
        assert raises(lambda: f.exact_division(f_cmp), ValueError)
        assert raises(lambda: f.exact_division("AAA"), TypeError)
        assert raises(lambda: f.exact_division(0), ZeroDivisionError)

        assert (f * g).exact_division(g) == f
        assert raises(lambda: f.exact_division(g), ValueError)

        # true div
        assert raises(lambda: f / "AAA", TypeError)
        assert raises(lambda: f / 0, ZeroDivisionError)
        assert raises(lambda: f_cmp / 2, ZeroDivisionError)
    
        assert (f + f) / 2  ==  f
        assert (f + f) / fmpz(2)  ==  f
        assert (f + f) / F_test(2)  ==  f

        # floor div
        assert raises(lambda: 1 // f_bad, ZeroDivisionError)
        assert raises(lambda: f // f_cmp, ValueError)
        assert raises(lambda: f // "AAA", TypeError)
        assert raises(lambda: "AAA" // f, TypeError)
        assert (f * g) // g == f
        assert (f + f) // 2  ==  f
        assert (f + f) // fmpz(2)  ==  f
        assert (f + f) // F_test(2)  ==  f
        assert 2 // R_test(2) == 1        
        assert (f + 1) // f == 1

        # pow
        assert raises(lambda: f**(-2), ValueError)
        assert f*f == f**2
        assert f*f == f**fmpz(2)

        # Shifts
        assert raises(lambda: R_test([1,2,3]).left_shift(-1), ValueError)
        assert raises(lambda: R_test([1,2,3]).right_shift(-1), ValueError)
        assert R_test([1,2,3]).left_shift(3) == R_test([0,0,0,1,2,3])
        assert R_test([1,2,3]).right_shift(1) == R_test([2,3])

        # Mod
        assert raises(lambda: f % f_bad, ValueError)
        assert raises(lambda: 123 % f_bad, ValueError)
        assert raises(lambda: f % "AAA", TypeError)
        assert raises(lambda: tuple() % f, TypeError), f'{"AAA" % f = }'

        assert f % 1 == 0
        assert R_test.one() % 1 == 0
        assert 100 % R_test.one() == 0
        assert (f*g + 1) % f == 1
        assert (f*g + g) % f == (g % f)
        assert f % R_test([0,1]) == f.constant_coefficient()

        # Evaluation
        h = R_test([0, 1])
        assert h(1) == 1
        assert h(-1) == R_test.modulus() - 1
        h = R_test([0, 0, 1])
        assert h(1) == h(-1)
        assert raises(lambda: h("AAA"), TypeError)
        assert f([-1,-2,-3]) == [f(x) for x in [-1, -2, -3]]

        # compose
        assert raises(lambda: h.compose("AAA"), TypeError)

        # Reverse
        assert raises(lambda: h.reverse(degree=-100), ValueError)
        assert R_test([-1,-2,-3]).reverse() == R_test([-3,-2,-1])

        # monic
        assert raises(lambda: f_bad.monic(), ValueError)
        assert R_test([1,2]).monic() == R_test([1 / F_test(2), 1])
        assert R_test([1,2]).monic(check=False) == R_test([1 / F_test(2), 1])

        # Square
        assert f*f == f**2 == f.square()

        # mulmod
        assert f.mulmod(f, g) == (f*f) % g
        assert raises(lambda: f.mulmod(f, "AAA"), TypeError)
        assert raises(lambda: f.mulmod("AAA", g), TypeError)

        # powmod
        assert f.powmod(2, g) == (f*f) % g
        assert raises(lambda: f.powmod(2, "AAA"), TypeError)

        # divmod
        S, T = f.divmod(g)
        assert S*g + T == f
        assert raises(lambda: f.divmod("AAA"), TypeError)
        assert raises(lambda: f_bad.divmod(f_bad), ValueError)

        # gcd
        assert raises(lambda: f_cmp.gcd(f_cmp), NotImplementedError)
        assert raises(lambda: f.gcd("f"), TypeError)

        # xgcd
        assert raises(lambda: (f_cmp).xgcd(f_cmp), ValueError)
        assert raises(lambda: (f).xgcd("f_cmp"), TypeError)

        # disc.
        assert raises(lambda: (f_cmp).discriminant(), NotImplementedError)

        # Radical
        assert raises(lambda: (f_cmp).radical(), NotImplementedError)

        # inverse_mod
        f_inv = f.inverse_mod(g)
        assert (f * f_inv) % g == 1
        assert raises(lambda: f.inverse_mod("AAA"), TypeError)
        assert raises(lambda: (f_cmp).inverse_mod(f_cmp), ValueError)

        f_inv = f.inverse_series_trunc(2)
        assert (f * f_inv) % R_test([0,0,1]) == 1
        assert raises(lambda: R_cmp([0,0,1]).inverse_series_trunc(2), ValueError)

        # Resultant
        f1 = R_test([-3, 1])
        f2 = R_test([-5, 1])
        assert f1.resultant(f2) == (3 - 5)
        assert raises(lambda: f.resultant("AAA"), TypeError)        

        # sqrt
        f1 = R_test.random_element(irreducible=True)
        assert raises(lambda: f1.sqrt(), ValueError)
        assert (f1*f1).sqrt() in [f1, -f1]

        # deflation
        f1 = R_test([1,0,2,0,3])
        assert raises(lambda: f1.deflate(100), ValueError)
        assert f1.deflate(2) == R_test([1,2,3])

        # factor
        ff = R_test([3,2,1]) * R_test([3,2,1]) * R_test([5,4,3])
        ff_rebuild = R_test.one()
        c, facs = ff.factor()
        ff_rebuild *= c
        for p, e in facs:
            assert p.is_irreducible()
            ff_rebuild *= p**e
        assert ff_rebuild == ff

        assert set(ff.factor()[1]) == set(ff.factor(algorithm="cantor_zassenhaus")[1])
        assert set(ff.factor()[1]) == set(ff.factor(algorithm="kaltofen_shoup")[1])
        assert set(ff.factor()[1]) == set(ff.factor(algorithm="berlekamp")[1])
        assert raises(lambda: R_test([0,0,1]).factor(algorithm="AAA"), ValueError)
        assert raises(lambda: R_test([0,0,1]).complex_roots(), DomainError)

        # composite moduli not supported
        assert raises(lambda: R_cmp([0,0,1]).factor(), NotImplementedError)
        assert raises(lambda: R_cmp([0,0,1]).factor_squarefree(), NotImplementedError)
        assert raises(lambda: R_cmp([0,0,1]).roots(), NotImplementedError)
        assert raises(lambda: R_cmp([0,0,1]).complex_roots(), DomainError)

        # minpoly
        assert raises(lambda: R_cmp.minpoly([1,2,3,4]), NotImplementedError)
        assert raises(lambda: R_test.minpoly(1), ValueError)
        assert raises(lambda: R_test.minpoly([1,2,3,"AAA"]), ValueError)

        # multipoint_evaluation
        assert raises(lambda: R_test([1,2,3]).multipoint_evaluate([1,2,3,"AAA"]), ValueError)
        assert raises(lambda: R_test([1,2,3]).multipoint_evaluate("AAA"), ValueError)

        f = R_test([1,2,3])
        l = [-1,-2,-3,-4,-5]
        assert [f(x) for x in l] == f.multipoint_evaluate(l)


def _all_polys():
    return [
        # (poly_type, scalar_type, is_field)
        (flint.fmpz_poly, flint.fmpz, False),
        (flint.fmpq_poly, flint.fmpq, True),
        (lambda *a: flint.nmod_poly(*a, 17), lambda x: flint.nmod(x, 17), True),
        (lambda *a: flint.fmpz_mod_poly(*a, flint.fmpz_mod_poly_ctx(163)), 
         lambda x: flint.fmpz_mod(x, flint.fmpz_mod_ctx(163)), 
         True),
        (lambda *a: flint.fmpz_mod_poly(*a, flint.fmpz_mod_poly_ctx(2**127 - 1)), 
         lambda x: flint.fmpz_mod(x, flint.fmpz_mod_ctx(2**127 - 1)), 
         True),
        (lambda *a: flint.fmpz_mod_poly(*a, flint.fmpz_mod_poly_ctx(2**255 - 19)), 
         lambda x: flint.fmpz_mod(x, flint.fmpz_mod_ctx(2**255 - 19)), 
         True),
    ]


def test_polys():
    for P, S, is_field in _all_polys():

        assert P([S(1)]) == P([1]) == P(P([1])) == P(1)

        assert raises(lambda: P([None]), TypeError)
        assert raises(lambda: P(object()), TypeError), f"{P(object()) = }"
        assert raises(lambda: P(None), TypeError)
        assert raises(lambda: P(None, None), TypeError)
        assert raises(lambda: P([1,2], None), TypeError)
        assert raises(lambda: P(1, None), TypeError)

        assert len(P([])) == P([]).length() == 0
        assert len(P([1])) == P([1]).length() == 1
        assert len(P([1,2])) == P([1,2]).length() == 2
        assert len(P([1,2,3])) == P([1,2,3]).length() == 3

        assert P([]).degree() == -1
        assert P([1]).degree() == 0
        assert P([1,2]).degree() == 1
        assert P([1,2,3]).degree() == 2

        assert (P([1]) == P([1])) is True
        assert (P([1]) != P([1])) is False
        assert (P([1]) == P([2])) is False
        assert (P([1]) != P([2])) is True

        assert (P([1]) == None) is False
        assert (P([1]) != None) is True
        assert (None == P([1])) is False
        assert (None != P([1])) is True

        assert raises(lambda: P([1]) < P([1]), TypeError)
        assert raises(lambda: P([1]) <= P([1]), TypeError)
        assert raises(lambda: P([1]) > P([1]), TypeError)
        assert raises(lambda: P([1]) >= P([1]), TypeError)
        assert raises(lambda: P([1]) < None, TypeError)
        assert raises(lambda: P([1]) <= None, TypeError)
        assert raises(lambda: P([1]) > None, TypeError)
        assert raises(lambda: P([1]) >= None, TypeError)
        assert raises(lambda: None < P([1]), TypeError)
        assert raises(lambda: None <= P([1]), TypeError)
        assert raises(lambda: None > P([1]), TypeError)
        assert raises(lambda: None >= P([1]), TypeError)

        assert P([1, 2, 3])[1] == S(2)
        assert P([1, 2, 3])[-1] == S(0)
        assert P([1, 2, 3])[3] == S(0)

        p = P([1, 2, 3])
        p[1] = S(4)
        assert p == P([1, 4, 3])

        def setbad(obj, i, val):
            obj[i] = val

        assert raises(lambda: setbad(p, 2, None), TypeError)
        assert raises(lambda: setbad(p, -1, 1), ValueError)

        for v in [], [1], [1, 2]:
            if P == flint.fmpz_poly:
                assert P(v).repr() == f'fmpz_poly({v!r})'
            elif P == flint.fmpq_poly:
                assert P(v).repr() == f'fmpq_poly({v!r})'
            elif P == flint.nmod_poly:
                assert P(v).repr() == f'nmod_poly({v!r}, 17)'

        assert repr(P([])) == '0'
        assert repr(P([1])) == '1'
        assert repr(P([1, 2])) == '2*x + 1'
        assert repr(P([1, 2, 3])) == '3*x^2 + 2*x + 1'

        p = P([1, 2, 3])
        assert p(0) == p(S(0)) == S(1) == 1
        assert p(1) == p(S(1)) == S(6) == 6
        assert p(p) == P([6, 16, 36, 36, 27])
        assert raises(lambda: p(None), TypeError)

        assert bool(P([])) is False
        assert bool(P([1])) is True

        assert +P([1, 2, 3]) == P([1, 2, 3])
        assert -P([1, 2, 3]) == P([-1, -2, -3])

        assert P([1, 2, 3]) + P([4, 5, 6]) == P([5, 7, 9])

        for T in [int, S, flint.fmpz]:
            assert P([1, 2, 3]) + T(1) == P([2, 2, 3])
            assert T(1) + P([1, 2, 3]) == P([2, 2, 3])

        assert raises(lambda: P([1, 2, 3]) + None, TypeError)
        assert raises(lambda: None + P([1, 2, 3]), TypeError)

        assert P([1, 2, 3]) - P([4, 5, 6]) == P([-3, -3, -3])

        for T in [int, S, flint.fmpz]:
            assert P([1, 2, 3]) - T(1) == P([0, 2, 3])
            assert T(1) - P([1, 2, 3]) == P([0, -2, -3])

        assert raises(lambda: P([1, 2, 3]) - None, TypeError)
        assert raises(lambda: None - P([1, 2, 3]), TypeError)

        assert P([1, 2, 3]) * P([4, 5, 6]) == P([4, 13, 28, 27, 18])

        for T in [int, S, flint.fmpz]:
            assert P([1, 2, 3]) * T(2) == P([2, 4, 6])
            assert T(2) * P([1, 2, 3]) == P([2, 4, 6])

        assert raises(lambda: P([1, 2, 3]) * None, TypeError)
        assert raises(lambda: None * P([1, 2, 3]), TypeError)

        assert P([1, 2, 1]) // P([1, 1]) == P([1, 1])
        assert P([1, 2, 1]) % P([1, 1]) == P([0])
        assert divmod(P([1, 2, 1]), P([1, 1])) == (P([1, 1]), P([0]))

        if is_field:
            assert P([1, 1]) // 2 == P([S(1)/2, S(1)/2])
            assert P([1, 1]) % 2 == P([0])
        else:
            assert P([1, 1]) // 2 == P([0, 0])
            assert P([1, 1]) % 2 == P([1, 1])

        assert 1 // P([1, 1]) == P([0])
        assert 1 % P([1, 1]) == P([1])
        assert divmod(1, P([1, 1])) == (P([0]), P([1]))

        assert raises(lambda: P([1, 2, 1]) // None, TypeError)
        assert raises(lambda: P([1, 2, 1]) % None, TypeError)
        assert raises(lambda: divmod(P([1, 2, 1]), None), TypeError)

        assert raises(lambda: None // P([1, 1]), TypeError)
        assert raises(lambda: None % P([1, 1]), TypeError)
        assert raises(lambda: divmod(None, P([1, 1])), TypeError)

        assert raises(lambda: P([1, 2, 1]) // 0, ZeroDivisionError)
        assert raises(lambda: P([1, 2, 1]) % 0, ZeroDivisionError)
        assert raises(lambda: divmod(P([1, 2, 1]), 0), ZeroDivisionError)

        assert raises(lambda: P([1, 2, 1]) // P([0]), ZeroDivisionError)
        assert raises(lambda: P([1, 2, 1]) % P([0]), ZeroDivisionError)
        assert raises(lambda: divmod(P([1, 2, 1]), P([0])), ZeroDivisionError)

        if is_field:
            assert P([2, 2]) / 2 == P([1, 1])
            assert P([1, 2]) / 2 == P([S(1)/2, 1])
            assert raises(lambda: P([1, 2]) / 0, ZeroDivisionError)
        else:
            assert raises(lambda: P([2, 2]) / 2, TypeError)

        assert raises(lambda: 1 / P([1, 1]), TypeError)
        assert raises(lambda: P([1, 2, 1]) / P([1, 1]), TypeError)
        assert raises(lambda: P([1, 2, 1]) / P([1, 2]), TypeError)

        assert P([1, 1]) ** 0 == P([1])
        assert P([1, 1]) ** 1 == P([1, 1])
        assert P([1, 1]) ** 2 == P([1, 2, 1])
        assert raises(lambda: P([1, 1]) ** -1, ValueError)
        assert raises(lambda: P([1, 1]) ** None, TypeError)
        
        # # XXX: Not sure what this should do in general:
        assert raises(lambda: pow(P([1, 1]), 2, 3), NotImplementedError)

        assert P([1, 2, 1]).gcd(P([1, 1])) == P([1, 1])
        assert raises(lambda: P([1, 2, 1]).gcd(None), TypeError)

        if is_field:
            p1 = P([1, 0, 1])
            p2 = P([2, 1])
            g, s, t = P([1]), P([1])/5, P([2, -1])/5
            assert p1.xgcd(p2) == (g, s, t)
            assert raises(lambda: p1.xgcd(None), TypeError)

        assert P([1, 2, 1]).factor() == (S(1), [(P([1, 1]), 2)])

        assert P([1, 2, 1]).sqrt() == P([1, 1])
        assert raises(lambda: P([1, 2, 2]).sqrt(), ValueError), f"{P}, {P([1, 2, 2]).sqrt()}"

        if P == flint.fmpq_poly:
            assert raises(lambda: P([1, 2, 1], 3).sqrt(), ValueError)
            assert P([1, 2, 1], 4).sqrt() == P([1, 1], 2)

        assert P([]).deflation() == (P([]), 1)
        assert P([1, 2]).deflation() == (P([1, 2]), 1)
        assert P([1, 0, 2]).deflation() == (P([1, 2]), 2)

        assert P([1, 2, 1]).derivative() == P([2, 2])

        if is_field:
            assert P([1, 2, 1]).integral() == P([0, 1, 1, S(1)/3])



all_tests = [
    test_pyflint,
    test_fmpz,
    test_fmpz_factor,
    test_fmpz_functions,
    test_fmpz_poly,
    test_fmpz_poly_factor,
    test_fmpz_poly_functions,
    test_fmpz_mat,
    test_fmpz_series,
    test_fmpq,
    test_fmpq_poly,
    test_fmpq_mat,
    test_fmpq_series,
    test_nmod,
    test_nmod_poly,
    test_nmod_mat,
    test_arb,
    test_fmpz_mod,
    test_fmpz_mod_dlog,
    test_fmpz_mod_poly,
    test_polys,
]
