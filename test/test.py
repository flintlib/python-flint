import sys
import flint
import operator
import doctest

if sys.version_info[0] >= 3:
    long = int

ctx = flint.ctx

def raises(f, exception):
    try:
        f()
    except exception:
        return True
    return False

def test_fmpz():
    assert flint.fmpz() == flint.fmpz(0)
    L = [0, 1, 2, 3, 2**31-1, 2**31, 2**63-1, 2**63, 2**64-1, 2**64]
    L += [-x for x in L]
    for s in L:
        for t in L:
            for ltype in (flint.fmpz, int, long):
                for rtype in (flint.fmpz, int, long):
                    assert ltype(s) + rtype(t) == s + t
                    assert ltype(s) - rtype(t) == s - t
                    assert ltype(s) * rtype(t) == s * t
                    if t == 0:
                        assert raises(lambda: ltype(s) // rtype(t), ZeroDivisionError)
                        assert raises(lambda: ltype(s) % rtype(t), ZeroDivisionError)
                        assert raises(lambda: divmod(ltype(s), rtype(t)), ZeroDivisionError)
                    else:
                        assert ltype(s) // rtype(t) == s // t
                        assert ltype(s) % rtype(t) == s % t
                        assert divmod(ltype(s), rtype(t)) == divmod(s, t)
                    assert (ltype(s) == rtype(t)) == (s == t)
                    assert (ltype(s) != rtype(t)) == (s != t)
                    assert (ltype(s) < rtype(t)) == (s < t)
                    assert (ltype(s) <= rtype(t)) == (s <= t)
                    assert (ltype(s) > rtype(t)) == (s > t)
                    assert (ltype(s) >= rtype(t)) == (s >= t)
                    if 0 <= t < 10:
                        assert (ltype(s) ** rtype(t)) == (s ** t)
    assert +flint.fmpz(0) == 0
    assert +flint.fmpz(1) == 1
    assert +flint.fmpz(-1) == -1
    assert -flint.fmpz(0) == 0
    assert -flint.fmpz(1) == -1
    assert -flint.fmpz(-1) == 1
    assert abs(flint.fmpz(0)) == 0
    assert abs(flint.fmpz(1)) == 1
    assert abs(flint.fmpz(-1)) == 1
    assert int(flint.fmpz(2)) == 2
    assert isinstance(int(flint.fmpz(2)), int)
    assert long(flint.fmpz(2)) == 2
    assert isinstance(long(flint.fmpz(2)), long)
    ctx.pretty = False
    assert repr(flint.fmpz(0)) == "fmpz(0)"
    assert repr(flint.fmpz(-27)) == "fmpz(-27)"
    ctx.pretty = True
    assert repr(flint.fmpz(0)) == "0"
    assert repr(flint.fmpz(-27)) == "-27"
    assert bool(flint.fmpz(0)) == False
    assert bool(flint.fmpz(1)) == True

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
    assert +Z([1,2]) == Z([1,2])
    assert -Z([1,2]) == Z([-1,-2])
    assert raises(lambda: Z([1,2,3]) ** -1, (OverflowError, ValueError))
    assert raises(lambda: Z([1,2,3]) ** Z([1,2]), TypeError)
    assert raises(lambda: Z([1,2]) // Z([]), ZeroDivisionError)
    assert raises(lambda: Z([]) // Z([]), ZeroDivisionError)
    assert raises(lambda: Z([1,2]) % Z([]), ZeroDivisionError)
    assert raises(lambda: divmod(Z([1,2]), Z([])), ZeroDivisionError)
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
    p = Z([3,4,5])
    assert p(2) == 31
    assert p(flint.fmpq(2,3)) == flint.fmpq(71,9)
    assert p(Z([1,-1])) == Z([12,-14,5])
    assert p(flint.fmpq_poly([2,3],5)) == flint.fmpq_poly([27,24,9],5)
    assert p(flint.arb("1.1")).overlaps(flint.arb("13.45"))
    assert p(flint.acb("1.1", "1.2")).overlaps(flint.acb("6.25", "18.00"))

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
    assert A*(B*C) == (A*B)*C
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
    assert raises(lambda: M(1,2) ** 3, ValueError)
    assert raises(lambda: M(1,1) ** M(1,1), TypeError)
    assert raises(lambda: 1 ** M(1,1), TypeError)
    assert raises(lambda: M([1]), TypeError)
    assert raises(lambda: M([[1],[2,3]]), ValueError)
    assert M([[1,2,3],[4,5,6]]) == M(2,3,[1,2,3,4,5,6])

def test_fmpq():
    Q = flint.fmpq
    assert Q() == Q(0)
    assert Q(0) != Q(1)
    assert Q(0) == 0
    assert 0 == Q(0)
    assert Q(2) != 1
    assert 1 != Q(2)
    assert Q(1,2) != 1
    assert Q(2,3) == Q(flint.fmpz(2),long(3))
    assert Q(-2,-4) == Q(1,2)
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
    assert operator.truediv(Q(2,3), 5) == Q(2,15)
    ctx.pretty = False
    assert repr(Q(-2,3)) == "fmpq(-2,3)"
    ctx.pretty = True
    assert str(Q(-2,3)) == "-2/3"
    assert Q(2,3).p == Q(2,3).numer() == 2
    assert Q(2,3).q == Q(2,3).denom() == 3
    assert +Q(5,7) == Q(5,7)
    assert -Q(5,7) == Q(-5,7)
    assert -Q(-5,7) == Q(5,7)
    assert abs(Q(5,7)) == Q(5,7)
    assert abs(-Q(5,7)) == Q(5,7)
    assert raises(lambda: Q(1,0), ZeroDivisionError)
    assert raises(lambda: Q(1,2) / Q(0), ZeroDivisionError)
    assert raises(lambda: Q(1,2) / 0, ZeroDivisionError)

def test_fmpq_poly():
    Q = flint.fmpq_poly
    Z = flint.fmpz_poly
    assert Q() == Q([]) == Q([0]) == Q([0,0])
    assert Q() != Q([1])
    assert Q([1]) == Q([1])
    assert bool(Q()) == False
    assert bool(Q([1])) == True
    assert Q(Q([1,2])) == Q([1,2])
    assert Q(Z([1,2])) == Q([1,2])
    assert Q([1,2]) + 3 == Q([4,2])
    assert 3 + Q([1,2]) == Q([4,2])
    assert Q([1,2]) - 3 == Q([-2,2])
    assert 3 - Q([1,2]) == Q([2,-2])
    assert -Q([1,2]) == Q([-1,-2])
    assert Q([flint.fmpq(1,2),1]) * 2 == Q([1,2])
    assert Q([1,2]) == Z([1,2])
    assert Z([1,2]) == Q([1,2])
    assert Q([1,2]) != Z([3,2])
    assert Z([1,2]) != Q([3,2])
    assert Q([1,2,3])*Q([1,2]) == Q([1,4,7,6])
    assert Q([1,2,3])*Z([1,2]) == Q([1,4,7,6])
    assert Q([1,2,3]) * 3 == Q([3,6,9])
    assert 3 * Q([1,2,3]) == Q([3,6,9])
    assert Q([1,2,3]) * flint.fmpq(2,3) == (Q([1,2,3]) * 2) / 3
    assert flint.fmpq(2,3) * Q([1,2,3]) == (Q([1,2,3]) * 2) / 3
    assert raises(lambda: Q([1,2]) / Q([1,2]), TypeError)
    assert Q([1,2,3]) / flint.fmpq(2,3) == Q([1,2,3]) * flint.fmpq(3,2)
    assert Q([1,2,3]) ** 2 == Q([1,2,3]) * Q([1,2,3])
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
    assert repr(Q([15,20,10]) / 25) == "fmpq_poly([3, 4, 2], 5)"
    ctx.pretty = True
    assert str(Q([3,4,2],5)) == "2/5*x^2 + 4/5*x + 3/5"
    a = Q([2,2,3],4)
    assert a[2] == flint.fmpq(3,4)
    a[2] = 4
    assert a == Q([1,1,8],2)
    p = Q([3,4,5],7)
    assert p(2) == flint.fmpq(31,7)
    assert p(flint.fmpq(2,3)) == flint.fmpq(71,63)
    assert p(Z([1,-1])) == Q([12,-14,5],7)
    assert p(flint.fmpq_poly([2,3],5)) == flint.fmpq_poly([27,24,9],7*5)

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
    assert Q(1,2,[3,4]) * 2 == Q(1,2,[6,8])
    assert Q(1,2,[3,4]) * flint.fmpq(1,3) == Q(1,2,[1,flint.fmpq(4,3)])
    assert Q(1,2,[3,4]) * flint.fmpq(5,3) == Q(1,2,[5,flint.fmpq(20,3)])
    assert 2 * Q(1,2,[3,4]) == Q(1,2,[6,8])
    assert flint.fmpq(1,3) * Q(1,2,[3,4]) == Q(1,2,[1,flint.fmpq(4,3)])
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
    assert Q([[1,2,3],[4,5,6]]) == Q(2,3,[1,2,3,4,5,6])

def test_nmod():
    G = flint.nmod
    assert G(0,2) == G(2,2) == G(-2,2)
    assert G(1,2) != G(0,2)
    assert G(0,2) != G(0,3)
    assert G(3,5) == G(8,5)
    #assert G(3,5) == 8        # do we want this?
    #assert 8 == G(3,5)
    assert G(3,5) != 7
    assert 7 != G(3,5)
    assert G(-3,5) == -G(3,5) == G(2,5)
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
    assert G(3,17) * flint.fmpq(11,5) == G(10,17)
    assert G(3,17) / flint.fmpq(11,5) == G(6,17)
    assert raises(lambda: G(2,5) / G(0,5), ZeroDivisionError)
    assert raises(lambda: G(2,5) / 0, ZeroDivisionError)
    assert raises(lambda: G(2,5) + G(2,7), ValueError)
    assert raises(lambda: G(2,5) - G(2,7), ValueError)
    assert raises(lambda: G(2,5) * G(2,7), ValueError)
    assert raises(lambda: G(2,5) / G(2,7), ValueError)
    assert G(3,17).modulus() == 17

def test_nmod_poly():
    P = flint.nmod_poly
    Z = flint.fmpz_poly
    assert P([],17) == P([0],17)
    assert P([1,2,3],17) == P([1,2,3],17)
    assert P([1,2,3],17) != P([1,2,3],15)
    assert P([1,2,3],17) != P([1,2,4],15)
    assert P(Z([1,2,3]),17) == P([1,2,3],17)
    assert P([1,2,flint.nmod(3,17)],17) == P([1,2,3],17)
    assert raises(lambda: P([1,2,flint.nmod(3,15)],17), ValueError)
    assert P([1,2,3],17).degree() == 2
    assert P([1,2,3],17).length() == 3
    assert P([1,2,3],17) + 2 == P([3,2,3],17)
    assert 2 + P([1,2,3],17) == P([3,2,3],17)
    assert P([1,2,3],17) + P([3,4,5],17) == P([4,6,8],17)
    assert P([1,2,3],17) + P([3,4,5],17) == P([4,6,8],17)
    assert P([1,2,3],17) - 2 == P([16,2,3],17)
    assert 2 - P([1,2,3],17) == -P([16,2,3],17)
    assert P([1,2,3],17) - P([3,4,6],17) == P([15,15,14],17)
    assert P([1,2,3],17) * 2 == P([2,4,6],17)
    assert 2 * P([1,2,3],17) == P([2,4,6],17)
    assert P([1,2,3],17) * P([1,2,3],17) == P([1,4,10,12,9], 17)
    assert P([1,2,3],17) * Z([1,2,3]) == P([1,4,10,12,9], 17)
    assert Z([1,2,3]) * P([1,2,3],17) == P([1,4,10,12,9], 17)
    assert P([1,2,3,4,5],17) % P([2,3,4],17) == P([12,12],17)
    assert P([1,2,3,4,5],17) // P([2,3,4],17) == P([3,16,14],17)
    assert P([1,2,3,4,5],17) ** 2 == P([1,2,3,4,5],17) * P([1,2,3,4,5],17)
    assert P([1,2,3],17) * flint.nmod(3,17) == P([3,6,9],17)
    assert str(P([1,2,3],17)) == "3*x^2 + 2*x + 1"
    p = P([3,4,5],17)
    assert p(14) == flint.nmod(2,17)
    assert p(P([1,2,3],17)) == P([12,11,11,9,11],17)

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
    assert M(1,2,[3,4],17) / 3 == M(1,2,[3,4],17) * (~G(3,17))
    assert M(2,2,[1,2,3,4], 17).inv().det() == ~(M(2,2,[1,2,3,4], 17).det())
    assert M(2,2,[1,2,3,4], 17).inv().inv() == M(2,2,[1,2,3,4], 17)
    assert M(2,2,[0,1,2,3],17) * M(2, 2, [2,3,4,5], 17) == M(2,2,[4,5,16,4],17)
    assert raises(lambda: M([1], 5), TypeError)
    assert raises(lambda: M([[1],[2,3]], 5), ValueError)
    assert M([[1,2,3],[4,5,6]], 5) == M(2,3,[1,2,3,4,5,6], 5)

def test_arb():
    A = flint.arb
    assert A(3) > A(2.5)
    assert A(3) >= A("2.5")
    assert A(3) < A((3,1))
    assert A(3) <= A("inf")
    assert A(3) == A(3)
    assert A(3) != A(2)
    assert not (A("1.1") == A("1.1"))

if __name__ == "__main__":
    sys.stdout.write("test_fmpz..."); test_fmpz(); print("OK")
    sys.stdout.write("test_fmpz_poly..."); test_fmpz_poly(); print("OK")
    sys.stdout.write("test_fmpz_mat..."); test_fmpz_mat(); print("OK")
    sys.stdout.write("test_fmpq..."); test_fmpq(); print("OK")
    sys.stdout.write("test_fmpq_poly..."); test_fmpq_poly(); print("OK")
    sys.stdout.write("test_fmpq_mat..."); test_fmpq_mat(); print("OK")
    sys.stdout.write("test_nmod..."); test_nmod(); print("OK")
    sys.stdout.write("test_nmod_poly..."); test_nmod_poly(); print("OK")
    sys.stdout.write("test_nmod_mat..."); test_nmod_mat(); print("OK")
    sys.stdout.write("test_arb.."); test_arb(); print("OK")
    sys.stdout.write("doctests...");
    fail, total = doctest.testmod(flint);
    if fail == 0:
        print("OK")
    else:
        raise AssertionError("%i of %i doctests failed" % (fail, total))
