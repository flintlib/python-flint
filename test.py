import sys
import flint

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
    assert repr(flint.fmpz(0)) == "fmpz(0)"
    assert repr(flint.fmpz(-27)) == "fmpz(-27)"

def test_fmpz_poly():
    Z = flint.fmpz_poly
    assert Z() == Z([])
    assert Z() == Z([0])
    assert Z() == Z([0,flint.fmpz(0),0])
    assert Z() == Z([0,0L,0])
    assert Z() != Z([1])
    assert Z([1]) == Z([1L])
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
    c = M(2,2,[1,2,3,4])
    assert c.det() == -2
    assert raises(lambda: a + c, ValueError)
    assert (a * 3).entries() == [3,6,9,12,15,18]
    assert (3 * a).entries() == [3,6,9,12,15,18]
    assert (a * 3L).entries() == [3,6,9,12,15,18]
    assert (3L * a).entries() == [3,6,9,12,15,18]
    assert (a * flint.fmpz(3)).entries() == [3,6,9,12,15,18]
    assert (flint.fmpz(3) * a).entries() == [3,6,9,12,15,18]
    assert M.randrank(5,7,3,10).rank() == 3
    A = M.randbits(5,3,2)
    B = M.randtest(3,7,3)
    C = M.randtest(7,2,4)
    assert A*(B*C) == (A*B)*C

if __name__ == "__main__":
    sys.stdout.write("test_fmpz..."); test_fmpz(); print("OK")
    sys.stdout.write("test_fmpz_poly..."); test_fmpz_poly(); print("OK")
    sys.stdout.write("test_fmpz_mat..."); test_fmpz_mat(); print("OK")

