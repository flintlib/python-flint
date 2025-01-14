from flint.flint_base.flint_base cimport flint_scalar
from flint.utils.typecheck cimport typecheck
from flint.types.fmpz cimport fmpz_set_any_ref
from flint.types.fmpz cimport fmpz
from flint.types.fmpz cimport any_as_fmpz

from flint.flintlib.types.flint cimport FMPZ_UNKNOWN, FMPZ_TMP, FMPZ_REF
from flint.flintlib.functions.fmpz cimport fmpz_set, fmpz_one
from flint.flintlib.functions.fmpz cimport fmpz_is_zero, fmpz_sgn
from flint.flintlib.functions.fmpz cimport fmpz_fdiv_q, fmpz_bits
from flint.flintlib.functions.fmpz cimport fmpz_cdiv_q
from flint.flintlib.functions.fmpz cimport fmpz_tdiv_q
from flint.flintlib.functions.fmpz cimport fmpz_clear
from flint.flintlib.functions.fmpq cimport *
from flint.flintlib.functions.bernoulli cimport *

cdef int fmpq_set_any_ref(fmpq_t x, obj):
    cdef int status
    fmpq_init(x)
    if typecheck(obj, fmpq):
        x[0] = (<fmpq>obj).val[0]
        return FMPZ_REF
    if typecheck(obj, fmpz):
        fmpz_set(fmpq_numref(x), (<fmpz>obj).val)
        fmpz_one(fmpq_denref(x))
        return FMPZ_TMP
    status = fmpz_set_any_ref(fmpq_numref(x), obj)
    if status != FMPZ_UNKNOWN:
        fmpz_one(fmpq_denref(x))
        return FMPZ_TMP
    fmpq_clear(x)
    return FMPZ_UNKNOWN

cdef any_as_fmpq(obj):
    cdef fmpq_t x
    cdef int status
    cdef fmpq q
    status = fmpq_set_any_ref(x, obj)
    if status == FMPZ_REF:
        q = fmpq.__new__(fmpq)
        fmpq_set(q.val, x)
        return q
    elif status == FMPZ_TMP:
        q = fmpq.__new__(fmpq)
        fmpq_clear(q.val)
        q.val[0] = x[0]
        return q
    else:
        return NotImplemented

cdef class fmpq(flint_scalar):
    """
    The fmpq type represents multiprecision rational numbers.

        >>> fmpq(1,7) + fmpq(50,51)
        401/357

    """

#    cdef fmpq_t val

    def __cinit__(self):
        fmpq_init(self.val)

    def __dealloc__(self):
        fmpq_clear(self.val)

    def __init__(self, *args):
        if not args:
            return  # zero
        elif len(args) == 2:
            p, q = args
        elif len(args) == 1:
            p = args[0]
            if typecheck(p, fmpq):
                fmpq_set(self.val, (<fmpq>p).val)
                return
            elif typecheck(p, str):
                if "/" in p:
                    p, q = p.split("/")
                    p = fmpz(p)
                    q = fmpz(q)
                else:
                    p = fmpz(p)
                    q = fmpz(1)
            else:
                p2 = any_as_fmpq(p)
                if p2 is NotImplemented:
                    raise TypeError("cannot create fmpq from object of type %s" % type(p))
                fmpq_set(self.val, (<fmpq>p2).val)
                return
        else:
            raise TypeError("fmpq() takes at most 2 arguments (%d given)" % len(args))

        p2 = any_as_fmpz(p)
        if p2 is NotImplemented:
            raise TypeError("cannot create fmpq from object of type %s" % type(p))
        q2 = any_as_fmpz(q)
        if q2 is NotImplemented:
            raise TypeError("cannot create fmpq from object of type %s" % type(q))
        if fmpz_is_zero((<fmpz>q2).val):
            raise ZeroDivisionError("cannot create rational number with zero denominator")

        fmpz_set(fmpq_numref(self.val), (<fmpz>p2).val)
        fmpz_set(fmpq_denref(self.val), (<fmpz>q2).val)
        fmpq_canonicalise(self.val)

    def __richcmp__(s, t, int op):
        cdef bint res
        s = any_as_fmpq(s)
        if s is NotImplemented:
            return s
        t = any_as_fmpq(t)
        if t is NotImplemented:
            return t
        if op == 2 or op == 3:
            res = fmpq_equal((<fmpq>s).val, (<fmpq>t).val)
            if op == 3:
                res = not res
            return res
        else:
            # todo: use fmpq_cmp when available
            if op == 0:
                res = (s-t).p < 0
            elif op == 1:
                res = (s-t).p <= 0
            elif op == 4:
                res = (s-t).p > 0
            elif op == 5:
                res = (s-t).p >= 0
            else:
                raise ValueError
            return res

    def numer(self):
        """
        Returns the numerator of *self* as an *fmpz*.
        """
        cdef fmpz x = fmpz.__new__(fmpz)
        fmpz_set(x.val, fmpq_numref(self.val))
        return x

    def denom(self):
        """
        Returns the denominator of *self* as an *fmpz*.
        """
        cdef fmpz x = fmpz.__new__(fmpz)
        fmpz_set(x.val, fmpq_denref(self.val))
        return x

    p = property(numer)
    q = property(denom)

    # These are the property names in the numeric tower.
    numerator = property(numer)
    denominator = property(denom)

    def __reduce__(self):
        return (fmpq, (int(self.p), int(self.q)))

    def repr(self):
        if self.q == 1:
            return "fmpq(%s)" % self.p
        else:
            return "fmpq(%s,%s)" % (self.p, self.q)

    def str(self, **kwargs):
        """
        Converts *self* to a string, forwarding optional keyword arguments
        to :meth:`.fmpz.str`.

            >>> fmpq.bernoulli(12).str()
            '-691/2730'
            >>> fmpq.bernoulli(100).str(base=2, condense=10)
            '-110001110{...257 digits...}0011011111/1000001000110010'
        """
        if self.q == 1:
            return self.p.str(**kwargs)
        else:
            return "%s/%s" % (self.p.str(**kwargs), self.q.str(**kwargs))

    def __int__(self):
        return int(self.trunc())

    def __floor__(self):
        return self.floor()

    def __ceil__(self):
        return self.ceil()

    def __trunc__(self):
        return self.trunc()

    def __bool__(self):
        return not fmpq_is_zero(self.val)

    def __round__(self, ndigits=None):
        return self.round(ndigits)

    def __pos__(self):
        return self

    def __neg__(self):
        cdef fmpq r = fmpq.__new__(fmpq)
        fmpq_neg(r.val, self.val)
        return r

    def __abs__(self):
        cdef fmpq r
        if fmpz_sgn(fmpq_numref(self.val)) >= 0:
            return self
        r = fmpq.__new__(fmpq)
        fmpq_neg(r.val, self.val)
        return r

    def __add__(s, t):
        cdef fmpq r
        t = any_as_fmpq(t)
        if t is NotImplemented:
            return t
        r = fmpq.__new__(fmpq)
        fmpq_add(r.val, (<fmpq>s).val, (<fmpq>t).val)
        return r

    def __radd__(s, t):
        cdef fmpq r
        t = any_as_fmpq(t)
        if t is NotImplemented:
            return t
        r = fmpq.__new__(fmpq)
        fmpq_add(r.val, (<fmpq>t).val, (<fmpq>s).val)
        return r

    def __sub__(s, t):
        cdef fmpq r
        t = any_as_fmpq(t)
        if t is NotImplemented:
            return t
        r = fmpq.__new__(fmpq)
        fmpq_sub(r.val, (<fmpq>s).val, (<fmpq>t).val)
        return r

    def __rsub__(s, t):
        cdef fmpq r
        t = any_as_fmpq(t)
        if t is NotImplemented:
            return t
        r = fmpq.__new__(fmpq)
        fmpq_sub(r.val, (<fmpq>t).val, (<fmpq>s).val)
        return r

    def __mul__(s, t):
        cdef fmpq r
        t = any_as_fmpq(t)
        if t is NotImplemented:
            return t
        r = fmpq.__new__(fmpq)
        fmpq_mul(r.val, (<fmpq>s).val, (<fmpq>t).val)
        return r

    def __rmul__(s, t):
        cdef fmpq r
        t = any_as_fmpq(t)
        if t is NotImplemented:
            return t
        r = fmpq.__new__(fmpq)
        fmpq_mul(r.val, (<fmpq>t).val, (<fmpq>s).val)
        return r

    @staticmethod
    def _div_(s, t):
        cdef fmpq r
        s = any_as_fmpq(s)
        if s is NotImplemented:
            return s
        t = any_as_fmpq(t)
        if t is NotImplemented:
            return t
        if fmpq_is_zero((<fmpq>t).val):
            raise ZeroDivisionError("fmpq division by zero")
        r = fmpq.__new__(fmpq)
        fmpq_div(r.val, (<fmpq>s).val, (<fmpq>t).val)
        return r

    def __truediv__(s, t):
        return fmpq._div_(s, t)

    def __rtruediv__(s, t):
        return fmpq._div_(t, s)

    def gcd(s, t):
        """GCD of two rational numbers.

            >>> fmpq(1,2).gcd(fmpq(3,4))
            1/4

        The GCD is defined as the GCD of the numerators divided by the LCM of
        the denominators. This is consistent with ``fmpz.gcd()`` but not with
        ``fmpq_poly.gcd()``.
        """
        cdef fmpq r
        t = any_as_fmpq(t)
        if t is NotImplemented:
            raise TypeError("fmpq expected")
        r = fmpq.__new__(fmpq)
        fmpq_gcd(r.val, (<fmpq>s).val, (<fmpq>t).val)
        return r

    def next(s, bint signed=True, bint minimal=True):
        """
        Returns the next rational number after *s* as ordered by
        minimal height (if *minimal* is True) or following the Calkin-Wilf
        sequence (if *minimal* is False). If *signed* is set to False,
        only the nonnegative rational numbers are considered.

            >>> fmpq(23456789,98765432).next()
            -23456789/98765432
            >>> fmpq(23456789,98765432).next(signed=False)
            98765432/23456789
            >>> fmpq(23456789,98765432).next(signed=False, minimal=False)
            98765432/75308643
            >>> a, b, c, d = [fmpq(0)], [fmpq(0)], [fmpq(0)], [fmpq(0)]
            >>> for i in range(20):
            ...     a.append(a[-1].next())
            ...     b.append(b[-1].next(signed=False))
            ...     c.append(c[-1].next(minimal=False))
            ...     d.append(d[-1].next(signed=False, minimal=False))
            ...
            >>> a
            [0, 1, -1, 1/2, -1/2, 2, -2, 1/3, -1/3, 3, -3, 2/3, -2/3, 3/2, -3/2, 1/4, -1/4, 4, -4, 3/4, -3/4]
            >>> b
            [0, 1, 1/2, 2, 1/3, 3, 2/3, 3/2, 1/4, 4, 3/4, 4/3, 1/5, 5, 2/5, 5/2, 3/5, 5/3, 4/5, 5/4, 1/6]
            >>> c
            [0, 1, -1, 1/2, -1/2, 2, -2, 1/3, -1/3, 3/2, -3/2, 2/3, -2/3, 3, -3, 1/4, -1/4, 4/3, -4/3, 3/5, -3/5]
            >>> d
            [0, 1, 1/2, 2, 1/3, 3/2, 2/3, 3, 1/4, 4/3, 3/5, 5/2, 2/5, 5/3, 3/4, 4, 1/5, 5/4, 4/7, 7/3, 3/8]

        """
        u = fmpq.__new__(fmpq)
        if signed:
            if minimal:
                fmpq_next_signed_minimal((<fmpq>u).val, (<fmpq>s).val)
            else:
                fmpq_next_signed_calkin_wilf((<fmpq>u).val, (<fmpq>s).val)
        else:
            if fmpz_sgn(fmpq_numref(s.val)) < 0:
                raise ValueError("s must be nonnegative")
            if minimal:
                fmpq_next_minimal((<fmpq>u).val, (<fmpq>s).val)
            else:
                fmpq_next_calkin_wilf((<fmpq>u).val, (<fmpq>s).val)
        return u

    @staticmethod
    def bernoulli(ulong n, bint cache=False):
        """
        Returns the Bernoulli number `B_n` as an *fmpq*.

            >>> [fmpq.bernoulli(n) for n in range(8)]
            [1, -1/2, 1/6, 0, -1/30, 0, 1/42, 0]
            >>> fmpq.bernoulli(50)
            495057205241079648212477525/66

        If *cache* is set to *True*, all the Bernoulli numbers up to *n* are
        computed and cached for fast subsequent retrieval. This feature should
        be used with caution if *n* is large. Calling
        :func:`ctx.cleanup()` frees cached Bernoulli numbers.
        """
        if cache:
            assert n <= 1000000
            bernoulli_cache_compute(n+1)
        u = fmpq.__new__(fmpq)
        bernoulli_fmpq_ui((<fmpq>u).val, n)
        return u

    @staticmethod
    def harmonic(ulong n):
        """
        Returns the harmonic number `H_n` as an *fmpq*.

            >>> [fmpq.harmonic(n) for n in range(6)]
            [0, 1, 3/2, 11/6, 25/12, 137/60]
            >>> fmpq.harmonic(50)
            13943237577224054960759/3099044504245996706400
        """
        cdef fmpq v = fmpq()
        fmpq_harmonic_ui(v.val, n)
        return v

    @staticmethod
    def dedekind_sum(n, k):
        """
        Dedekind sum.

            >>> fmpq.dedekind_sum(10, 3)
            1/18
        """
        cdef fmpz nv, kv
        cdef fmpq v
        nv = fmpz(n)
        kv = fmpz(k)
        v = fmpq()
        fmpq_dedekind_sum(v.val, nv.val, kv.val)
        return v

    def floor(self):
        """
        Floor function.

            >>> fmpq(3,2).floor()
            1
        """
        cdef fmpz r = fmpz.__new__(fmpz)
        fmpz_fdiv_q(r.val, fmpq_numref(self.val), fmpq_denref(self.val))
        return r

    def ceil(self):
        """
        Ceiling function.

            >>> fmpq(3,2).ceil()
            2
        """
        cdef fmpz r = fmpz.__new__(fmpz)
        fmpz_cdiv_q(r.val, fmpq_numref(self.val), fmpq_denref(self.val))
        return r

    def trunc(self):
        """
        Truncation function.

            >>> fmpq(3,2).trunc()
            1
            >>> fmpq(-3,2).trunc()
            -1
        """
        cdef fmpz r = fmpz.__new__(fmpz)
        fmpz_tdiv_q(r.val, fmpq_numref(self.val), fmpq_denref(self.val))
        return r

    def round(self, ndigits=None):
        """
        Rounding function.

            >>> fmpq(3,2).round()
            2
            >>> fmpq(-3,2).round()
            -2
        """
        from fractions import Fraction
        fself = Fraction(int(self.p), int(self.q))
        if ndigits is not None:
            fround = round(fself, ndigits)
            return fmpq(fround.numerator, fround.denominator)
        else:
            fround = round(fself)
            return fmpz(fround)

    def __hash__(self):
        import sys
        from fractions import Fraction
        if sys.version_info < (3, 12):
            return hash(Fraction(int(self.p), int(self.q), _normalize=False))
        else:
            return hash(Fraction._from_coprime_ints(int(self.p), int(self.q)))

    def height_bits(self, bint signed=False):
        """
        Returns the bit length of the maximum of the numerator and denominator.
        With signed=True, returns the negative value if the number is
        negative.

            >>> fmpq(1001,5).height_bits()
            10
            >>> fmpq(-5,1001).height_bits(signed=True)
            -10
        """
        cdef long b1, b2
        b1 = fmpz_bits(fmpq_numref(self.val))
        b2 = fmpz_bits(fmpq_denref(self.val))
        if signed and fmpz_sgn(fmpq_numref(self.val)) < 0:
            return -max(b1, b2)
        else:
            return max(b1, b2)

    def __pow__(self, n, z):
        cdef fmpz_struct nval[1]
        cdef int ntype = FMPZ_UNKNOWN
        cdef fmpq v
        cdef int success

        assert z is None

        ntype = fmpz_set_any_ref(nval, n)
        if ntype == FMPZ_UNKNOWN:
            return NotImplemented

        if fmpq_is_zero((<fmpq>self).val) and fmpz_sgn(nval) == -1:
            if ntype == FMPZ_TMP:
                fmpz_clear(nval)
            raise ZeroDivisionError

        v = fmpq.__new__(fmpq)
        success = fmpq_pow_fmpz(v.val, (<fmpq>self).val, nval)

        if ntype == FMPZ_TMP:
            fmpz_clear(nval)

        if success:
            return v
        else:
            raise OverflowError("fmpq_pow_fmpz(): exponent too large")

    def sqrt(self):
        """
        Return exact rational square root of self or raise an error.

            >>> fmpq(9, 4).sqrt()
            3/2
            >>> fmpq(8).sqrt()
            Traceback (most recent call last):
                ...
            flint.utils.flint_exceptions.DomainError: not a square number

        """
        return fmpq(self.numer().sqrt(), self.denom().sqrt())
