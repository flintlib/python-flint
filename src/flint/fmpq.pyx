cdef any_as_fmpq(obj):
    if typecheck(obj, fmpq):
        return obj
    z = any_as_fmpz(obj)
    if z is NotImplemented:
        return z
    q = fmpq.__new__(fmpq)
    fmpz_set(fmpq_numref((<fmpq>q).val), (<fmpz>z).val)
    fmpz_one(fmpq_denref((<fmpq>q).val))
    return q

cdef class fmpq(flint_scalar):
    """
    The fmpq type represents multiprecision rational numbers.

        >>> fmpq(1,7) + fmpq(50,51)
        401/357

    """

    cdef fmpq_t val

    def __cinit__(self):
        fmpq_init(self.val)

    def __dealloc__(self):
        fmpq_clear(self.val)

    def __init__(self, p=None, q=None):
        cdef long x
        if q is None:
            if p is None:
                return # zero
            elif typecheck(p, fmpq):
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
                p = any_as_fmpq(p)
                if p is NotImplemented:
                    raise ValueError("cannot create fmpq from object of type %s" % type(p))
                fmpq_set(self.val, (<fmpq>p).val)
                return
        p = any_as_fmpz(p)
        if p is NotImplemented:
            raise ValueError("cannot create fmpq from object of type %s" % type(p))
        q = any_as_fmpz(q)
        if q is NotImplemented:
            raise ValueError("cannot create fmpq from object of type %s" % type(q))
        if fmpz_is_zero((<fmpz>q).val):
            raise ZeroDivisionError("cannot create rational number with zero denominator")
        fmpz_set(fmpq_numref(self.val), (<fmpz>p).val)
        fmpz_set(fmpq_denref(self.val), (<fmpz>q).val)
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
            if op == 0: res = (s-t).p < 0
            elif op == 1: res = (s-t).p <= 0
            elif op == 4: res = (s-t).p > 0
            elif op == 5: res = (s-t).p >= 0
            else: raise ValueError
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

    def __nonzero__(self):
        return not fmpq_is_zero(self.val)

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
        s = any_as_fmpq(s)
        if s is NotImplemented:
            return s
        t = any_as_fmpq(t)
        if t is NotImplemented:
            return t
        r = fmpq.__new__(fmpq)
        fmpq_add(r.val, (<fmpq>s).val, (<fmpq>t).val)
        return r

    def __sub__(s, t):
        cdef fmpq r
        s = any_as_fmpq(s)
        if s is NotImplemented:
            return s
        t = any_as_fmpq(t)
        if t is NotImplemented:
            return t
        r = fmpq.__new__(fmpq)
        fmpq_sub(r.val, (<fmpq>s).val, (<fmpq>t).val)
        return r

    def __mul__(s, t):
        cdef fmpq r
        s = any_as_fmpq(s)
        if s is NotImplemented:
            return s
        t = any_as_fmpq(t)
        if t is NotImplemented:
            return t
        r = fmpq.__new__(fmpq)
        fmpq_mul(r.val, (<fmpq>s).val, (<fmpq>t).val)
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

    def __div__(s, t):
        return fmpq._div_(s, t)

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

    def __hash__(self):
        from fractions import Fraction
        return hash(Fraction(int(self.p), int(self.q), _normalize=False))

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
        cdef fmpq v
        cdef long e
        assert z is None
        e = n
        if type(self) is fmpq:
            v = fmpq.__new__(fmpq)
            if e >= 0:
                fmpz_pow_ui(fmpq_numref(v.val), fmpq_numref((<fmpq>self).val), e)
                fmpz_pow_ui(fmpq_denref(v.val), fmpq_denref((<fmpq>self).val), e)
            else:
                if fmpq_is_zero((<fmpq>self).val):
                    raise ZeroDivisionError
                fmpz_pow_ui(fmpq_denref(v.val), fmpq_numref((<fmpq>self).val), -e)
                fmpz_pow_ui(fmpq_numref(v.val), fmpq_denref((<fmpq>self).val), -e)
            return v
        return NotImplemented

