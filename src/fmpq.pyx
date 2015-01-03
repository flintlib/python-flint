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
        fmpq(401,357)

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
            raise NotImplementedError("fmpq comparisons")

    def numer(self):
        cdef fmpz x = fmpz.__new__(fmpz)
        fmpz_set(x.val, fmpq_numref(self.val))
        return x

    def denom(self):
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

    def str(self):
        if self.q == 1:
            return str(self.p)
        else:
            return "%s/%s" % (self.p, self.q)

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

    def __div__(s, t):
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

    # __truediv__ = __div__ doesn't seem to work?
    def __truediv__(s, t):
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

    def next(s, bint signed=True, bint minimal=True):
        """
        Returns the next rational number after *s* as ordered by
        minimal height (if *minimal* is True) or following the Calkin-Wilf
        sequence (if *minimal* is False). If *signed* is set to False,
        only the nonnegative rational numbers are considered.

            >>> fmpq(23456789,98765432).next()
            fmpq(-23456789,98765432)
            >>> fmpq(23456789,98765432).next(signed=False)
            fmpq(98765432,23456789)
            >>> fmpq(23456789,98765432).next(signed=False, minimal=False)
            fmpq(98765432,75308643)
            >>> a = fmpq(0)
            >>> for i in range(20):
            ...     a = a.next(); print a,
            ... 
            1 -1 1/2 -1/2 2 -2 1/3 -1/3 3 -3 2/3 -2/3 3/2 -3/2 1/4 -1/4 4 -4 3/4 -3/4
            >>> a = fmpq(0)
            >>> for i in range(20):
            ...     a = a.next(signed=False); print a,
            ... 
            1 1/2 2 1/3 3 2/3 3/2 1/4 4 3/4 4/3 1/5 5 2/5 5/2 3/5 5/3 4/5 5/4 1/6
            >>> a = fmpq(0)
            >>> for i in range(20):
            ...     a = a.next(minimal=False); print a,
            ... 
            1 -1 1/2 -1/2 2 -2 1/3 -1/3 3/2 -3/2 2/3 -2/3 3 -3 1/4 -1/4 4/3 -4/3 3/5 -3/5
            >>> a = fmpq(0)
            >>> for i in range(20):
            ...     a = a.next(signed=False, minimal=False); print a,
            ... 
            1 1/2 2 1/3 3/2 2/3 3 1/4 4/3 3/5 5/2 2/5 5/3 3/4 4 1/5 5/4 4/7 7/3 3/8
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
    def bernoulli_ui(ulong n):
        """
        Returns the nth Bernoulli number B_n as an fmpq.

            >>> [fmpq.bernoulli_ui(n) for n in range(8)]
            [fmpq(1), fmpq(-1,2), fmpq(1,6), fmpq(0), fmpq(-1,30), fmpq(0), fmpq(1,42), fmpq(0)]
            >>> fmpq.bernoulli_ui(50)
            fmpq(495057205241079648212477525,66)
        """
        u = fmpq.__new__(fmpq)
        bernoulli_fmpq_ui((<fmpq>u).val, n)
        return u

    @staticmethod
    def harmonic_ui(ulong n):
        """
        Returns the harmonic number H_n as an fmpq.

            >>> [fmpq.harmonic_ui(n) for n in range(6)]
            [fmpq(0), fmpq(1), fmpq(3,2), fmpq(11,6), fmpq(25,12), fmpq(137,60)]
            >>> fmpq.harmonic_ui(50)
            fmpq(13943237577224054960759,3099044504245996706400)
        """
        cdef fmpq v = fmpq()
        fmpq_harmonic_ui(v.val, n)
        return v

