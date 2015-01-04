cdef inline int fmpz_set_pylong(fmpz_t x, obj):
    cdef int overflow
    cdef long longval
    longval = PyLong_AsLongAndOverflow(<PyObject*>obj, &overflow)
    if overflow:
        s = "%x" % obj
        fmpz_set_str(x, s, 16)
    else:
        fmpz_set_si(x, longval)

cdef inline int fmpz_set_python(fmpz_t x, obj):
    if PyInt_Check(<PyObject*>obj):
        fmpz_set_si(x, PyInt_AS_LONG(<PyObject*>obj))
        return 1
    if PyLong_Check(<PyObject*>obj):
        fmpz_set_pylong(x, obj)
        return 1
    return 0

cdef fmpz_get_intlong(fmpz_t x):
    """
    Convert fmpz_t to a Python int or long.
    """
    cdef char * s
    if COEFF_IS_MPZ(x[0]):
        s = fmpz_get_str(NULL, 16, x)
        v = int(s, 16)
        libc.stdlib.free(s)
        return v
    else:
        return <long>x[0]

cdef int fmpz_set_any_ref(fmpz_t x, obj):
    if typecheck(obj, fmpz):
        x[0] = (<fmpz>obj).val[0]
        return FMPZ_REF
    if PyInt_Check(<PyObject*>obj):
        fmpz_init(x)
        fmpz_set_si(x, PyInt_AS_LONG(<PyObject*>obj))
        return FMPZ_TMP
    if PyLong_Check(<PyObject*>obj):
        fmpz_init(x)
        s = "%x" % obj
        fmpz_set_str(x, s, 16)
        return FMPZ_TMP
    return FMPZ_UNKNOWN

cdef any_as_fmpz(obj):
    cdef fmpz_struct x[1]
    cdef bint xtype
    cdef fmpz v
    xtype = fmpz_set_any_ref(x, obj)
    if xtype == FMPZ_REF:
        v = fmpz.__new__(fmpz)
        fmpz_set(v.val, x)
        return v
    elif xtype == FMPZ_TMP:
        v = fmpz.__new__(fmpz)
        fmpz_clear(v.val)
        v.val[0] = x[0]
        return v
    else:
        return NotImplemented

cdef class fmpz(flint_scalar):
    """
    The fmpz type represents multiprecision integers.

        >>> fmpz(3) ** 25
        fmpz(847288609443)

    """

    cdef fmpz_t val

    def __cinit__(self):
        fmpz_init(self.val)

    def __dealloc__(self):
        fmpz_clear(self.val)

    def __init__(self, val=None):
        cdef long x
        if val is not None:
            if typecheck(val, fmpz):
                fmpz_set(self.val, (<fmpz>val).val)
            else:
                if fmpz_set_any_ref(self.val, val) == FMPZ_UNKNOWN: # XXX
                    if typecheck(val, str):
                        if fmpz_set_str(self.val, val, 10) != 0:
                            raise ValueError("invalid string for fmpz")
                        return
                    raise TypeError("cannot create fmpz from type %s" % type(val))

    # XXX: improve!
    def __int__(self):
        return fmpz_get_intlong(self.val)

    def __long__(self):
        return long(fmpz_get_intlong(self.val))

    def __index__(self):
        return fmpz_get_intlong(self.val)

    def __richcmp__(s, t, int op):
        cdef bint res = 0
        cdef long tl
        cdef fmpz_struct tval[1]
        cdef fmpz_struct * sval
        cdef int ttype
        sval = &((<fmpz>s).val[0])
        if PyInt_Check(<PyObject*>t):
            tl = PyInt_AS_LONG(<PyObject*>t)
            if   op == 2: res = fmpz_cmp_si(sval, tl) == 0
            elif op == 3: res = fmpz_cmp_si(sval, tl) != 0
            elif op == 0: res = fmpz_cmp_si(sval, tl) < 0
            elif op == 1: res = fmpz_cmp_si(sval, tl) <= 0
            elif op == 4: res = fmpz_cmp_si(sval, tl) > 0
            elif op == 5: res = fmpz_cmp_si(sval, tl) >= 0
        else:
            ttype = fmpz_set_any_ref(tval, t)
            if ttype != FMPZ_UNKNOWN:
                if   op == 2: res = fmpz_equal(sval, tval)
                elif op == 3: res = not fmpz_equal(sval, tval)
                elif op == 0: res = fmpz_cmp(sval, tval) < 0
                elif op == 1: res = fmpz_cmp(sval, tval) <= 0
                elif op == 4: res = fmpz_cmp(sval, tval) > 0
                elif op == 5: res = fmpz_cmp(sval, tval) >= 0
            if ttype == FMPZ_TMP:
                fmpz_clear(tval)
            if ttype == FMPZ_UNKNOWN:
                return NotImplemented
        return res

    def str(self):
        cdef char * s = fmpz_get_str(NULL, 10, self.val)
        try:
            res = str(s)
        finally:
            libc.stdlib.free(s)
        return res

    def repr(self):
        return "fmpz(%s)" % self.str()

    def __nonzero__(self):
        return not fmpz_is_zero(self.val)

    def __pos__(self):
        return self

    def __neg__(self):
        cdef fmpz res = fmpz.__new__(fmpz)
        fmpz_neg(res.val, self.val)
        return res

    def __abs__(self):
        if fmpz_sgn(self.val) >= 0:
            return self
        return -self

    def __add__(s, t):
        cdef fmpz_struct sval[1]
        cdef fmpz_struct tval[1]
        cdef int stype = FMPZ_UNKNOWN
        cdef int ttype = FMPZ_UNKNOWN
        u = NotImplemented
        stype = fmpz_set_any_ref(sval, s)
        if stype != FMPZ_UNKNOWN:
            ttype = fmpz_set_any_ref(tval, t)
            if ttype != FMPZ_UNKNOWN:
                u = fmpz.__new__(fmpz)
                fmpz_add((<fmpz>u).val, sval, tval)
        if stype == FMPZ_TMP: fmpz_clear(sval)
        if ttype == FMPZ_TMP: fmpz_clear(tval)
        return u

    def __sub__(s, t):
        cdef fmpz_struct sval[1]
        cdef fmpz_struct tval[1]
        cdef int stype = FMPZ_UNKNOWN
        cdef int ttype = FMPZ_UNKNOWN
        u = NotImplemented
        stype = fmpz_set_any_ref(sval, s)
        if stype != FMPZ_UNKNOWN:
            ttype = fmpz_set_any_ref(tval, t)
            if ttype != FMPZ_UNKNOWN:
                u = fmpz.__new__(fmpz)
                fmpz_sub((<fmpz>u).val, sval, tval)
        if stype == FMPZ_TMP: fmpz_clear(sval)
        if ttype == FMPZ_TMP: fmpz_clear(tval)
        return u

    def __mul__(s, t):
        cdef fmpz_struct sval[1]
        cdef fmpz_struct tval[1]
        cdef int stype = FMPZ_UNKNOWN
        cdef int ttype = FMPZ_UNKNOWN
        u = NotImplemented
        stype = fmpz_set_any_ref(sval, s)
        if stype != FMPZ_UNKNOWN:
            ttype = fmpz_set_any_ref(tval, t)
            if ttype != FMPZ_UNKNOWN:
                u = fmpz.__new__(fmpz)
                fmpz_mul((<fmpz>u).val, sval, tval)
        if stype == FMPZ_TMP: fmpz_clear(sval)
        if ttype == FMPZ_TMP: fmpz_clear(tval)
        return u

    def __floordiv__(s, t):
        cdef fmpz_struct sval[1]
        cdef fmpz_struct tval[1]
        cdef int stype = FMPZ_UNKNOWN
        cdef int ttype = FMPZ_UNKNOWN
        u = NotImplemented
        ttype = fmpz_set_any_ref(tval, t)
        if ttype != FMPZ_UNKNOWN:
            if fmpz_is_zero(tval):
                if ttype == FMPZ_TMP:
                    fmpz_clear(tval)
                raise ZeroDivisionError("fmpz division by zero")
            stype = fmpz_set_any_ref(sval, s)
            if stype != FMPZ_UNKNOWN:
                u = fmpz.__new__(fmpz)
                fmpz_fdiv_q((<fmpz>u).val, sval, tval)
        if stype == FMPZ_TMP: fmpz_clear(sval)
        if ttype == FMPZ_TMP: fmpz_clear(tval)
        return u

    def __mod__(s, t):
        cdef fmpz_struct sval[1]
        cdef fmpz_struct tval[1]
        cdef int stype = FMPZ_UNKNOWN
        cdef int ttype = FMPZ_UNKNOWN
        u = NotImplemented
        ttype = fmpz_set_any_ref(tval, t)
        if ttype != FMPZ_UNKNOWN:
            if fmpz_is_zero(tval):
                if ttype == FMPZ_TMP:
                    fmpz_clear(tval)
                raise ZeroDivisionError("fmpz division by zero")
            stype = fmpz_set_any_ref(sval, s)
            if stype != FMPZ_UNKNOWN:
                u = fmpz.__new__(fmpz)
                fmpz_fdiv_r((<fmpz>u).val, sval, tval)
        if stype == FMPZ_TMP: fmpz_clear(sval)
        if ttype == FMPZ_TMP: fmpz_clear(tval)
        return u

    def __divmod__(s, t):
        cdef fmpz_struct sval[1]
        cdef fmpz_struct tval[1]
        cdef int stype = FMPZ_UNKNOWN
        cdef int ttype = FMPZ_UNKNOWN
        u = NotImplemented
        ttype = fmpz_set_any_ref(tval, t)
        if ttype != FMPZ_UNKNOWN:
            if fmpz_is_zero(tval):
                if ttype == FMPZ_TMP:
                    fmpz_clear(tval)
                raise ZeroDivisionError("fmpz division by zero")
            stype = fmpz_set_any_ref(sval, s)
            if stype != FMPZ_UNKNOWN:
                u1 = fmpz.__new__(fmpz)
                u2 = fmpz.__new__(fmpz)
                fmpz_fdiv_qr((<fmpz>u1).val, (<fmpz>u2).val, sval, tval)
                u = u1, u2
        if stype == FMPZ_TMP: fmpz_clear(sval)
        if ttype == FMPZ_TMP: fmpz_clear(tval)
        return u

    def __pow__(s, t, m):
        cdef fmpz_struct sval[1]
        cdef int stype = FMPZ_UNKNOWN
        cdef ulong exp
        u = NotImplemented
        if m is not None:
            raise NotImplementedError("modular exponentiation")
        stype = fmpz_set_any_ref(sval, s)
        if stype != FMPZ_UNKNOWN:
            c = t
            u = fmpz.__new__(fmpz)
            fmpz_pow_ui((<fmpz>u).val, sval, c)
        if stype == FMPZ_TMP: fmpz_clear(sval)
        return u

    def gcd(self, other):
        """
        Returns the greatest common divisor of self and other.

            >>> fmpz(30).gcd(45)
            fmpz(15)
        """
        cdef fmpz_struct tval[1]
        cdef int ttype = FMPZ_UNKNOWN
        ttype = fmpz_set_any_ref(tval, other)
        if ttype == FMPZ_UNKNOWN:
            raise TypeError("input must be an integer")
        u = fmpz.__new__(fmpz)
        fmpz_gcd((<fmpz>u).val, self.val, tval)
        if ttype == FMPZ_TMP:
            fmpz_clear(tval)
        return u

    def factor(self):
        """
        Factors self into pseudoprimes, returning a list of
        (prime, exp) pairs. The sign is ignored.

            >>> fmpz(5040).factor()
            [(fmpz(2), 4), (fmpz(3), 2), (fmpz(5), 1), (fmpz(7), 1)]

        Warning: this is a preliminary implementation and should only
        be used for integers that fit in a single word, or which are
        larger but only have very small prime factors. It will
        crash if trying to factor a large integer with large
        prime factors.
        """
        cdef fmpz_factor_t fac
        cdef int i
        fmpz_factor_init(fac)
        fmpz_factor(fac, self.val)
        res = [0] * fac.num
        for 0 <= i < fac.num:
            u = fmpz.__new__(fmpz)
            fmpz_set((<fmpz>u).val, &fac.p[i])
            exp = <long> fac.exp[i]
            res[i] = (u, exp)
        fmpz_factor_clear(fac)
        return res

    def number_of_partitions(n):
        """
        Returns `p(n)`, the number of partitions of `n`, as an *fmpz*.

            >>> [int(fmpz(n).number_of_partitions()) for n in range(8)]
            [1, 1, 2, 3, 5, 7, 11, 15]
            >>> fmpz(100).number_of_partitions()
            fmpz(190569292)
            >>> len(str(fmpz(10**9).number_of_partitions()))
            35219

        Warning: the partition function grows rapidly.
        On a 32-bit system, `n` must not be larger than about `10^{16}`.
        On a 64-bit system, `n` must not be larger than about `10^{20}`.
        For large `n`, this function benefits from setting ``ctx.threads = 2``
        on multicore systems.
        """
        cdef fmpz v = fmpz()
        partitions_fmpz_fmpz(v.val, n.val, 0)
        return v

    def moebius_mu(s):
        """
        Returns the Moebius function `\mu(n)` as an *fmpz*.

            >>> [int(fmpz(n).moebius_mu()) for n in range(10)]
            [0, 1, -1, -1, 0, -1, 1, -1, 0, 0]
        """
        cdef fmpz v = fmpz()
        fmpz_set_si(v.val, fmpz_moebius_mu(s.val))
        return v


