cdef inline int fmpz_set_pylong(fmpz_t x, obj):
    cdef int overflow
    cdef slong longval
    longval = pylong_as_slong(<PyObject*>obj, &overflow)
    if overflow:
        s = "%x" % obj
        fmpz_set_str(x, chars_from_str(s), 16)
    else:
        fmpz_set_si(x, longval)

cdef inline int fmpz_set_python(fmpz_t x, obj):
    if PY_MAJOR_VERSION < 3 and PyInt_Check(<PyObject*>obj):
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
        v = int(str_from_chars(s), 16)
        libc.stdlib.free(s)
        return v
    else:
        return <slong>x[0]

cdef int fmpz_set_any_ref(fmpz_t x, obj):
    if typecheck(obj, fmpz):
        x[0] = (<fmpz>obj).val[0]
        return FMPZ_REF
    if PY_MAJOR_VERSION < 3 and PyInt_Check(<PyObject*>obj):
        fmpz_init(x)
        fmpz_set_si(x, PyInt_AS_LONG(<PyObject*>obj))
        return FMPZ_TMP
    if PyLong_Check(<PyObject*>obj):
        fmpz_init(x)
        fmpz_set_pylong(x, obj)
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
    The *fmpz* type represents an arbitrary-size integer.

        >>> fmpz(3) ** 25
        847288609443

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
                        if fmpz_set_str(self.val, chars_from_str(val), 10) != 0:
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
        if PY_MAJOR_VERSION < 3 and PyInt_Check(<PyObject*>t):
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

    def bit_length(self):
        return fmpz_bits(self.val)

    def str(self, int base=10, long condense=0):
        """
        Converts *self* to a string, optionally in a non-decimal base
        and optionally showing only leading and trailing digits.

            >>> (fmpz(3) ** 100).str()
            '515377520732011331036461129765621272702107522001'
            >>> (fmpz(3) ** 100).str(base=3, condense=10)
            '1000000000{...81 digits...}0000000000'
        """
        assert 2 <= base <= 36
        cdef char * s = fmpz_get_str(NULL, base, self.val)
        try:
            res = str_from_chars(s)
        finally:
            libc.stdlib.free(s)
        if condense > 0:
            res = _str_trunc(res, condense)
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
            15
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

    def lcm(self, other):
        """
        Returns the greatest common divisor of self and other.

            >>> fmpz(30).gcd(45)
            15
        """
        cdef fmpz_struct tval[1]
        cdef int ttype = FMPZ_UNKNOWN
        ttype = fmpz_set_any_ref(tval, other)
        if ttype == FMPZ_UNKNOWN:
            raise TypeError("input must be an integer")
        u = fmpz.__new__(fmpz)
        fmpz_lcm((<fmpz>u).val, self.val, tval)
        if ttype == FMPZ_TMP:
            fmpz_clear(tval)
        return u

    def factor(self, trial_limit=None):
        """
        Factors self into prime numbers, returning a list of
        (prime, exp) pairs. The sign is ignored.

            >>> fmpz(5040).factor()
            [(2, 4), (3, 2), (5, 1), (7, 1)]
            >>> fmpz(10**35 + 1).factor()
            [(11, 1), (9091, 1), (909091, 1), (4147571, 1), (265212793249617641, 1)]
            >>> len(fmpz.fac_ui(1000).factor())
            168

        Warning: factoring large integers can be slow unless all
        prime factors are small.

        If *trial_limit* is set, perform trial division with at most
        this many primes, returning an incomplete factorization in which
        the largest factor may be composite. Factors largers than the trial
        division limit may still be found if it is cheap to do so, but no
        expensive algorithms will be run.

            >>> fmpz(2**128+10).factor()
            [(2, 1), (7, 1), (23, 1), (677, 1), (2957, 1), (1042733, 1), (506256324715258822390969, 1)]
            >>> fmpz(2**128+10).factor(trial_limit=200)
            [(2, 1), (7, 1), (23, 1), (677, 1), (1560971251139657345905734136865089, 1)]
        """
        cdef fmpz_factor_t fac
        cdef fmpz_t tmp
        cdef int i
        fmpz_factor_init(fac)
        if trial_limit is not None:
            fmpz_factor_trial_range(fac, self.val, 0, trial_limit)
            fmpz_init(tmp)
            fmpz_factor_expand(tmp, fac)
            fmpz_divexact(tmp, self.val, tmp)
            fmpz_abs(tmp, tmp)
            if not fmpz_is_one(tmp):
                _fmpz_factor_append(fac, tmp, 1)
            fmpz_clear(tmp)
        else:
            fmpz_factor(fac, self.val)
        res = [0] * fac.num
        for 0 <= i < fac.num:
            u = fmpz.__new__(fmpz)
            fmpz_set((<fmpz>u).val, &fac.p[i])
            exp = <long> fac.exp[i]
            res[i] = (u, exp)
        fmpz_factor_clear(fac)
        return res

    def factor_smooth(self, bits=15, int proved=-1):
        cdef fmpz_factor_t fac
        cdef int i
        fmpz_factor_init(fac)
        fmpz_factor_smooth(fac, self.val, bits, proved)
        res = [0] * fac.num
        for 0 <= i < fac.num:
            u = fmpz.__new__(fmpz)
            fmpz_set((<fmpz>u).val, &fac.p[i])
            exp = <long> fac.exp[i]
            res[i] = (u, exp)
        fmpz_factor_clear(fac)
        return res

    def is_prime(self):
        return fmpz_is_prime(self.val)

    def is_probable_prime(self):
        return fmpz_is_probabprime(self.val)

    def is_perfect_power(self):
        cdef int k
        cdef fmpz v = fmpz()
        k = fmpz_is_perfect_power(v.val, self.val)
        return k != 0

    def partitions_p(n):
        r"""
        Returns `p(n)`, the number of partitions of `n`, as an *fmpz*.

            >>> [fmpz(n).partitions_p() for n in range(8)]
            [1, 1, 2, 3, 5, 7, 11, 15]
            >>> fmpz(100).partitions_p()
            190569292
            >>> len(str(fmpz(10**9).partitions_p()))
            35219

        The partition function grows rapidly.
        On a 32-bit system, `n` must not be larger than about `10^{16}`.
        On a 64-bit system, `n` must not be larger than about `10^{20}`.
        For large `n`, this function benefits from setting ``ctx.threads = 2``
        on multicore systems.
        """
        cdef fmpz v = fmpz()
        partitions_fmpz_fmpz(v.val, n.val, 0)
        return v

    def moebius_mu(n):
        r"""
        Returns the Moebius function `\mu(n)` as an *fmpz*.

            >>> [fmpz(n).moebius_mu() for n in range(10)]
            [0, 1, -1, -1, 0, -1, 1, -1, 0, 0]
        """
        cdef fmpz v = fmpz()
        fmpz_set_si(v.val, fmpz_moebius_mu(n.val))
        return v

    @classmethod
    def fac_ui(cls, ulong n):
        r"""
        Returns the factorial `n!` as an *fmpz*.

            >>> fmpz.fac_ui(10)
            3628800
        """
        u = fmpz.__new__(fmpz)
        fmpz_fac_ui((<fmpz>u).val, n)
        return u

    @classmethod
    def primorial_ui(cls, ulong n):
        r"""
        Returns the product of all primes less than or equal to *n*
        (*n* primorial) as an *fmpz*.

            >>> fmpz.primorial_ui(10)
            210
        """
        u = fmpz.__new__(fmpz)
        fmpz_primorial((<fmpz>u).val, n)
        return u

    @classmethod
    def fib_ui(cls, ulong n):
        r"""
        Returns the Fibonacci number `F_n` as an *fmpz*.

            >>> fmpz.fib_ui(10)
            55
        """
        u = fmpz.__new__(fmpz)
        fmpz_fib_ui((<fmpz>u).val, n)
        return u

    def rising(s, ulong n):
        r"""
        Returns the rising factorial `s (s+1) \cdots (s+n-1)` as an *fmpz*.

            >>> fmpz(10).rising(5)
            240240
        """
        u = fmpz.__new__(fmpz)
        fmpz_rfac_ui((<fmpz>u).val, (<fmpz>s).val, n)
        return u

    @classmethod
    def bin_uiui(cls, ulong n, ulong k):
        r"""
        Returns the binomial coefficient `{n \choose k}` as an *fmpz*.
        """
        u = fmpz.__new__(fmpz)
        fmpz_bin_uiui((<fmpz>u).val, n, k)
        return u

    @classmethod
    def bell_number(cls, ulong n):
        r"""
        Returns the Bell number `B_n` as an *fmpz*.

            >>> fmpz.bell_number(10)
            115975
        """
        u = fmpz.__new__(fmpz)
        arith_bell_number((<fmpz>u).val, n)
        return u

    @classmethod
    def euler_number(cls, ulong n):
        r"""
        Returns the Euler number `E_n` as an *fmpz*.

            >>> fmpz.euler_number(10)
            -50521
        """
        u = fmpz.__new__(fmpz)
        arith_euler_number((<fmpz>u).val, n)
        return u

    @classmethod
    def stirling_s1(cls, ulong n, ulong k):
        r"""
        Returns the Stirling number of the first kind `S_1(n,k)` as an *fmpz*.

            >>> fmpz.stirling_s1(10,5)
            -269325
        """
        u = fmpz.__new__(fmpz)
        arith_stirling_number_1((<fmpz>u).val, n, k)
        return u

    @classmethod
    def stirling_s2(cls, ulong n, ulong k):
        r"""
        Returns the Stirling number of the second kind `S_2(n,k)` as an *fmpz*.

            >>> fmpz.stirling_s2(10,5)
            42525
        """
        u = fmpz.__new__(fmpz)
        arith_stirling_number_2((<fmpz>u).val, n, k)
        return u

    def divisor_sigma(n, k):
        r"""
        Returns the divisor sum `\sigma_k(n)` as an *fmpz*.

            >>> fmpz(60).divisor_sigma(0)
            12
            >>> fmpz(60).divisor_sigma(1)
            168
            >>> fmpz(60).divisor_sigma(10)
            605263138639095300
        """
        cdef fmpz v = fmpz()
        arith_divisor_sigma(v.val, k, n.val)
        return v

    def euler_phi(n):
        r"""
        Returns the Euler totient function `\varphi(n)` as an *fmpz*.

            >>> fmpz(60).euler_phi()
            16
            >>> fmpz(3**10).euler_phi()
            39366
        """
        cdef fmpz v = fmpz()
        arith_euler_phi(v.val, n.val)
        return v

    def __hash__(self):
        return hash(int(self))

    def height_bits(self, bint signed=False):
        if signed and fmpz_sgn(self.val) < 0:
            return -self.bit_length()
        else:
            return self.bit_length()

    def isqrt(self):
        cdef fmpz v
        if fmpz_sgn(self.val) < 0:
            raise ValueError("integer square root of a negative number")
        v = fmpz()
        fmpz_sqrt(v.val, self.val)
        return v

    def sqrtrem(self):
        cdef fmpz u, v
        if fmpz_sgn(self.val) < 0:
            raise ValueError("integer square root of a negative number")
        u = fmpz()
        v = fmpz()
        fmpz_sqrtrem(u.val, v.val, self.val)
        return u, v

    # warning: m should be prime!
    def sqrtmod(self, m):
        cdef fmpz v
        v = fmpz()
        m = fmpz(m)
        if not fmpz_sqrtmod(v.val, self.val, (<fmpz>m).val):
            raise ValueError("unable to compute modular square root")
        return v

    def root(self, long n):
        cdef fmpz v
        if fmpz_sgn(self.val) < 0:
            raise ValueError("integer root of a negative number")
        if n <= 0:
            raise ValueError("n >= 1 is required")
        v = fmpz()
        fmpz_root(v.val, self.val, n)
        return v

    def jacobi(self, other):
        cdef fmpz_struct tval[1]
        cdef int ttype = FMPZ_UNKNOWN
        ttype = fmpz_set_any_ref(tval, other)
        if ttype == FMPZ_UNKNOWN:
            raise TypeError("input must be an integer")
        u = fmpz.__new__(fmpz)
        v = fmpz_jacobi(self.val, tval)
        if ttype == FMPZ_TMP:
            fmpz_clear(tval)
        return fmpz(v)

