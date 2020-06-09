cdef _str_trunc(s, trunc=0):
    if trunc > 0 and len(s) > 3 * trunc:
        left = right = trunc
        omitted = len(s) - left - right
        return s[:left] + ("{...%s digits...}" % omitted) + s[-right:]
    return s

cdef arb_from_str(str s):
    s = s.strip()
    if ("/" in s) and ("+/-" not in s):
        return arb(fmpq(s))
    s = s.replace("±", "+/-")
    a = arb.__new__(arb)
    if arb_set_str((<arb>a).val, chars_from_str(s), getprec()) == 0:
        return a
    else:
        raise ValueError("invalid string for arb()")

cdef arb_set_mpmath_mpf(arb_t x, obj):
    sgn, man, exp, bc = obj

    if not man:
        if not exp:
            arb_zero(x)
        else:
            arb_indeterminate(x)
    else:
        man = fmpz(long(man))
        exp = fmpz(exp)

        arb_set_fmpz(x, (<fmpz>man).val)
        arb_mul_2exp_fmpz(x, x, (<fmpz>exp).val)

        if sgn:
            arb_neg(x, x)

cdef int arb_set_python(arb_t x, obj, bint allow_conversion) except -1:
    """
    Sets an arb_t given any Python object. If allow_conversion is set,
    conversions for the arb() constructor are done (from tuples, strings, etc.)
    """
    cdef fmpz_t t

    if typecheck(obj, arb):
        arb_set(x, (<arb>obj).val)
        return 1

    if typecheck(obj, arf):
        arb_set_arf(x, (<arf>obj).val)
        return 1

    if typecheck(obj, fmpz):
        arb_set_fmpz(x, (<fmpz>obj).val)
        return 1

    if typecheck(obj, fmpq):
        arb_set_fmpq(x, (<fmpq>obj).val, getprec())
        return 1

    if PY_MAJOR_VERSION < 3 and PyInt_Check(<PyObject*>obj):
        arb_set_si(x, PyInt_AS_LONG(<PyObject*>obj))
        return 1

    if PyLong_Check(<PyObject*>obj):
        fmpz_init(t)
        fmpz_set_pylong(t, obj)
        arb_set_fmpz(x, t)
        fmpz_clear(t)
        return 1

    if typecheck(obj, float):
        arf_set_d(arb_midref(x), PyFloat_AS_DOUBLE(<PyObject*>obj))
        mag_zero(arb_radref(x))
        return 1

    if hasattr(obj, "_mpf_"):
        arb_set_mpmath_mpf(x, obj._mpf_)
        return 1

    if allow_conversion and typecheck(obj, tuple):
        v = arf(obj)
        arb_set_arf(x, (<arf>v).val)
        return 1

    if allow_conversion and typecheck(obj, str):
        obj = arb_from_str(obj)
        arb_set(x, (<arb>obj).val)
        return 1

    return 0

cdef inline int arb_set_any_ref(arb_t x, obj):
    if typecheck(obj, arb):
        x[0] = (<arb>obj).val[0]
        return FMPZ_REF
    arb_init(x)
    if arb_set_python(x, obj, 0):
        return FMPZ_TMP
    return FMPZ_UNKNOWN

cdef any_as_arb(x):
    cdef arb t
    if typecheck(x, arb):
        return x
    t = arb()
    if arb_set_python(t.val, x, 0) == 0:
        raise TypeError("cannot create arb from type %s" % type(x))
    return t

cdef any_as_arb_or_notimplemented(x):
    cdef arb t
    if typecheck(x, arb):
        return x
    t = arb()
    if arb_set_python(t.val, x, 0) == 0:
        return NotImplemented
    return t

cdef _arb_div_(s, t):
    cdef arb_struct sval[1]
    cdef arb_struct tval[1]
    cdef int stype, ttype
    stype = arb_set_any_ref(sval, s)
    if stype == FMPZ_UNKNOWN:
        return NotImplemented
    ttype = arb_set_any_ref(tval, t)
    if ttype == FMPZ_UNKNOWN:
        return NotImplemented
    u = arb.__new__(arb)
    arb_div((<arb>u).val, sval, tval, getprec())
    if stype == FMPZ_TMP: arb_clear(sval)
    if ttype == FMPZ_TMP: arb_clear(tval)
    return u

cdef class arb(flint_scalar):
    ur"""
    Represents a real number `x` by a midpoint `m` and a radius `r`
    such that `x \in [m \pm r] = [m-r, m+r]`.
    The midpoint and radius are both floating-point numbers. The radius
    uses a fixed, implementation-defined precision (30 bits).
    The precision used for midpoints is controlled by :attr:`ctx.prec` (bits)
    or equivalently :attr:`ctx.dps` (digits).

    The constructor accepts a midpoint *mid* and a radius *rad*, either of
    which defaults to zero if omitted. The arguments can be tuples
    `(a, b)` representing exact floating-point data `a 2^b`, integers,
    floating-point numbers, rational strings, or decimal strings.
    If the radius is nonzero, it might be rounded up to a slightly larger
    value than the exact value passed by the user.

        >>> arb(10.25)
        10.2500000000000
        >>> print(1 / arb(4))  # exact
        0.250000000000000
        >>> print(1 / arb(3))  # approximate
        [0.333333333333333 +/- 3.71e-16]
        >>> print(arb("3.0"))
        3.00000000000000
        >>> print(arb("0.1"))
        [0.100000000000000 +/- 2.23e-17]
        >>> print(arb("1/10"))
        [0.100000000000000 +/- 2.23e-17]
        >>> print(arb("3.14159 +/- 0.00001"))
        [3.1416 +/- 2.01e-5]
        >>> ctx.dps = 50
        >>> print(arb("1/3"))
        [0.33333333333333333333333333333333333333333333333333 +/- 3.78e-51]
        >>> ctx.default()

    Converting to or from decimal results in some loss of accuracy.
    See :meth:`.arb.str` for details.
    """

    cdef arb_t val

    def __cinit__(self):
        arb_init(self.val)

    def __dealloc__(self):
        arb_clear(self.val)

    def __init__(self, mid=None, rad=None):
        if mid is not None:
            if arb_set_python(self.val, mid, 1) == 0:
                raise TypeError("cannot create arb from type %s" % type(mid))
        if rad is not None:
            rad = arb(rad)
            arb_add_error(self.val, (<arb>rad).val)
            #rad = arf(rad)
            #arb_add_error_arf(self.val, (<arf>rad).val)

    cpdef bint is_zero(self):
        return arb_is_zero(self.val)

    cpdef bint is_finite(self):
        return arb_is_finite(self.val)

    cpdef bint is_nan(self):
        return arf_is_nan(arb_midref(self.val))

    cpdef bint is_exact(self):
        return arb_is_exact(self.val)

    def man_exp(self):
        """
        Decomposes *self* into an integer mantissa and an exponent,
        returning an *fmpz* pair. Requires that *self* is exact
        and finite.

            >>> arb("1.1").mid().man_exp()
            (4953959590107545, -52)
            >>> arb("1.1").rad().man_exp()
            (1, -52)
            >>> arb(0).man_exp()
            (0, 0)
            >>> arb("1.1").man_exp()
            Traceback (most recent call last):
              ...
            ValueError: man_exp requires an exact, finite value
            >>> arb("+inf").man_exp()
            Traceback (most recent call last):
              ...
            ValueError: man_exp requires an exact, finite value
        """
        cdef fmpz man, exp
        if not self.is_finite() or not self.is_exact():
            raise ValueError("man_exp requires an exact, finite value")
        man = fmpz()
        exp = fmpz()
        arf_get_fmpz_2exp(man.val, exp.val, arb_midref(self.val))
        return man, exp

    def fmpq(self):
        cdef fmpq res
        if not self.is_finite() or not self.is_exact():
            raise ValueError("fmpq requires an exact, finite value")
        res = fmpq()
        arf_get_fmpq(res.val, arb_midref(self.val))
        return res

    def fmpz(self):
        cdef fmpz res
        if not self.is_integer():
            raise ValueError("fmpz requires an exact, integer value")
        res = fmpz()
        arf_get_fmpz(res.val, arb_midref(self.val), ARF_RND_DOWN)
        return res

    cpdef bint is_integer(self):
        return arb_is_int(self.val)

    def mid(self):
        """
        Returns the midpoint of *self* as an exact *arb*:

            >>> arb("1 +/- 0.3").mid()
            1.00000000000000
        """
        cdef arb x = arb()
        arf_set(arb_midref(x.val), arb_midref(self.val))
        return x

    def rad(self):
        """
        Returns the radius of *self* as an exact *arb*:

            >>> print(arb("1 +/- 0.3").rad().str(5, radius=False))
            0.30000
        """
        cdef arb x = arb()
        arf_set_mag(arb_midref(x.val), arb_radref(self.val))
        return x

    def abs_lower(self):
        """
        Lower bound for the absolute value of *self*.
        The output is an *arb* holding an exact floating-point number
        that has been rounded down to the current precision.

            >>> print(arb("-5 +/- 2").abs_lower().str(5, radius=False))
            3.0000
        """
        cdef arb x = arb()
        arb_get_abs_lbound_arf(arb_midref(x.val), self.val, getprec())
        return x

    def abs_upper(self):
        """
        Upper bound for the absolute value of *self*.
        The output is an *arb* holding an exact floating-point number
        that has been rounded up to the current precision.

            >>> print(arb("-5 +/- 2").abs_upper().str(5, radius=False))
            7.0000
        """
        cdef arb x = arb()
        arb_get_abs_ubound_arf(arb_midref(x.val), self.val, getprec())
        return x

    def lower(self):
        """
        Lower bound for *self* (towards `-\infty`).
        The output is an *arb* holding an exact floating-point number
        that has been rounded down to the current precision.

            >>> print(arb("-5 +/- 2").lower().str(5, radius=False))
            -7.0000
        """
        cdef arb x = arb()
        arb_get_lbound_arf(arb_midref(x.val), self.val, getprec())
        return x

    def upper(self):
        """
        Upper bound for *self* (towards `+\infty`).
        The output is an *arb* holding an exact floating-point number
        that has been rounded up to the current precision.

            >>> print(arb("-5 +/- 2").upper().str(5, radius=False))
            -3.0000
        """
        cdef arb x = arb()
        arb_get_ubound_arf(arb_midref(x.val), self.val, getprec())
        return x

    def mid_rad_10exp(self, long n=0):
        """
        Returns an *fmpz* triple (*mid*, *rad*, *exp*) where the larger of *mid*
        and *rad* has *n* digits plus a few digits (*n* defaults to the current
        precision), such that *self* is contained in
        `[\operatorname{mid} \pm \operatorname{rad}] 10^{\operatorname{exp}}`.

            >>> (arb(1) / 3).mid_rad_10exp(10)
            (333333333333333, 2, -15)
            >>> (arb(1) / 3).mid_rad_10exp(20)
            (3333333333333333148296162, 555111516, -25)
            >>> arb(0,1e-100).mid_rad_10exp(10)
            (0, 100000000376832, -114)
            >>> arb(-1,1e100).mid_rad_10exp()
            (0, 10000000083585662976, 81)

        """
        cdef fmpz mid, rad, exp
        if n <= 0:
            n = ctx.dps
        mid = fmpz()
        rad = fmpz()
        exp = fmpz()
        arb_get_fmpz_mid_rad_10exp(mid.val, rad.val, exp.val, self.val, n)
        return mid, rad, exp

    @property
    def _mpf_(self):
        try:
            import mpmath
            mpmath_mpz = mpmath.libmp.MPZ
        except ImportError:
            mpmath_mpz = long
        if not self.is_finite():
            return (0, mpmath_mpz(0), -123, -1)
        man, exp = self.mid().man_exp()
        man = mpmath_mpz(long(man))
        if man < 0:
            return (1, -man, long(exp), man.bit_length())
        else:
            return (0, man, long(exp), man.bit_length())

    def repr(self):
        mid = self.mid()
        rad = self.rad()
        if rad.is_zero():
            return "arb(%s)" % mid._repr_str()
        else:
            return "arb(%s, %s)" % (mid._repr_str(), rad._repr_str())

    def str(self, long n=0, bint radius=True, bint more=False, long condense=0):
        ur"""
        Produces a human-readable decimal representation of self, with
        up to *n* printed digits (which defaults to the current precision)
        for the midpoint. The output can be parsed by the *arb* constructor.

        Since the internal representation is binary, conversion
        to decimal (and back from decimal) is generally inexact.
        Binary-decimal-binary roundtrips may result in significantly
        larger intervals, and should therefore be done sparingly.

            >>> print(arb.pi().str())
            [3.14159265358979 +/- 3.34e-15]
            >>> print(arb.pi().str(5))
            [3.1416 +/- 7.35e-6]
            >>> print(arb.pi().str(5, radius=False))
            3.1416

        By default, the output is truncated so that all displayed digits
        are guaranteed to be correct, up to adding or subtracting 1 in the
        last displayed digit (as a special case, if the output ends with a
        string of 0s, the correct decimal expansion to infinite precision
        could have a string of 9s).

            >>> print((arb(1) - arb("1e-10")).str(5))
            [1.0000 +/- 4e-10]
            >>> print((arb(1) - arb("1e-10")).str(10))
            [0.9999999999 +/- 3e-15]

        To force more digits, set *more* to *True*.

            >>> print(arb("0.1").str(30))
            [0.100000000000000 +/- 2.23e-17]
            >>> print(arb("0.1").str(30, more=True))
            [0.0999999999999999916733273153113 +/- 1.39e-17]

        Note that setting *more* to *True* results in a smaller printed radius,
        since there is less error from the conversion back to decimal.

            >>> x = arb.pi().sin()
            >>> print(x.str())
            [+/- 3.46e-16]
            >>> print(x.str(more=True))
            [1.22460635382238e-16 +/- 2.23e-16]

        The error indicated in the output may be much larger than the actual
        error in the internal representation of *self*. For example, if *self*
        is accurate to 1000 digits and printing is done at 10-digit precision,
        the output might only appear as being accurate to 10 digits. It is
        even possible for *self* to be exact and have an inexact decimal
        representation.

        The *condense* option can be passed to show only leading and trailing
        digits of the fractional, integer and exponent parts of the output.

            >>> ctx.dps = 1000
            >>> print(arb.pi().str(condense=10))
            [3.1415926535{...979 digits...}9216420199 +/- 9.28e-1001]
            >>> print(arb.fac_ui(300).str(condense=10))
            3060575122{...595 digits...}0000000000.0000000000{...365 digits...}0000000000
            >>> print(arb(10**100).exp().str(condense=5))
            [1.53837{...989 digits...}96534e+43429{...90 digits...}17483 +/- 4.84e+43429{...90 digits...}16483]
            >>> ctx.default()

        """
        cdef ulong flags
        cdef char * s
        flags = 0
        if not radius:
            flags |= ARB_STR_NO_RADIUS
        if more:
            flags |= ARB_STR_MORE
        if condense > 0:
            flags |= ARB_STR_CONDENSE * condense
        if n <= 0:
            n = ctx.dps
        s = arb_get_str(self.val, n, flags)
        try:
            res = str_from_chars(s)
        finally:
            libc.stdlib.free(s)
        if ctx.unicode:
            res = res.replace("+/-", "±")
        return res

    def __float__(self):
        return arf_get_d(arb_midref(self.val), ARF_RND_NEAR)

    def __richcmp__(s, t, int op):
        cdef bint res
        cdef arb_struct sval[1]
        cdef arb_struct tval[1]
        cdef int stype, ttype
        stype = arb_set_any_ref(sval, s)
        if stype == FMPZ_UNKNOWN:
            return NotImplemented
        ttype = arb_set_any_ref(tval, t)
        if ttype == FMPZ_UNKNOWN:
            return NotImplemented
        res = 0
        if   op == 2: res = arb_eq(sval, tval)
        elif op == 3: res = arb_ne(sval, tval)
        elif op == 0: res = arb_lt(sval, tval)
        elif op == 1: res = arb_le(sval, tval)
        elif op == 4: res = arb_gt(sval, tval)
        elif op == 5: res = arb_ge(sval, tval)
        if stype == FMPZ_TMP: arb_clear(sval)
        if ttype == FMPZ_TMP: arb_clear(tval)
        return res

    def __contains__(self, other):
        other = any_as_arb(other)
        return arb_contains(self.val, (<arb>other).val)

    def contains(self, other):
        other = any_as_arb(other)
        return bool(arb_contains(self.val, (<arb>other).val))

    def contains_interior(self, other):
        other = any_as_arb(other)
        return bool(arb_contains_interior(self.val, (<arb>other).val))

    def overlaps(self, other):
        other = any_as_arb(other)
        return bool(arb_overlaps((<arb>self).val, (<arb>other).val))

    def contains_integer(self):
        return bool(arb_contains_int(self.val))

    @property
    def real(self):
        return self

    @property
    def imag(self):
        return arb()

    def __pos__(self):
        res = arb.__new__(arb)
        arb_set_round((<arb>res).val, (<arb>self).val, getprec())
        return res

    def __neg__(self):
        res = arb.__new__(arb)
        arb_neg_round((<arb>res).val, (<arb>self).val, getprec())
        return res

    def neg(self, bint exact=False):
        res = arb.__new__(arb)
        if exact:
            arb_set((<arb>res).val, (<arb>self).val)
        else:
            arb_set_round((<arb>res).val, (<arb>self).val, getprec())
        return res

    def __abs__(self):
        res = arb.__new__(arb)
        arb_abs((<arb>res).val, (<arb>self).val)
        arb_set_round((<arb>res).val, (<arb>res).val, getprec())
        return res

    def sgn(self):
        """
        Sign function, returning an *arb*.

            >>> arb(-3).sgn()
            -1.00000000000000
            >>> arb(0).sgn()
            0
            >>> arb("0 +/- 1").sgn()
            [+/- 1.01]
        """
        res = arb.__new__(arb)
        arb_sgn((<arb>res).val, (<arb>self).val)
        return res

    def __add__(s, t):
        cdef arb_struct sval[1]
        cdef arb_struct tval[1]
        cdef int stype, ttype
        stype = arb_set_any_ref(sval, s)
        if stype == FMPZ_UNKNOWN:
            return NotImplemented
        ttype = arb_set_any_ref(tval, t)
        if ttype == FMPZ_UNKNOWN:
            return NotImplemented
        u = arb.__new__(arb)
        arb_add((<arb>u).val, sval, tval, getprec())
        if stype == FMPZ_TMP: arb_clear(sval)
        if ttype == FMPZ_TMP: arb_clear(tval)
        return u

    def __sub__(s, t):
        cdef arb_struct sval[1]
        cdef arb_struct tval[1]
        cdef int stype, ttype
        stype = arb_set_any_ref(sval, s)
        if stype == FMPZ_UNKNOWN:
            return NotImplemented
        ttype = arb_set_any_ref(tval, t)
        if ttype == FMPZ_UNKNOWN:
            return NotImplemented
        u = arb.__new__(arb)
        arb_sub((<arb>u).val, sval, tval, getprec())
        if stype == FMPZ_TMP: arb_clear(sval)
        if ttype == FMPZ_TMP: arb_clear(tval)
        return u

    def __mul__(s, t):
        cdef arb_struct sval[1]
        cdef arb_struct tval[1]
        cdef int stype, ttype
        stype = arb_set_any_ref(sval, s)
        if stype == FMPZ_UNKNOWN:
            return NotImplemented
        ttype = arb_set_any_ref(tval, t)
        if ttype == FMPZ_UNKNOWN:
            return NotImplemented
        u = arb.__new__(arb)
        arb_mul((<arb>u).val, sval, tval, getprec())
        if stype == FMPZ_TMP: arb_clear(sval)
        if ttype == FMPZ_TMP: arb_clear(tval)
        return u

    def __truediv__(s, t):
        return _arb_div_(s, t)

    def __div__(s, t):
        return _arb_div_(s, t)

    def __pow__(s, t, modulus):
        cdef arb_struct sval[1]
        cdef arb_struct tval[1]
        cdef int stype, ttype
        if modulus is not None:
            raise TypeError("three-argument pow() not supported by arb type")
        stype = arb_set_any_ref(sval, s)
        if stype == FMPZ_UNKNOWN:
            return NotImplemented
        ttype = arb_set_any_ref(tval, t)
        if ttype == FMPZ_UNKNOWN:
            return NotImplemented
        u = arb.__new__(arb)
        arb_pow((<arb>u).val, sval, tval, getprec())
        if stype == FMPZ_TMP: arb_clear(sval)
        if ttype == FMPZ_TMP: arb_clear(tval)
        return u

    def floor(s):
        ur"""
        Floor function `\lfloor s \rfloor`.

            >>> print(arb.pi().floor())
            3.00000000000000
            >>> print((arb.pi() - arb.pi()).floor().str(more=True))
            [-0.500000000000000 +/- 0.501]
        """
        u = arb.__new__(arb)
        arb_floor((<arb>u).val, (<arb>s).val, getprec())
        return u

    def ceil(s):
        ur"""
        Ceiling function `\lceil s \rceil`.

            >>> print(arb.pi().ceil())
            4.00000000000000
            >>> print((arb.pi() - arb.pi()).ceil().str(more=True))
            [0.500000000000000 +/- 0.501]
        """
        u = arb.__new__(arb)
        arb_ceil((<arb>u).val, (<arb>s).val, getprec())
        return u

    def sqrt(s):
        r"""
        Square root `\sqrt{s}`.

            >>> showgood(lambda: arb(3).sqrt(), dps=25)
            1.732050807568877293527446
            >>> showgood(lambda: arb(0).sqrt(), dps=25)
            0
            >>> showgood(lambda: arb(-1).sqrt(), dps=25)
            Traceback (most recent call last):
              ...
            ValueError: no convergence (maxprec=960, try higher maxprec)

        This function is undefined for negative input.
        Use :meth:`.acb.sqrt` for the complex extension.
        """
        u = arb.__new__(arb)
        arb_sqrt((<arb>u).val, (<arb>s).val, getprec())
        return u

    def rsqrt(s):
        r"""
        Reciprocal square root `1/\sqrt{s}`.

            >>> showgood(lambda: arb(3).rsqrt(), dps=25)
            0.5773502691896257645091488
            >>> showgood(lambda: arb(0).rsqrt(), dps=25)
            Traceback (most recent call last):
              ...
            ValueError: no convergence (maxprec=960, try higher maxprec)
            >>> showgood(lambda: arb(-1).rsqrt(), dps=25)
            Traceback (most recent call last):
              ...
            ValueError: no convergence (maxprec=960, try higher maxprec)

        This function is undefined for negative input.
        Use :meth:`.acb.rsqrt` for the complex extension.
        """
        u = arb.__new__(arb)
        arb_rsqrt((<arb>u).val, (<arb>s).val, getprec())
        return u

    def exp(s):
        r"""
        Exponential function `\exp(s)`.

            >>> showgood(lambda: arb(1).exp(), dps=25)
            2.718281828459045235360287
        """
        u = arb.__new__(arb)
        arb_exp((<arb>u).val, (<arb>s).val, getprec())
        return u

    def expm1(s):
        r"""
        Exponential function `\exp(s) - 1`, computed accurately for small *s*.

            >>> showgood(lambda: (arb(10) ** -8).expm1(), dps=25)
            1.000000005000000016666667e-8
        """
        u = arb.__new__(arb)
        arb_expm1((<arb>u).val, (<arb>s).val, getprec())
        return u

    def log(s):
        r"""
        Natural logarithm `\log(s)`.

            >>> showgood(lambda: arb(2).log(), dps=25)
            0.6931471805599453094172321
            >>> showgood(lambda: arb(100).exp().log(), dps=25)
            100.0000000000000000000000
            >>> showgood(lambda: arb(-1).sqrt(), dps=25)
            Traceback (most recent call last):
              ...
            ValueError: no convergence (maxprec=960, try higher maxprec)

        This function is undefined for negative input.
        Use :meth:`.acb.log` for the complex extension.
        """
        u = arb.__new__(arb)
        arb_log((<arb>u).val, (<arb>s).val, getprec())
        return u

    def log1p(s):
        r"""
        Natural logarithm `\log(1+s)`, computed accurately for small *s*.

            >>> showgood(lambda: acb(1).log1p(), dps=25)
            0.6931471805599453094172321
            >>> showgood(lambda: arb("1e-100000000000000000").log1p(), dps=25)
            1.000000000000000000000000e-100000000000000000

        This function is undefined for `s \le -1`.
        Use :meth:`.acb.log1p` for the complex extension.
        """
        u = arb.__new__(arb)
        arb_log1p((<arb>u).val, (<arb>s).val, getprec())
        return u

    def sin(s):
        r"""
        Sine function `\sin(s)`.

            >>> showgood(lambda: arb(1).sin(), dps=25)
            0.8414709848078965066525023
        """
        u = arb.__new__(arb)
        arb_sin((<arb>u).val, (<arb>s).val, getprec())
        return u

    def cos(s):
        r"""
        Cosine function `\cos(s)`.

            >>> showgood(lambda: arb(1).cos(), dps=25)
            0.5403023058681397174009366
        """
        u = arb.__new__(arb)
        arb_cos((<arb>u).val, (<arb>s).val, getprec())
        return u

    def sin_cos(s):
        r"""
        Computes `\sin(s)` and `\cos(s)` simultaneously.

            >>> showgood(lambda: arb(1).sin_cos(), dps=25)
            (0.8414709848078965066525023, 0.5403023058681397174009366)
        """
        u = arb.__new__(arb)
        v = arb.__new__(arb)
        arb_sin_cos((<arb>u).val, (<arb>v).val, (<arb>s).val, getprec())
        return u, v

    def sin_pi(s):
        r"""
        Sine function `\sin(\pi s)`.

            >>> showgood(lambda: arb(0.75).sin_pi(), dps=25)
            0.7071067811865475244008444
            >>> showgood(lambda: arb(1).sin_pi(), dps=25)
            0
        """
        u = arb.__new__(arb)
        arb_sin_pi((<arb>u).val, (<arb>s).val, getprec())
        return u

    def cos_pi(s):
        r"""
        Cosine function `\cos(\pi s)`.

            >>> showgood(lambda: arb(0.75).cos_pi(), dps=25)
            -0.7071067811865475244008444
            >>> showgood(lambda: arb(0.5).cos_pi(), dps=25)
            0
        """
        u = arb.__new__(arb)
        arb_cos_pi((<arb>u).val, (<arb>s).val, getprec())
        return u

    def sin_cos_pi(s):
        r"""
        Computes `\sin(\pi s)` and `\cos(\pi s)` simultaneously.

            >>> showgood(lambda: arb(0.75).sin_cos_pi(), dps=25)
            (0.7071067811865475244008444, -0.7071067811865475244008444)
        """
        u = arb.__new__(arb)
        v = arb.__new__(arb)
        arb_sin_cos_pi((<arb>u).val, (<arb>v).val, (<arb>s).val, getprec())
        return u, v

    def tan(s):
        r"""
        Tangent function `\tan(s)`.

            >>> showgood(lambda: arb(1).tan(), dps=25)
            1.557407724654902230506975
        """
        u = arb.__new__(arb)
        arb_tan((<arb>u).val, (<arb>s).val, getprec())
        return u

    def cot(s):
        r"""
        Cotangent function `\cot(s)`.

            >>> showgood(lambda: arb(1).cot(), dps=25)
            0.6420926159343307030064200
        """
        u = arb.__new__(arb)
        arb_cot((<arb>u).val, (<arb>s).val, getprec())
        return u

    def tan_pi(s):
        r"""
        Tangent function `\tan(\pi s)`.

            >>> showgood(lambda: arb(0.125).tan_pi(), dps=25)
            0.4142135623730950488016887
            >>> showgood(lambda: arb(1).tan_pi(), dps=25)
            0
        """
        u = arb.__new__(arb)
        arb_tan_pi((<arb>u).val, (<arb>s).val, getprec())
        return u

    def cot_pi(s):
        r"""
        Cotangent function `\cot(\pi s)`.

            >>> showgood(lambda: arb(0.125).cot_pi(), dps=25)
            2.414213562373095048801689
            >>> showgood(lambda: arb(0.5).cot_pi(), dps=25)
            0
        """
        u = arb.__new__(arb)
        arb_cot_pi((<arb>u).val, (<arb>s).val, getprec())
        return u

    @staticmethod
    def sin_pi_fmpq(fmpq s):
        r"""
        Returns the algebraic sine value `\sin(\pi s)`.

            >>> showgood(lambda: arb.sin_pi_fmpq(fmpq(3,4)), dps=25)
            0.7071067811865475244008444
        """
        u = arb.__new__(arb)
        arb_sin_pi_fmpq((<arb>u).val, (<fmpq>s).val, getprec())
        return u

    @staticmethod
    def cos_pi_fmpq(fmpq s):
        r"""
        Returns the algebraic cosine value `\cos(\pi s)`.

            >>> showgood(lambda: arb.cos_pi_fmpq(fmpq(3,4)), dps=25)
            -0.7071067811865475244008444
        """
        u = arb.__new__(arb)
        arb_cos_pi_fmpq((<arb>u).val, (<fmpq>s).val, getprec())
        return u

    @staticmethod
    def sin_cos_pi_fmpq(fmpq s):
        r"""
        Computes `\sin(\pi s)` and `\cos(\pi s)` simultaneously.

            >>> showgood(lambda: arb.sin_cos_pi_fmpq(fmpq(3,4)), dps=25)
            (0.7071067811865475244008444, -0.7071067811865475244008444)
        """
        u = arb.__new__(arb)
        v = arb.__new__(arb)
        arb_sin_cos_pi_fmpq((<arb>u).val, (<arb>v).val, (<fmpq>s).val, getprec())
        return u, v

    def sec(s):
        """
        Secant function `\operatorname{sec}(s)`.

            >>> showgood(lambda: arb(1).sec(), dps=25)
            1.850815717680925617911753
        """
        u = arb.__new__(arb)
        arb_sec((<arb>u).val, (<arb>s).val, getprec())
        return u

    def csc(s):
        """
        Cosecant function `\operatorname{csc}(s)`.

            >>> showgood(lambda: arb(1).csc(), dps=25)
            1.188395105778121216261599
        """
        u = arb.__new__(arb)
        arb_csc((<arb>u).val, (<arb>s).val, getprec())
        return u

    def sinc(s):
        r"""
        Sinc function, `\operatorname{sinc}(x) = \sin(x)/x`.

            >>> showgood(lambda: arb(3).sinc(), dps=25)
            0.04704000268662240736691493
        """
        u = arb.__new__(arb)
        arb_sinc((<arb>u).val, (<arb>s).val, getprec())
        return u

    def sinc_pi(s):
        r"""
        Normalized sinc function, `\operatorname{sinc}(\pi x) = \sin(\pi x)/(\pi x)`.

            >>> showgood(lambda: arb(1.5).sinc_pi(), dps=25)
            -0.2122065907891937810251784
            >>> showgood(lambda: arb(2).sinc_pi(), dps=25)
            0
        """
        u = arb.__new__(arb)
        arb_sinc_pi((<arb>u).val, (<arb>s).val, getprec())
        return u

    def atan(s):
        r"""
        Inverse tangent `\operatorname{atan}(s)`.

            >>> showgood(lambda: arb(1).atan(), dps=25)
            0.7853981633974483096156608
        """
        u = arb.__new__(arb)
        arb_atan((<arb>u).val, (<arb>s).val, getprec())
        return u

    @staticmethod
    def atan2(s, t):
        r"""
        Two-argument inverse tangent `\operatorname{atan2}(s,t)`.

            >>> showgood(lambda: arb.atan2(-10,-5), dps=25)
            -2.034443935795702735445578
        """
        s = any_as_arb(s)
        t = any_as_arb(t)
        u = arb.__new__(arb)
        arb_atan2((<arb>u).val, (<arb>s).val, (<arb>t).val, getprec())
        return u

    def acos(s):
        r"""
        Inverse cosine `\operatorname{acos}(s)`.

            >>> showgood(lambda: arb(0).acos(), dps=25)
            1.570796326794896619231322
        """
        u = arb.__new__(arb)
        arb_acos((<arb>u).val, (<arb>s).val, getprec())
        return u

    def asin(s):
        r"""
        Inverse sine `\operatorname{asin}(s)`.

            >>> showgood(lambda: arb(1).asin(), dps=25)
            1.570796326794896619231322
        """
        u = arb.__new__(arb)
        arb_asin((<arb>u).val, (<arb>s).val, getprec())
        return u

    def atanh(s):
        r"""
        Inverse hyperbolic tangent `\operatorname{atanh}(s)`.

            >>> showgood(lambda: arb("0.99").atanh(), dps=25)
            2.646652412362246197705061
        """
        u = arb.__new__(arb)
        arb_atanh((<arb>u).val, (<arb>s).val, getprec())
        return u

    def asinh(s):
        r"""
        Inverse hyperbolic sine `\operatorname{asinh}(s)`.

            >>> showgood(lambda: arb(1000).asinh(), dps=25)
            7.600902709541988611523290
        """
        u = arb.__new__(arb)
        arb_asinh((<arb>u).val, (<arb>s).val, getprec())
        return u

    def acosh(s):
        r"""
        Inverse hyperbolic cosine `\operatorname{acosh}(s)`.

            >>> showgood(lambda: arb(2).acosh(), dps=25)
            1.316957896924816708625046
        """
        u = arb.__new__(arb)
        arb_acosh((<arb>u).val, (<arb>s).val, getprec())
        return u

    def sinh(s):
        r"""
        Hyperbolic sine `\sinh(s)`.

            >>> showgood(lambda: arb(1).sinh(), dps=25)
            1.175201193643801456882382
        """
        u = arb.__new__(arb)
        arb_sinh((<arb>u).val, (<arb>s).val, getprec())
        return u

    def cosh(s):
        r"""
        Hyperbolic cosine `\cosh(s)`.

            >>> showgood(lambda: arb(1).cosh(), dps=25)
            1.543080634815243778477906
        """
        u = arb.__new__(arb)
        arb_cosh((<arb>u).val, (<arb>s).val, getprec())
        return u

    def sinh_cosh(s):
        r"""
        Computes `\sinh(s)` and `\cosh(s)` simultaneously.

            >>> showgood(lambda: arb(1).sinh_cosh(), dps=25)
            (1.175201193643801456882382, 1.543080634815243778477906)
        """
        u = arb.__new__(arb)
        v = arb.__new__(arb)
        arb_sinh_cosh((<arb>u).val, (<arb>v).val, (<arb>s).val, getprec())
        return u, v

    def tanh(s):
        r"""
        Hyperbolic tangent `\tanh(s)`.

            >>> showgood(lambda: arb(1).tanh(), dps=25)
            0.7615941559557648881194583
        """
        u = arb.__new__(arb)
        arb_tanh((<arb>u).val, (<arb>s).val, getprec())
        return u

    def coth(s):
        r"""
        Hyperbolic cotangent `\coth(s)`.

            >>> showgood(lambda: arb(1).coth(), dps=25)
            1.313035285499331303636161
        """
        u = arb.__new__(arb)
        arb_coth((<arb>u).val, (<arb>s).val, getprec())
        return u

    def sech(s):
        r"""
        Hyperbolic secant `\operatorname{sech}(s)`.

            >>> showgood(lambda: arb(1).sech(), dps=25)
            0.6480542736638853995749774
        """
        u = arb.__new__(arb)
        arb_sech((<arb>u).val, (<arb>s).val, getprec())
        return u

    def csch(s):
        r"""
        Hyperbolic cosecant `\operatorname{csch}(s)`.

            >>> showgood(lambda: arb(1).csch(), dps=25)
            0.8509181282393215451338428
        """
        u = arb.__new__(arb)
        arb_csch((<arb>u).val, (<arb>s).val, getprec())
        return u

    def gamma(s):
        """
        Gamma function `\Gamma(s)`.

            >>> showgood(lambda: arb(10).gamma(), dps=25)
            362880.0000000000000000000
            >>> showgood(lambda: arb(-2.5).gamma(), dps=25)
            -0.9453087204829418812256893
            >>> showgood(lambda: (arb.pi() ** 10).gamma(), dps=25)
            1.705646271897306403570389e+424898
            >>> showgood(lambda: arb(0).gamma(), dps=25)  # pole
            Traceback (most recent call last):
              ...
            ValueError: no convergence (maxprec=960, try higher maxprec)
        """
        u = arb.__new__(arb)
        arb_gamma((<arb>u).val, (<arb>s).val, getprec())
        return u

    @staticmethod
    def gamma_fmpq(fmpq s):
        """
        Computes the gamma function `\Gamma(s)` of a given *fmpq* *s*,
        exploiting the fact that *s* is an exact rational number to
        improve performance.

            >>> showgood(lambda: arb.gamma_fmpq(fmpq(1,4)), dps=25)
            3.625609908221908311930685
        """
        u = arb.__new__(arb)
        arb_gamma_fmpq((<arb>u).val, (<fmpq>s).val, getprec())
        return u

    def rgamma(s):
        ur"""
        Reciprocal gamma function `1/\Gamma(s)`, avoiding
        division by zero at the poles of the gamma function.

            >>> showgood(lambda: arb(1.5).rgamma(), dps=25)
            1.128379167095512573896159
            >>> print(arb(0).rgamma())
            0
            >>> print(arb(-1).rgamma())
            0
            >>> print(arb(-3,1e-10).rgamma())
            [+/- 6.01e-10]
        """
        u = arb.__new__(arb)
        arb_rgamma((<arb>u).val, (<arb>s).val, getprec())
        return u

    def lgamma(s):
        """
        Logarithmic gamma function `\log \Gamma(s)`.

            >>> showgood(lambda: arb(100).lgamma(), dps=25)
            359.1342053695753987760440

        This function is undefined for negative `s`. Use :meth:`.acb.lgamma`
        for the complex extension.
        """
        u = arb.__new__(arb)
        arb_lgamma((<arb>u).val, (<arb>s).val, getprec())
        return u

    def digamma(s):
        """
        Digamma function `\psi(s)`.

            >>> showgood(lambda: arb(1).digamma(), dps=25)
            -0.5772156649015328606065121
        """
        u = arb.__new__(arb)
        arb_digamma((<arb>u).val, (<arb>s).val, getprec())
        return u

    def rising(s, n):
        """
        Rising factorial `(s)_n`.

            >>> showgood(lambda: arb.pi().rising(0), dps=25)
            1.000000000000000000000000
            >>> showgood(lambda: arb.pi().rising(10), dps=25)
            299606572.3661012684972888
            >>> showgood(lambda: arb.pi().rising(0.5), dps=25)
            1.703592785410167015590330
        """
        u = arb.__new__(arb)
        n = any_as_arb(n)
        arb_rising((<arb>u).val, (<arb>s).val, (<arb>n).val, getprec())
        return u

    @staticmethod
    def rising_fmpq_ui(fmpq s, ulong n):
        """
        Computes the rising factorial `(s)_n` where *s* is a rational
        number and *n* is an unsigned integer. The current implementation
        does not use the gamma function, so *n* should be moderate.

            >>> showgood(lambda: arb.rising_fmpq_ui(fmpq(-1,3), 100), dps=25)
            -4.960517984074284420131903e+154
        """
        u = arb.__new__(arb)
        arb_rising_fmpq_ui((<arb>u).val, (<fmpq>s).val, n, getprec())
        return u

    def rising2(s, ulong n):
        """
        Computes the rising factorial `(s)_n` where *n* is an unsigned
        integer, along with the first derivative with respect to `(s)_n`.
        The current implementation does not use the gamma function,
        so *n* should be moderate.

            >>> u, v = arb(3).rising2(5)
            >>> print(u); print(v)
            2520.00000000000
            2754.00000000000
        """
        u = arb.__new__(arb)
        v = arb.__new__(arb)
        arb_rising2_ui((<arb>u).val, (<arb>v).val, (<arb>s).val, n, getprec())
        return u, v

    def zeta(s, a=None):
        """
        Riemann zeta function `\zeta(s)` or the Hurwitz
        zeta function `\zeta(s,a)` if a second parameter is passed.

            >>> showgood(lambda: arb(4.25).zeta(), dps=25)
            1.066954190711214532450370
            >>> showgood(lambda: arb(4.25).zeta(2.75), dps=25)
            0.01991885526414599096374229

        This function is undefined for some
        combinations of `s, a`. Use :meth:`.acb.zeta` for the complex extension.
        """
        u = arb.__new__(arb)
        if a is None:
            arb_zeta((<arb>u).val, (<arb>s).val, getprec())
        else:
            a = any_as_arb(a)
            arb_hurwitz_zeta((<arb>u).val, (<arb>s).val, (<arb>a).val, getprec())
        return u

    def agm(s, t=1):
        """
        The arithmetic-geometric mean `M(s,t)`, or `M(s) = M(s,1)`
        if no extra parameter is passed.

            >>> showgood(lambda: arb(2).sqrt().agm(), dps=25)
            1.198140234735592207439922
            >>> showgood(lambda: arb(2).agm(100), dps=25)
            29.64467336236643624632443
            >>> showgood(lambda: arb(0).agm(0), dps=25)
            0
            >>> showgood(arb(0).agm, dps=25)
            0
            >>> showgood(arb(-2).agm, dps=25)   # not defined for negatives
            Traceback (most recent call last):
              ...
            ValueError: no convergence (maxprec=960, try higher maxprec)

        This function is undefined for negative input.
        Use :meth:`.acb.agm` for the complex extension.
        """
        t = any_as_arb(t)
        u = arb.__new__(arb)
        arb_agm((<arb>u).val, (<arb>s).val, (<arb>t).val, getprec())
        return u

    @staticmethod
    def bernoulli(n):
        """
        Computes the Bernoulli number `B_n` as an *arb*,
        where *n* is a given integer.

            >>> showgood(lambda: arb.bernoulli(1), dps=25)
            -0.5000000000000000000000000
            >>> showgood(lambda: arb.bernoulli(2), dps=25)
            0.1666666666666666666666667
            >>> showgood(lambda: arb.bernoulli(10**7), dps=25)
            -4.983176441416329828329241e+57675260
            >>> showgood(lambda: arb.bernoulli(10**30), dps=25)
            -5.048207707665410387507573e+28767525649738633122783863898083
        """
        u = arb.__new__(arb)
        n = fmpz(n)
        arb_bernoulli_fmpz((<arb>u).val, (<fmpz>n).val, getprec())
        return u

    @staticmethod
    def bell_number(n):
        r"""
        Computes the Bell number `B_n` as an *arb*, where *n* is a given integer.

            >>> showgood(lambda: arb.bell_number(1), dps=25)
            1.000000000000000000000000
            >>> showgood(lambda: arb.bell_number(10), dps=25)
            115975.0000000000000000000
            >>> showgood(lambda: arb.bell_number(10**20), dps=25)
            5.382701131762816107395343e+1794956117137290721328
        """
        u = arb.__new__(arb)
        n = fmpz(n)
        arb_bell_fmpz((<arb>u).val, (<fmpz>n).val, getprec())
        return u

    @staticmethod
    def partitions_p(n):
        """
        Number of partitions of the integer *n*, evaluated as an *arb*.

            >>> showgood(lambda: arb.partitions_p(10), dps=25)
            42.00000000000000000000000
            >>> showgood(lambda: arb.partitions_p(100), dps=25)
            190569292.0000000000000000
            >>> showgood(lambda: arb.partitions_p(10**50), dps=25)
            3.285979358867807890529967e+11140086280105007830283557
        """
        u = arb.__new__(arb)
        n = fmpz(n)
        arb_partitions_fmpz((<arb>u).val, (<fmpz>n).val, getprec())
        return u

    def bernoulli_poly(s, ulong n):
        """
        Returns the value of the Bernoulli polynomial `B_n(s)`.

            >>> showgood(lambda: arb("3.141").bernoulli_poly(5), dps=25)
            113.5165117028492010000000
        """
        u = arb.__new__(arb)
        arb_bernoulli_poly_ui((<arb>u).val, n, s.val, getprec())
        return u

    def fac(s):
        r"""
        Factorial, equivalent to `s! = \Gamma(s+1)`.

            >>> showgood(lambda: arb(5).fac(), dps=25)
            120.0000000000000000000000
            >>> showgood(lambda: arb("0.1").fac(), dps=25)
            0.9513507698668731836292487
        """
        u = arb.__new__(arb)
        arb_add_ui((<arb>u).val, s.val, 1, getprec())
        arb_gamma((<arb>u).val, (<arb>u).val, getprec())
        return u

    @staticmethod
    def fac_ui(ulong n):
        """
        Factorial `n!`, given an integer.

            >>> print(arb.fac_ui(10))
            3628800.00000000
            >>> showgood(lambda: arb.fac_ui(10**9).log(), dps=25)
            19723265848.22698260792313
        """
        u = arb.__new__(arb)
        arb_fac_ui((<arb>u).val, n, getprec())
        return u

    def bin(s, ulong k):
        """
        Binomial coefficient `{s \choose k}`. Currently *k* is limited
        to an integer; this restriction will be removed in the future
        by using the gamma function.

            >>> print(arb(10).bin(5))
            252.000000000000
            >>> showgood(lambda: arb.pi().bin(100), dps=25)
            5.478392395095119521549286e-9
        """
        u = arb.__new__(arb)
        arb_bin_ui((<arb>u).val, (<arb>s).val, k, getprec())
        return u

    @staticmethod
    def bin_uiui(ulong n, ulong k):
        """
        Binomial coefficient `{n \choose k}`.

            >>> print(arb.bin_uiui(10, 5))
            252.000000000000
        """
        u = arb.__new__(arb)
        arb_bin_uiui((<arb>u).val, n, k, getprec())
        return u

    @staticmethod
    def fib(n):
        """
        Computes the Fibonacci number `F_n` as an *arb*,
        where *n* is a given integer.

            >>> print(arb.fib(10))
            55.0000000000000
            >>> showgood(lambda: arb.fib(10**100).log(), dps=25)
            4.812118250596034474977589e+99
        """
        u = arb.__new__(arb)
        n = fmpz(n)
        arb_fib_fmpz((<arb>u).val, (<fmpz>n).val, getprec())
        return u

    def polylog(self, s):
        """
        Polylogarithm `\operatorname{Li}_s(z)` where
        the argument *z* is given by *self* and the order *s* is given
        as an extra parameter.

            >>> showgood(lambda: arb(-1).polylog(2), dps=25)
            -0.8224670334241132182362076
            >>> showgood(lambda: arb(-3).polylog(1.75), dps=25)
            -1.813689945878060161612620
            >>> showgood(lambda: arb(-2.5).polylog(4.75), dps=25)
            -2.322090601785704585092044

        This function is undefined for some
        combinations of `s, z`. Use :meth:`.acb.polylog` for the complex extension.
        """
        u = arb.__new__(arb)
        s = any_as_arb(s)
        arb_polylog((<arb>u).val, (<arb>s).val, (<arb>self).val, getprec())
        return u

    def airy_ai(s, int derivative=0):
        r"""
        Airy function `\operatorname{Ai}(s)`, or
        `\operatorname{Ai}'(s)` if *derivative* is 1.

            >>> showgood(lambda: arb(-1).airy_ai(), dps=25)
            0.5355608832923521187995166
            >>> showgood(lambda: arb(-1).airy_ai(derivative=1), dps=25)
            -0.01016056711664520939504547
        """
        u = arb.__new__(arb)
        if derivative == 0:
            arb_hypgeom_airy((<arb>u).val, NULL, NULL, NULL, (<arb>s).val, getprec())
        elif derivative == 1:
            arb_hypgeom_airy(NULL, (<arb>u).val, NULL, NULL, (<arb>s).val, getprec())
        else:
            raise ValueError("derivative must be 0 or 1")
        return u

    def airy_bi(s, int derivative=0):
        r"""
        Airy function `\operatorname{Bi}(s)`, or
        `\operatorname{Bi}'(s)` if *derivative* is 1.

            >>> showgood(lambda: arb(-1).airy_bi(), dps=25)
            0.1039973894969446118886900
            >>> showgood(lambda: arb(-1).airy_bi(derivative=1), dps=25)
            0.5923756264227923508167792
        """
        u = arb.__new__(arb)
        if derivative == 0:
            arb_hypgeom_airy(NULL, NULL, (<arb>u).val, NULL, (<arb>s).val, getprec())
        elif derivative == 1:
            arb_hypgeom_airy(NULL, NULL, NULL, (<arb>u).val, (<arb>s).val, getprec())
        else:
            raise ValueError("derivative must be 0 or 1")
        return u

    def airy(s):
        r"""
        Computes the Airy function values `\operatorname{Ai}(s)`,
        `\operatorname{Ai}'(s)`, `\operatorname{Bi}(s)`,
        `\operatorname{Bi}'(s)` simultaneously, returning a tuple.

            >>> showgood(lambda: arb(-1).airy(), dps=10)
            (0.5355608833, -0.01016056712, 0.1039973895, 0.5923756264)
        """
        u = arb.__new__(arb)
        v = arb.__new__(arb)
        w = arb.__new__(arb)
        z = arb.__new__(arb)
        arb_hypgeom_airy((<arb>u).val, (<arb>v).val,
                        (<arb>w).val, (<arb>z).val, (<arb>s).val, getprec())
        return u, v, w, z

    @staticmethod
    def airy_ai_zero(n, int derivative=0):
        r"""
        For positive integer *n*, returns the zero `a_n` of the
        Airy function `\operatorname{Ai}(s)`, or the corresponding
        zero `a'_n` of `\operatorname{Ai}'(s)` if *derivative* is 1.

            >>> showgood(lambda: arb.airy_ai_zero(1), dps=25)
            -2.338107410459767038489197
            >>> showgood(lambda: arb.airy_ai_zero(1000), dps=25)
            -281.0315196125215528353364
            >>> showgood(lambda: arb.airy_ai_zero(1, derivative=1), dps=25)
            -1.018792971647471089017325
        """
        n = fmpz(n)
        if n <= 0:
            raise ValueError("index must be >= 1")
        u = arb.__new__(arb)
        if derivative == 0:
            arb_hypgeom_airy_zero((<arb>u).val, NULL, NULL, NULL, (<fmpz>n).val, getprec())
        elif derivative == 1:
            arb_hypgeom_airy_zero(NULL, (<arb>u).val, NULL, NULL, (<fmpz>n).val, getprec())
        else:
            raise ValueError("derivative must be 0 or 1")
        return u

    @staticmethod
    def airy_bi_zero(n, int derivative=0):
        r"""
        For positive integer *n*, returns the zero `b_n` of the
        Airy function `\operatorname{Bi}(s)`, or the corresponding
        zero `b'_n` of `\operatorname{Bi}'(s)` if *derivative* is 1.

            >>> showgood(lambda: arb.airy_bi_zero(1), dps=25)
            -1.173713222709127924919980
            >>> showgood(lambda: arb.airy_bi_zero(1000), dps=25)
            -280.9378112034152401578834
            >>> showgood(lambda: arb.airy_bi_zero(1, derivative=1), dps=25)
            -2.294439682614123246622459
        """
        n = fmpz(n)
        if n <= 0:
            raise ValueError("index must be >= 1")
        u = arb.__new__(arb)
        if derivative == 0:
            arb_hypgeom_airy_zero(NULL, NULL, (<arb>u).val, NULL, (<fmpz>n).val, getprec())
        elif derivative == 1:
            arb_hypgeom_airy_zero(NULL, NULL, NULL, (<arb>u).val, (<fmpz>n).val, getprec())
        else:
            raise ValueError("derivative must be 0 or 1")
        return u

    def chebyshev_t(s, n):
        r"""
        Chebyshev function of the first kind `T_n(s)`.

            >>> showgood(lambda: (arb(1)/3).chebyshev_t(3), dps=25)
            -0.8518518518518518518518519
        """
        v = arb.__new__(arb)
        n = any_as_arb(n)
        arb_hypgeom_chebyshev_t((<arb>v).val, (<arb>n).val, (<arb>s).val, getprec())
        return v

    def chebyshev_u(s, n):
        r"""
        Chebyshev function of the second kind `U_n(s)`.

            >>> showgood(lambda: (arb(1)/3).chebyshev_u(3), dps=25)
            -1.037037037037037037037037
        """
        v = arb.__new__(arb)
        n = any_as_arb(n)
        arb_hypgeom_chebyshev_u((<arb>v).val, (<arb>n).val, (<arb>s).val, getprec())
        return v

    def jacobi_p(s, n, a, b):
        r"""
        Jacobi polynomial (or Jacobi function) `P_n^{a,b}(s)`.

            >>> showgood(lambda: (arb(1)/3).jacobi_p(5, 0.25, 0.5), dps=25)
            0.4131944444444444444444444
        """
        v = arb.__new__(arb)
        n = any_as_arb(n)
        a = any_as_arb(a)
        b = any_as_arb(b)
        arb_hypgeom_jacobi_p((<arb>v).val, (<arb>n).val, (<arb>a).val, (<arb>b).val, (<arb>s).val, getprec())
        return v

    def gegenbauer_c(s, n, m):
        r"""
        Gegenbauer function `C_n^{m}(s)`.

            >>> showgood(lambda: (arb(1)/3).gegenbauer_c(5, 0.25), dps=25)
            0.1321855709876543209876543
        """
        v = arb.__new__(arb)
        n = any_as_arb(n)
        m = any_as_arb(m)
        arb_hypgeom_gegenbauer_c((<arb>v).val, (<arb>n).val, (<arb>m).val, (<arb>s).val, getprec())
        return v

    def laguerre_l(s, n, m=0):
        r"""
        Laguerre function `L_n^{m}(s)`.

            >>> showgood(lambda: (arb(1)/3).laguerre_l(5, 0.25), dps=25)
            0.03871323490012002743484225
        """
        v = arb.__new__(arb)
        n = any_as_arb(n)
        m = any_as_arb(m)
        arb_hypgeom_laguerre_l((<arb>v).val, (<arb>n).val, (<arb>m).val, (<arb>s).val, getprec())
        return v

    def hermite_h(s, n):
        r"""
        Hermite function `H_n(s)`.

            >>> showgood(lambda: (arb(1)/3).hermite_h(5), dps=25)
            34.20576131687242798353909
        """
        v = arb.__new__(arb)
        n = any_as_arb(n)
        arb_hypgeom_hermite_h((<arb>v).val, (<arb>n).val, (<arb>s).val, getprec())
        return v

    def legendre_p(s, n, m=0, int type=2):
        r"""
        Legendre function of the first kind `P_n^m(z)`.

            >>> showgood(lambda: (arb(1)/3).legendre_p(5), dps=25)
            0.3333333333333333333333333
            >>> showgood(lambda: (arb(1)/3).legendre_p(5, 1.5), dps=25)
            -2.372124991643971726805456
            >>> showgood(lambda: (arb(3)).legendre_p(5, 1.5, type=3), dps=25)
            17099.70021476473458984981

        The optional parameter *type* can be 2 or 3, and selects between
        two different branch cut conventions (see *Mathematica* and *mpmath*).
        """
        v = arb.__new__(arb)
        n = any_as_arb(n)
        m = any_as_arb(m)
        if type != 2 and type != 3:
            raise ValueError("type must be 2 or 3")
        type -= 2
        arb_hypgeom_legendre_p((<arb>v).val, (<arb>n).val, (<arb>m).val, (<arb>s).val, type, getprec())
        return v

    def legendre_q(s, n, m=0, int type=2):
        r"""
        Legendre function of the second kind `Q_n^m(z)`.

            >>> showgood(lambda: (arb(1)/3).legendre_q(5), dps=25)
            0.1655245300933242182362054
            >>> showgood(lambda: (arb(1)/3).legendre_q(5, 1.5), dps=25)
            -6.059967350218583975575616

        The optional parameter *type* can be 2 or 3, and selects between
        two different branch cut conventions (see *Mathematica* and *mpmath*).
        """
        v = arb.__new__(arb)
        n = any_as_arb(n)
        m = any_as_arb(m)
        if type != 2 and type != 3:
            raise ValueError("type must be 2 or 3")
        type -= 2
        arb_hypgeom_legendre_q((<arb>v).val, (<arb>n).val, (<arb>m).val, (<arb>s).val, type, getprec())
        return v

    @staticmethod
    def legendre_p_root(ulong n, ulong k, bint weight=False):
        r"""
        Returns the index-*k* zero of the Legendre polynomial `P_n(x)`.
        The zeros are indexed in decreasing order.

        If *weight* is True, returns a tuple (*x*, *w*) where *x* is
        the zero and *w* is the corresponding weight for Gauss-Legendre
        quadrature on `(-1,1)`.

            >>> for k in range(5):
            ...     showgood(lambda: arb.legendre_p_root(5,k), dps=25)
            ...
            0.9061798459386639927976269
            0.5384693101056830910363144
            0
            -0.5384693101056830910363144
            -0.9061798459386639927976269

            >>> for k in range(3):
            ...     showgood(lambda: arb.legendre_p_root(3,k,weight=True), dps=15)
            ...
            (0.774596669241483, 0.555555555555556)
            (0, 0.888888888888889)
            (-0.774596669241483, 0.555555555555556)

        """
        cdef arb x, w
        if k >= n:
            raise ValueError("require k < n")
        if weight:
            x = arb.__new__(arb)
            w = arb.__new__(arb)
            arb_hypgeom_legendre_p_ui_root(x.val, w.val, n, k, getprec())
            return x, w
        else:
            x = arb.__new__(arb)
            arb_hypgeom_legendre_p_ui_root(x.val, NULL, n, k, getprec())
            return x

    def erf(s):
        r"""
        Error function `\operatorname{erf}(s)`.

            >>> showgood(lambda: arb(3).erf(), dps=25)
            0.9999779095030014145586272
        """
        u = arb.__new__(arb)
        arb_hypgeom_erf((<arb>u).val, (<arb>s).val, getprec())
        return u

    def erfc(s):
        r"""
        Complementary error function `\operatorname{erfc}(s)`.

            >>> showgood(lambda: arb(3).erfc(), dps=25)
            2.209049699858544137277613e-5
        """
        u = arb.__new__(arb)
        arb_hypgeom_erfc((<arb>u).val, (<arb>s).val, getprec())
        return u

    def erfi(s):
        r"""
        Imaginary error function `\operatorname{erfi}(s)`.

            >>> showgood(lambda: arb(3).erfi(), dps=25)
            1629.994622601565651061648
        """
        u = arb.__new__(arb)
        arb_hypgeom_erfi((<arb>u).val, (<arb>s).val, getprec())
        return u

    def fresnel_s(s, bint normalized=True):
        r"""
        Fresnel sine integral `S(s)`, optionally not normalized.

            >>> showgood(lambda: arb(3).fresnel_s(), dps=25)
            0.4963129989673750360976123
            >>> showgood(lambda: arb(3).fresnel_s(normalized=False), dps=25)
            0.7735625268937690171497722
        """
        u = arb.__new__(arb)
        arb_hypgeom_fresnel((<arb>u).val, NULL, (<arb>s).val, normalized, getprec())
        return u

    def fresnel_c(s, bint normalized=True):
        r"""
        Fresnel cosine integral `C(s)`, optionally not normalized.

            >>> showgood(lambda: arb(3).fresnel_c(), dps=25)
            0.6057207892976856295561611
            >>> showgood(lambda: arb(3).fresnel_c(normalized=False), dps=25)
            0.7028635577302687301744099
        """
        u = arb.__new__(arb)
        arb_hypgeom_fresnel(NULL, (<arb>u).val, (<arb>s).val, normalized, getprec())
        return u

    def ei(s):
        r"""
        Exponential integral `\operatorname{Ei}(s)`.

            >>> showgood(lambda: arb(3).ei(), dps=25)
            9.933832570625416558008336
        """
        u = arb.__new__(arb)
        arb_hypgeom_ei((<arb>u).val, (<arb>s).val, getprec())
        return u

    def si(s):
        r"""
        Sine integral `\operatorname{Si}(s)`.

            >>> showgood(lambda: arb(3).si(), dps=25)
            1.848652527999468256397730
        """
        u = arb.__new__(arb)
        arb_hypgeom_si((<arb>u).val, (<arb>s).val, getprec())
        return u

    def ci(s):
        r"""
        Cosine integral `\operatorname{Ci}(s)`.

            >>> showgood(lambda: arb(3).ci(), dps=25)
            0.1196297860080003276264723
        """
        u = arb.__new__(arb)
        arb_hypgeom_ci((<arb>u).val, (<arb>s).val, getprec())
        return u

    def shi(s):
        r"""
        Hyperbolic sine integral `\operatorname{Shi}(s)`.

            >>> showgood(lambda: arb(3).shi(), dps=25)
            4.973440475859806797710418
        """
        u = arb.__new__(arb)
        arb_hypgeom_shi((<arb>u).val, (<arb>s).val, getprec())
        return u

    def chi(s):
        r"""
        Hyperbolic cosine integral `\operatorname{Chi}(s)`.

            >>> showgood(lambda: arb(3).chi(), dps=25)
            4.960392094765609760297918
        """
        u = arb.__new__(arb)
        arb_hypgeom_chi((<arb>u).val, (<arb>s).val, getprec())
        return u

    def li(s, bint offset=False):
        r"""
        Logarithmic integral `\operatorname{li}(s)`, optionally
        the offset logarithmic integral
        `\operatorname{Li}(s) = \operatorname{li}(s) - \operatorname{li}(2)`.

            >>> showgood(lambda: arb(10).li(), dps=25)
            6.165599504787297937522982
            >>> showgood(lambda: arb(10).li(offset=True), dps=25)
            5.120435724669805152678393
        """
        u = arb.__new__(arb)
        arb_hypgeom_li((<arb>u).val, (<arb>s).val, offset, getprec())
        return u

    def bessel_j(self, n):
        r"""
        Bessel function `J_n(z)`, where the argument *z* is given by
        *self* and the order *n* is passed as an extra parameter.

            >>> showgood(lambda: arb(5).bessel_j(1), dps=25)
            -0.3275791375914652220377343
        """
        n = any_as_arb(n)
        u = arb.__new__(arb)
        arb_hypgeom_bessel_j((<arb>u).val, (<arb>n).val, (<arb>self).val, getprec())
        return u

    def bessel_y(self, n):
        r"""
        Bessel function `Y_n(z)`, where the argument *z* is given by
        *self* and the order *n* is passed as an extra parameter.

            >>> showgood(lambda: arb(5).bessel_y(1), dps=25)
            0.1478631433912268448010507
        """
        n = any_as_arb(n)
        u = arb.__new__(arb)
        arb_hypgeom_bessel_y((<arb>u).val, (<arb>n).val, (<arb>self).val, getprec())
        return u

    def bessel_k(self, n, bint scaled=False):
        r"""
        Bessel function `K_n(z)`, where the argument *z* is given by
        *self* and the order *n* is passed as an extra parameter.
        Optionally a scaled Bessel function can be computed.

            >>> showgood(lambda: arb(5).bessel_k(1), dps=25)
            0.004044613445452164208365022
            >>> showgood(lambda: arb(5).bessel_k(1, scaled=True), dps=25)
            0.6002738587883125829360457
        """
        n = any_as_arb(n)
        u = arb.__new__(arb)
        if scaled:
            arb_hypgeom_bessel_k_scaled((<arb>u).val, (<arb>n).val, (<arb>self).val, getprec())
        else:
            arb_hypgeom_bessel_k((<arb>u).val, (<arb>n).val, (<arb>self).val, getprec())
        return u

    def bessel_i(self, n, bint scaled=False):
        r"""
        Bessel function `I_n(z)`, where the argument *z* is given by
        *self* and the order *n* is passed as an extra parameter.
        Optionally a scaled Bessel function can be computed.

            >>> showgood(lambda: arb(5).bessel_i(1), dps=25)
            24.33564214245052719914305
            >>> showgood(lambda: arb(5).bessel_i(1, scaled=True), dps=25)
            0.1639722669445423569261229
        """
        n = any_as_arb(n)
        u = arb.__new__(arb)
        if scaled:
            arb_hypgeom_bessel_i_scaled((<arb>u).val, (<arb>n).val, (<arb>self).val, getprec())
        else:
            arb_hypgeom_bessel_i((<arb>u).val, (<arb>n).val, (<arb>self).val, getprec())
        return u

    def gamma_upper(self, s, int regularized=0):
        r"""
        Upper incomplete gamma function `\Gamma(s,z)`. The argument *z*
        is given by *self* and the order *s* is passed as an extra parameter.
        Optionally the regularized version `Q(s,z)` can be computed.

            >>> showgood(lambda: arb(5).gamma_upper(2.5), dps=25)
            0.1000132513171574223384170
            >>> showgood(lambda: arb(5).gamma_upper(2.5, regularized=1), dps=25)
            0.07523524614651217872207687
        """
        s = any_as_arb(s)
        u = arb.__new__(arb)
        arb_hypgeom_gamma_upper((<arb>u).val, (<arb>s).val, (<arb>self).val, regularized, getprec())
        return u

    def gamma_lower(self, s, int regularized=0):
        r"""
        Lower incomplete gamma function `\gamma(s,z)`. The argument *z*
        is given by *self* and the order *s* is passed as an extra parameter.
        Optionally the regularized versions `P(s,z)` and
        `\gamma^{*}(s,z) = z^{-s} P(s,z)` can be computed.

            >>> showgood(lambda: arb(5).gamma_lower(2.5), dps=25)
            1.229327136861979598135209
            >>> showgood(lambda: arb(5).gamma_lower(2.5, regularized=1), dps=25)
            0.9247647538534878212779231
            >>> showgood(lambda: arb(5).gamma_lower(2.5, regularized=2), dps=25)
            0.01654269482249807489997922
        """
        s = any_as_arb(s)
        u = arb.__new__(arb)
        arb_hypgeom_gamma_lower((<arb>u).val, (<arb>s).val, (<arb>self).val, regularized, getprec())
        return u

    def beta_lower(self, a, b, int regularized=0):
        r"""
        Lower incomplete beta function `B(a,b;z)`. The argument *z*
        is given by *self* and the parameters *a* and *b* are passed
        as extra function arguments.
        Optionally the regularized version of this function can be computed.

            >>> showgood(lambda: arb("0.9").beta_lower(1, 2.5), dps=25)
            0.3987350889359326482672004
            >>> showgood(lambda: arb("0.9").beta_lower(1, 2.5, regularized=True), dps=25)
            0.9968377223398316206680011
        """
        a = any_as_arb(a)
        b = any_as_arb(b)
        u = arb.__new__(arb)
        arb_hypgeom_beta_lower((<arb>u).val, (<arb>a).val, (<arb>b).val, (<arb>self).val, regularized, getprec())
        return u

    def expint(self, s):
        r"""
        Generalized exponential integral `E_s(z)`. The argument *z*
        is given by *self* and the order *s* is passed as an extra parameter.

            >>> showgood(lambda: arb(5).expint(2), dps=25)
            0.0009964690427088381099832386
        """
        s = any_as_arb(s)
        u = arb.__new__(arb)
        arb_hypgeom_expint((<arb>u).val, (<arb>s).val, (<arb>self).val, getprec())
        return u

    def hypgeom(self, a, b, bint regularized=False):
        r"""
        Generalized hypergeometric function `{}_pF_q(a;b;z)`.
        The argument *z* is given by *self* and *a* and *b*
        are additional lists of complex numbers defining the parameters.
        Optionally the regularized hypergeometric function can be
        computed.

            >>> showgood(lambda: arb(5).hypgeom([1,2,3],[5,4.5,6]), dps=25)
            1.301849229178968153998454
            >>> showgood(lambda: arb(5).hypgeom([1,2,3],[5,4.5,6],regularized=True), dps=25)
            3.886189282817193519132054e-5
        """
        cdef long i, p, q, prec
        cdef arb_ptr aa, bb
        a = [any_as_arb(t) for t in a]
        b = [any_as_arb(t) for t in b]
        p = len(a)
        q = len(b)
        aa = <arb_ptr>libc.stdlib.malloc(p * cython.sizeof(arb_struct))
        bb = <arb_ptr>libc.stdlib.malloc(q * cython.sizeof(arb_struct))
        for i in range(p):
            aa[i] = (<arb>(a[i])).val[0]
        for i in range(q):
            bb[i] = (<arb>(b[i])).val[0]
        u = arb.__new__(arb)
        arb_hypgeom_pfq((<arb>u).val, aa, p, bb, q, (<arb>self).val, regularized, getprec())
        libc.stdlib.free(aa)
        libc.stdlib.free(bb)
        return u

    def hypgeom_u(self, a, b):
        r"""
        The hypergeometric function `U(a,b,z)` where the argument *z*
        is given by *self*

            >>> showgood(lambda: arb(5).hypgeom_u(1, 2.5), dps=25)
            0.2184157028890778783289036
        """
        a = any_as_arb(a)
        b = any_as_arb(b)
        u = arb.__new__(arb)
        arb_hypgeom_u((<arb>u).val, (<arb>a).val, (<arb>b).val, (<arb>self).val, getprec())
        return u

    def hypgeom_1f1(self, a, b, bint regularized=False):
        r"""
        The hypergeometric function `{}_1F_1(a,b,z)` where the argument *z*
        is given by *self*
        Optionally the regularized version of this function can be computed.

            >>> showgood(lambda: arb(-5).hypgeom_1f1(1, 2.5), dps=25)
            0.2652884733179767675077230
            >>> showgood(lambda: arb(-5).hypgeom_1f1(1, 2.5, regularized=True), dps=25)
            0.1995639910417191573048664
        """
        a = any_as_arb(a)
        b = any_as_arb(b)
        u = arb.__new__(arb)
        arb_hypgeom_1f1((<arb>u).val, (<arb>a).val, (<arb>b).val, (<arb>self).val, regularized, getprec())
        return u

    def hypgeom_0f1(self, a, bint regularized=False):
        r"""
        The hypergeometric function `{}_0F_1(a,z)` where the argument *z*
        is given by *self*
        Optionally the regularized version of this function can be computed.

            >>> showgood(lambda: arb(-5).hypgeom_0f1(2.5), dps=25)
            0.003114611044402738470826907
            >>> showgood(lambda: arb(-5).hypgeom_0f1(2.5, regularized=True), dps=25)
            0.002342974810739764377177885
        """
        a = any_as_arb(a)
        u = arb.__new__(arb)
        arb_hypgeom_0f1((<arb>u).val, (<arb>a).val, (<arb>self).val, regularized, getprec())
        return u

    def hypgeom_2f1(self, a, b, c, bint regularized=False, bint ab=False, bint ac=False, bc=False, abc=False):
        r"""
        The hypergeometric function `{}_2F_1(a,b,c,z)` where the argument *z*
        is given by *self*
        Optionally the regularized version of this function can be computed.

            >>> showgood(lambda: arb(-5).hypgeom_2f1(1,2,3), dps=25)
            0.2566592424617555999350018
            >>> showgood(lambda: arb(-5).hypgeom_2f1(1,2,3,regularized=True), dps=25)
            0.1283296212308777999675009

        The flags *ab*, *ac*, *bc*, *abc* can be used to specify whether the
        parameter differences `a-b`, `a-c`, `b-c` and `a+b-c` represent
        exact integers, even if the input intervals are inexact.
        If the parameters are exact, these flags are not needed.

            >>> showgood(lambda: arb("9/10").hypgeom_2f1(arb(2).sqrt(), 0.5, arb(2).sqrt()+1.5, abc=True), dps=25)
            1.447530478120770807945697
            >>> showgood(lambda: arb("9/10").hypgeom_2f1(arb(2).sqrt(), 0.5, arb(2).sqrt()+1.5), dps=25)
            Traceback (most recent call last):
              ...
            ValueError: no convergence (maxprec=960, try higher maxprec)
        """
        cdef int flags
        a = any_as_arb(a)
        b = any_as_arb(b)
        c = any_as_arb(c)
        u = arb.__new__(arb)
        flags = 0
        if regularized: flags |= 1
        if ab: flags |= 2
        if ac: flags |= 4
        if bc: flags |= 8
        if abc: flags |= 16
        arb_hypgeom_2f1((<arb>u).val, (<arb>a).val, (<arb>b).val, (<arb>c).val,
            (<arb>self).val, flags, getprec())
        return u

    @staticmethod
    def pi():
        """
        Returns the constant `\pi` as an *arb*.

            >>> showgood(arb.pi, dps=25)
            3.141592653589793238462643
        """
        u = arb.__new__(arb)
        arb_const_pi((<arb>u).val, getprec())
        return u

    @staticmethod
    def const_sqrt_pi():
        """
        The constant `\sqrt{\pi}`.

            >>> showgood(arb.const_sqrt_pi, dps=25)
            1.772453850905516027298167
        """
        u = arb.__new__(arb)
        arb_const_sqrt_pi((<arb>u).val, getprec())
        return u

    @staticmethod
    def const_log2():
        """
        The constant `\log(2)`.

            >>> showgood(arb.const_log2, dps=25)
            0.6931471805599453094172321
        """
        u = arb.__new__(arb)
        arb_const_log2((<arb>u).val, getprec())
        return u

    @staticmethod
    def const_log10():
        """
        The constant `\log(10)`.

            >>> showgood(arb.const_log10, dps=25)
            2.302585092994045684017991
        """
        u = arb.__new__(arb)
        arb_const_log10((<arb>u).val, getprec())
        return u

    @staticmethod
    def const_euler():
        """
        Euler's constant `\gamma`.

            >>> showgood(arb.const_euler, dps=25)
            0.5772156649015328606065121
        """
        u = arb.__new__(arb)
        arb_const_euler((<arb>u).val, getprec())
        return u

    @staticmethod
    def const_catalan():
        """
        Catalan's constant.

            >>> showgood(arb.const_catalan, dps=25)
            0.9159655941772190150546035
        """
        u = arb.__new__(arb)
        arb_const_catalan((<arb>u).val, getprec())
        return u

    @staticmethod
    def const_e():
        """
        The constant `e`.

            >>> showgood(arb.const_e, dps=25)
            2.718281828459045235360287
        """
        u = arb.__new__(arb)
        arb_const_e((<arb>u).val, getprec())
        return u

    @staticmethod
    def const_khinchin():
        """
        Khinchin's constant `K`.

            >>> showgood(arb.const_khinchin, dps=25)
            2.685452001065306445309715
        """
        u = arb.__new__(arb)
        arb_const_khinchin((<arb>u).val, getprec())
        return u

    @staticmethod
    def const_glaisher():
        """
        Glaisher's constant.

            >>> showgood(arb.const_glaisher, dps=25)
            1.282427129100622636875343
        """
        u = arb.__new__(arb)
        arb_const_glaisher((<arb>u).val, getprec())
        return u

    @staticmethod
    def pos_inf():
        u = arb.__new__(arb)
        arb_pos_inf((<arb>u).val)
        return u

    @staticmethod
    def neg_inf():
        u = arb.__new__(arb)
        arb_neg_inf((<arb>u).val)
        return u

    @staticmethod
    def nan():
        u = arb.__new__(arb)
        arb_indeterminate((<arb>u).val)
        return u

    def unique_fmpz(self):
        r"""
        If *self* represents exactly one integer, returns this value
        as an *fmpz*; otherwise returns *None*.

            >>> arb("5 +/- 0.1").unique_fmpz()
            5
            >>> arb("5 +/- 0.9").unique_fmpz()
            5
            >>> arb("5.1 +/- 0.9").unique_fmpz()
            >>>

        """
        u = fmpz.__new__(fmpz)
        if arb_get_unique_fmpz((<fmpz>u).val, self.val):
            return u
        else:
            return None

    def rel_accuracy_bits(self):
        return arb_rel_accuracy_bits(self.val)

    def lambertw(s, int branch=0):
        r"""
        Lambert *W* function, `W_k(s)`. Either the principal
        branch (`k = 0`) or the `k = -1` branch can be computed
        by specifying *branch*.

            >>> showgood(lambda: arb(10).lambertw(), dps=25)
            1.745528002740699383074301
            >>> showgood(lambda: arb("-0.1").lambertw(-1), dps=25)
            -3.577152063957297218409392
        """
        cdef int flags
        if branch == 0:
            flags = 0
        elif branch == -1:
            flags = 1
        else:
            raise ValueError("invalid branch")
        u = arb.__new__(arb)
        arb_lambertw((<arb>u).val, (<arb>s).val, flags, getprec())
        return u

    def nonnegative_part(self):
        res = arb.__new__(arb)
        arb_set_round((<arb>res).val, (<arb>self).val, getprec())
        arb_nonnegative_part((<arb>res).val, (<arb>res).val)
        return res

    def union(s, t):
        """
        Returns a ball containing the union of *s* and *t*.

            >>> x = arb(3).union(5); x.lower(); x.upper()
            [2.99999999813735 +/- 4.86e-15]
            [5.00000000186265 +/- 4.86e-15]
        """
        v = arb.__new__(arb)
        t = any_as_arb(t)
        arb_union((<arb>v).val, (<arb>s).val, (<arb>t).val, getprec())
        return v

    def intersection(s, t):
        """
        Returns a ball containing the intersection of *s* and *t*.
        If *s* and *t* are non-overlapping, raises ValueError.

            >>> arb("10 +/- 8.001").intersection(arb("0 +/- 2.001"))
            [2.00 +/- 1.01e-3]
            >>> arb(2).intersection(3)
            Traceback (most recent call last):
              ...
            ValueError: empty intersection
        """
        v = arb.__new__(arb)
        t = any_as_arb(t)
        if arb_intersection((<arb>v).val, (<arb>s).val, (<arb>t).val, getprec()):
            return v
        raise ValueError("empty intersection")

    def min(s, t):
        """
        Minimum value of *s* and *t*.

            >>> arb(2).min(3)
            2.00000000000000
            >>> arb(2).min(arb("3 +/- 1.1"))
            [2e+0 +/- 0.101]
        """
        v = arb.__new__(arb)
        t = any_as_arb(t)
        arb_min((<arb>v).val, (<arb>s).val, (<arb>t).val, getprec())
        return v

    def max(s, t):
        """
        Maximum value of *s* and *t*.

            >>> arb(2).max(arb("3 +/- 1.1"))
            [+/- 4.11]
            >>> arb(4).max(arb("3 +/- 1.1"))
            [4e+0 +/- 0.101]
        """
        v = arb.__new__(arb)
        t = any_as_arb(t)
        arb_max((<arb>v).val, (<arb>s).val, (<arb>t).val, getprec())
        return v

    def root(s, ulong n):
        """
        Principal *n*-th root of *s*.

            >>> showgood(lambda: arb(10).root(3), dps=25)
            2.154434690031883721759294
        """
        v = arb.__new__(arb)
        arb_root_ui((<arb>v).val, (<arb>s).val, n, getprec())
        return v

    @staticmethod
    def gram_point(n):
        """
        Returns the *n*-th Gram point.

            >>> showgood(lambda: arb.gram_point(-1), dps=25)
            9.666908056130192141261536
            >>> showgood(lambda: arb.gram_point(0), dps=25)
            17.84559954041086081682634
            >>> showgood(lambda: arb.gram_point(10**30), dps=25)
            9.829776286927442475869051e+28
        """
        n = fmpz(n)
        v = arb.__new__(arb)
        acb_dirichlet_gram_point((<arb>v).val, (<fmpz>n).val, NULL, NULL, getprec())
        return v

    def zeta_nzeros(x):
        """
        Number of zeros of the Riemann zeta function with positive
        imaginary part between 0 and *x*.

            >>> arb("-5").zeta_nzeros()
            0
            >>> arb("14").zeta_nzeros()
            0
            >>> arb("15").zeta_nzeros()
            1.00000000000000
            >>> arb("14.1 +/- 0.1").zeta_nzeros()
            [+/- 1.01]
            >>> arb("100").zeta_nzeros()
            29.0000000000000
            >>> arb("1e6").zeta_nzeros()
            1747146.00000000

        """
        v = arb.__new__(arb)
        acb_dirichlet_zeta_nzeros((<arb>v).val, (<arb>x).val, getprec())
        return v

    def backlund_s(x):
        """
        Backlund *S* function related to the Riemann zeta function.

            >>> showgood(lambda: arb(123).backlund_s(), dps=25)
            0.4757920863536796196115749
        """
        v = arb.__new__(arb)
        acb_dirichlet_backlund_s((<arb>v).val, (<arb>x).val, getprec())
        return v

    def coulomb(self, l, eta):
        r"""
        Computes the Coulomb wave functions `F_{\ell}(\eta,z)`,
        `G_{\ell}(\eta,z)`, where *z* is given by *self*.
        Both function values are computed simultaneously and a tuple
        is returned.

            >>> showgood(lambda: arb(1).coulomb(0.5, 0.25), dps=10)
            (0.4283180781, 1.218454487)
        """
        l = any_as_arb(l)
        eta = any_as_arb(eta)
        F = arb.__new__(arb)
        G = arb.__new__(arb)
        arb_hypgeom_coulomb((<arb>F).val, (<arb>G).val,
                        (<arb>l).val, (<arb>eta).val, (<arb>self).val, getprec())
        return F, G

    def coulomb_f(self, l, eta):
        r"""
        Regular Coulomb wave function `F_{\ell}(\eta,z)` where
        *z* is given by *self*.

            >>> showgood(lambda: arb(1).coulomb_f(0.5, 0.25), dps=25)
            0.4283180781043541845555944
        """
        l = any_as_arb(l)
        eta = any_as_arb(eta)
        F = arb.__new__(arb)
        arb_hypgeom_coulomb((<arb>F).val, NULL,
                        (<arb>l).val, (<arb>eta).val, (<arb>self).val, getprec())
        return F

    def coulomb_g(self, l, eta):
        r"""
        Irregular Coulomb wave function `G_{\ell}(\eta,z)` where
        *z* is given by *self*.

            >>> showgood(lambda: arb(1).coulomb_g(0.5, 0.25), dps=25)
            1.218454487206367973745641
        """
        l = any_as_arb(l)
        eta = any_as_arb(eta)
        G = arb.__new__(arb)
        arb_hypgeom_coulomb(NULL, (<arb>G).val,
                        (<arb>l).val, (<arb>eta).val, (<arb>self).val, getprec())
        return G

