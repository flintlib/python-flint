def _str_trunc(s, trunc=0):
    if trunc > 0 and len(s) > 3 * trunc:
        left = right = trunc
        omitted = len(s) - left - right
        return s[:left] + ("{...%s digits...}" % omitted) + s[-right:]
    return s

def arb_from_str(str s):
    s = s.strip()
    if ("/" in s) and ("+/-" not in s):
        return arb(fmpq(s))
    s = s.replace("±", "+/-")
    a = arb.__new__(arb)
    if arb_set_str((<arb>a).val, s, getprec()) == 0:
        return a
    else:
        raise ValueError("invalid string for arb()")

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

    if PyInt_Check(<PyObject*>obj):
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

def any_as_arb(x):
    if typecheck(x, arb):
        return x
    return arb(x)

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
        arb((41, -2))
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

    def mid(self):
        cdef arf x = arf()
        arf_set(x.val, arb_midref(self.val))
        return x

    def rad(self):
        cdef arf x = arf()
        arf_set_mag(x.val, arb_radref(self.val))
        return x

    def mid_rad_10exp(self, long n=0):
        """
        Returns an *fmpz* triple (*mid*, *rad*, *exp*) where the larger of *mid*
        and *rad* has *n* digits plus a few digits (*n* defaults to the current
        precision), such that *self* is contained in
        `[\operatorname{mid} \pm \operatorname{rad}] 10^{\operatorname{exp}}`.

            >>> (arb(1) / 3).mid_rad_10exp(10)
            (fmpz(333333333333333), fmpz(2), fmpz(-15))
            >>> (arb(1) / 3).mid_rad_10exp(20)
            (fmpz(3333333333333333148296162), fmpz(555111516), fmpz(-25))
            >>> arb(0,1e-100).mid_rad_10exp(10)
            (fmpz(0), fmpz(100000000507904), fmpz(-114))
            >>> arb(-1,1e100).mid_rad_10exp()
            (fmpz(0), fmpz(10000000169485008897), fmpz(81))

        """
        cdef fmpz mid, rad, exp
        if n <= 0:
            n = ctx.dps
        mid = fmpz()
        rad = fmpz()
        exp = fmpz()
        arb_get_fmpz_mid_rad_10exp(mid.val, rad.val, exp.val, self.val, n)
        return mid, rad, exp

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

            >>> print arb.pi().str()
            [3.14159265358979 +/- 3.57e-15]
            >>> print arb.pi().str(5)
            [3.1416 +/- 7.35e-6]
            >>> print arb.pi().str(5, radius=False)
            3.1416

        By default, the output is truncated so that all displayed digits
        are guaranteed to be correct, up to adding or subtracting 1 in the
        last displayed digit (as a special case, if the output ends with a
        string of 0s, the correct decimal expansion to infinite precision
        could have a string of 9s).

            >>> print (arb(1) - arb("1e-10")).str(5)
            [1.0000 +/- 4e-10]
            >>> print (arb(1) - arb("1e-10")).str(10)
            [0.9999999999 +/- 3e-15]

        To force more digits, set *more* to *True*.

            >>> print arb("0.1").str(30)
            [0.100000000000000 +/- 2.23e-17]
            >>> print arb("0.1").str(30, more=True)
            [0.0999999999999999916733273153113 +/- 1.39e-17]

        Note that setting *more* to *True* results in a smaller printed radius,
        since there is less error from the conversion back to decimal.

            >>> x = arb.pi().sin()
            >>> print(x.str())
            [+/- 5.68e-16]
            >>> print(x.str(more=True))
            [1.22460635382238e-16 +/- 4.45e-16]

        The error indicated in the output may be much larger than the actual
        error in the internal representation of *self*. For example, if *self*
        is accurate to 1000 digits and printing is done at 10-digit precision,
        the output might only appear as being accurate to 10 digits. It is
        even possible for *self* to be exact and have an inexact decimal
        representation.

        The *condense* option can be passed to show only leading and trailing
        digits of the fractional, integer and exponent parts of the output.

            >>> ctx.dps = 1000
            >>> print arb.pi().str(condense=10)
            [3.1415926535{...979 digits...}9216420199 +/- 1.17e-1000]
            >>> print arb.fac_ui(300).str(condense=10)
            3060575122{...595 digits...}0000000000.0000000000{...365 digits...}0000000000
            >>> print arb(10**100).exp().str(condense=5)
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
            res = str(s)
        finally:
            libc.stdlib.free(s)
        if ctx.unicode:
            res = res.replace("+/-", "±")
        return res

    def __float__(self):
        return arf_get_d(arb_midref(self.val), ARF_RND_DOWN)

    def __richcmp__(s, t, int op):
        cdef bint res
        d = s - t
        if not typecheck(s, arb):
            return NotImplemented
        res = 0
        if   op == 2: res = arb_is_zero((<arb>d).val)
        elif op == 3: res = arb_is_nonzero((<arb>d).val)
        elif op == 0: res = arb_is_negative((<arb>d).val)
        elif op == 1: res = arb_is_nonpositive((<arb>d).val)
        elif op == 4: res = arb_is_positive((<arb>d).val)
        elif op == 5: res = arb_is_nonnegative((<arb>d).val)
        return res

    def __contains__(self, other):
        other = any_as_arb(other)
        return arb_contains(self.val, (<arb>other).val)

    def overlaps(self, other):
        other = any_as_arb(other)
        return bool(arb_overlaps((<arb>self).val, (<arb>other).val))

    def __pos__(self):
        res = arb.__new__(arb)
        arb_set_round((<arb>res).val, (<arb>self).val, getprec())
        return res

    def __neg__(self):
        res = arb.__new__(arb)
        arb_neg_round((<arb>res).val, (<arb>self).val, getprec())
        return res

    def __abs__(self):
        res = arb.__new__(arb)
        arb_abs((<arb>res).val, (<arb>self).val)
        arb_set_round((<arb>res).val, (<arb>res).val, getprec())
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

    def __div__(s, t):
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
        Computes the floor function `\lfloor s \rfloor`.

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
        Computes the ceiling function `\lceil s \rceil`.

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
        Computes the square root `\sqrt{s}`.

            >>> showgood(lambda: arb(3).sqrt(), dps=25)
            1.732050807568877293527446
            >>> showgood(lambda: arb(0).sqrt(), dps=25)
            0
            >>> showgood(lambda: arb(-1).sqrt(), dps=25)
            Traceback (most recent call last):
              ...
            ValueError: no convergence

        This function is undefined for negative input.
        Use :meth:`.acb.sqrt` for the complex extension.
        """
        u = arb.__new__(arb)
        arb_sqrt((<arb>u).val, (<arb>s).val, getprec())
        return u

    def rsqrt(s):
        r"""
        Computes the reciprocal square root `1/\sqrt{s}`.

            >>> showgood(lambda: arb(3).rsqrt(), dps=25)
            0.5773502691896257645091488
            >>> showgood(lambda: arb(0).rsqrt(), dps=25)
            Traceback (most recent call last):
              ...
            ValueError: no convergence
            >>> showgood(lambda: arb(-1).rsqrt(), dps=25)
            Traceback (most recent call last):
              ...
            ValueError: no convergence

        This function is undefined for negative input.
        Use :meth:`.acb.rsqrt` for the complex extension.
        """
        u = arb.__new__(arb)
        arb_rsqrt((<arb>u).val, (<arb>s).val, getprec())
        return u

    def exp(s):
        r"""
        Computes the exponential function `\exp(s)`.

            >>> showgood(lambda: arb(1).exp(), dps=25)
            2.718281828459045235360287
        """
        u = arb.__new__(arb)
        arb_exp((<arb>u).val, (<arb>s).val, getprec())
        return u

    def expm1(s):
        r"""
        Computes `\exp(s) - 1`, accurately for small *s*.

            >>> showgood(lambda: (arb(10) ** -8).expm1(), dps=25)
            1.000000005000000016666667e-8
        """
        u = arb.__new__(arb)
        arb_expm1((<arb>u).val, (<arb>s).val, getprec())
        return u

    def log(s):
        r"""
        Computes the natural logarithm `\log(s)`.

            >>> showgood(lambda: arb(2).log(), dps=25)
            0.6931471805599453094172321
            >>> showgood(lambda: arb(100).exp().log(), dps=25)
            100.0000000000000000000000
            >>> showgood(lambda: arb(-1).sqrt(), dps=25)
            Traceback (most recent call last):
              ...
            ValueError: no convergence

        This function is undefined for negative input.
        Use :meth:`.acb.log` for the complex extension.
        """
        u = arb.__new__(arb)
        arb_log((<arb>u).val, (<arb>s).val, getprec())
        return u

    def log1p(s):
        r"""
        Computes `\log(1+s)`, accurately for small *s*.

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
        Computes the sine `\sin(s)`.

            >>> showgood(lambda: arb(1).sin(), dps=25)
            0.8414709848078965066525023
        """
        u = arb.__new__(arb)
        arb_sin((<arb>u).val, (<arb>s).val, getprec())
        return u

    def cos(s):
        r"""
        Computes the cosine `\cos(s)`.

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
        Computes the sine `\sin(\pi s)`.

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
        Computes the cosine `\cos(\pi s)`.

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
        Computes the tangent `\tan(s)`.

            >>> showgood(lambda: arb(1).tan(), dps=25)
            1.557407724654902230506975
        """
        u = arb.__new__(arb)
        arb_tan((<arb>u).val, (<arb>s).val, getprec())
        return u

    def cot(s):
        r"""
        Computes the cotangent `\cot(s)`.

            >>> showgood(lambda: arb(1).cot(), dps=25)
            0.6420926159343307030064200
        """
        u = arb.__new__(arb)
        arb_cot((<arb>u).val, (<arb>s).val, getprec())
        return u

    def tan_pi(s):
        r"""
        Computes the tangent `\tan(\pi s)`.

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
        Computes the cotangent `\cot(\pi s)`.

            >>> showgood(lambda: arb(0.125).cot_pi(), dps=25)
            2.414213562373095048801689
            >>> showgood(lambda: arb(0.5).cot_pi(), dps=25)
            0
        """
        u = arb.__new__(arb)
        arb_cot_pi((<arb>u).val, (<arb>s).val, getprec())
        return u

    @classmethod
    def sin_pi_fmpq(cls, fmpq s):
        r"""
        Computes the algebraic sine value `\sin(\pi s)`.

            >>> showgood(lambda: arb.sin_pi_fmpq(fmpq(3,4)), dps=25)
            0.7071067811865475244008444
        """
        u = arb.__new__(arb)
        arb_sin_pi_fmpq((<arb>u).val, (<fmpq>s).val, getprec())
        return u

    @classmethod
    def cos_pi_fmpq(cls, fmpq s):
        r"""
        Computes the algebraic cosine value `\cos(\pi s)`.

            >>> showgood(lambda: arb.cos_pi_fmpq(fmpq(3,4)), dps=25)
            -0.7071067811865475244008444
        """
        u = arb.__new__(arb)
        arb_cos_pi_fmpq((<arb>u).val, (<fmpq>s).val, getprec())
        return u

    @classmethod
    def sin_cos_pi_fmpq(cls, fmpq s):
        r"""
        Computes `\sin(\pi s)` and `\cos(\pi s)` simultaneously.

            >>> showgood(lambda: arb.sin_cos_pi_fmpq(fmpq(3,4)), dps=25)
            (0.7071067811865475244008444, -0.7071067811865475244008444)
        """
        u = arb.__new__(arb)
        v = arb.__new__(arb)
        arb_sin_cos_pi_fmpq((<arb>u).val, (<arb>v).val, (<fmpq>s).val, getprec())
        return u, v

    def atan(s):
        r"""
        Computes the inverse tangent `\operatorname{atan}(s)`.

            >>> showgood(lambda: arb(1).atan(), dps=25)
            0.7853981633974483096156608
        """
        u = arb.__new__(arb)
        arb_atan((<arb>u).val, (<arb>s).val, getprec())
        return u

    @classmethod
    def atan2(cls, s, t):
        r"""
        Computes the two-argument inverse tangent `\operatorname{atan2}(s,t)`.

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
        Computes the inverse cosine `\operatorname{acos}(s)`.

            >>> showgood(lambda: arb(0).acos(), dps=25)
            1.570796326794896619231322
        """
        u = arb.__new__(arb)
        arb_acos((<arb>u).val, (<arb>s).val, getprec())
        return u

    def asin(s):
        r"""
        Computes the inverse sine `\operatorname{asin}(s)`.

            >>> showgood(lambda: arb(1).asin(), dps=25)
            1.570796326794896619231322
        """
        u = arb.__new__(arb)
        arb_asin((<arb>u).val, (<arb>s).val, getprec())
        return u

    def atanh(s):
        u = arb.__new__(arb)
        arb_atanh((<arb>u).val, (<arb>s).val, getprec())
        return u

    def asinh(s):
        u = arb.__new__(arb)
        arb_asinh((<arb>u).val, (<arb>s).val, getprec())
        return u

    def acosh(s):
        u = arb.__new__(arb)
        arb_acosh((<arb>u).val, (<arb>s).val, getprec())
        return u

    def sinh(s):
        r"""
        Computes the hyperbolic sine `\sinh(s)`.

            >>> showgood(lambda: arb(1).sinh(), dps=25)
            1.175201193643801456882382
        """
        u = arb.__new__(arb)
        arb_sinh((<arb>u).val, (<arb>s).val, getprec())
        return u

    def cosh(s):
        r"""
        Computes the hyperbolic cosine `\cosh(s)`.

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
        Computes the hyperbolic tangent `\tanh(s)`.

            >>> showgood(lambda: arb(1).tanh(), dps=25)
            0.7615941559557648881194583
        """
        u = arb.__new__(arb)
        arb_tanh((<arb>u).val, (<arb>s).val, getprec())
        return u

    def coth(s):
        r"""
        Computes the hyperbolic cotangent `\coth(s)`.

            >>> showgood(lambda: arb(1).coth(), dps=25)
            1.313035285499331303636161
        """
        u = arb.__new__(arb)
        arb_coth((<arb>u).val, (<arb>s).val, getprec())
        return u

    def gamma(s):
        """
        Computes the gamma function `\Gamma(s)`.

            >>> showgood(lambda: arb(10).gamma(), dps=25)
            362880.0000000000000000000
            >>> showgood(lambda: arb(-2.5).gamma(), dps=25)
            -0.9453087204829418812256893
            >>> showgood(lambda: (arb.pi() ** 10).gamma(), dps=25)
            1.705646271897306403570389e+424898
            >>> showgood(lambda: arb(0).gamma(), dps=25)  # pole
            Traceback (most recent call last):
              ...
            ValueError: no convergence
        """
        u = arb.__new__(arb)
        arb_gamma((<arb>u).val, (<arb>s).val, getprec())
        return u

    @classmethod
    def gamma_fmpq(cls, fmpq s):
        """
        Computes the gamma function `\Gamma(s)`, exploiting the fact that
        *s* is a rational number to improve performance.

            >>> showgood(lambda: arb.gamma_fmpq(fmpq(1,4)), dps=25)
            3.625609908221908311930685
        """
        u = arb.__new__(arb)
        arb_gamma_fmpq((<arb>u).val, (<fmpq>s).val, getprec())
        return u

    def rgamma(s):
        ur"""
        Computes the reciprocal gamma function `1/\Gamma(s)`, avoiding
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
        Computes the logarithmic gamma function `\log \Gamma(s)`.

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
        Computes the digamma function `\psi(s)`.

            >>> showgood(lambda: arb(1).digamma(), dps=25)
            -0.5772156649015328606065121
        """
        u = arb.__new__(arb)
        arb_digamma((<arb>u).val, (<arb>s).val, getprec())
        return u

    def rising_ui(s, ulong n):
        """
        Computes the rising factorial `(s)_n` where *n* is an unsigned
        integer. The current implementation does not use the gamma function,
        so *n* should be moderate.

            >>> showgood(lambda: arb.pi().rising_ui(0), dps=25)
            1.000000000000000000000000
            >>> showgood(lambda: arb.pi().rising_ui(10), dps=25)
            299606572.3661012684972888
        """
        u = arb.__new__(arb)
        arb_rising_ui((<arb>u).val, (<arb>s).val, n, getprec())
        return u

    @classmethod
    def rising_fmpq_ui(cls, fmpq s, ulong n):
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

    def rising2_ui(s, ulong n):
        """
        Computes the rising factorial `(s)_n` where *n* is an unsigned
        integer, along with the first derivative with respect to `(s)_n`.
        The current implementation does not use the gamma function,
        so *n* should be moderate.

            >>> u, v = arb(3).rising2_ui(5)
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
        Computes the Riemann zeta function `\zeta(s)` or the Hurwitz
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
        Computes the arithmetic-geometric mean `M(s,t)`, or `M(s) = M(s,1)`
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
            ValueError: no convergence

        This function is undefined for negative input.
        Use :meth:`.acb.agm` for the complex extension.
        """
        t = any_as_arb(t)
        u = arb.__new__(arb)
        arb_agm((<arb>u).val, (<arb>s).val, (<arb>t).val, getprec())
        return u

    @classmethod
    def zeta_ui(cls, ulong n):
        """
        Computes the Riemann zeta function value `\zeta(n)`. This is
        faster than :meth:`.arb.zeta`.

            >>> showgood(lambda: arb.zeta_ui(2), dps=25)
            1.644934066848226436472415
            >>> showgood(lambda: arb.zeta_ui(500) - 1, dps=25)
            3.054936363499604682051979e-151
            >>> showgood(lambda: arb.zeta_ui(1))  # pole
            Traceback (most recent call last):
              ...
            ValueError: no convergence
        """
        u = arb.__new__(arb)
        arb_zeta_ui((<arb>u).val, n, getprec())
        return u

    @classmethod
    def bernoulli_ui(cls, ulong n):
        """
        Computes the Bernoulli number `B_n`.

            >>> showgood(lambda: arb.bernoulli_ui(1), dps=25)
            -0.5000000000000000000000000
            >>> showgood(lambda: arb.bernoulli_ui(2), dps=25)
            0.1666666666666666666666667
            >>> showgood(lambda: arb.bernoulli_ui(10**7), dps=25)
            -4.983176441416329828329241e+57675260
        """
        u = arb.__new__(arb)
        arb_bernoulli_ui((<arb>u).val, n, getprec())
        return u

    @classmethod
    def fac_ui(cls, ulong n):
        """
        Computes the factorial `n!`.

            >>> print(arb.fac_ui(10))
            3628800.00000000
            >>> showgood(lambda: arb.fac_ui(10**9).log(), dps=25)
            19723265848.22698260792313
        """
        u = arb.__new__(arb)
        arb_fac_ui((<arb>u).val, n, getprec())
        return u

    def bin_ui(s, ulong k):
        """
        Computes the binomial coefficient `{s \choose k}`.
        The current implementation does not use the gamma function,
        so *k* should be moderate.

            >>> print(arb(10).bin_ui(5))
            252.000000000000
            >>> showgood(lambda: arb.pi().bin_ui(100), dps=25)
            5.478392395095119521549286e-9
        """
        u = arb.__new__(arb)
        arb_bin_ui((<arb>u).val, (<arb>s).val, k, getprec())
        return u

    @classmethod
    def bin_uiui(cls, ulong n, ulong k):
        """
        Computes the binomial coefficient `{n \choose k}`.
        The current implementation does not use the gamma function,
        so *k* should be moderate.

            >>> print(arb.bin_uiui(10, 5))
            252.000000000000
        """
        u = arb.__new__(arb)
        arb_bin_uiui((<arb>u).val, n, k, getprec())
        return u

    @classmethod
    def fib_ui(cls, ulong n):
        """
        Computes the Fibonacci number `F_n`.

            >>> print(arb.fib_ui(10))
            55.0000000000000
        """
        u = arb.__new__(arb)
        arb_fib_ui((<arb>u).val, n, getprec())
        return u

    @classmethod
    def fib_fmpz(cls, fmpz n):
        """
        Computes the Fibonacci number `F_n`, where *n* may be a bignum.

            >>> print(arb.fib_fmpz(fmpz(10)))
            55.0000000000000
            >>> showgood(lambda: arb.fib_fmpz(fmpz(10**100)).log(), dps=25)
            4.812118250596034474977589e+99
        """
        u = arb.__new__(arb)
        arb_fib_fmpz((<arb>u).val, (<fmpz>n).val, getprec())
        return u

    @classmethod
    def polylog(cls, arb s, arb z):
        """
        Computes the polylogarithm `\operatorname{Li}_s(z)`.

            >>> showgood(lambda: arb.polylog(arb(2), arb(-1)), dps=25)
            -0.8224670334241132182362076
            >>> showgood(lambda: arb.polylog(arb(1.75), arb(-3)), dps=25)
            -1.813689945878060161612620
            >>> showgood(lambda: arb.polylog(arb(4.75), arb(-2.5)), dps=25)
            -2.322090601785704585092044

        This function is undefined for some
        combinations of `s, z`. Use :meth:`.acb.polylog` for the complex extension.
        """
        u = arb.__new__(arb)
        arb_polylog((<arb>u).val, (<arb>s).val, (<arb>z).val, getprec())
        return u

    def chebyshev_t_ui(s, ulong n, bint pair=False):
        """
        Computes the Chebyshev polynomial of the first kind `T_n(s)`.
        If *pair* is True, returns the pair `(T_n(s), T_{n-1}(s))`.

            >>> showgood(lambda: (arb(1)/3).chebyshev_t_ui(3), dps=25)
            -0.8518518518518518518518519
            >>> showgood(lambda: (arb(1)/3).chebyshev_t_ui(4), dps=25)
            0.2098765432098765432098765
            >>> showgood(lambda: (arb(1)/3).chebyshev_t_ui(4, pair=True), dps=25)
            (0.2098765432098765432098765, -0.8518518518518518518518519)
        """
        if pair:
            u = arb.__new__(arb)
            v = arb.__new__(arb)
            arb_chebyshev_t2_ui((<arb>u).val, (<arb>v).val, n, (<arb>s).val, getprec())
            return u, v
        else:
            u = arb.__new__(arb)
            arb_chebyshev_t_ui((<arb>u).val, n, (<arb>s).val, getprec())
            return u

    def chebyshev_u_ui(s, ulong n, bint pair=False):
        """
        Computes the Chebyshev polynomial of the second kind `U_n(s)`.
        If *pair* is True, returns the pair `(U_n(s), U_{n-1}(s))`.

            >>> showgood(lambda: (arb(1)/3).chebyshev_u_ui(3), dps=25)
            -1.037037037037037037037037
            >>> showgood(lambda: (arb(1)/3).chebyshev_u_ui(4), dps=25)
            -0.1358024691358024691358025
            >>> showgood(lambda: (arb(1)/3).chebyshev_u_ui(4, pair=True), dps=25)
            (-0.1358024691358024691358025, -1.037037037037037037037037)
        """
        if pair:
            u = arb.__new__(arb)
            v = arb.__new__(arb)
            arb_chebyshev_u2_ui((<arb>u).val, (<arb>v).val, n, (<arb>s).val, getprec())
            return u, v
        else:
            u = arb.__new__(arb)
            arb_chebyshev_u_ui((<arb>u).val, n, (<arb>s).val, getprec())
            return u

    @classmethod
    def const_pi(cls):
        """
        Computes the constant `\pi`.

            >>> showgood(arb.const_pi, dps=25)
            3.141592653589793238462643
            >>> showgood(arb.pi, dps=25)    # alias
            3.141592653589793238462643
        """
        u = arb.__new__(arb)
        arb_const_pi((<arb>u).val, getprec())
        return u

    pi = const_pi

    @classmethod
    def const_sqrt_pi(cls):
        """
        Computes the constant `\sqrt{\pi}`.

            >>> showgood(arb.const_sqrt_pi, dps=25)
            1.772453850905516027298167
        """
        u = arb.__new__(arb)
        arb_const_sqrt_pi((<arb>u).val, getprec())
        return u

    @classmethod
    def const_log2(cls):
        """
        Computes the constant `\log(2)`.

            >>> showgood(arb.const_log2, dps=25)
            0.6931471805599453094172321
        """
        u = arb.__new__(arb)
        arb_const_log2((<arb>u).val, getprec())
        return u

    @classmethod
    def const_log10(cls):
        """
        Computes the constant `\log(10)`.

            >>> showgood(arb.const_log10, dps=25)
            2.302585092994045684017991
        """
        u = arb.__new__(arb)
        arb_const_log10((<arb>u).val, getprec())
        return u

    @classmethod
    def const_euler(cls):
        """
        Computes Euler's constant `\gamma`.

            >>> showgood(arb.const_euler, dps=25)
            0.5772156649015328606065121
        """
        u = arb.__new__(arb)
        arb_const_euler((<arb>u).val, getprec())
        return u

    @classmethod
    def const_catalan(cls):
        """
        Computes Catalan's constant.

            >>> showgood(arb.const_catalan, dps=25)
            0.9159655941772190150546035
        """
        u = arb.__new__(arb)
        arb_const_catalan((<arb>u).val, getprec())
        return u

    @classmethod
    def const_e(cls):
        """
        Computes the constant `e`.

            >>> showgood(arb.const_e, dps=25)
            2.718281828459045235360287
        """
        u = arb.__new__(arb)
        arb_const_e((<arb>u).val, getprec())
        return u

    @classmethod
    def const_khinchin(cls):
        """
        Computes Khinchin's constant `K`.

            >>> showgood(arb.const_khinchin, dps=25)
            2.685452001065306445309715
        """
        u = arb.__new__(arb)
        arb_const_khinchin((<arb>u).val, getprec())
        return u

    @classmethod
    def const_glaisher(cls):
        """
        Computes Glaisher's constant.

            >>> showgood(arb.const_glaisher, dps=25)
            1.282427129100622636875343
        """
        u = arb.__new__(arb)
        arb_const_glaisher((<arb>u).val, getprec())
        return u

    @classmethod
    def pos_inf(cls):
        u = arb.__new__(arb)
        arb_pos_inf((<arb>u).val)
        return u

    @classmethod
    def neg_inf(cls):
        u = arb.__new__(arb)
        arb_neg_inf((<arb>u).val)
        return u

    @classmethod
    def nan(cls):
        u = arb.__new__(arb)
        arb_indeterminate((<arb>u).val)
        return u

    def unique_fmpz(self):
        u = fmpz.__new__(fmpz)
        if arb_get_unique_fmpz((<fmpz>u).val, self.val):
            return u
        else:
            return None

    def rel_accuracy_bits(self):
        return arb_rel_accuracy_bits(self.val)

