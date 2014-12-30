cdef int arb_set_python(arb_t x, obj, bint allow_conversion):
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


cdef class arb:
    r"""
    Represents a real number `x` by a midpoint `m` and a radius `r`
    such that `x \in [m \pm r] = [m-r, m+r]`.
    The midpoint and radius are both floating-point numbers. The radius
    uses a fixed, implementation-defined precision (30 bits).
    The precision used for midpoints is controlled by :attr:`ctx.prec` (bits)
    or equivalently :attr:`ctx.dps` (digits).

    The constructor accepts a midpoint *mid* and a radius *rad*, either of
    which defaults to zero if omitted. The arguments can be tuples
    `(a, b)` representing exact floating-point data `a 2^b`.
    If the radius is nonzero, it might be rounded up to a slightly larger
    value than the exact value passed by the user.

        >>> arb(10.25)
        arb((41, -2))
        >>> print(1 / arb(4))  # exact
        0.25
        >>> print(1 / arb(3))  # approximate
        [0.333333333333333 ± 5.55e-17]
        >>> ctx.dps = 50
        >>> print(1 / arb(3))
        [0.33333333333333333333333333333333333333333333333333 ± 6.68e-52]
        >>> ctx.default()

    .. warning::

        Decimal printing introduces rounding error in the last displayed digit.
        This error is not added to the printed radius. For example, `2^{-1000}`
        is represented exactly internally and therefore prints without a radius,
        but the printed decimal value is not exact:

            >>> print(arb((1,-1000)))
            9.33263618503219e-302

        A workaround for rigorous decimal conversion is to compute
        `\lfloor x 10^n \rfloor`, `\lceil x 10^n \rceil`:

            >>> x = arb((1,-1000)) * 10**310
            >>> print(x.floor()); print(x.ceil()); print(x.floor().rad()); print(x.ceil().rad())
            933263618.0
            933263619.0
            0.0
            0.0

        More robust radix conversion functions will be added in the future.
        When in doubt, the power-of-two representation shown by :func:`repr`
        is authoritative.

    """

    cdef arb_t val

    def __cinit__(self):
        arb_init(self.val)

    def __dealloc__(self):
        arb_clear(self.val)

    def __init__(self, mid=None, rad=None):
        if mid is not None:
            if not arb_set_python(self.val, mid, 1):
                raise TypeError("cannot create arb from type %s" % type(mid))
        if rad is not None:
            rad = arf(rad)
            arb_add_error_arf(self.val, (<arf>rad).val)

    def is_zero(self):
        return arb_is_zero(self.val)

    def mid(self):
        cdef arf x = arf()
        arf_set(x.val, arb_midref(self.val))
        return x

    def rad(self):
        cdef arf x = arf()
        arf_set_mag(x.val, arb_radref(self.val))
        return x

    def __repr__(self):
        if ctx.pretty:
            return str(self)
        mid = self.mid()
        rad = self.rad()
        if rad.is_zero():
            return "arb(%s)" % mid._repr_str()
        else:
            return "arb(%s, %s)" % (mid._repr_str(), rad._repr_str())

    def __str__(self):
        mid = self.mid()
        rad = self.rad()
        if rad.is_zero():
            return mid._dec_str()
        else:
            return "[%s ± %s]" % (mid._dec_str(), rad._dec_str(3))

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
        r"""
        Computes the floor function `\lfloor s \rfloor`.

            >>> print(arb.pi().floor())
            3.0
            >>> print((arb.pi() - arb.pi()).floor())
            [-0.5 ± 0.5]
        """
        u = arb.__new__(arb)
        arb_floor((<arb>u).val, (<arb>s).val, getprec())
        return u

    def ceil(s):
        r"""
        Computes the ceiling function `\lceil s \rceil`.

            >>> print(arb.pi().ceil())
            4.0
            >>> print((arb.pi() - arb.pi()).ceil())
            [0.5 ± 0.5]
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
            0.0
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
            100.0
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
            0.0
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
            0.0
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
            0.64209261593433070300642
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
            0.0
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
            0.0
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
            362880.0
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
        """
        Computes the reciprocal gamma function `1/\Gamma(s)`, avoiding
        division by zero at the poles of the gamma function.

            >>> showgood(lambda: arb(1.5).rgamma(), dps=25)
            1.128379167095512573896159
            >>> print(arb(0).rgamma())
            0.0
            >>> print(arb(-1).rgamma())
            0.0
            >>> print(arb(-3,1e-10).rgamma())
            [0.0 ± 6.0e-10]
        """
        u = arb.__new__(arb)
        arb_rgamma((<arb>u).val, (<arb>s).val, getprec())
        return u

    def lgamma(s):
        """
        Computes the logarithmic gamma function `\log \Gamma(s)`.

            >>> showgood(lambda: arb(100).lgamma(), dps=25)
            359.134205369575398776044

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
            1.0
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
            2520.0
            2754.0
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
            1.06695419071121453245037
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
            0.0
            >>> showgood(arb(0).agm, dps=25)
            0.0
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
            -0.5
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
            3628800.0
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
            252.0
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
            252.0
        """
        u = arb.__new__(arb)
        arb_bin_uiui((<arb>u).val, n, k, getprec())
        return u

    @classmethod
    def fib_ui(cls, ulong n):
        """
        Computes the Fibonacci number `F_n`.

            >>> print(arb.fib_ui(10))
            55.0
        """
        u = arb.__new__(arb)
        arb_fib_ui((<arb>u).val, n, getprec())
        return u

    @classmethod
    def fib_fmpz(cls, fmpz n):
        """
        Computes the Fibonacci number `F_n`, where *n* may be a bignum.

            >>> print(arb.fib_fmpz(fmpz(10)))
            55.0
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
            -1.81368994587806016161262
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

