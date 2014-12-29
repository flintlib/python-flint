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
        u = arb.__new__(arb)
        arb_floor((<arb>u).val, (<arb>s).val, getprec())
        return u

    def ceil(s):
        u = arb.__new__(arb)
        arb_ceil((<arb>u).val, (<arb>s).val, getprec())
        return u

    def sqrt(s):
        u = arb.__new__(arb)
        arb_sqrt((<arb>u).val, (<arb>s).val, getprec())
        return u

    def rsqrt(s):
        u = arb.__new__(arb)
        arb_rsqrt((<arb>u).val, (<arb>s).val, getprec())
        return u

    def exp(s):
        u = arb.__new__(arb)
        arb_exp((<arb>u).val, (<arb>s).val, getprec())
        return u

    def expm1(s):
        u = arb.__new__(arb)
        arb_expm1((<arb>u).val, (<arb>s).val, getprec())
        return u

    def log(s):
        u = arb.__new__(arb)
        arb_log((<arb>u).val, (<arb>s).val, getprec())
        return u

    def sin(s):
        u = arb.__new__(arb)
        arb_sin((<arb>u).val, (<arb>s).val, getprec())
        return u

    def cos(s):
        u = arb.__new__(arb)
        arb_cos((<arb>u).val, (<arb>s).val, getprec())
        return u

    def sin_cos(s):
        u = arb.__new__(arb)
        v = arb.__new__(arb)
        arb_sin_cos((<arb>u).val, (<arb>v).val, (<arb>s).val, getprec())
        return u, v

    def sin_pi(s):
        u = arb.__new__(arb)
        arb_sin_pi((<arb>u).val, (<arb>s).val, getprec())
        return u

    def cos_pi(s):
        u = arb.__new__(arb)
        arb_cos_pi((<arb>u).val, (<arb>s).val, getprec())
        return u

    def sin_cos_pi(s):
        u = arb.__new__(arb)
        v = arb.__new__(arb)
        arb_sin_cos_pi((<arb>u).val, (<arb>v).val, (<arb>s).val, getprec())
        return u, v

    def tan(s):
        u = arb.__new__(arb)
        arb_tan((<arb>u).val, (<arb>s).val, getprec())
        return u

    def cot(s):
        u = arb.__new__(arb)
        arb_cot((<arb>u).val, (<arb>s).val, getprec())
        return u

    def tan_pi(s):
        u = arb.__new__(arb)
        arb_tan_pi((<arb>u).val, (<arb>s).val, getprec())
        return u

    def cot_pi(s):
        u = arb.__new__(arb)
        arb_cot_pi((<arb>u).val, (<arb>s).val, getprec())
        return u

    @classmethod
    def sin_pi_fmpq(cls, fmpq s):
        u = arb.__new__(arb)
        arb_sin_pi_fmpq((<arb>u).val, (<fmpq>s).val, getprec())
        return u

    @classmethod
    def cos_pi_fmpq(cls, fmpq s):
        u = arb.__new__(arb)
        arb_cos_pi_fmpq((<arb>u).val, (<fmpq>s).val, getprec())
        return u

    @classmethod
    def sin_cos_pi_fmpq(cls, fmpq s):
        u = arb.__new__(arb)
        v = arb.__new__(arb)
        arb_sin_cos_pi_fmpq((<arb>u).val, (<arb>v).val, (<fmpq>s).val, getprec())
        return u, v

    def atan(s):
        u = arb.__new__(arb)
        arb_atan((<arb>u).val, (<arb>s).val, getprec())
        return u

    def atan2(s, t):
        t = any_as_arb(t)
        u = arb.__new__(arb)
        arb_atan2((<arb>u).val, (<arb>s).val, (<arb>t).val, getprec())
        return u

    def acos(s):
        u = arb.__new__(arb)
        arb_acos((<arb>u).val, (<arb>s).val, getprec())
        return u

    def asin(s):
        u = arb.__new__(arb)
        arb_asin((<arb>u).val, (<arb>s).val, getprec())
        return u

    def sinh(s):
        u = arb.__new__(arb)
        arb_sinh((<arb>u).val, (<arb>s).val, getprec())
        return u

    def cosh(s):
        u = arb.__new__(arb)
        arb_sinh((<arb>u).val, (<arb>s).val, getprec())
        return u

    def sinh_cosh(s):
        u = arb.__new__(arb)
        v = arb.__new__(arb)
        arb_sinh_cosh((<arb>u).val, (<arb>v).val, (<arb>s).val, getprec())
        return u, v

    def tanh(s):
        u = arb.__new__(arb)
        arb_tanh((<arb>u).val, (<arb>s).val, getprec())
        return u

    def coth(s):
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
            >>> showgood(lambda: (arb.const_pi() ** 10).gamma(), dps=25)
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
        u = arb.__new__(arb)
        arb_gamma_fmpq((<arb>u).val, (<fmpq>s).val, getprec())
        return u

    def rgamma(s):
        u = arb.__new__(arb)
        arb_rgamma((<arb>u).val, (<arb>s).val, getprec())
        return u

    def lgamma(s):
        u = arb.__new__(arb)
        arb_lgamma((<arb>u).val, (<arb>s).val, getprec())
        return u

    def digamma(s):
        u = arb.__new__(arb)
        arb_digamma((<arb>u).val, (<arb>s).val, getprec())
        return u

    def rising_ui(s, ulong n):
        """
        Computes the rising factorial rf(s,n) where n is an unsigned
        integer. The current implementation does not use the gamma function,
        so n should be moderate.

            >>> showgood(lambda: arb.const_pi().rising_ui(0), dps=25)
            1.0
            >>> showgood(lambda: arb.const_pi().rising_ui(10), dps=25)
            299606572.3661012684972888
        """
        u = arb.__new__(arb)
        arb_rising_ui((<arb>u).val, (<arb>s).val, n, getprec())
        return u

    @classmethod
    def rising_fmpq_ui(cls, fmpq s, ulong n):
        """
        Computes the rising factorial rf(s,n) where s is an fmpq and n is an
        unsigned integer. The current implementation does not use the
        gamma function, so n should be moderate.

            >>> showgood(lambda: arb.rising_fmpq_ui(fmpq(-1,3), 100), dps=25)
            -4.960517984074284420131903e+154
        """
        u = arb.__new__(arb)
        arb_rising_fmpq_ui((<arb>u).val, (<fmpq>s).val, n, getprec())
        return u

    def rising2_ui(s, ulong n):
        """
        Computes the rising factorial rf(s,n) where n is an unsigned
        integer, along with the first derivative with respect to s.
        The current implementation does not use the gamma function,
        so n should be moderate.

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
        u = arb.__new__(arb)
        arb_fac_ui((<arb>u).val, n, getprec())
        return u

    def bin_ui(s, ulong k):
        u = arb.__new__(arb)
        arb_bin_ui((<arb>u).val, (<arb>s).val, k, getprec())
        return u

    @classmethod
    def bin_uiui(cls, ulong n, ulong k):
        u = arb.__new__(arb)
        arb_bin_uiui((<arb>u).val, n, k, getprec())
        return u

    @classmethod
    def fib_ui(cls, ulong n):
        u = arb.__new__(arb)
        arb_fib_ui((<arb>u).val, n, getprec())
        return u

    @classmethod
    def fib_fmpz(cls, fmpz n):
        u = arb.__new__(arb)
        arb_fib_fmpz((<arb>u).val, (<fmpz>n).val, getprec())
        return u

    @classmethod
    def polylog(cls, arb s, arb z):
        u = arb.__new__(arb)
        arb_polylog((<arb>u).val, (<arb>s).val, (<arb>z).val, getprec())
        return u

    @classmethod
    def polylog_si(cls, long s, arb z):
        u = arb.__new__(arb)
        arb_polylog_si((<arb>u).val, s, (<arb>z).val, getprec())
        return u

    def dilog(s):
        u = arb.__new__(arb)
        arb_polylog_si((<arb>u).val, 2, (<arb>s).val, getprec())
        return u

    def chebyshev_t_ui(s, ulong n):
        u = arb.__new__(arb)
        arb_chebyshev_t_ui((<arb>u).val, n, (<arb>s).val, getprec())
        return u

    def chebyshev_t2_ui(s, ulong n):
        u = arb.__new__(arb)
        v = arb.__new__(arb)
        arb_chebyshev_t2_ui((<arb>u).val, (<arb>v).val, n, (<arb>s).val, getprec())
        return u

    def chebyshev_u_ui(s, ulong n):
        u = arb.__new__(arb)
        arb_chebyshev_u_ui((<arb>u).val, n, (<arb>s).val, getprec())
        return u

    def chebyshev_u2_ui(s, ulong n):
        u = arb.__new__(arb)
        v = arb.__new__(arb)
        arb_chebyshev_u2_ui((<arb>u).val, (<arb>v).val, n, (<arb>s).val, getprec())
        return u

    @classmethod
    def atan(cls, arb b, arb a):
        u = arb.__new__(arb)
        arb_atan2((<arb>u).val, (<arb>b).val, (<arb>a).val, getprec())
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

