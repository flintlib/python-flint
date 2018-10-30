cdef int acb_set_python(acb_t x, obj, bint allow_conversion):
    cdef double re, im

    if typecheck(obj, acb):
        acb_set(x, (<acb>obj).val)
        return 1

    if arb_set_python(acb_realref(x), obj, allow_conversion):
        arb_zero(acb_imagref(x))
        return 1

    if PyComplex_Check(<PyObject*>obj):
        re = PyComplex_RealAsDouble(<PyObject*>obj)
        im = PyComplex_ImagAsDouble(<PyObject*>obj)
        arf_set_d(arb_midref(acb_realref(x)), re)
        arf_set_d(arb_midref(acb_imagref(x)), im)
        mag_zero(arb_radref(acb_realref(x)))
        mag_zero(arb_radref(acb_imagref(x)))
        return 1

    if hasattr(obj, "_mpc_"):
        xre, xim = obj._mpc_
        arb_set_mpmath_mpf(acb_realref(x), xre)
        arb_set_mpmath_mpf(acb_imagref(x), xim)
        return 1

    return 0

cdef inline int acb_set_any_ref(acb_t x, obj):
    if typecheck(obj, acb):  # should be exact check?
        x[0] = (<acb>obj).val[0]
        return FMPZ_REF

    acb_init(x)
    if acb_set_python(x, obj, 0):
        return FMPZ_TMP

    return FMPZ_UNKNOWN

cdef any_as_acb(x):
    cdef acb t
    if typecheck(x, acb):
        return x
    t = acb()
    if acb_set_python(t.val, x, 0) == 0:
        raise TypeError("cannot create acb from type %s" % type(x))
    return t

cdef any_as_acb_or_notimplemented(x):
    cdef acb t
    if typecheck(x, acb):
        return x
    t = acb()
    if acb_set_python(t.val, x, 0) == 0:
        return NotImplemented
    return t

"""
cdef any_as_arb_or_acb(x):
    if typecheck(x, arb) or typecheck(x, acb):
        return x
    try:
        return arb(x)
    except (TypeError, ValueError):
        return acb(x)
"""

# Copied with modifications from sage/rings/complex_arb.pyx
@cython.internal
cdef class IntegrationContext:
    cdef object f
    cdef object exn_type
    cdef object exn_obj
    cdef object exn_tb

cdef int acb_calc_func_callback(acb_ptr out, const acb_t inp, void * param, long order, long prec):
    cdef IntegrationContext ictx
    cdef acb x
    try:
        ictx = <IntegrationContext>param
        if ictx.exn_type is not None or order >= 2:
            acb_indeterminate(out)
            return 0
        x = acb.__new__(acb)
        acb_set(x.val, inp)
        try:
            y = ictx.f(x, (order == 1))
            if not typecheck(y, acb):
                raise TypeError("integrand must return an acb")
            acb_set(out, (<acb> y).val)
        except:
            import sys
            ictx.exn_type, ictx.exn_obj, ictx.exn_tb = sys.exc_info()
            acb_indeterminate(out)
        return 0
    finally:
        pass


cdef class acb(flint_scalar):

    cdef acb_t val

    def __cinit__(self):
        acb_init(self.val)

    def __dealloc__(self):
        acb_clear(self.val)

    def __init__(self, real=None, imag=None):
        if real is not None:
            if not acb_set_python(self.val, real, 1):
                raise TypeError("cannot create acb from type %s" % type(real))
        if imag is not None:
            if not arb_is_zero(acb_imagref(self.val)):
                raise ValueError("must create acb from one complex number or two real numbers")
            if not arb_set_python(acb_imagref(self.val), imag, 1):
                raise TypeError("cannot create arb from type %s" % type(imag))

    @property
    def real(self):
        cdef arb re = arb()
        arb_set(re.val, acb_realref(self.val))
        return re

    @property
    def imag(self):
        cdef arb im = arb()
        arb_set(im.val, acb_imagref(self.val))
        return im

    @property
    def _mpc_(self):
        return (self.real._mpf_, self.imag._mpf_)

    def __richcmp__(s, t, int op):
        cdef acb_struct sval[1]
        cdef acb_struct tval[1]
        cdef bint res
        cdef int stype, ttype
        if not (op == 2 or op == 3):
            raise ValueError("comparing complex numbers")
        stype = acb_set_any_ref(sval, s)
        if stype == FMPZ_UNKNOWN:
            return NotImplemented
        ttype = acb_set_any_ref(tval, t)
        if ttype == FMPZ_UNKNOWN:
            return NotImplemented
        if op == 2:
            res = acb_eq(sval, tval)
        else:
            res = acb_ne(sval, tval)
        if stype == FMPZ_TMP: acb_clear(sval)
        if ttype == FMPZ_TMP: acb_clear(tval)
        return res

    def __contains__(self, other):
        other = any_as_acb(other)
        return acb_contains(self.val, (<acb>other).val)

    def contains(self, other):
        other = any_as_acb(other)
        return bool(acb_contains(self.val, (<acb>other).val))

    def overlaps(self, other):
        other = any_as_acb(other)
        return bool(acb_overlaps((<acb>self).val, (<acb>other).val))

    def mid(self):
        """
        Returns an exact *acb* representing the midpoint of *self*:

            >>> acb("1 +/- 0.3", "2 +/- 0.4").mid()
            1.00000000000000 + 2.00000000000000j
        """
        cdef acb u = acb()
        arf_set(arb_midref(acb_realref(u.val)), arb_midref(acb_realref(self.val)))
        arf_set(arb_midref(acb_imagref(u.val)), arb_midref(acb_imagref(self.val)))
        return u

    def rad(self):
        """
        Returns an upper bound for the radius (magnitude of error) of self as an *arb*.

            >>> print(acb("1 +/- 0.3", "2 +/- 0.4").rad().str(5, radius=False))
            0.50000
        """
        cdef arb u = arb()
        mag_hypot(arb_radref(u.val), arb_radref(acb_realref(self.val)), arb_radref(acb_imagref(self.val)))
        arf_set_mag(arb_midref(u.val), arb_radref(u.val))
        mag_zero(arb_radref(u.val))
        return u

    def complex_rad(self):
        """
        Returns an *acb* representing the radii of the real and imaginary parts of *self*
        together a single complex number.

            >>> print(acb("1 +/- 0.3", "2 +/- 0.4").complex_rad().str(5, radius=False))
            0.30000 + 0.40000j
        """
        cdef acb u = acb()
        arf_set_mag(arb_midref(acb_realref(u.val)), arb_radref(acb_realref(self.val)))
        arf_set_mag(arb_midref(acb_imagref(u.val)), arb_radref(acb_imagref(self.val)))
        return u

    def repr(self):
        real = self.real
        imag = self.imag
        if imag.is_zero():
            return "acb(%s)" % real.repr()
        else:
            return "acb(%s, %s)" % (real.repr(), imag.repr())

    def str(self, *args, **kwargs):
        real = self.real
        imag = self.imag
        if imag.is_zero():
            return real.str(*args, **kwargs)
        elif real.is_zero():
            return imag.str(*args, **kwargs) + "j"
        else:
            re = real.str(*args, **kwargs)
            im = imag.str(*args, **kwargs)
            if im.startswith("-"):
                return "%s - %sj" % (re, im[1:])
            else:
                return "%s + %sj" % (re, im)

    def __complex__(self):
        return complex(arf_get_d(arb_midref(acb_realref(self.val)), ARF_RND_NEAR),
                       arf_get_d(arb_midref(acb_imagref(self.val)), ARF_RND_NEAR))

    def __pos__(self):
        res = acb.__new__(acb)
        acb_set_round((<acb>res).val, (<acb>self).val, getprec())
        return res

    def __neg__(self):
        res = acb.__new__(acb)
        acb_neg_round((<acb>res).val, (<acb>self).val, getprec())
        return res

    def __abs__(self):
        res = arb.__new__(arb)
        acb_abs((<arb>res).val, (<acb>self).val, getprec())
        return res

    def abs_lower(self):
        """
        Lower bound for the absolute value of *self*.
        The output is an *arb* holding an exact floating-point number
        that has been rounded down to the current precision.

            >>> print(acb(3, "-5 +/- 2").abs_lower().str(5, radius=False))
            4.2426
        """
        cdef arb x = arb()
        acb_get_abs_lbound_arf(arb_midref(x.val), self.val, getprec())
        return x

    def abs_upper(self):
        """
        Upper bound for the absolute value of *self*.
        The output is an *arb* holding an exact floating-point number
        that has been rounded up to the current precision.

            >>> print(acb(3, "-5 +/- 2").abs_upper().str(5, radius=False))
            7.6158
        """
        cdef arb x = arb()
        acb_get_abs_ubound_arf(arb_midref(x.val), self.val, getprec())
        return x

    def csgn(self):
        res = arb.__new__(arb)
        acb_csgn((<arb>res).val, (<acb>self).val)
        return res

    def sgn(self):
        res = acb.__new__(acb)
        acb_sgn((<acb>res).val, (<acb>self).val, getprec())
        return res

    def arg(self):
        res = arb.__new__(arb)
        acb_arg((<arb>res).val, (<acb>self).val, getprec())
        return res

    def __add__(s, t):
        cdef acb_struct sval[1]
        cdef acb_struct tval[1]
        cdef int stype, ttype
        stype = acb_set_any_ref(sval, s)
        if stype == FMPZ_UNKNOWN:
            return NotImplemented
        ttype = acb_set_any_ref(tval, t)
        if ttype == FMPZ_UNKNOWN:
            return NotImplemented
        u = acb.__new__(acb)
        acb_add((<acb>u).val, sval, tval, getprec())
        if stype == FMPZ_TMP: acb_clear(sval)
        if ttype == FMPZ_TMP: acb_clear(tval)
        return u

    def __sub__(s, t):
        cdef acb_struct sval[1]
        cdef acb_struct tval[1]
        cdef int stype, ttype
        stype = acb_set_any_ref(sval, s)
        if stype == FMPZ_UNKNOWN:
            return NotImplemented
        ttype = acb_set_any_ref(tval, t)
        if ttype == FMPZ_UNKNOWN:
            return NotImplemented
        u = acb.__new__(acb)
        acb_sub((<acb>u).val, sval, tval, getprec())
        if stype == FMPZ_TMP: acb_clear(sval)
        if ttype == FMPZ_TMP: acb_clear(tval)
        return u

    def __mul__(s, t):
        cdef acb_struct sval[1]
        cdef acb_struct tval[1]
        cdef int stype, ttype
        stype = acb_set_any_ref(sval, s)
        if stype == FMPZ_UNKNOWN:
            return NotImplemented
        ttype = acb_set_any_ref(tval, t)
        if ttype == FMPZ_UNKNOWN:
            return NotImplemented
        u = acb.__new__(acb)
        acb_mul((<acb>u).val, sval, tval, getprec())
        if stype == FMPZ_TMP: acb_clear(sval)
        if ttype == FMPZ_TMP: acb_clear(tval)
        return u

    @staticmethod
    cdef _div_(s, t):
        cdef acb_struct sval[1]
        cdef acb_struct tval[1]
        cdef int stype, ttype
        stype = acb_set_any_ref(sval, s)
        if stype == FMPZ_UNKNOWN:
            return NotImplemented
        ttype = acb_set_any_ref(tval, t)
        if ttype == FMPZ_UNKNOWN:
            return NotImplemented
        u = acb.__new__(acb)
        acb_div((<acb>u).val, sval, tval, getprec())
        if stype == FMPZ_TMP: acb_clear(sval)
        if ttype == FMPZ_TMP: acb_clear(tval)
        return u

    def __truediv__(s, t):
        return acb._div_(s, t)

    def __div__(s, t):
        return acb._div_(s, t)

    def __pow__(s, t, u):
        cdef acb_struct sval[1]
        cdef acb_struct tval[1]
        cdef int stype, ttype
        if u is not None:
            raise ValueError("modular exponentiation of complex number")
        stype = acb_set_any_ref(sval, s)
        if stype == FMPZ_UNKNOWN:
            return NotImplemented
        ttype = acb_set_any_ref(tval, t)
        if ttype == FMPZ_UNKNOWN:
            return NotImplemented
        u = acb.__new__(acb)
        acb_pow((<acb>u).val, sval, tval, getprec())
        if stype == FMPZ_TMP: acb_clear(sval)
        if ttype == FMPZ_TMP: acb_clear(tval)
        return u

    def pow(s, t, bint analytic=False):
        t = any_as_acb(t)
        u = acb.__new__(acb)
        acb_pow_analytic((<acb>u).val, (<acb>s).val, (<acb>t).val, analytic, getprec())
        return u

    def log(s, bint analytic=False):
        r"""
        Computes the natural logarithm `\log(s)`.

            >>> showgood(lambda: acb(1,2).log(), dps=25)
            0.8047189562170501873003797 + 1.107148717794090503017065j
            >>> showgood(lambda: acb(-5).log(), dps=25)
            1.609437912434100374600759 + 3.141592653589793238462643j
        """
        u = acb.__new__(acb)
        acb_log_analytic((<acb>u).val, (<acb>s).val, analytic, getprec())
        return u

    def log1p(s):
        r"""
        Computes `\log(1+s)`, accurately for small *s*.

            >>> showgood(lambda: acb(1,2).log1p(), dps=25)
            1.039720770839917964125848 + 0.7853981633974483096156608j
            >>> showgood(lambda: acb(0,"1e-100000000000000000").log1p(), dps=25)
            5.000000000000000000000000e-200000000000000001 + 1.000000000000000000000000e-100000000000000000j
        """
        u = acb.__new__(acb)
        acb_log1p((<acb>u).val, (<acb>s).val, getprec())
        return u

    def asin(s):
        u = acb.__new__(acb)
        acb_asin((<acb>u).val, (<acb>s).val, getprec())
        return u

    def acos(s):
        u = acb.__new__(acb)
        acb_acos((<acb>u).val, (<acb>s).val, getprec())
        return u

    def atan(s):
        r"""
        Computes the inverse tangent `\operatorname{atan}(s)`.

            >>> showgood(lambda: acb(1,2).atan(), dps=25)
            1.338972522294493561124194 + 0.4023594781085250936501898j
        """
        u = acb.__new__(acb)
        acb_atan((<acb>u).val, (<acb>s).val, getprec())
        return u

    def asinh(s):
        u = acb.__new__(acb)
        acb_asinh((<acb>u).val, (<acb>s).val, getprec())
        return u

    def acosh(s):
        u = acb.__new__(acb)
        acb_acosh((<acb>u).val, (<acb>s).val, getprec())
        return u

    def atanh(s):
        u = acb.__new__(acb)
        acb_atanh((<acb>u).val, (<acb>s).val, getprec())
        return u

    def agm(s, t=None):
        """
        Computes the arithmetic-geometric mean `M(s,t)`, or `M(s) = M(s,1)`
        if no extra parameter is passed.

            >>> showgood(lambda: acb(2).agm(), dps=25)
            1.456791031046906869186432
            >>> showgood(lambda: acb(1,1).agm(), dps=25)
            1.049160528732780220531827 + 0.4781557460881612293261882j
            >>> showgood(lambda: (acb(-95,-65)/100).agm(acb(684,747)/1000), dps=25)
            -0.3711072435676023931065922 + 0.3199561471173686568561674j

        """
        if t is None:
            u = acb.__new__(acb)
            acb_agm1((<acb>u).val, (<acb>s).val, getprec())
            return u
        else:
            return (s / t).agm() * t

    def gamma(s):
        """
        Computes the gamma function `\Gamma(s)`.

            >>> showgood(lambda: acb(1,2).gamma(), dps=25)
            0.1519040026700361374481610 + 0.01980488016185498197191013j
        """
        u = acb.__new__(acb)
        acb_gamma((<acb>u).val, (<acb>s).val, getprec())
        return u

    def rgamma(s):
        """
        Computes the reciprocal gamma function `1/\Gamma(s)`, avoiding
        division by zero at the poles of the gamma function.

            >>> showgood(lambda: acb(1,2).rgamma(), dps=25)
            6.473073626019134501563613 - 0.8439438407732021454882999j
            >>> print(acb(0).rgamma())
            0
            >>> print(acb(-1).rgamma())
            0

        """
        u = acb.__new__(acb)
        acb_rgamma((<acb>u).val, (<acb>s).val, getprec())
        return u

    def lgamma(s):
        """
        Computes the logarithmic gamma function `\log \Gamma(s)`.
        The function is defined to be continuous away from the
        negative half-axis and thus differs from `\log(\Gamma(s))` in general.

            >>> showgood(lambda: acb(1,2).lgamma(), dps=25)
            -1.876078786430929341229996 + 0.1296463163097883113837075j
            >>> showgood(lambda: (acb(0,10).lgamma() - acb(0,10).gamma().log()).imag / arb.pi(), dps=25)
            4.000000000000000000000000
        """
        u = acb.__new__(acb)
        acb_lgamma((<acb>u).val, (<acb>s).val, getprec())
        return u

    def digamma(s):
        """
        Computes the digamma function `\psi(s)`.

            >>> showgood(lambda: acb(1,2).digamma(), dps=25)
            0.7145915153739775266568699 + 1.320807282642230228386088j
        """
        u = acb.__new__(acb)
        acb_digamma((<acb>u).val, (<acb>s).val, getprec())
        return u

    def zeta(s, a=None):
        """
        The Riemann zeta function `\zeta(s)`, or the Hurwitz
        zeta function `\zeta(s,a)` if a second parameter is passed.

            >>> showgood(lambda: acb(0.5,1000).zeta(), dps=25)
            0.3563343671943960550744025 + 0.9319978312329936651150604j
            >>> showgood(lambda: acb(1,2).zeta(acb(2,3)), dps=25)
            -2.953059572088556722876240 + 3.410962524512050603254574j
        """
        if a is None:
            u = acb.__new__(acb)
            acb_zeta((<acb>u).val, (<acb>s).val, getprec())
            return u
        else:
            a = any_as_acb(a)
            u = acb.__new__(acb)
            acb_hurwitz_zeta((<acb>u).val, (<acb>s).val, (<acb>a).val, getprec())
            return u

    @classmethod
    def pi(cls):
        """
        Computes the constant `\pi`.

            >>> showgood(lambda: acb.pi(), dps=25)
            3.141592653589793238462643
        """
        u = acb.__new__(acb)
        acb_const_pi((<acb>u).val, getprec())
        return u

    def sqrt(s, bint analytic=False):
        r"""
        Computes the square root `\sqrt{s}`.

            >>> showgood(lambda: acb(1,2).sqrt(), dps=25)
            1.272019649514068964252422 + 0.7861513777574232860695586j
        """
        u = acb.__new__(acb)
        acb_sqrt_analytic((<acb>u).val, (<acb>s).val, analytic, getprec())
        return u

    def rsqrt(s, bint analytic=False):
        r"""
        Computes the reciprocal square root `1/\sqrt{s}`.

            >>> showgood(lambda: acb(1,2).rsqrt(), dps=25)
            0.5688644810057831072783079 - 0.3515775842541429284870573j
        """
        u = acb.__new__(acb)
        acb_rsqrt_analytic((<acb>u).val, (<acb>s).val, analytic, getprec())
        return u

    def exp(s):
        r"""
        Computes the exponential function `\exp(s)`.

            >>> showgood(lambda: acb(1,2).exp(), dps=25)
            -1.131204383756813638431255 + 2.471726672004818927616931j
        """
        u = acb.__new__(acb)
        acb_exp((<acb>u).val, (<acb>s).val, getprec())
        return u

    def exp_pi_i(s):
        r"""
        Computes the exponential function of modified argument `\exp(\pi i s)`.

            >>> showgood(lambda: acb(1,2).exp_pi_i(), dps=25)
            -0.001867442731707988814430213
            >>> showgood(lambda: acb(1.5,2.5).exp_pi_i(), dps=25)
            -0.0003882032039267662472325299j
            >>> showgood(lambda: acb(1.25,2.25).exp_pi_i(), dps=25)
            -0.0006020578259597635239581705 - 0.0006020578259597635239581705j
        """
        u = acb.__new__(acb)
        acb_exp_pi_i((<acb>u).val, (<acb>s).val, getprec())
        return u

    def expm1(s):
        u = acb.__new__(acb)
        acb_expm1((<acb>u).val, (<acb>s).val, getprec())
        return u

    def sin(s):
        r"""
        Computes the sine `\sin(s)`.

            >>> showgood(lambda: acb(1,2).sin(), dps=25)
            3.165778513216168146740735 + 1.959601041421605897070352j
        """
        u = acb.__new__(acb)
        acb_sin((<acb>u).val, (<acb>s).val, getprec())
        return u

    def cos(s):
        r"""
        Computes the cosine `\cos(s)`.

            >>> showgood(lambda: acb(1,2).cos(), dps=25)
            2.032723007019665529436343 - 3.051897799151800057512116j
        """
        u = acb.__new__(acb)
        acb_cos((<acb>u).val, (<acb>s).val, getprec())
        return u

    def sin_cos(s):
        r"""
        Computes `\sin(s)` and `\cos(s)` simultaneously.

            >>> showgood(lambda: acb(1,2).sin_cos(), dps=15)
            (3.16577851321617 + 1.95960104142161j, 2.03272300701967 - 3.05189779915180j)
        """
        u = acb.__new__(acb)
        v = acb.__new__(acb)
        acb_sin_cos((<acb>u).val, (<acb>v).val, (<acb>s).val, getprec())
        return u, v

    def tan(s):
        r"""
        Computes the tangent `\tan(s)`.

            >>> showgood(lambda: acb(1,2).tan(), dps=25)
            0.03381282607989669028437056 + 1.014793616146633568117054j
        """
        u = acb.__new__(acb)
        acb_tan((<acb>u).val, (<acb>s).val, getprec())
        return u

    def cot(s):
        r"""
        Computes the cotangent `\cot(s)`.

            >>> showgood(lambda: acb(1,2).cot(), dps=25)
            0.03279775553375259406276455 - 0.9843292264581910294718882j
        """
        u = acb.__new__(acb)
        acb_cot((<acb>u).val, (<acb>s).val, getprec())
        return u

    def sin_pi(s):
        r"""
        Computes the sine `\sin(\pi s)`.

            >>> showgood(lambda: acb(1,2).sin_pi(), dps=25)
            -267.7448940410165142571174j
        """
        u = acb.__new__(acb)
        acb_sin_pi((<acb>u).val, (<acb>s).val, getprec())
        return u

    def cos_pi(s):
        r"""
        Computes the cosine `\cos(\pi s)`.

            >>> showgood(lambda: acb(1,2).cos_pi(), dps=25)
            -267.7467614837482222459319
        """
        u = acb.__new__(acb)
        acb_cos_pi((<acb>u).val, (<acb>s).val, getprec())
        return u

    def sin_cos_pi(s):
        r"""
        Computes `\sin(\pi s)` and `\cos(\pi s)` simultaneously.

            >>> showgood(lambda: acb(1,2).sin_cos_pi(), dps=25)
            (-267.7448940410165142571174j, -267.7467614837482222459319)
        """
        u = acb.__new__(acb)
        v = acb.__new__(acb)
        acb_sin_cos_pi((<acb>u).val, (<acb>v).val, (<acb>s).val, getprec())
        return u, v

    def tan_pi(s):
        r"""
        Computes the tangent `\tan(\pi s)`.

            >>> showgood(lambda: acb(1,2).tan_pi(), dps=25)
            0.9999930253396106106051072j
        """
        u = acb.__new__(acb)
        acb_tan_pi((<acb>u).val, (<acb>s).val, getprec())
        return u

    def cot_pi(s):
        r"""
        Computes the cotangent `\cot(\pi s)`.

            >>> showgood(lambda: acb(1,2).cot_pi(), dps=25)
            -1.000006974709035616233122j
        """
        u = acb.__new__(acb)
        acb_cot_pi((<acb>u).val, (<acb>s).val, getprec())
        return u

    def sec(s):
        u = acb.__new__(acb)
        acb_sec((<acb>u).val, (<acb>s).val, getprec())
        return u

    def csc(s):
        u = acb.__new__(acb)
        acb_csc((<acb>u).val, (<acb>s).val, getprec())
        return u

    def sinh(s):
        u = acb.__new__(acb)
        acb_sinh((<acb>u).val, (<acb>s).val, getprec())
        return u

    def cosh(s):
        u = acb.__new__(acb)
        acb_cosh((<acb>u).val, (<acb>s).val, getprec())
        return u

    def sinh_cosh(s):
        u = acb.__new__(acb)
        v = acb.__new__(acb)
        acb_sinh_cosh((<acb>u).val, (<acb>v).val, (<acb>s).val, getprec())
        return u, v

    def tanh(s):
        u = acb.__new__(acb)
        acb_tanh((<acb>u).val, (<acb>s).val, getprec())
        return u

    def coth(s):
        u = acb.__new__(acb)
        acb_coth((<acb>u).val, (<acb>s).val, getprec())
        return u

    def sech(s):
        u = acb.__new__(acb)
        acb_sech((<acb>u).val, (<acb>s).val, getprec())
        return u

    def csch(s):
        u = acb.__new__(acb)
        acb_csch((<acb>u).val, (<acb>s).val, getprec())
        return u

    def sinc(s):
        u = acb.__new__(acb)
        acb_sinc((<acb>u).val, (<acb>s).val, getprec())
        return u

    def sinc_pi(s):
        u = acb.__new__(acb)
        acb_sinc_pi((<acb>u).val, (<acb>s).val, getprec())
        return u

    def rising(s, n):
        """
        Computes the rising factorial `(s)_n`.

            >>> showgood(lambda: acb(1,2).rising(5), dps=25)
            -540.0000000000000000000000 - 100.0000000000000000000000j
            >>> showgood(lambda: acb(1,2).rising(2+3j), dps=25)
            0.3898076751098812033498554 - 0.01296289149274537721607465j
        """
        u = acb.__new__(acb)
        n = any_as_acb(n)
        acb_rising((<acb>u).val, (<acb>s).val, (<acb>n).val, getprec())
        return u

    def rising2(s, ulong n):
        """
        Computes the rising factorial `(s)_n` where *n* is an unsigned
        integer, along with the first derivative with respect to `(s)_n`.
        The current implementation does not use the gamma function,
        so *n* should be moderate.

            >>> showgood(lambda: acb(1,2).rising2(5), dps=25)
            (-540.0000000000000000000000 - 100.0000000000000000000000j, -666.0000000000000000000000 + 420.0000000000000000000000j)
        """
        u = acb.__new__(acb)
        v = acb.__new__(acb)
        acb_rising2_ui((<acb>u).val, (<acb>v).val, (<acb>s).val, n, getprec())
        return u, v

    @classmethod
    def polylog(cls, s, z):
        """
        Computes the polylogarithm `\operatorname{Li}_s(z)`.

            >>> showgood(lambda: acb.polylog(acb(2), acb(3)), dps=25)
            2.320180423313098396406194 - 3.451392295223202661433821j
            >>> showgood(lambda: acb.polylog(acb(1,2), acb(2,3)), dps=25)
            -6.854984251829306907116775 + 7.375884252099214498443165j
        """
        u = acb.__new__(acb)
        s = any_as_acb(s)
        z = any_as_acb(z)
        acb_polylog((<acb>u).val, (<acb>s).val, (<acb>z).val, getprec())
        return u

    @classmethod
    def polygamma(cls, s, z):
        u = acb.__new__(acb)
        s = any_as_acb(s)
        z = any_as_acb(z)
        acb_polygamma((<acb>u).val, (<acb>s).val, (<acb>z).val, getprec())
        return u

    def log_barnes_g(s):
        u = acb.__new__(acb)
        acb_log_barnes_g((<acb>u).val, (<acb>s).val, getprec())
        return u

    def barnes_g(s):
        u = acb.__new__(acb)
        acb_barnes_g((<acb>u).val, (<acb>s).val, getprec())
        return u

    @classmethod
    def hypgeom(cls, a, b, z, bint regularized=False, long n=-1):
        r"""
        Computes the generalized hypergeometric function `{}_pF_q(a;b;z)`
        given lists of complex numbers `a` and `b` and a complex number `z`.

        The optional parameter *n*, if nonnegative, controls the number
        of terms to add in the hypergeometric series. This is just a tuning
        parameter: a rigorous error bound is computed regardless of *n*.

            >>> showgood(lambda: acb.hypgeom([1+1j, 2-2j], [3, fmpq(1,3)], acb.pi()), dps=25)  # 2F2
            144.9760711583421645394627 - 51.06535684838559608699106j

        """
        cdef long i, p, q, prec
        cdef acb_ptr aa, bb
        a = [any_as_acb(t) for t in a]
        b = [any_as_acb(t) for t in b]
        if n != -1:
            b += [acb(1)]
        z = acb(z)
        p = len(a)
        q = len(b)
        aa = <acb_ptr>libc.stdlib.malloc(p * cython.sizeof(acb_struct))
        bb = <acb_ptr>libc.stdlib.malloc(q * cython.sizeof(acb_struct))
        for i in range(p):
            aa[i] = (<acb>(a[i])).val[0]
        for i in range(q):
            bb[i] = (<acb>(b[i])).val[0]
        u = acb.__new__(acb)
        if n == -1:
            acb_hypgeom_pfq((<acb>u).val, aa, p, bb, q, (<acb>z).val, regularized, getprec())
        else:
            if regularized:
                raise NotImplementedError
            acb_hypgeom_pfq_direct((<acb>u).val, aa, p, bb, q, (<acb>z).val, n, getprec())
        libc.stdlib.free(aa)
        libc.stdlib.free(bb)
        return u

    @classmethod
    def hypgeom_u(cls, a, b, z, long n=-1, bint asymp=False):
        r"""
        Computes Tricomi's confluent hypergeometric function `U(a,b,z)`
        given complex numbers `a`, `b`, `z`.

        If *asymp* is set to *True* the asymptotic series is forced.
        If `|z|` is small, the attainable accuracy is then limited.
        The optional parameter *n*, if nonnegative, controls the number
        of terms to add in the asymptotic series. This is just a tuning
        parameter: a rigorous error bound is computed regardless of *n*.

            >>> showgood(lambda: acb.hypgeom_u(1+1j, 2+3j, 400+500j), dps=25)
            0.001836433961463105354717547 - 0.003358699641979853540147122j
            >>> showgood(lambda: acb.hypgeom_u(1+1j, 2+3j, -30), dps=25)
            0.7808944974399200669072087 - 0.2674783064947089569672470j
            >>> print(acb.hypgeom_u(1+1j, 2+3j, -30, n=0, asymp=True))
            [+/- 3.41] + [+/- 3.41]j
            >>> print(acb.hypgeom_u(1+1j, 2+3j, -30, n=30, asymp=True))
            [0.7808944974 +/- 7.46e-11] + [-0.2674783065 +/- 3.31e-11]j
            >>> print(acb.hypgeom_u(1+1j, 2+3j, -30, n=60, asymp=True))
            [0.78089 +/- 8.14e-6] + [-0.26748 +/- 5.34e-6]j
        """
        a = any_as_acb(a)
        b = any_as_acb(b)
        z = any_as_acb(z)
        if asymp:
            t = z**(-a)
            u = acb.__new__(acb)
            acb_hypgeom_u_asymp((<acb>u).val, (<acb>a).val, (<acb>b).val, (<acb>z).val, n, getprec())
            return t * u
        else:
            u = acb.__new__(acb)
            acb_hypgeom_u((<acb>u).val, (<acb>a).val, (<acb>b).val, (<acb>z).val, getprec())
            return u

    @classmethod
    def hypgeom_1f1(cls, a, b, z, bint regularized=False):
        r"""
        Computes Kummer's confluent hypergeometric function `{}_1F_1(a,b,z)`
        given complex numbers `a`, `b`, `z`. Optionally, computes
        the regularized version.

            >>> showgood(lambda: acb.hypgeom_1f1(2+3j, 3+4j, 40000+50000j), dps=25)
            3.730925582634533963357515e+17366 + 3.199717318207534638202987e+17367j
            >>> showgood(lambda: acb.hypgeom_1f1(2+3j, 3+4j, 40000+50000j) / acb(3+4j).gamma(), dps=25)
            -1.846160890579724375436801e+17368 + 2.721369772032882467996588e+17367j
            >>> showgood(lambda: acb.hypgeom_1f1(5, -3, 10, regularized=True), dps=25)
            832600407043.6938843410086
            >>> showgood(lambda: acb.hypgeom_1f1(-5,-6,10), dps=25)
            403.7777777777777777777778
            >>> showgood(lambda: acb.hypgeom_1f1(-5,-5,10,regularized=True), dps=25)
            0
            >>> showgood(lambda: acb.hypgeom_1f1(-5,-6,10,regularized=True), dps=25)
            0
            >>> showgood(lambda: acb.hypgeom_1f1(-5,-4,10,regularized=True), dps=25)
            -100000.0000000000000000000

        """
        a = any_as_acb(a)
        b = any_as_acb(b)
        z = any_as_acb(z)
        u = acb.__new__(acb)
        acb_hypgeom_1f1((<acb>u).val, (<acb>a).val, (<acb>b).val, (<acb>z).val, regularized, getprec())
        return u

    @classmethod
    def bessel_j(cls, a, z):
        r"""
        Computes the Bessel function of the first kind `J_a(z)`.

            >>> showgood(lambda: acb.bessel_j(1+2j, 2+3j), dps=25)
            0.5041904509946947234759103 - 0.1765180072689450645147231j
        """
        a = any_as_acb(a)
        z = any_as_acb(z)
        u = acb.__new__(acb)
        acb_hypgeom_bessel_j((<acb>u).val, (<acb>a).val, (<acb>z).val, getprec())
        return u

    @classmethod
    def bessel_k(cls, a, z, bint scaled=False):
        r"""
        Computes the modified Bessel function of the second kind `K_a(z)`.

            >>> showgood(lambda: acb.bessel_k(1+2j, 2+3j), dps=25)
            -0.09884736370006798963642719 - 0.02870616366668971734065520j
            >>> showgood(lambda: acb.bessel_k(3, 2), dps=25)
            0.6473853909486341531592356

        """
        a = any_as_acb(a)
        z = any_as_acb(z)
        u = acb.__new__(acb)
        if scaled:
            acb_hypgeom_bessel_k_scaled((<acb>u).val, (<acb>a).val, (<acb>z).val, getprec())
        else:
            acb_hypgeom_bessel_k((<acb>u).val, (<acb>a).val, (<acb>z).val, getprec())
        return u

    @classmethod
    def bessel_i(cls, a, z, bint scaled=False):
        r"""
        Computes the modified Bessel function of the first kind `I_a(z)`.

        """
        a = any_as_acb(a)
        z = any_as_acb(z)
        u = acb.__new__(acb)
        if scaled:
            acb_hypgeom_bessel_i_scaled((<acb>u).val, (<acb>a).val, (<acb>z).val, getprec())
        else:
            acb_hypgeom_bessel_i((<acb>u).val, (<acb>a).val, (<acb>z).val, getprec())
        return u

    @classmethod
    def bessel_y(cls, a, z):
        r"""
        Computes the Bessel function of the second kind `Y_a(z)`.

        """
        a = any_as_acb(a)
        z = any_as_acb(z)
        u = acb.__new__(acb)
        acb_hypgeom_bessel_y((<acb>u).val, (<acb>a).val, (<acb>z).val, getprec())
        return u

    def erf(s):
        r"""
        Computes the error function `\operatorname{erf}(s)`.

            >>> showgood(lambda: acb(2+3j).erf() - 1, dps=25)
            -21.82946142761456838910309 + 8.687318271470163144428079j
            >>> showgood(lambda: acb("77.7").erf() - 1, dps=25, maxdps=10000)
            -7.929310690520378873143053e-2625
        """
        u = acb.__new__(acb)
        acb_hypgeom_erf((<acb>u).val, (<acb>s).val, getprec())
        return u

    @classmethod
    def modular_theta(cls, z, tau):
        r"""
        Computes the Jacobi theta functions `\theta_1(z,\tau)`,
        `\theta_2(z,\tau)`, `\theta_3(z,\tau)`, `\theta_4(z,\tau)`,
        returning all four values as a tuple.
        We define the theta functions with a factor `\pi` for *z* included;
        for example
        `\theta_3(z,\tau) = 1 + 2 \sum_{n=1}^{\infty} q^{n^2} \cos(2n\pi z)`.

            >>> for i in range(4):
            ...     showgood(lambda: acb.modular_theta(1+1j, 1.25+3j)[i], dps=25)
            ... 
            1.820235910124989594900076 - 1.216251950154477951760042j
            -1.220790267576967690128359 - 1.827055516791154669091679j
            0.9694430387796704100046143 - 0.03055696120816803328582847j
            1.030556961196006476576271 + 0.03055696120816803328582847j
        """
        z = any_as_acb(z)
        tau = any_as_acb(tau)
        t1 = acb.__new__(acb)
        t2 = acb.__new__(acb)
        t3 = acb.__new__(acb)
        t4 = acb.__new__(acb)
        acb_modular_theta((<acb>t1).val, (<acb>t2).val,
                          (<acb>t3).val, (<acb>t4).val,
                          (<acb>z).val, (<acb>tau).val, getprec())
        return (t1, t2, t3, t4)

    def modular_eta(tau):
        r"""
        Computes the Dedekind eta function `\eta(\tau)`.

            >>> showgood(lambda: acb(1+1j).modular_eta(), dps=25)
            0.7420487758365647263392722 + 0.1988313702299107190516142j
        """
        u = acb.__new__(acb)
        acb_modular_eta((<acb>u).val, (<acb>tau).val, getprec())
        return u

    def modular_j(tau):
        r"""
        Computes the modular *j*-invariant `j(\tau)`.

            >>> showgood(lambda: (1 + acb(-163).sqrt()/2).modular_j(), dps=25)
            262537412640769488.0000000
            >>> showgood(lambda: acb(3.25+100j).modular_j() - 744, dps=25, maxdps=2000)
            -3.817428033843319485125773e-539 - 7.503618895582604309279721e+272j
        """
        u = acb.__new__(acb)
        acb_modular_j((<acb>u).val, (<acb>tau).val, getprec())
        return u

    def modular_lambda(tau):
        r"""
        Computes the modular lambda function `\lambda(\tau)`.

            >>> showgood(lambda: acb(0.25+5j).modular_lambda(), dps=25)
            1.704995415668039343330405e-6 + 1.704992508662079437786598e-6j
        """
        u = acb.__new__(acb)
        acb_modular_lambda((<acb>u).val, (<acb>tau).val, getprec())
        return u

    def modular_delta(tau):
        r"""
        Computes the modular discriminant `\Delta(\tau)`.

            >>> showgood(lambda: acb(0.25+5j).modular_delta(), dps=25)
            1.237896015010281684100435e-26 + 2.271101068324093838679275e-14j
        """
        u = acb.__new__(acb)
        acb_modular_delta((<acb>u).val, (<acb>tau).val, getprec())
        return u

    def elliptic_k(m):
        """
        Computes the complete elliptic integral of the first kind `K(m)`.

            >>> showgood(lambda: 2 * acb(0).elliptic_k(), dps=25)
            3.141592653589793238462643
            >>> showgood(lambda: acb(100+50j).elliptic_k(), dps=25)
            0.2052037361984861505113972 + 0.3158446040520529200980591j
            >>> showgood(lambda: (1 - acb("1e-100")).elliptic_k(), dps=25)
            116.5155490108221748197340

        """
        u = acb.__new__(acb)
        acb_modular_elliptic_k((<acb>u).val, (<acb>m).val, getprec())
        return u

    def elliptic_e(m):
        """
        Computes the complete elliptic integral of the second kind `E(m)`.

            >>> showgood(lambda: 2 * acb(0).elliptic_e(), dps=25)
            3.141592653589793238462643
            >>> showgood(lambda: acb(100+50j).elliptic_e(), dps=25)
            2.537235535915583230880553 - 10.10997759792420947130321j
            >>> showgood(lambda: (1 - acb("1e-100")).elliptic_e(), dps=25)
            1.000000000000000000000000

        """
        u = acb.__new__(acb)
        acb_modular_elliptic_e((<acb>u).val, (<acb>m).val, getprec())
        return u

    def unique_fmpz(self):
        u = fmpz.__new__(fmpz)
        if acb_get_unique_fmpz((<fmpz>u).val, self.val):
            return u
        else:
            return None

    @classmethod
    def gamma_upper(cls, s, z, int regularized=0):
        r"""
        Computes the upper incomplete gamma function `\Gamma(s,z)`.

            >>> showgood(lambda: acb.gamma_upper(1+2j, 2+3j), dps=25)
            0.02614303924198793235765248 - 0.0007536537278463329391666679j
        """
        s = any_as_acb(s)
        z = any_as_acb(z)
        u = acb.__new__(acb)
        acb_hypgeom_gamma_upper((<acb>u).val, (<acb>s).val, (<acb>z).val, regularized, getprec())
        return u

    @classmethod
    def gamma_lower(cls, s, z, int regularized=0):
        r"""
        Computes the lower incomplete gamma function `\gamma(s,z)`.
        """
        s = any_as_acb(s)
        z = any_as_acb(z)
        u = acb.__new__(acb)
        acb_hypgeom_gamma_lower((<acb>u).val, (<acb>s).val, (<acb>z).val, regularized, getprec())
        return u

    @classmethod
    def beta_lower(cls, a, b, z, int regularized=0):
        r"""
        Computes the lower incomplete beta function `B(a,b;z)`.
        """
        a = any_as_acb(a)
        b = any_as_acb(b)
        z = any_as_acb(z)
        u = acb.__new__(acb)
        acb_hypgeom_beta_lower((<acb>u).val, (<acb>a).val, (<acb>b).val, (<acb>z).val, regularized, getprec())
        return u

    @classmethod
    def expint(cls, s, z):
        r"""
        Computes the exponential integral `E_s(z)`.

            >>> showgood(lambda: acb.expint(1+2j, 2+3j), dps=25)
            -0.01442661495527080336037156 + 0.01942348372986687164972794j
        """
        s = any_as_acb(s)
        z = any_as_acb(z)
        u = acb.__new__(acb)
        acb_hypgeom_expint((<acb>u).val, (<acb>s).val, (<acb>z).val, getprec())
        return u

    def erfc(s):
        r"""
        Computes the complementary error function `\operatorname{erfc}(s)`.

            >>> showgood(lambda: acb("77.7").erfc(), dps=25)
            7.929310690520378873143053e-2625
        """
        u = acb.__new__(acb)
        acb_hypgeom_erfc((<acb>u).val, (<acb>s).val, getprec())
        return u

    def erfi(s):
        r"""
        Computes the imaginary error function `\operatorname{erfi}(s)`.

            >>> showgood(lambda: acb(10).erfi(), dps=25)
            1.524307422708669699360547e+42
        """
        u = acb.__new__(acb)
        acb_hypgeom_erfi((<acb>u).val, (<acb>s).val, getprec())
        return u

    def rel_accuracy_bits(self):
        return acb_rel_accuracy_bits(self.val)

    def ei(s):
        r"""
        Computes the exponential integral `\operatorname{Ei}(s)`.

            >>> showgood(lambda: acb(10).ei(), dps=25)
            2492.228976241877759138440
        """
        u = acb.__new__(acb)
        acb_hypgeom_ei((<acb>u).val, (<acb>s).val, getprec())
        return u

    def si(s):
        r"""
        Computes the sine integral `\operatorname{Si}(s)`.

            >>> showgood(lambda: acb(10).si(), dps=25)
            1.658347594218874049330972
        """
        u = acb.__new__(acb)
        acb_hypgeom_si((<acb>u).val, (<acb>s).val, getprec())
        return u

    def ci(s):
        r"""
        Computes the cosine integral `\operatorname{Ci}(s)`.

            >>> showgood(lambda: acb(10).ci(), dps=25)
            -0.04545643300445537263453283
        """
        u = acb.__new__(acb)
        acb_hypgeom_ci((<acb>u).val, (<acb>s).val, getprec())
        return u

    def shi(s):
        r"""
        Computes the hyperbolic sine integral `\operatorname{Shi}(s)`.

            >>> showgood(lambda: acb(10).shi(), dps=25)
            1246.114490199423344411882
        """
        u = acb.__new__(acb)
        acb_hypgeom_shi((<acb>u).val, (<acb>s).val, getprec())
        return u

    def chi(s):
        r"""
        Computes the hyperbolic cosine integral `\operatorname{Chi}(s)`.

            >>> showgood(lambda: acb(10).chi(), dps=25)
            1246.114486042454414726558
        """
        u = acb.__new__(acb)
        acb_hypgeom_chi((<acb>u).val, (<acb>s).val, getprec())
        return u

    def li(s, bint offset=False):
        r"""
        Computes the logarithmic integral `\operatorname{li}(s)`, optionally
        the offset logarithmic integral
        `\operatorname{Li}(s) = \operatorname{li}(s) - \operatorname{li}(2)`.

            >>> showgood(lambda: acb(10).li(), dps=25)
            6.165599504787297937522982
            >>> showgood(lambda: acb(10).li(offset=True), dps=25)
            5.120435724669805152678393
        """
        u = acb.__new__(acb)
        acb_hypgeom_li((<acb>u).val, (<acb>s).val, offset, getprec())
        return u

    @classmethod
    def hypgeom_2f1(cls, a, b, c, z, bint regularized=False, bint ab=False, bint ac=False, bc=False, abc=False):
        r"""
        """
        cdef int flags
        a = any_as_acb(a)
        b = any_as_acb(b)
        c = any_as_acb(c)
        z = any_as_acb(z)
        u = acb.__new__(acb)
        flags = 0
        if regularized: flags |= 1
        if ab: flags |= 2
        if ac: flags |= 4
        if bc: flags |= 8
        if abc: flags |= 16
        acb_hypgeom_2f1((<acb>u).val, (<acb>a).val, (<acb>b).val, (<acb>c).val,
            (<acb>z).val, flags, getprec())
        return u

    def chebyshev_t(s, n):
        r"""
        Chebyshev function of the first kind `T_n(s)`.

            >>> showgood(lambda: (acb(1)/3).chebyshev_t(3), dps=25)
            -0.8518518518518518518518519
        """
        v = acb.__new__(acb)
        n = any_as_acb(n)
        acb_hypgeom_chebyshev_t((<acb>v).val, (<acb>n).val, (<acb>s).val, getprec())
        return v

    def chebyshev_u(s, n):
        r"""
        Chebyshev function of the second kind `U_n(s)`.

            >>> showgood(lambda: (acb(1)/3).chebyshev_u(3), dps=25)
            -1.037037037037037037037037
        """
        v = acb.__new__(acb)
        n = any_as_acb(n)
        acb_hypgeom_chebyshev_u((<acb>v).val, (<acb>n).val, (<acb>s).val, getprec())
        return v

    def jacobi_p(s, n, a, b):
        r"""
        Jacobi polynomial (or Jacobi function) `P_n^{a,b}(s)`.

            >>> showgood(lambda: (acb(1)/3).jacobi_p(5, 0.25, 0.5), dps=25)
            0.4131944444444444444444444
        """
        v = acb.__new__(acb)
        n = any_as_acb(n)
        a = any_as_acb(a)
        b = any_as_acb(b)
        acb_hypgeom_jacobi_p((<acb>v).val, (<acb>n).val, (<acb>a).val, (<acb>b).val, (<acb>s).val, getprec())
        return v

    def gegenbauer_c(s, n, m):
        r"""
        Gegenbauer function `C_n^{m}(s)`.

            >>> showgood(lambda: (acb(1)/3).gegenbauer_c(5, 0.25), dps=25)
            0.1321855709876543209876543
        """
        v = acb.__new__(acb)
        n = any_as_acb(n)
        m = any_as_acb(m)
        acb_hypgeom_gegenbauer_c((<acb>v).val, (<acb>n).val, (<acb>m).val, (<acb>s).val, getprec())
        return v

    def laguerre_l(s, n, m=0):
        r"""
        Laguerre function `L_n^{m}(s)`.

            >>> showgood(lambda: (acb(1)/3).laguerre_l(5, 0.25), dps=25)
            0.03871323490012002743484225
        """
        v = acb.__new__(acb)
        n = any_as_acb(n)
        m = any_as_acb(m)
        acb_hypgeom_laguerre_l((<acb>v).val, (<acb>n).val, (<acb>m).val, (<acb>s).val, getprec())
        return v

    def hermite_h(s, n):
        r"""
        Hermite function `H_n(s)`.

            >>> showgood(lambda: (acb(1)/3).hermite_h(5), dps=25)
            34.20576131687242798353909
        """
        v = acb.__new__(acb)
        n = any_as_acb(n)
        acb_hypgeom_hermite_h((<acb>v).val, (<acb>n).val, (<acb>s).val, getprec())
        return v

    def legendre_p(s, n, m=0, int type=2):
        r"""
        Legendre function of the first kind `P_n^m(z)`.

            >>> showgood(lambda: (acb(1)/3).legendre_p(5), dps=25)
            0.3333333333333333333333333
            >>> showgood(lambda: (acb(1)/3).legendre_p(5, 1.5), dps=25)
            -2.372124991643971726805456
            >>> showgood(lambda: (acb(3)).legendre_p(5, 1.5, type=2), dps=25)
            -12091.31397811720689120900 - 12091.31397811720689120900j
            >>> showgood(lambda: (acb(3)).legendre_p(5, 1.5, type=3), dps=25)
            17099.70021476473458984981

        The optional parameter *type* can be 2 or 3, and selects between
        two different branch cut conventions (see *Mathematica* and *mpmath*).
        """
        v = acb.__new__(acb)
        n = any_as_acb(n)
        m = any_as_acb(m)
        if type != 2 and type != 3:
            raise ValueError("type must be 2 or 3")
        type -= 2
        acb_hypgeom_legendre_p((<acb>v).val, (<acb>n).val, (<acb>m).val, (<acb>s).val, type, getprec())
        return v

    def legendre_q(s, n, m=0, int type=2):
        r"""
        Legendre function of the second kind `Q_n^m(z)`.

            >>> showgood(lambda: (acb(1)/3).legendre_q(5), dps=25)
            0.1655245300933242182362054
            >>> showgood(lambda: (acb(1)/3).legendre_q(5, 1.5), dps=25)
            -6.059967350218583975575616
            >>> showgood(lambda: (acb(3)).legendre_q(5, 1.5, type=2), dps=25)
            -18992.99179585607569083751 + 18992.99179585607569083751j
            >>> showgood(lambda: (acb(3)).legendre_q(5, 1.5, type=3), dps=25)
            -0.0003010942389043591251820434j

        The optional parameter *type* can be 2 or 3, and selects between
        two different branch cut conventions (see *Mathematica* and *mpmath*).
        """
        v = acb.__new__(acb)
        n = any_as_acb(n)
        m = any_as_acb(m)
        if type != 2 and type != 3:
            raise ValueError("type must be 2 or 3")
        type -= 2
        acb_hypgeom_legendre_q((<acb>v).val, (<acb>n).val, (<acb>m).val, (<acb>s).val, type, getprec())
        return v

    @classmethod
    def hypgeom_0f1(cls, a, z, bint regularized=False):
        r"""
        """
        a = any_as_acb(a)
        z = any_as_acb(z)
        u = acb.__new__(acb)
        acb_hypgeom_0f1((<acb>u).val, (<acb>a).val, (<acb>z).val, regularized, getprec())
        return u

    def dirichlet_eta(s):
        u = acb.__new__(acb)
        acb_dirichlet_eta((<acb>u).val, (<acb>s).val, getprec())
        return u

    def airy_ai(s, int derivative=0):
        r"""
        Airy function `\operatorname{Ai}(s)`, or
        `\operatorname{Ai}'(s)` if *derivative* is 1.

            >>> showgood(lambda: acb(-1+1j).airy_ai(), dps=25)
            0.8221174265552725939610246 - 0.1199663426644243438939006j
            >>> showgood(lambda: acb(-1+1j).airy_ai(derivative=1), dps=25)
            -0.3790604792268334962367164 - 0.6045001308622460716372591j
        """
        u = acb.__new__(acb)
        if derivative == 0:
            acb_hypgeom_airy((<acb>u).val, NULL, NULL, NULL, (<acb>s).val, getprec())
        elif derivative == 1:
            acb_hypgeom_airy(NULL, (<acb>u).val, NULL, NULL, (<acb>s).val, getprec())
        else:
            raise ValueError("derivative must be 0 or 1")
        return u

    def airy_bi(s, int derivative=0):
        r"""
        Airy function `\operatorname{Bi}(s)`, or
        `\operatorname{Bi}'(s)` if *derivative* is 1.

            >>> showgood(lambda: acb(-1+1j).airy_bi(), dps=25)
            0.2142904015348735739780868 + 0.6739169237227052095951775j
            >>> showgood(lambda: acb(-1+1j).airy_bi(derivative=1), dps=25)
            0.8344734885227826369040951 - 0.3465260632668285285537957j
        """
        u = acb.__new__(acb)
        if derivative == 0:
            acb_hypgeom_airy(NULL, NULL, (<acb>u).val, NULL, (<acb>s).val, getprec())
        elif derivative == 1:
            acb_hypgeom_airy(NULL, NULL, NULL, (<acb>u).val, (<acb>s).val, getprec())
        else:
            raise ValueError("derivative must be 0 or 1")
        return u

    def airy(s):
        r"""
        Computes the Airy function values `\operatorname{Ai}(s)`,
        `\operatorname{Ai}'(s)`, `\operatorname{Bi}(s)`,
        `\operatorname{Bi}'(s)` simultaneously, returning a tuple.

            >>> showgood(lambda: acb(-1+1j).airy(), dps=5)
            (0.82212 - 0.11997j, -0.37906 - 0.60450j, 0.21429 + 0.67392j, 0.83447 - 0.34653j)
        """
        u = acb.__new__(acb)
        v = acb.__new__(acb)
        w = acb.__new__(acb)
        z = acb.__new__(acb)
        acb_hypgeom_airy((<acb>u).val, (<acb>v).val,
                        (<acb>w).val, (<acb>z).val, (<acb>s).val, getprec())
        return u, v, w, z

    def fresnel_s(s, bint normalized=True):
        u = acb.__new__(acb)
        acb_hypgeom_fresnel((<acb>u).val, NULL, (<acb>s).val, normalized, getprec())
        return u

    def fresnel_c(s, bint normalized=True):
        u = acb.__new__(acb)
        acb_hypgeom_fresnel(NULL, (<acb>u).val, (<acb>s).val, normalized, getprec())
        return u

    def log_sin_pi(s):
        u = acb.__new__(acb)
        acb_log_sin_pi((<acb>u).val, (<acb>s).val, getprec())
        return u

    @classmethod
    def elliptic_rf(cls, x, y, z):
        r"""
        """
        x = any_as_acb(x)
        y = any_as_acb(y)
        z = any_as_acb(z)
        u = acb.__new__(acb)
        acb_elliptic_rf((<acb>u).val, (<acb>x).val, (<acb>y).val, (<acb>z).val, 0, getprec())
        return u

    @classmethod
    def elliptic_rc(cls, x, y):
        r"""
        """
        x = any_as_acb(x)
        y = any_as_acb(y)
        u = acb.__new__(acb)
        acb_elliptic_rf((<acb>u).val, (<acb>x).val, (<acb>y).val, (<acb>y).val, 0, getprec())
        return u

    @classmethod
    def elliptic_rj(cls, x, y, z, p):
        r"""
        """
        x = any_as_acb(x)
        y = any_as_acb(y)
        z = any_as_acb(z)
        p = any_as_acb(p)
        u = acb.__new__(acb)
        acb_elliptic_rj((<acb>u).val, (<acb>x).val, (<acb>y).val, (<acb>z).val, (<acb>p).val, 0, getprec())
        return u

    @classmethod
    def elliptic_rd(cls, x, y, z):
        r"""
        """
        x = any_as_acb(x)
        y = any_as_acb(y)
        z = any_as_acb(z)
        u = acb.__new__(acb)
        acb_elliptic_rj((<acb>u).val, (<acb>x).val, (<acb>y).val, (<acb>z).val, (<acb>z).val, 0, getprec())
        return u

    @classmethod
    def elliptic_rg(cls, x, y, z):
        r"""
        """
        x = any_as_acb(x)
        y = any_as_acb(y)
        z = any_as_acb(z)
        u = acb.__new__(acb)
        acb_elliptic_rg((<acb>u).val, (<acb>x).val, (<acb>y).val, (<acb>z).val, 0, getprec())
        return u

    @classmethod
    def elliptic_f(cls, phi, m, bint pi=False):
        r"""
        """
        phi = any_as_acb(phi)
        m = any_as_acb(m)
        u = acb.__new__(acb)
        acb_elliptic_f((<acb>u).val, (<acb>phi).val, (<acb>m).val, pi, getprec())
        return u

    @classmethod
    def elliptic_e_inc(cls, phi, m, bint pi=False):
        r"""
        """
        phi = any_as_acb(phi)
        m = any_as_acb(m)
        u = acb.__new__(acb)
        acb_elliptic_e_inc((<acb>u).val, (<acb>phi).val, (<acb>m).val, pi, getprec())
        return u

    @classmethod
    def elliptic_pi(cls, n, m):
        r"""
        """
        n = any_as_acb(n)
        m = any_as_acb(m)
        u = acb.__new__(acb)
        acb_elliptic_pi((<acb>u).val, (<acb>n).val, (<acb>m).val, getprec())
        return u

    @classmethod
    def elliptic_pi_inc(cls, n, phi, m, bint pi=False):
        r"""
        """
        n = any_as_acb(n)
        phi = any_as_acb(phi)
        m = any_as_acb(m)
        u = acb.__new__(acb)
        acb_elliptic_pi_inc((<acb>u).val, (<acb>n).val, (<acb>phi).val, (<acb>m).val, pi, getprec())
        return u

    @classmethod
    def elliptic_p(cls, z, tau):
        r"""
        """
        z = any_as_acb(z)
        tau = any_as_acb(tau)
        u = acb.__new__(acb)
        acb_elliptic_p((<acb>u).val, (<acb>z).val, (<acb>tau).val, getprec())
        return u

    @classmethod
    def elliptic_zeta(cls, z, tau):
        r"""
        """
        z = any_as_acb(z)
        tau = any_as_acb(tau)
        u = acb.__new__(acb)
        acb_elliptic_zeta((<acb>u).val, (<acb>z).val, (<acb>tau).val, getprec())
        return u

    @classmethod
    def elliptic_sigma(cls, z, tau):
        r"""
        """
        z = any_as_acb(z)
        tau = any_as_acb(tau)
        u = acb.__new__(acb)
        acb_elliptic_sigma((<acb>u).val, (<acb>z).val, (<acb>tau).val, getprec())
        return u

    @classmethod
    def elliptic_inv_p(cls, z, tau):
        r"""
        """
        z = any_as_acb(z)
        tau = any_as_acb(tau)
        u = acb.__new__(acb)
        acb_elliptic_inv_p((<acb>u).val, (<acb>z).val, (<acb>tau).val, getprec())
        return u

    @classmethod
    def elliptic_roots(cls, tau):
        r"""
        """
        tau = any_as_acb(tau)
        e1 = acb.__new__(acb)
        e2 = acb.__new__(acb)
        e3 = acb.__new__(acb)
        acb_elliptic_roots((<acb>e1).val, (<acb>e2).val, (<acb>e3).val, (<acb>tau).val, getprec())
        return (e1, e2, e3)

    @classmethod
    def elliptic_invariants(cls, tau):
        r"""
        """
        tau = any_as_acb(tau)
        g1 = acb.__new__(acb)
        g2 = acb.__new__(acb)
        acb_elliptic_invariants((<acb>g1).val, (<acb>g2).val, (<acb>tau).val, getprec())
        return (g1, g2)

    def lambertw(s, branch=0, bint left=False, bint middle=False):
        r"""
        Lambert *W* function, `W_k(s)` where *k* is given by *branch*.

            >>> showgood(lambda: acb(1).lambertw(), dps=25)
            0.5671432904097838729999687
            >>> showgood(lambda: acb(1).lambertw(-1), dps=25)
            -1.533913319793574507919741 - 4.375185153061898385470907j
            >>> showgood(lambda: acb(1).lambertw(-2), dps=25)
            -2.401585104868002884174140 - 10.77629951611507089849710j

        The branch cuts follow Corless et al. by default. This function
        allows selecting alternative branch cuts in order to support
        continuous analytic continuation on complex intervals.

        If *left* is set, computes `W_{\mathrm{left}|k}(z)`
        which corresponds to `W_k(z)` in the upper
        half plane and `W_{k+1}(z)` in the lower half plane, connected continuously
        to the left of the branch points.
        In other words, the branch cut on `(-\infty,0)` is rotated counterclockwise
        to `(0,+\infty)`.
        (For `k = -1` and `k = 0`, there is also a branch cut on `(-1/e,0)`,
        continuous from below instead of from above to maintain counterclockwise
        continuity.)

        If *middle* is set, computes
        `W_{\mathrm{middle}}(z)` which corresponds to
        `W_{-1}(z)` in the upper half plane and `W_{1}(z)` in the lower half
        plane, connected continuously through `(-1/e,0)` with branch cuts
        on `(-\infty,-1/e)` and `(0,+\infty)`. `W_{\mathrm{middle}}(z)` extends the
        real analytic function `W_{-1}(x)` defined on `(-1/e,0)` to a complex
        analytic function, whereas the standard branch `W_{-1}(z)` has a branch
        cut along the real segment.

            >>> acb(-5,"+/- 1e-20").lambertw()
            [0.844844605432170 +/- 6.18e-16] + [+/- 1.98]j
            >>> acb(-5,"+/- 1e-20").lambertw(left=True)
            [0.844844605432170 +/- 7.85e-16] + [1.97500875488903 +/- 4.05e-15]j
            >>> acb(-5,"+/- 1e-20").lambertw(-1, left=True)
            [0.844844605432170 +/- 7.85e-16] + [-1.97500875488903 +/- 4.05e-15]j
            >>> acb(-0.25,"+/- 1e-20").lambertw(middle=True)
            [-2.15329236411035 +/- 1.48e-15] + [+/- 4.87e-16]j

        """
        cdef int flags
        u = acb.__new__(acb)
        flags = 0
        if middle:
            flags = 4
            branch = -1
        elif left:
            flags = 2
        k = any_as_fmpz(branch)
        acb_lambertw((<acb>u).val, (<acb>s).val, (<fmpz>k).val, flags, getprec())
        return u

    def real_abs(s, bint analytic=False):
        r"""
        Absolute value of a real variable, extended to a piecewise complex
        analytic function. This function is useful for integration.

            >>> f = lambda z, a: (z**3).real_abs(analytic=a)
            >>> showgood(lambda: acb.integral(f, -1, 1), dps=25)
            0.5000000000000000000000000

        """
        u = acb.__new__(acb)
        acb_real_abs((<acb>u).val, (<acb>s).val, analytic, getprec())
        return u

    def real_sgn(s, bint analytic=False):
        r"""
        Sign function of a real variable, extended to a piecewise complex
        analytic function. This function is useful for integration.
        """
        u = acb.__new__(acb)
        acb_real_sgn((<acb>u).val, (<acb>s).val, analytic, getprec())
        return u

    def real_heaviside(s, bint analytic=False):
        r"""
        Heaviside step function of a real variable, extended to a piecewise complex
        analytic function. This function is useful for integration.
        """
        u = acb.__new__(acb)
        acb_real_heaviside((<acb>u).val, (<acb>s).val, analytic, getprec())
        return u

    def real_floor(s, bint analytic=False):
        r"""
        Floor function of a real variable, extended to a piecewise complex
        analytic function. This function is useful for integration.
        """
        u = acb.__new__(acb)
        acb_real_floor((<acb>u).val, (<acb>s).val, analytic, getprec())
        return u

    def real_ceil(s, bint analytic=False):
        r"""
        Ceiling function of a real variable, extended to a piecewise complex
        analytic function. This function is useful for integration.
        """
        u = acb.__new__(acb)
        acb_real_ceil((<acb>u).val, (<acb>s).val, analytic, getprec())
        return u

    def real_max(s, t, bint analytic=False):
        r"""
        Maximum value of two real variables, extended to a piecewise complex
        analytic function. This function is useful for integration.

            >>> f = lambda x, a: x.sin().real_max(x.cos(), analytic=a)
            >>> showgood(lambda: acb.integral(f, 0, acb.pi()))
            2.41421356237310

        """
        u = acb.__new__(acb)
        acb_real_max((<acb>u).val, (<acb>s).val, (<acb>t).val, analytic, getprec())
        return u

    def real_min(s, t, bint analytic=False):
        r"""
        Minimum value of two real variables, extended to a piecewise complex
        analytic function. This function is useful for integration.

            >>> f = lambda x, a: x.sin().real_min(x.cos(), analytic=a)
            >>> showgood(lambda: acb.integral(f, 0, acb.pi()))
            -0.414213562373095

        """
        u = acb.__new__(acb)
        acb_real_min((<acb>u).val, (<acb>s).val, (<acb>t).val, analytic, getprec())
        return u

    def real_sqrt(s, bint analytic=False):
        r"""
        Square root of a real variable assumed to be nonnegative, extended
        to a piecewise complex analytic function. This function is
        useful for integration.

            >>> acb.integral(lambda x, a: x.sqrt(analytic=a), 0, 1)
            [0.6666666666667 +/- 4.19e-14] + [+/- 1.12e-16]j
            >>> acb.integral(lambda x, a: x.real_sqrt(analytic=a), 0, 1)
            [0.6666666666667 +/- 4.19e-14]
        """
        u = acb.__new__(acb)
        acb_real_sqrtpos((<acb>u).val, (<acb>s).val, analytic, getprec())
        return u

    @classmethod
    def stieltjes(cls, n, a=1):
        r"""
        Generalized Stieltjes constant `\gamma_n(a)`.

            >>> showgood(lambda: acb.stieltjes(1), dps=25)
            -0.07281584548367672486058638
            >>> showgood(lambda: acb.stieltjes(10**10), dps=25)
            7.588362123713105194822403e+12397849705
            >>> showgood(lambda: acb.stieltjes(1000, 2+3j), dps=25)
            -1.206122870741999199264747e+494 - 1.389205283963836265123849e+494j
        """
        cdef acb u
        n = fmpz(n)
        if n < 0:
            raise ValueError("n must be nonnegative")
        a = acb(a)
        u = acb.__new__(acb)
        acb_dirichlet_stieltjes(u.val, (<fmpz>n).val, (<acb>a).val, getprec())
        return u

    @classmethod
    def integral(cls, func, a, b, params=None,
            rel_tol=None, abs_tol=None,
            deg_limit=None, eval_limit=None, depth_limit=None,
            use_heap=None, verbose=None):

        cdef IntegrationContext ictx = IntegrationContext()
        cdef acb_calc_integrate_opt_t arb_opts
        cdef long cgoal, prec
        cdef mag_t ctol
        cdef arb tmp
        cdef acb ca, cb, res

        ca = any_as_acb(a)
        cb = any_as_acb(b)

        mag_init(ctol)

        ictx.f = func
        ictx.exn_type = None

        acb_calc_integrate_opt_init(arb_opts)
        if deg_limit is not None:
            arb_opts.deg_limit = deg_limit
        if eval_limit is not None:
            arb_opts.eval_limit = eval_limit
        if depth_limit is not None:
            arb_opts.depth_limit = depth_limit
        if use_heap is not None:
            arb_opts.use_heap = use_heap
        if verbose is not None:
            arb_opts.verbose = verbose

        prec = ctx.prec

        if rel_tol is None:
            cgoal = prec
        else:
            tmp = any_as_arb(rel_tol)
            cgoal = arf_abs_bound_lt_2exp_si(arb_midref(tmp.val))
            cgoal = -cgoal

        if abs_tol is None:
            mag_set_ui_2exp_si(ctol, 1, -prec)
        else:
            tmp = any_as_arb(abs_tol)
            arb_get_mag(ctol, tmp.val)

        res = acb.__new__(acb)

        try:
            acb_calc_integrate(
                    res.val,
                    <acb_calc_func_t> acb_calc_func_callback,
                    <void *> ictx,
                    ca.val, cb.val,
                    cgoal, ctol, arb_opts, prec)
        finally:
            mag_clear(ctol)

        if ictx.exn_type is not None:
            raise ictx.exn_type, ictx.exn_obj, ictx.exn_tb

        return res

