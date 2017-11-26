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

def any_as_acb(x):
    if typecheck(x, acb):
        return x
    return acb(x)

def any_as_arb_or_acb(x):
    if typecheck(x, arb) or typecheck(x, acb):
        return x
    try:
        return arb(x)
    except (TypeError, ValueError):
        return acb(x)

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

    def mid(self):
        cdef acb u = acb()
        arf_set(arb_midref(acb_realref(u.val)), arb_midref(acb_realref(self.val)))
        arf_set(arb_midref(acb_imagref(u.val)), arb_midref(acb_imagref(self.val)))
        return u

    def rad(self):
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
        return complex(arf_get_d(arb_midref(acb_realref(self.val)), ARF_RND_DOWN),
                       arf_get_d(arb_midref(acb_imagref(self.val)), ARF_RND_DOWN))

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

    def __div__(s, t):
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

    def log(s):
        r"""
        Computes the natural logarithm `\log(s)`.

            >>> showgood(lambda: acb(1,2).log(), dps=25)
            0.8047189562170501873003797 + 1.107148717794090503017065j
            >>> showgood(lambda: acb(-5).log(), dps=25)
            1.609437912434100374600759 + 3.141592653589793238462643j
        """
        u = acb.__new__(acb)
        acb_log((<acb>u).val, (<acb>s).val, getprec())
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
        Computes the Riemann zeta function `\zeta(s)` or the Hurwitz
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
    def const_pi(cls):
        """
        Computes the constant `\pi`.

            >>> showgood(lambda: acb.const_pi(), dps=25)
            3.141592653589793238462643
            >>> showgood(lambda: acb.pi(), dps=25)  # alias
            3.141592653589793238462643
        """
        u = acb.__new__(acb)
        acb_const_pi((<acb>u).val, getprec())
        return u

    pi = const_pi

    def sqrt(s):
        r"""
        Computes the square root `\sqrt{s}`.

            >>> showgood(lambda: acb(1,2).sqrt(), dps=25)
            1.272019649514068964252422 + 0.7861513777574232860695586j
        """
        u = acb.__new__(acb)
        acb_sqrt((<acb>u).val, (<acb>s).val, getprec())
        return u

    def rsqrt(s):
        r"""
        Computes the reciprocal square root `1/\sqrt{s}`.

            >>> showgood(lambda: acb(1,2).rsqrt(), dps=25)
            0.5688644810057831072783079 - 0.3515775842541429284870573j
        """
        u = acb.__new__(acb)
        acb_rsqrt((<acb>u).val, (<acb>s).val, getprec())
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

    def sinc(s):
        u = acb.__new__(acb)
        acb_sinc((<acb>u).val, (<acb>s).val, getprec())
        return u

    def sinc_pi(s):
        u = acb.__new__(acb)
        acb_sinc_pi((<acb>u).val, (<acb>s).val, getprec())
        return u

    def rising_ui(s, ulong n):
        """
        Computes the rising factorial `(s)_n` where *n* is an unsigned
        integer. The current implementation does not use the gamma function,
        so *n* should be moderate.

            >>> showgood(lambda: acb(1,2).rising_ui(5), dps=25)
            -540.0000000000000000000000 - 100.0000000000000000000000j
        """
        u = acb.__new__(acb)
        acb_rising_ui((<acb>u).val, (<acb>s).val, n, getprec())
        return u

    def rising2_ui(s, ulong n):
        """
        Computes the rising factorial `(s)_n` where *n* is an unsigned
        integer, along with the first derivative with respect to `(s)_n`.
        The current implementation does not use the gamma function,
        so *n* should be moderate.

            >>> showgood(lambda: acb(1,2).rising2_ui(5), dps=25)
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

        Currently, the implementation only uses direct summation
        of the hypergeometric series. No analytic continuation is performed,
        and no asymptotic expansions are used.
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
    def hypgeom_m(cls, a, b, z, bint regularized=False):
        r"""
        Computes Kummer's confluent hypergeometric function `M(a,b,z)`
        given complex numbers `a`, `b`, `z`. Optionally, computes
        the regularized version.

            >>> showgood(lambda: acb.hypgeom_m(2+3j, 3+4j, 40000+50000j), dps=25)
            3.730925582634533963357515e+17366 + 3.199717318207534638202987e+17367j
            >>> showgood(lambda: acb.hypgeom_m(2+3j, 3+4j, 40000+50000j) / acb(3+4j).gamma(), dps=25)
            -1.846160890579724375436801e+17368 + 2.721369772032882467996588e+17367j
            >>> showgood(lambda: acb.hypgeom_m(5, -3, 10, regularized=True), dps=25)
            832600407043.6938843410086
            >>> showgood(lambda: acb.hypgeom_m(-5,-6,10), dps=25)
            403.7777777777777777777778
            >>> showgood(lambda: acb.hypgeom_m(-5,-5,10,regularized=True), dps=25)
            0
            >>> showgood(lambda: acb.hypgeom_m(-5,-6,10,regularized=True), dps=25)
            0
            >>> showgood(lambda: acb.hypgeom_m(-5,-4,10,regularized=True), dps=25)
            -100000.0000000000000000000

        """
        a = any_as_acb(a)
        b = any_as_acb(b)
        z = any_as_acb(z)
        u = acb.__new__(acb)
        acb_hypgeom_m((<acb>u).val, (<acb>a).val, (<acb>b).val, (<acb>z).val, regularized, getprec())
        return u

    hypgeom_1f1 = hypgeom_m

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
    def bessel_k(cls, a, z):
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
        acb_hypgeom_bessel_k((<acb>u).val, (<acb>a).val, (<acb>z).val, getprec())
        return u

    @classmethod
    def bessel_i(cls, a, z):
        r"""
        Computes the modified Bessel function of the first kind `I_a(z)`.

        """
        a = any_as_acb(a)
        z = any_as_acb(z)
        u = acb.__new__(acb)
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

    @classmethod
    def legendre_p(cls, n, m, z, int type_=0):
        r"""
        """
        n = any_as_acb(n)
        m = any_as_acb(m)
        z = any_as_acb(z)
        u = acb.__new__(acb)
        acb_hypgeom_legendre_p((<acb>u).val, (<acb>n).val, (<acb>m).val,
            (<acb>z).val, type_, getprec())
        return u

    @classmethod
    def legendre_q(cls, n, m, z, int type_=0):
        r"""
        """
        n = any_as_acb(n)
        m = any_as_acb(m)
        z = any_as_acb(z)
        u = acb.__new__(acb)
        acb_hypgeom_legendre_q((<acb>u).val, (<acb>n).val, (<acb>m).val,
            (<acb>z).val, type_, getprec())
        return u

    @classmethod
    def spherical_y(cls, n, m, theta, phi):
        r"""
        """
        theta = any_as_acb(theta)
        phi = any_as_acb(phi)
        u = acb.__new__(acb)
        acb_hypgeom_spherical_y((<acb>u).val, n, m,
            (<acb>theta).val, (<acb>phi).val, getprec())
        return u

    @classmethod
    def jacobi_p(cls, n, a, b, z):
        r"""
        """
        n = any_as_acb(n)
        a = any_as_acb(a)
        b = any_as_acb(b)
        z = any_as_acb(z)
        u = acb.__new__(acb)
        acb_hypgeom_jacobi_p((<acb>u).val, (<acb>n).val, (<acb>a).val, (<acb>b).val,
            (<acb>z).val, getprec())
        return u

    @classmethod
    def gegenbauer_c(cls, n, m, z):
        r"""
        """
        n = any_as_acb(n)
        m = any_as_acb(m)
        z = any_as_acb(z)
        u = acb.__new__(acb)
        acb_hypgeom_gegenbauer_c((<acb>u).val, (<acb>n).val, (<acb>m).val,
            (<acb>z).val, getprec())
        return u

    @classmethod
    def laguerre_l(cls, n, m, z):
        r"""
        """
        n = any_as_acb(n)
        m = any_as_acb(m)
        z = any_as_acb(z)
        u = acb.__new__(acb)
        acb_hypgeom_laguerre_l((<acb>u).val, (<acb>n).val, (<acb>m).val,
            (<acb>z).val, getprec())
        return u

    @classmethod
    def hermite_h(cls, n, z):
        r"""
        """
        n = any_as_acb(n)
        z = any_as_acb(z)
        u = acb.__new__(acb)
        acb_hypgeom_hermite_h((<acb>u).val, (<acb>n).val, (<acb>z).val, getprec())
        return u

    @classmethod
    def chebyshev_t(cls, n, z):
        r"""
        """
        n = any_as_acb(n)
        z = any_as_acb(z)
        u = acb.__new__(acb)
        acb_hypgeom_chebyshev_t((<acb>u).val, (<acb>n).val, (<acb>z).val, getprec())
        return u

    @classmethod
    def chebyshev_u(cls, n, z):
        r"""
        """
        n = any_as_acb(n)
        z = any_as_acb(z)
        u = acb.__new__(acb)
        acb_hypgeom_chebyshev_u((<acb>u).val, (<acb>n).val, (<acb>z).val, getprec())
        return u

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

    def airy_ai(s):
        u = acb.__new__(acb)
        acb_hypgeom_airy((<acb>u).val, NULL, NULL, NULL, (<acb>s).val, getprec())
        return u

    def airy_ai_prime(s):
        u = acb.__new__(acb)
        acb_hypgeom_airy(NULL, (<acb>u).val, NULL, NULL, (<acb>s).val, getprec())
        return u

    def airy_bi(s):
        u = acb.__new__(acb)
        acb_hypgeom_airy(NULL, NULL, (<acb>u).val, NULL, (<acb>s).val, getprec())
        return u

    def airy_bi_prime(s):
        u = acb.__new__(acb)
        acb_hypgeom_airy(NULL, NULL, NULL, (<acb>u).val, (<acb>s).val, getprec())
        return u

    def airy(s):
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

    def lambertw(s, branch=0, int flags=0):
        u = acb.__new__(acb)
        k = any_as_fmpz(branch)
        acb_lambertw((<acb>u).val, (<acb>s).val, (<fmpz>k).val, flags, getprec())
        return u

