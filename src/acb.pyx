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

cdef class acb:

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

    def __repr__(self):
        if ctx.pretty:
            return str(self)
        real = self.real
        imag = self.imag
        if imag.is_zero():
            return "acb(%s)" % repr(real)
        else:
            return "acb(%s, %s)" % (repr(real), repr(imag))

    def __str__(self):
        real = self.real
        imag = self.imag
        if imag.is_zero():
            return str(real)
        elif real.is_zero():
            return str(imag) + "j"
        else:
            if arf_sgn(arb_midref((<arb>imag).val)) < 0:
                return "%s - %sj" % (str(real), str(imag)[1:])
            else:
                return "%s + %sj" % (str(real), str(imag))

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
            0.151904002670036137448161 + 0.01980488016185498197191013j
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
            0.0
            >>> print(acb(-1).rgamma())
            0.0

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
            4.0
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
            -2.95305957208855672287624 + 3.410962524512050603254574j
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

            >>> showgood(lambda: acb.const_pi().real, dps=25)
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
            (3.16577851321617 + 1.95960104142161j, 2.03272300701967 - 3.0518977991518j)
        """
        u = acb.__new__(acb)
        v = acb.__new__(acb)
        acb_sin_cos((<acb>u).val, (<acb>v).val, (<acb>s).val, getprec())
        return u, v

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

    def rising_ui(s, ulong n):
        """
        Computes the rising factorial `(s)_n` where *n* is an unsigned
        integer. The current implementation does not use the gamma function,
        so *n* should be moderate.

            >>> showgood(lambda: acb(1,2).rising_ui(5), dps=25)
            -540.0 - 100.0j
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
            (-540.0 - 100.0j, -666.0 + 420.0j)
        """
        u = acb.__new__(acb)
        v = acb.__new__(acb)
        acb_rising2_ui((<acb>u).val, (<acb>v).val, (<acb>s).val, n, getprec())
        return u, v

    @classmethod
    def polylog(cls, acb s, acb z):
        """
        Computes the polylogarithm `\operatorname{Li}_s(z)`.

            >>> showgood(lambda: acb.polylog(acb(2), acb(3)), dps=25)
            2.320180423313098396406194 - 3.451392295223202661433821j
            >>> showgood(lambda: acb.polylog(acb(1,2), acb(2,3)), dps=25)
            -6.854984251829306907116775 + 7.375884252099214498443165j
        """
        u = acb.__new__(acb)
        acb_polylog((<acb>u).val, (<acb>s).val, (<acb>z).val, getprec())
        return u

