cdef arb_series_coerce_operands(x, y):
    if typecheck(x, arb_series):
        if isinstance(y, (int, long, float, fmpz, fmpz_poly, fmpz_series, fmpq, fmpq_poly, fmpq_series, arb, arb_poly)):
            return x, arb_series(y)
        if isinstance(y, (complex, acb, acb_poly, acb_series)):
            return acb_series(x), acb_series(y)
    else:
        if isinstance(x, (int, long, float, fmpz, fmpz_poly, fmpz_series, fmpq, fmpq_poly, fmpq_series, arb, arb_poly)):
            return arb_series(x), y
        if isinstance(x, (complex, acb, acb_poly, acb_series)):
            return acb_series(x), acb_series(y)
    return NotImplemented, NotImplemented

cdef class arb_series(flint_series):

    cdef arb_poly_t val
    cdef long prec

    def __cinit__(self):
        arb_poly_init(self.val)
        self.prec = 0

    def __dealloc__(self):
        arb_poly_clear(self.val)

    def __init__(self, val=None, prec=None):
        if prec is None:
            self.prec = getcap()
        else:
            self.prec = prec
        if self.prec < 0:
            self.prec = -1
        if val is not None:
            if typecheck(val, arb_series):
                arb_poly_set(self.val, (<arb_series>val).val)
                self.prec = min((<fmpz_series>val).prec, getcap())
            elif typecheck(val, fmpz_series):
                arb_poly_set_fmpz_poly(self.val, (<fmpz_series>val).val, getprec())
                self.prec = min((<fmpz_series>val).prec, getcap())
            elif typecheck(val, fmpz_poly):
                arb_poly_set_fmpz_poly(self.val, (<fmpz_poly>val).val, getprec())
            elif typecheck(val, arb_poly):
                arb_poly_set(self.val, (<arb_poly>val).val)
            elif typecheck(val, list):
                arb_poly_set_list(self.val, val, getprec())
            else:
                arb_poly_set_list(self.val, [val], getprec())
        arb_poly_truncate(self.val, max(0, self.prec))

    def __len__(self):
        return arb_poly_length(self.val)

    cpdef long length(self):
        return arb_poly_length(self.val)

    def __getitem__(self, long i):
        cdef arb x
        x = arb()
        if i < 0:
            return x
        arb_poly_get_coeff_arb(x.val, self.val, i)
        return x

    def __setitem__(self, long i, x):
        if i < 0:
            raise ValueError("cannot assign to index < 0 of polynomial")
        if not typecheck(x, arb):
            x = arb(x)
        arb_poly_set_coeff_arb(self.val, i, (<arb>x).val)

    def repr(self, **kwargs):
        return "arb_series([%s], prec=%s)" % (", ".join(map(str, self)), self.prec)

    def str(self, **kwargs):
        if self.prec > 0:
            s = arb_poly(list(self)).str(ascending=True)
            return s + (" + O(x^%s)" % self.prec)
        elif self.prec == 0:
            return "O(x^0)"
        else:
            return "(invalid power series)"

    def __pos__(self):
        return self

    def __neg__(s):
        cdef long cap
        u = arb_series.__new__(arb_series)
        cap = getcap()
        cap = min(cap, (<arb_series>s).prec)
        if cap > 0:
            arb_poly_neg((<arb_series>u).val, (<arb_series>s).val)
            arb_poly_truncate((<arb_series>u).val, cap)
        (<arb_series>u).prec = cap
        return u

    def __add__(s, t):
        cdef long cap
        if type(s) is type(t):
            u = arb_series.__new__(arb_series)
            cap = getcap()
            cap = min(cap, (<arb_series>s).prec)
            cap = min(cap, (<arb_series>t).prec)
            if cap > 0:
                arb_poly_add((<arb_series>u).val, (<arb_series>s).val, (<arb_series>t).val, getprec())
                arb_poly_truncate((<arb_series>u).val, cap)
            (<arb_series>u).prec = cap
            return u
        s, t = arb_series_coerce_operands(s, t)
        if s is NotImplemented:
            return s
        return s + t

    def __sub__(s, t):
        cdef long cap
        if type(s) is type(t):
            u = arb_series.__new__(arb_series)
            cap = getcap()
            cap = min(cap, (<arb_series>s).prec)
            cap = min(cap, (<arb_series>t).prec)
            if cap > 0:
                arb_poly_sub((<arb_series>u).val, (<arb_series>s).val, (<arb_series>t).val, getprec())
                arb_poly_truncate((<arb_series>u).val, cap)
            (<arb_series>u).prec = cap
            return u
        s, t = arb_series_coerce_operands(s, t)
        if s is NotImplemented:
            return s
        return s - t

    def __mul__(s, t):
        cdef long cap
        if type(s) is type(t):
            u = arb_series.__new__(arb_series)
            cap = getcap()
            cap = min(cap, (<arb_series>s).prec)
            cap = min(cap, (<arb_series>t).prec)
            if cap > 0:
                arb_poly_mullow((<arb_series>u).val, (<arb_series>s).val, (<arb_series>t).val, cap, getprec())
            (<arb_series>u).prec = cap
            return u
        s, t = arb_series_coerce_operands(s, t)
        if s is NotImplemented:
            return s
        return s * t

    cpdef valuation(self):
        cdef long i
        if self.length() == 0:
            return -1
        i = 0
        while arb_is_zero(&(self.val.coeffs[i])):
            i += 1
        return i

    def __div__(s, t):
        cdef long cap, sval, tval
        cdef arb_poly_t stmp, ttmp
        if type(s) is type(t):
            cap = getcap()
            cap = min(cap, (<arb_series>s).prec)
            cap = min(cap, (<arb_series>t).prec)

            if (<arb_series>t).length() == 0:
                raise ZeroDivisionError("power series division")

            u = arb_series.__new__(arb_series)

            if (<arb_series>s).length() == 0:
                u.cap = cap
                return u

            sval = (<arb_series>s).valuation()
            tval = (<arb_series>t).valuation()

            if sval < tval:
                raise ValueError("quotient would not be a power series")

            if not arb_is_nonzero(&((<arb_series>t).val.coeffs[tval])):
                raise ValueError("leading term in denominator is not nonzero")

            if tval == 0:
                arb_poly_div_series((<arb_series>u).val, (<arb_series>s).val, (<arb_series>t).val, cap, getprec())
            else:
                arb_poly_init(stmp)
                arb_poly_init(ttmp)
                arb_poly_shift_right(stmp, (<arb_series>s).val, tval)
                arb_poly_shift_right(ttmp, (<arb_series>t).val, tval)
                cap -= tval
                arb_poly_div_series((<arb_series>u).val, stmp, ttmp, cap, getprec())
                arb_poly_clear(stmp)
                arb_poly_clear(ttmp)

            (<arb_series>u).prec = cap
            return u

        s, t = arb_series_coerce_operands(s, t)
        if s is NotImplemented:
            return s
        return s / t

    def __truediv__(s, t):
        return arb_series.__div__(s, t)

    def __pow__(s, t, mod):
        cdef long cap
        if mod is not None:
            raise NotImplementedError("modular exponentiation")
        if type(s) is type(t):
            u = arb_series.__new__(arb_series)
            cap = getcap()
            cap = min(cap, (<arb_series>s).prec)
            cap = min(cap, (<arb_series>t).prec)
            if cap > 0:
                arb_poly_pow_series((<arb_series>u).val, (<arb_series>s).val, (<arb_series>t).val, cap, getprec())
            (<arb_series>u).prec = cap
            return u
        s, t = arb_series_coerce_operands(s, t)
        if s is NotImplemented:
            return s
        return s ** t

    def __call__(s, t):
        cdef long cap
        if typecheck(t, arb_series):
            u = arb_series.__new__(arb_series)
            if (<arb_series>t).valuation() < 1:
                raise ValueError("power series composition with nonzero constant term")
            cap = getcap()
            cap = min(cap, (<arb_series>s).prec)
            cap = min(cap, (<arb_series>t).prec)
            arb_poly_compose_series((<arb_series>u).val, (<arb_series>s).val, (<arb_series>t).val, cap, getprec())
            (<arb_series>u).prec = cap
            return u
        raise TypeError("cannot call arb_series with input of type %s", type(t))

    def reversion(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<arb_series>s).prec)
        if s.length() < 2 or (not arb_is_zero(&s.val.coeffs[0])) or \
            (not arb_is_nonzero(&s.val.coeffs[1])):
                raise ValueError("power series reversion requires valuation 1")
        u = arb_series.__new__(arb_series)
        arb_poly_revert_series((<arb_series>u).val, (<arb_series>s).val, cap, getprec())
        (<arb_series>u).prec = cap
        return u

    def inv(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<arb_series>s).prec)
        if s.length() == 0:
            raise ZeroDivisionError
        u = arb_series.__new__(arb_series)
        arb_poly_inv_series((<arb_series>u).val, (<arb_series>s).val, cap, getprec())
        (<arb_series>u).prec = cap
        return u

    def derivative(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<arb_series>s).prec - 1)
        u = arb_series.__new__(arb_series)
        arb_poly_derivative((<arb_series>u).val, (<arb_series>s).val, getprec())
        arb_poly_truncate((<arb_series>u).val, max(0, cap))
        (<arb_series>u).prec = cap
        return u

    def integral(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<arb_series>s).prec + 1)
        u = arb_series.__new__(arb_series)
        arb_poly_integral((<arb_series>u).val, (<arb_series>s).val, getprec())
        arb_poly_truncate((<arb_series>u).val, max(0, cap))
        (<arb_series>u).prec = cap
        return u

    def sqrt(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<arb_series>s).prec)
        u = arb_series.__new__(arb_series)
        arb_poly_sqrt_series((<arb_series>u).val, (<arb_series>s).val, cap, getprec())
        (<arb_series>u).prec = cap
        return u

    def rsqrt(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<arb_series>s).prec)
        u = arb_series.__new__(arb_series)
        arb_poly_rsqrt_series((<arb_series>u).val, (<arb_series>s).val, cap, getprec())
        (<arb_series>u).prec = cap
        return u

    def exp(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<arb_series>s).prec)
        u = arb_series.__new__(arb_series)
        arb_poly_exp_series((<arb_series>u).val, (<arb_series>s).val, cap, getprec())
        (<arb_series>u).prec = cap
        return u

    def log(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<arb_series>s).prec)
        u = arb_series.__new__(arb_series)
        arb_poly_log_series((<arb_series>u).val, (<arb_series>s).val, cap, getprec())
        (<arb_series>u).prec = cap
        return u

    def atan(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<arb_series>s).prec)
        u = arb_series.__new__(arb_series)
        arb_poly_atan_series((<arb_series>u).val, (<arb_series>s).val, cap, getprec())
        (<arb_series>u).prec = cap
        return u

    def asin(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<arb_series>s).prec)
        u = arb_series.__new__(arb_series)
        arb_poly_asin_series((<arb_series>u).val, (<arb_series>s).val, cap, getprec())
        (<arb_series>u).prec = cap
        return u

    def acos(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<arb_series>s).prec)
        u = arb_series.__new__(arb_series)
        arb_poly_acos_series((<arb_series>u).val, (<arb_series>s).val, cap, getprec())
        (<arb_series>u).prec = cap
        return u

    def sin(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<arb_series>s).prec)
        u = arb_series.__new__(arb_series)
        arb_poly_sin_series((<arb_series>u).val, (<arb_series>s).val, cap, getprec())
        (<arb_series>u).prec = cap
        return u

    def cos(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<arb_series>s).prec)
        u = arb_series.__new__(arb_series)
        arb_poly_cos_series((<arb_series>u).val, (<arb_series>s).val, cap, getprec())
        (<arb_series>u).prec = cap
        return u

    def sin_cos(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<arb_series>s).prec)
        u = arb_series.__new__(arb_series)
        v = arb_series.__new__(arb_series)
        arb_poly_sin_cos_series((<arb_series>u).val, (<arb_series>v).val, (<arb_series>s).val, cap, getprec())
        (<arb_series>u).prec = cap
        (<arb_series>v).prec = cap
        return u, v

    def sin_pi(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<arb_series>s).prec)
        u = arb_series.__new__(arb_series)
        arb_poly_sin_pi_series((<arb_series>u).val, (<arb_series>s).val, cap, getprec())
        (<arb_series>u).prec = cap
        return u

    def cos_pi(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<arb_series>s).prec)
        u = arb_series.__new__(arb_series)
        arb_poly_cos_pi_series((<arb_series>u).val, (<arb_series>s).val, cap, getprec())
        (<arb_series>u).prec = cap
        return u

    def sin_cos_pi(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<arb_series>s).prec)
        u = arb_series.__new__(arb_series)
        v = arb_series.__new__(arb_series)
        arb_poly_sin_cos_pi_series((<arb_series>u).val, (<arb_series>v).val, (<arb_series>s).val, cap, getprec())
        (<arb_series>u).prec = cap
        (<arb_series>v).prec = cap
        return u, v

    def cot_pi(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<arb_series>s).prec)
        u = arb_series.__new__(arb_series)
        arb_poly_cot_pi_series((<arb_series>u).val, (<arb_series>s).val, cap, getprec())
        (<arb_series>u).prec = cap
        return u

    def tan(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<arb_series>s).prec)
        u = arb_series.__new__(arb_series)
        arb_poly_tan_series((<arb_series>u).val, (<arb_series>s).val, cap, getprec())
        (<arb_series>u).prec = cap
        return u

    def gamma(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<arb_series>s).prec)
        u = arb_series.__new__(arb_series)
        arb_poly_gamma_series((<arb_series>u).val, (<arb_series>s).val, cap, getprec())
        (<arb_series>u).prec = cap
        return u

    def rgamma(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<arb_series>s).prec)
        u = arb_series.__new__(arb_series)
        arb_poly_rgamma_series((<arb_series>u).val, (<arb_series>s).val, cap, getprec())
        (<arb_series>u).prec = cap
        return u

    def lgamma(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<arb_series>s).prec)
        u = arb_series.__new__(arb_series)
        arb_poly_lgamma_series((<arb_series>u).val, (<arb_series>s).val, cap, getprec())
        (<arb_series>u).prec = cap
        return u

    def rising_ui(s, ulong n):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<arb_series>s).prec)
        u = arb_series.__new__(arb_series)
        arb_poly_rising_ui_series((<arb_series>u).val, (<arb_series>s).val, n, cap, getprec())
        (<arb_series>u).prec = cap
        return u

    def zeta(s, a=1, bint deflate=0):
        cdef long cap
        a = arb(a)
        cap = getcap()
        cap = min(cap, (<arb_series>s).prec)
        u = arb_series.__new__(arb_series)
        arb_poly_zeta_series((<arb_series>u).val, (<arb_series>s).val, (<arb>a).val, deflate, cap, getprec())
        (<arb_series>u).prec = cap
        return u

    def riemann_siegel_theta(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<arb_series>s).prec)
        u = arb_series.__new__(arb_series)
        arb_poly_riemann_siegel_theta_series((<arb_series>u).val, (<arb_series>s).val, cap, getprec())
        (<arb_series>u).prec = cap
        return u

    def riemann_siegel_z(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<arb_series>s).prec)
        u = arb_series.__new__(arb_series)
        arb_poly_riemann_siegel_z_series((<arb_series>u).val, (<arb_series>s).val, cap, getprec())
        (<arb_series>u).prec = cap
        return u

    def erf(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<arb_series>s).prec)
        u = arb_series.__new__(arb_series)
        arb_hypgeom_erf_series((<arb_series>u).val, (<arb_series>s).val, cap, getprec())
        (<arb_series>u).prec = cap
        return u

    def erfc(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<arb_series>s).prec)
        u = arb_series.__new__(arb_series)
        arb_hypgeom_erfc_series((<arb_series>u).val, (<arb_series>s).val, cap, getprec())
        (<arb_series>u).prec = cap
        return u

    def erfi(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<arb_series>s).prec)
        u = arb_series.__new__(arb_series)
        arb_hypgeom_erfi_series((<arb_series>u).val, (<arb_series>s).val, cap, getprec())
        (<arb_series>u).prec = cap
        return u

    def fresnel(s, bint normalized=True):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<arb_series>s).prec)
        u = arb_series.__new__(arb_series)
        v = arb_series.__new__(arb_series)
        arb_hypgeom_fresnel_series((<arb_series>u).val, (<arb_series>v).val, (<arb_series>s).val, normalized, cap, getprec())
        (<arb_series>u).prec = cap
        (<arb_series>v).prec = cap
        return u, v

    def fresnel_s(s, bint normalized=True):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<arb_series>s).prec)
        u = arb_series.__new__(arb_series)
        arb_hypgeom_fresnel_series((<arb_series>u).val, NULL, (<arb_series>s).val, normalized, cap, getprec())
        (<arb_series>u).prec = cap
        return u

    def fresnel_c(s, bint normalized=True):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<arb_series>s).prec)
        u = arb_series.__new__(arb_series)
        arb_hypgeom_fresnel_series(NULL, (<arb_series>u).val, (<arb_series>s).val, normalized, cap, getprec())
        (<arb_series>u).prec = cap
        return u

    def ei(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<arb_series>s).prec)
        u = arb_series.__new__(arb_series)
        arb_hypgeom_ei_series((<arb_series>u).val, (<arb_series>s).val, cap, getprec())
        (<arb_series>u).prec = cap
        return u

    def si(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<arb_series>s).prec)
        u = arb_series.__new__(arb_series)
        arb_hypgeom_si_series((<arb_series>u).val, (<arb_series>s).val, cap, getprec())
        (<arb_series>u).prec = cap
        return u

    def ci(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<arb_series>s).prec)
        u = arb_series.__new__(arb_series)
        arb_hypgeom_ci_series((<arb_series>u).val, (<arb_series>s).val, cap, getprec())
        (<arb_series>u).prec = cap
        return u

    def shi(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<arb_series>s).prec)
        u = arb_series.__new__(arb_series)
        arb_hypgeom_shi_series((<arb_series>u).val, (<arb_series>s).val, cap, getprec())
        (<arb_series>u).prec = cap
        return u

    def chi(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<arb_series>s).prec)
        u = arb_series.__new__(arb_series)
        arb_hypgeom_chi_series((<arb_series>u).val, (<arb_series>s).val, cap, getprec())
        (<arb_series>u).prec = cap
        return u

    def li(s, bint offset=False):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<arb_series>s).prec)
        u = arb_series.__new__(arb_series)
        arb_hypgeom_li_series((<arb_series>u).val, (<arb_series>s).val, offset, cap, getprec())
        (<arb_series>u).prec = cap
        return u

    def airy_ai(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<arb_series>s).prec)
        u = arb_series.__new__(arb_series)
        arb_hypgeom_airy_series((<arb_series>u).val, NULL, NULL, NULL, (<arb_series>s).val, cap, getprec())
        (<arb_series>u).prec = cap
        return u

    def airy_ai_prime(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<arb_series>s).prec)
        u = arb_series.__new__(arb_series)
        arb_hypgeom_airy_series(NULL, (<arb_series>u).val, NULL, NULL, (<arb_series>s).val, cap, getprec())
        (<arb_series>u).prec = cap
        return u

    def airy_bi(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<arb_series>s).prec)
        u = arb_series.__new__(arb_series)
        arb_hypgeom_airy_series(NULL, NULL, (<arb_series>u).val, NULL, (<arb_series>s).val, cap, getprec())
        (<arb_series>u).prec = cap
        return u

    def airy_bi_prime(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<arb_series>s).prec)
        u = arb_series.__new__(arb_series)
        arb_hypgeom_airy_series(NULL, NULL, NULL, (<arb_series>u).val,(<arb_series>s).val, cap, getprec())
        (<arb_series>u).prec = cap
        return u

    def airy(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<arb_series>s).prec)
        u = arb_series.__new__(arb_series)
        v = arb_series.__new__(arb_series)
        w = arb_series.__new__(arb_series)
        z = arb_series.__new__(arb_series)
        arb_hypgeom_airy_series((<arb_series>u).val, (<arb_series>v).val, (<arb_series>w).val, (<arb_series>z).val, (<arb_series>s).val, cap, getprec())
        (<arb_series>u).prec = cap
        (<arb_series>v).prec = cap
        (<arb_series>w).prec = cap
        (<arb_series>z).prec = cap
        return u, v, w, z

    @classmethod
    def gamma_upper(cls, s, z, int regularized=0):
        cdef long cap
        s = arb(s)
        z = arb_series(z)
        cap = getcap()
        cap = min(cap, (<arb_series>z).prec)
        u = arb_series.__new__(arb_series)
        arb_hypgeom_gamma_upper_series((<arb_series>u).val, (<arb>s).val, (<arb_series>z).val, regularized, cap, getprec())
        (<arb_series>u).prec = cap
        return u

    @classmethod
    def gamma_lower(cls, s, z, int regularized=0):
        cdef long cap
        s = arb(s)
        z = arb_series(z)
        cap = getcap()
        cap = min(cap, (<arb_series>z).prec)
        u = arb_series.__new__(arb_series)
        arb_hypgeom_gamma_lower_series((<arb_series>u).val, (<arb>s).val, (<arb_series>z).val, regularized, cap, getprec())
        (<arb_series>u).prec = cap
        return u

    @classmethod
    def beta_lower(cls, a, b, z, int regularized=0):
        cdef long cap
        a = arb(a)
        b = arb(b)
        z = arb_series(z)
        cap = getcap()
        cap = min(cap, (<arb_series>z).prec)
        u = arb_series.__new__(arb_series)
        arb_hypgeom_beta_lower_series((<arb_series>u).val, (<arb>a).val, (<arb>b).val, (<arb_series>z).val, regularized, cap, getprec())
        (<arb_series>u).prec = cap
        return u
