cdef acb_series_coerce_operands(x, y):
    if typecheck(x, acb_series):
        if isinstance(y, (int, long, float, complex, fmpz, fmpz_poly, fmpz_series, fmpq, fmpq_poly, fmpq_series, arb, arb_poly, arb_series, acb, acb_poly)):
            return x, acb_series(y)
    else:
        if isinstance(x, (int, long, float, complex, fmpz, fmpz_poly, fmpz_series, fmpq, fmpq_poly, fmpq_series, arb, arb_poly, arb_series, acb, acb_poly)):
            return acb_series(x), y
    return NotImplemented, NotImplemented

cdef class acb_series(flint_series):

    cdef acb_poly_t val
    cdef long prec

    def __cinit__(self):
        acb_poly_init(self.val)
        self.prec = 0

    def __dealloc__(self):
        acb_poly_clear(self.val)

    def __init__(self, val=None, prec=None):
        if prec is None:
            self.prec = getcap()
        else:
            self.prec = prec
        if self.prec < 0:
            self.prec = -1
        if val is not None:
            if typecheck(val, acb_series):
                acb_poly_set(self.val, (<acb_series>val).val)
                self.prec = min((<fmpz_series>val).prec, getcap())
            elif typecheck(val, fmpz_series):
                acb_poly_set_fmpz_poly(self.val, (<fmpz_series>val).val, getprec())
                self.prec = min((<fmpz_series>val).prec, getcap())
            elif typecheck(val, fmpz_poly):
                acb_poly_set_fmpz_poly(self.val, (<fmpz_poly>val).val, getprec())
            elif typecheck(val, acb_poly):
                acb_poly_set(self.val, (<acb_poly>val).val)
            elif typecheck(val, list):
                acb_poly_set_list(self.val, val, getprec())
            else:
                acb_poly_set_list(self.val, [val], getprec())
        acb_poly_truncate(self.val, max(0, self.prec))

    def __len__(self):
        return acb_poly_length(self.val)

    cpdef long length(self):
        return acb_poly_length(self.val)

    def __getitem__(self, long i):
        cdef acb x
        x = acb()
        if i < 0:
            return x
        acb_poly_get_coeff_acb(x.val, self.val, i)
        return x

    def __setitem__(self, long i, x):
        if i < 0:
            raise ValueError("cannot assign to index < 0 of polynomial")
        if not typecheck(x, acb):
            x = acb(x)
        acb_poly_set_coeff_acb(self.val, i, (<acb>x).val)

    def repr(self, **kwargs):
        return "acb_series([%s], prec=%s)" % (", ".join(map(str, self)), self.prec)

    def str(self, **kwargs):
        if self.prec > 0:
            s = acb_poly(list(self)).str(ascending=True)
            return s + (" + O(x^%s)" % self.prec)
        elif self.prec == 0:
            return "O(x^0)"
        else:
            return "(invalid power series)"

    def __pos__(self):
        return self

    def __neg__(s):
        cdef long cap
        u = acb_series.__new__(acb_series)
        cap = getcap()
        cap = min(cap, (<acb_series>s).prec)
        if cap > 0:
            acb_poly_neg((<acb_series>u).val, (<acb_series>s).val)
            acb_poly_truncate((<acb_series>u).val, cap)
        (<acb_series>u).prec = cap
        return u

    def __add__(s, t):
        cdef long cap
        if type(s) is type(t):
            u = acb_series.__new__(acb_series)
            cap = getcap()
            cap = min(cap, (<acb_series>s).prec)
            cap = min(cap, (<acb_series>t).prec)
            if cap > 0:
                acb_poly_add((<acb_series>u).val, (<acb_series>s).val, (<acb_series>t).val, getprec())
                acb_poly_truncate((<acb_series>u).val, cap)
            (<acb_series>u).prec = cap
            return u
        s, t = acb_series_coerce_operands(s, t)
        if s is NotImplemented:
            return s
        return s + t

    def __sub__(s, t):
        cdef long cap
        if type(s) is type(t):
            u = acb_series.__new__(acb_series)
            cap = getcap()
            cap = min(cap, (<acb_series>s).prec)
            cap = min(cap, (<acb_series>t).prec)
            if cap > 0:
                acb_poly_sub((<acb_series>u).val, (<acb_series>s).val, (<acb_series>t).val, getprec())
                acb_poly_truncate((<acb_series>u).val, cap)
            (<acb_series>u).prec = cap
            return u
        s, t = acb_series_coerce_operands(s, t)
        if s is NotImplemented:
            return s
        return s - t

    def __mul__(s, t):
        cdef long cap
        if type(s) is type(t):
            u = acb_series.__new__(acb_series)
            cap = getcap()
            cap = min(cap, (<acb_series>s).prec)
            cap = min(cap, (<acb_series>t).prec)
            if cap > 0:
                acb_poly_mullow((<acb_series>u).val, (<acb_series>s).val, (<acb_series>t).val, cap, getprec())
            (<acb_series>u).prec = cap
            return u
        s, t = acb_series_coerce_operands(s, t)
        if s is NotImplemented:
            return s
        return s * t

    cpdef valuation(self):
        cdef long i
        if self.length() == 0:
            return -1
        i = 0
        while acb_is_zero(&(self.val.coeffs[i])):
            i += 1
        return i

    def __div__(s, t):
        cdef long cap, sval, tval
        cdef acb_poly_t stmp, ttmp
        if type(s) is type(t):
            cap = getcap()
            cap = min(cap, (<acb_series>s).prec)
            cap = min(cap, (<acb_series>t).prec)

            if (<acb_series>t).length() == 0:
                raise ZeroDivisionError("power series division")

            u = acb_series.__new__(acb_series)

            if (<acb_series>s).length() == 0:
                u.cap = cap
                return u

            sval = (<acb_series>s).valuation()
            tval = (<acb_series>t).valuation()

            if sval < tval:
                raise ValueError("quotient would not be a power series")

            if acb_contains_zero(&((<acb_series>t).val.coeffs[tval])):
                raise ValueError("leading term in denominator is not nonzero")

            if tval == 0:
                acb_poly_div_series((<acb_series>u).val, (<acb_series>s).val, (<acb_series>t).val, cap, getprec())
            else:
                acb_poly_init(stmp)
                acb_poly_init(ttmp)
                acb_poly_shift_right(stmp, (<acb_series>s).val, tval)
                acb_poly_shift_right(ttmp, (<acb_series>t).val, tval)
                cap -= tval
                acb_poly_div_series((<acb_series>u).val, stmp, ttmp, cap, getprec())
                acb_poly_clear(stmp)
                acb_poly_clear(ttmp)

            (<acb_series>u).prec = cap
            return u

        s, t = acb_series_coerce_operands(s, t)
        if s is NotImplemented:
            return s
        return s / t

    def __truediv__(s, t):
        return acb_series.__div__(s, t)

    def __pow__(s, t, mod):
        cdef long cap
        if mod is not None:
            raise NotImplementedError("modular exponentiation")
        if type(s) is type(t):
            u = acb_series.__new__(acb_series)
            cap = getcap()
            cap = min(cap, (<acb_series>s).prec)
            cap = min(cap, (<acb_series>t).prec)
            if cap > 0:
                acb_poly_pow_series((<acb_series>u).val, (<acb_series>s).val, (<acb_series>t).val, cap, getprec())
            (<acb_series>u).prec = cap
            return u
        s, t = acb_series_coerce_operands(s, t)
        if s is NotImplemented:
            return s
        return s ** t

    def __call__(s, t):
        cdef long cap
        if typecheck(t, acb_series):
            u = acb_series.__new__(acb_series)
            if (<acb_series>t).valuation() < 1:
                raise ValueError("power series composition with nonzero constant term")
            cap = getcap()
            cap = min(cap, (<acb_series>s).prec)
            cap = min(cap, (<acb_series>t).prec)
            acb_poly_compose_series((<acb_series>u).val, (<acb_series>s).val, (<acb_series>t).val, cap, getprec())
            (<acb_series>u).prec = cap
            return u
        raise TypeError("cannot call acb_series with input of type %s", type(t))

    def reversion(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<acb_series>s).prec)
        if s.length() < 2 or (not acb_is_zero(&s.val.coeffs[0])) or \
            (acb_contains_zero(&s.val.coeffs[1])):
                raise ValueError("power series reversion requires valuation 1")
        u = acb_series.__new__(acb_series)
        acb_poly_revert_series((<acb_series>u).val, (<acb_series>s).val, cap, getprec())
        (<acb_series>u).prec = cap
        return u

    def inv(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<acb_series>s).prec)
        if s.length() == 0:
            raise ZeroDivisionError
        u = acb_series.__new__(acb_series)
        acb_poly_inv_series((<acb_series>u).val, (<acb_series>s).val, cap, getprec())
        (<acb_series>u).prec = cap
        return u

    def derivative(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<acb_series>s).prec - 1)
        u = acb_series.__new__(acb_series)
        acb_poly_derivative((<acb_series>u).val, (<acb_series>s).val, getprec())
        acb_poly_truncate((<acb_series>u).val, max(0, cap))
        (<acb_series>u).prec = cap
        return u

    def integral(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<acb_series>s).prec + 1)
        u = acb_series.__new__(acb_series)
        acb_poly_integral((<acb_series>u).val, (<acb_series>s).val, getprec())
        acb_poly_truncate((<acb_series>u).val, max(0, cap))
        (<acb_series>u).prec = cap
        return u

    def sqrt(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<acb_series>s).prec)
        u = acb_series.__new__(acb_series)
        acb_poly_sqrt_series((<acb_series>u).val, (<acb_series>s).val, cap, getprec())
        (<acb_series>u).prec = cap
        return u

    def rsqrt(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<acb_series>s).prec)
        u = acb_series.__new__(acb_series)
        acb_poly_rsqrt_series((<acb_series>u).val, (<acb_series>s).val, cap, getprec())
        (<acb_series>u).prec = cap
        return u

    def exp(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<acb_series>s).prec)
        u = acb_series.__new__(acb_series)
        acb_poly_exp_series((<acb_series>u).val, (<acb_series>s).val, cap, getprec())
        (<acb_series>u).prec = cap
        return u

    def log(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<acb_series>s).prec)
        u = acb_series.__new__(acb_series)
        acb_poly_log_series((<acb_series>u).val, (<acb_series>s).val, cap, getprec())
        (<acb_series>u).prec = cap
        return u

    def atan(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<acb_series>s).prec)
        u = acb_series.__new__(acb_series)
        acb_poly_atan_series((<acb_series>u).val, (<acb_series>s).val, cap, getprec())
        (<acb_series>u).prec = cap
        return u

    def sin(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<acb_series>s).prec)
        u = acb_series.__new__(acb_series)
        acb_poly_sin_series((<acb_series>u).val, (<acb_series>s).val, cap, getprec())
        (<acb_series>u).prec = cap
        return u

    def cos(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<acb_series>s).prec)
        u = acb_series.__new__(acb_series)
        acb_poly_cos_series((<acb_series>u).val, (<acb_series>s).val, cap, getprec())
        (<acb_series>u).prec = cap
        return u

    def sin_cos(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<acb_series>s).prec)
        u = acb_series.__new__(acb_series)
        v = acb_series.__new__(acb_series)
        acb_poly_sin_cos_series((<acb_series>u).val, (<acb_series>v).val, (<acb_series>s).val, cap, getprec())
        (<acb_series>u).prec = cap
        (<acb_series>v).prec = cap
        return u, v

    def tan(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<acb_series>s).prec)
        u = acb_series.__new__(acb_series)
        acb_poly_tan_series((<acb_series>u).val, (<acb_series>s).val, cap, getprec())
        (<acb_series>u).prec = cap
        return u

    def gamma(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<acb_series>s).prec)
        u = acb_series.__new__(acb_series)
        acb_poly_gamma_series((<acb_series>u).val, (<acb_series>s).val, cap, getprec())
        (<acb_series>u).prec = cap
        return u

    def rgamma(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<acb_series>s).prec)
        u = acb_series.__new__(acb_series)
        acb_poly_rgamma_series((<acb_series>u).val, (<acb_series>s).val, cap, getprec())
        (<acb_series>u).prec = cap
        return u

    def lgamma(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<acb_series>s).prec)
        u = acb_series.__new__(acb_series)
        acb_poly_lgamma_series((<acb_series>u).val, (<acb_series>s).val, cap, getprec())
        (<acb_series>u).prec = cap
        return u

    def rising_ui(s, ulong n):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<acb_series>s).prec)
        u = acb_series.__new__(acb_series)
        acb_poly_rising_ui_series((<acb_series>u).val, (<acb_series>s).val, n, cap, getprec())
        (<acb_series>u).prec = cap
        return u

    def zeta(s, a=1, bint deflate=0):
        cdef long cap
        a = acb(a)
        cap = getcap()
        cap = min(cap, (<acb_series>s).prec)
        u = acb_series.__new__(acb_series)
        acb_poly_zeta_series((<acb_series>u).val, (<acb_series>s).val, (<acb>a).val, deflate, cap, getprec())
        (<acb_series>u).prec = cap
        return u

    @classmethod
    def polylog(cls, s, z):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<acb_series>s).prec)
        s = acb_series(s)
        z = acb(z)
        u = acb_series.__new__(acb_series)
        acb_poly_polylog_series((<acb_series>u).val, (<acb_series>s).val, (<acb>z).val, cap, getprec())
        (<acb_series>u).prec = cap
        return u

    def agm(s, t=None):
        cdef long cap
        if t is None:
            cap = getcap()
            cap = min(cap, (<acb_series>s).prec)
            u = acb_series.__new__(acb_series)
            acb_poly_agm1_series((<acb_series>u).val, (<acb_series>s).val, cap, getprec())
            (<acb_series>u).prec = cap
            return u
        else:
            return (s / t).agm() * t

    def elliptic_k(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<acb_series>s).prec)
        u = acb_series.__new__(acb_series)
        acb_poly_elliptic_k_series((<acb_series>u).val, (<acb_series>s).val, cap, getprec())
        (<acb_series>u).prec = cap
        return u

    def elliptic_p(s, tau):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<acb_series>s).prec)
        tau = acb(tau)
        u = acb_series.__new__(acb_series)
        acb_poly_elliptic_p_series((<acb_series>u).val, (<acb_series>s).val, (<acb>tau).val, cap, getprec())
        (<acb_series>u).prec = cap
        return u

    def erf(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<acb_series>s).prec)
        u = acb_series.__new__(acb_series)
        acb_poly_erf_series((<acb_series>u).val, (<acb_series>s).val, cap, getprec())
        (<acb_series>u).prec = cap
        return u

    @classmethod
    def gamma_upper(cls, s, z):
        cdef long cap
        s = acb(s)
        z = acb_series(z)
        cap = getcap()
        cap = min(cap, (<acb_series>z).prec)
        u = acb_series.__new__(acb_series)
        acb_poly_gamma_upper_series((<acb_series>u).val, (<acb>s).val, (<acb_series>z).val, cap, getprec())
        (<acb_series>u).prec = cap
        return u

