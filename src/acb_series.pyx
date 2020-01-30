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
                self.prec = min((<acb_series>val).prec, getcap())
            elif typecheck(val, arb_series):
                acb_poly_set_arb_poly(self.val, (<arb_series>val).val)
                self.prec = min((<arb_series>val).prec, getcap())
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

    @staticmethod
    def _truediv_(s, t):
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
        return acb_series._div_(s, t)

    def __div__(s, t):
        return acb_series._div_(s, t)

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

    def sin_pi(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<acb_series>s).prec)
        u = acb_series.__new__(acb_series)
        acb_poly_sin_pi_series((<acb_series>u).val, (<acb_series>s).val, cap, getprec())
        (<acb_series>u).prec = cap
        return u

    def cos_pi(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<acb_series>s).prec)
        u = acb_series.__new__(acb_series)
        acb_poly_cos_pi_series((<acb_series>u).val, (<acb_series>s).val, cap, getprec())
        (<acb_series>u).prec = cap
        return u

    def sin_cos_pi(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<acb_series>s).prec)
        u = acb_series.__new__(acb_series)
        v = acb_series.__new__(acb_series)
        acb_poly_sin_cos_pi_series((<acb_series>u).val, (<acb_series>v).val, (<acb_series>s).val, cap, getprec())
        (<acb_series>u).prec = cap
        (<acb_series>v).prec = cap
        return u, v

    def cot_pi(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<acb_series>s).prec)
        u = acb_series.__new__(acb_series)
        acb_poly_cot_pi_series((<acb_series>u).val, (<acb_series>s).val, cap, getprec())
        (<acb_series>u).prec = cap
        return u

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

    def rising(s, ulong n):
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

    def dirichlet_l(s, chi, bint deflate=0):
        cdef long cap
        cdef dirichlet_char cchar
        if isinstance(chi, dirichlet_char):
            cchar = chi
        else:
            cchar = dirichlet_char(chi[0], chi[1])
        cap = getcap()
        cap = min(cap, (<acb_series>s).prec)
        u = acb_series.__new__(acb_series)
        acb_dirichlet_l_series((<acb_series>u).val, (<acb_series>s).val, cchar.G.val, cchar.val, deflate, cap, getprec())
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
        acb_hypgeom_erf_series((<acb_series>u).val, (<acb_series>s).val, cap, getprec())
        (<acb_series>u).prec = cap
        return u

    def erfc(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<acb_series>s).prec)
        u = acb_series.__new__(acb_series)
        acb_hypgeom_erfc_series((<acb_series>u).val, (<acb_series>s).val, cap, getprec())
        (<acb_series>u).prec = cap
        return u

    def erfi(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<acb_series>s).prec)
        u = acb_series.__new__(acb_series)
        acb_hypgeom_erfi_series((<acb_series>u).val, (<acb_series>s).val, cap, getprec())
        (<acb_series>u).prec = cap
        return u

    @classmethod
    def gamma_upper(cls, s, z, int regularized=0):
        cdef long cap
        s = acb(s)
        z = acb_series(z)
        cap = getcap()
        cap = min(cap, (<acb_series>z).prec)
        u = acb_series.__new__(acb_series)
        acb_hypgeom_gamma_upper_series((<acb_series>u).val, (<acb>s).val, (<acb_series>z).val, regularized, cap, getprec())
        (<acb_series>u).prec = cap
        return u

    @classmethod
    def gamma_lower(cls, s, z, int regularized=0):
        cdef long cap
        s = acb(s)
        z = acb_series(z)
        cap = getcap()
        cap = min(cap, (<acb_series>z).prec)
        u = acb_series.__new__(acb_series)
        acb_hypgeom_gamma_lower_series((<acb_series>u).val, (<acb>s).val, (<acb_series>z).val, regularized, cap, getprec())
        (<acb_series>u).prec = cap
        return u

    @classmethod
    def beta_lower(cls, a, b, z, int regularized=0):
        cdef long cap
        a = acb(a)
        b = acb(b)
        z = acb_series(z)
        cap = getcap()
        cap = min(cap, (<acb_series>z).prec)
        u = acb_series.__new__(acb_series)
        acb_hypgeom_beta_lower_series((<acb_series>u).val, (<acb>a).val, (<acb>b).val, (<acb_series>z).val, regularized, cap, getprec())
        (<acb_series>u).prec = cap
        return u

    @classmethod
    def hypgeom(cls, a, b, z, long n=-1, bint regularized=False):
        r"""
        Computes the generalized hypergeometric function `{}_pF_q(a;b;z)`
        given lists of power series `a` and `b` and a power series `z`.

        The optional parameter *n*, if nonnegative, controls the number
        of terms to add in the hypergeometric series. This is just a tuning
        parameter: a rigorous error bound is computed regardless of *n*.
        """
        cdef long i, p, q, prec, cap
        cdef acb_poly_struct * aa
        cdef acb_poly_struct * bb
        a = [acb_series(t) for t in a]
        b = [acb_series(t) for t in b] + [acb_series(1)]  # todo: remove from a if there
        z = acb_series(z)
        p = len(a)
        q = len(b)
        aa = <acb_poly_struct *>libc.stdlib.malloc(p * cython.sizeof(acb_poly_struct))
        bb = <acb_poly_struct *>libc.stdlib.malloc(q * cython.sizeof(acb_poly_struct))
        cap = getcap()
        cap = min(cap, (<acb_series>z).prec)
        for i in range(p):
            cap = min(cap, (<acb_series>(a[i])).prec)
            aa[i] = (<acb_series>(a[i])).val[0]
        for i in range(q):
            cap = min(cap, (<acb_series>(b[i])).prec)
            bb[i] = (<acb_series>(b[i])).val[0]
        u = acb_series.__new__(acb_series)
        acb_hypgeom_pfq_series_direct((<acb_series>u).val, aa, p, bb, q, (<acb_series>z).val, regularized, n, cap, getprec())
        libc.stdlib.free(aa)
        libc.stdlib.free(bb)
        (<acb_series>u).prec = cap
        return u    

    def airy_ai(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<acb_series>s).prec)
        u = acb_series.__new__(acb_series)
        acb_hypgeom_airy_series((<acb_series>u).val, NULL, NULL, NULL, (<acb_series>s).val, cap, getprec())
        (<acb_series>u).prec = cap
        return u

    def airy_ai_prime(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<acb_series>s).prec)
        u = acb_series.__new__(acb_series)
        acb_hypgeom_airy_series(NULL, (<acb_series>u).val, NULL, NULL, (<acb_series>s).val, cap, getprec())
        (<acb_series>u).prec = cap
        return u

    def airy_bi(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<acb_series>s).prec)
        u = acb_series.__new__(acb_series)
        acb_hypgeom_airy_series(NULL, NULL, (<acb_series>u).val, NULL, (<acb_series>s).val, cap, getprec())
        (<acb_series>u).prec = cap
        return u

    def airy_bi_prime(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<acb_series>s).prec)
        u = acb_series.__new__(acb_series)
        acb_hypgeom_airy_series(NULL, NULL, NULL, (<acb_series>u).val,(<acb_series>s).val, cap, getprec())
        (<acb_series>u).prec = cap
        return u

    def airy(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<acb_series>s).prec)
        u = acb_series.__new__(acb_series)
        v = acb_series.__new__(acb_series)
        w = acb_series.__new__(acb_series)
        z = acb_series.__new__(acb_series)
        acb_hypgeom_airy_series((<acb_series>u).val, (<acb_series>v).val, (<acb_series>w).val, (<acb_series>z).val, (<acb_series>s).val, cap, getprec())
        (<acb_series>u).prec = cap
        (<acb_series>v).prec = cap
        (<acb_series>w).prec = cap
        (<acb_series>z).prec = cap
        return u, v, w, z

    def modular_theta(self, tau):
        cdef long cap
        tau = acb(tau)
        cap = getcap()
        cap = min(cap, (<acb_series>self).prec)
        t1 = acb_series.__new__(acb_series)
        t2 = acb_series.__new__(acb_series)
        t3 = acb_series.__new__(acb_series)
        t4 = acb_series.__new__(acb_series)
        acb_modular_theta_series((<acb_series>t1).val, (<acb_series>t2).val, (<acb_series>t3).val, (<acb_series>t4).val, (<acb_series>self).val, (<acb>tau).val, cap, getprec())
        (<acb_series>t1).prec = cap
        (<acb_series>t2).prec = cap
        (<acb_series>t3).prec = cap
        (<acb_series>t4).prec = cap
        return t1, t2, t3, t4

    def coulomb(self, l, eta):
        cdef long cap
        l = acb(l)
        eta = acb(eta)
        cap = getcap()
        cap = min(cap, (<acb_series>self).prec)
        F = acb_series.__new__(acb_series)
        G = acb_series.__new__(acb_series)
        Hpos = acb_series.__new__(acb_series)
        Hneg = acb_series.__new__(acb_series)
        acb_hypgeom_coulomb_series((<acb_series>F).val, (<acb_series>G).val, (<acb_series>Hpos).val, (<acb_series>Hneg).val,
            (<acb>l).val, (<acb>eta).val, (<acb_series>self).val, cap, getprec())
        (<acb_series>F).prec = cap
        (<acb_series>G).prec = cap
        (<acb_series>Hpos).prec = cap
        (<acb_series>Hneg).prec = cap
        return F, G, Hpos, Hneg

    def coulomb_f(self, l, eta):
        cdef long cap
        l = acb(l)
        eta = acb(eta)
        cap = getcap()
        cap = min(cap, (<acb_series>self).prec)
        F = acb_series.__new__(acb_series)
        acb_hypgeom_coulomb_series((<acb_series>F).val, NULL, NULL, NULL,
            (<acb>l).val, (<acb>eta).val, (<acb_series>self).val, cap, getprec())
        (<acb_series>F).prec = cap
        return F

    def coulomb_g(self, l, eta):
        cdef long cap
        l = acb(l)
        eta = acb(eta)
        cap = getcap()
        cap = min(cap, (<acb_series>self).prec)
        G = acb_series.__new__(acb_series)
        acb_hypgeom_coulomb_series(NULL, (<acb_series>G).val, NULL, NULL,
            (<acb>l).val, (<acb>eta).val, (<acb_series>self).val, cap, getprec())
        (<acb_series>G).prec = cap
        return G

    def fresnel(s, bint normalized=True):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<acb_series>s).prec)
        u = acb_series.__new__(acb_series)
        v = acb_series.__new__(acb_series)
        acb_hypgeom_fresnel_series((<acb_series>u).val, (<acb_series>v).val, (<acb_series>s).val, normalized, cap, getprec())
        (<acb_series>u).prec = cap
        (<acb_series>v).prec = cap
        return u, v

    def fresnel_s(s, bint normalized=True):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<acb_series>s).prec)
        u = acb_series.__new__(acb_series)
        acb_hypgeom_fresnel_series((<acb_series>u).val, NULL, (<acb_series>s).val, normalized, cap, getprec())
        (<acb_series>u).prec = cap
        return u

    def fresnel_c(s, bint normalized=True):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<acb_series>s).prec)
        u = acb_series.__new__(acb_series)
        acb_hypgeom_fresnel_series(NULL, (<acb_series>u).val, (<acb_series>s).val, normalized, cap, getprec())
        (<acb_series>u).prec = cap
        return u

    def ei(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<acb_series>s).prec)
        u = acb_series.__new__(acb_series)
        acb_hypgeom_ei_series((<acb_series>u).val, (<acb_series>s).val, cap, getprec())
        (<acb_series>u).prec = cap
        return u

    def si(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<acb_series>s).prec)
        u = acb_series.__new__(acb_series)
        acb_hypgeom_si_series((<acb_series>u).val, (<acb_series>s).val, cap, getprec())
        (<acb_series>u).prec = cap
        return u

    def ci(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<acb_series>s).prec)
        u = acb_series.__new__(acb_series)
        acb_hypgeom_ci_series((<acb_series>u).val, (<acb_series>s).val, cap, getprec())
        (<acb_series>u).prec = cap
        return u

    def shi(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<acb_series>s).prec)
        u = acb_series.__new__(acb_series)
        acb_hypgeom_shi_series((<acb_series>u).val, (<acb_series>s).val, cap, getprec())
        (<acb_series>u).prec = cap
        return u

    def chi(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<acb_series>s).prec)
        u = acb_series.__new__(acb_series)
        acb_hypgeom_chi_series((<acb_series>u).val, (<acb_series>s).val, cap, getprec())
        (<acb_series>u).prec = cap
        return u

    def li(s, bint offset=False):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<acb_series>s).prec)
        u = acb_series.__new__(acb_series)
        acb_hypgeom_li_series((<acb_series>u).val, (<acb_series>s).val, offset, cap, getprec())
        (<acb_series>u).prec = cap
        return u

    def lambertw(s, branch=0):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<acb_series>s).prec)
        k = any_as_fmpz(branch)
        u = acb_series.__new__(acb_series)
        acb_poly_lambertw_series((<acb_series>u).val, (<acb_series>s).val, (<fmpz>k).val, 0, cap, getprec())
        (<acb_series>u).prec = cap
        return u

