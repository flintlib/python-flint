cdef fmpq_series_coerce_operands(x, y):
    if typecheck(x, fmpq_series):
        if isinstance(y, (int, long, fmpz, fmpz_poly, fmpz_series, fmpq, fmpq_poly)):
            return x, fmpq_series(y)
        #if isinstance(y, (nmod, nmod_poly, nmod_series)):
        #    return nmod_series(x), nmod_series(y)
        if isinstance(y, (float, arb, arb_poly, arb_series)):
            return arb_series(x), arb_series(y)
        if isinstance(y, (complex, acb, acb_poly, acb_series)):
            return acb_series(x), acb_series(y)
    else:
        if isinstance(x,(int, long, fmpz, fmpz_poly, fmpz_series, fmpq, fmpq_poly)):
            return fmpq_series(x), y
        #if isinstance(x, (nmod, nmod_poly, nmod_series)):
        #    return nmod_series(x), nmod_series(y)
        if isinstance(x, (float, arb, arb_poly, arb_series)):
            return arb_series(x), arb_series(y)
        if isinstance(x, (complex, acb, acb_poly, acb_series)):
            return acb_series(x), acb_series(y)
    return NotImplemented, NotImplemented

cdef class fmpq_series(flint_series):

    cdef fmpq_poly_t val
    cdef long prec

    def __cinit__(self):
        fmpq_poly_init(self.val)
        self.prec = 0

    def __dealloc__(self):
        fmpq_poly_clear(self.val)

    def __init__(self, val=None, den=None, prec=None):
        if prec is None:
            self.prec = getcap()
        else:
            self.prec = prec
        if self.prec < 0:
            self.prec = -1
        if val is not None:
            if typecheck(val, fmpq_series):
                fmpq_poly_set(self.val, (<fmpq_series>val).val)
                self.prec = min((<fmpz_series>val).prec, getcap())
            elif typecheck(val, fmpz_series):
                fmpq_poly_set_fmpz_poly(self.val, (<fmpz_series>val).val)
                self.prec = min((<fmpz_series>val).prec, getcap())
            elif typecheck(val, fmpz_poly):
                fmpq_poly_set_fmpz_poly(self.val, (<fmpz_poly>val).val)
            elif typecheck(val, fmpq_poly):
                fmpq_poly_set(self.val, (<fmpq_poly>val).val)
            elif typecheck(val, list):
                fmpq_poly_set_list(self.val, val)
            else:
                fmpq_poly_set_list(self.val, [val])
        fmpq_poly_truncate(self.val, max(0, self.prec))
        if den is not None:
            den = any_as_fmpz(den)
            if den is NotImplemented:
                raise TypeError("denominator must be an integer, got %s", type(den))
            if fmpz_is_zero((<fmpz>den).val):
                raise ZeroDivisionError("cannot create fmpq_series with zero denominator")
            fmpq_poly_scalar_div_fmpz(self.val, self.val, (<fmpz>den).val)

    def __len__(self):
        return fmpq_poly_length(self.val)

    cpdef long length(self):
        return fmpq_poly_length(self.val)

    def numer(self):
        cdef fmpz_series x = fmpz_series.__new__(fmpz_series)
        fmpq_poly_get_numerator(x.val, self.val)
        x.prec = self.prec
        return x

    def denom(self):
        cdef fmpz x = fmpz.__new__(fmpz)
        fmpz_set(x.val, fmpq_poly_denref(self.val))
        return x

    p = property(numer)
    q = property(denom)

    def __getitem__(self, long i):
        cdef fmpq x
        x = fmpq()
        if i < 0:
            return x
        fmpq_poly_get_coeff_fmpq(x.val, self.val, i)
        return x

    def __setitem__(self, long i, x):
        if i < 0:
            raise ValueError("cannot assign to index < 0 of polynomial")
        if not typecheck(x, fmpq):
            x = fmpq(x)
        fmpq_poly_set_coeff_fmpq(self.val, i, (<fmpq>x).val)

    def repr(self, **kwargs):
        return "fmpq_series([%s], %s, prec=%s)" % (", ".join(map(str, self.numer())), str(self.denom()), self.prec)

    def str(self, **kwargs):
        if self.prec > 0:
            s = fmpq_poly(list(self)).str(ascending=True)
            return s + (" + O(x^%s)" % self.prec)
        elif self.prec == 0:
            return "O(x^0)"
        else:
            return "(invalid power series)"

    def __pos__(self):
        return self

    def __neg__(s):
        cdef long cap
        u = fmpq_series.__new__(fmpq_series)
        cap = getcap()
        cap = min(cap, (<fmpq_series>s).prec)
        if cap > 0:
            fmpq_poly_neg((<fmpq_series>u).val, (<fmpq_series>s).val)
            fmpq_poly_truncate((<fmpq_series>u).val, cap)
        (<fmpq_series>u).prec = cap
        return u

    def __add__(s, t):
        cdef long cap
        if type(s) is type(t):
            u = fmpq_series.__new__(fmpq_series)
            cap = getcap()
            cap = min(cap, (<fmpq_series>s).prec)
            cap = min(cap, (<fmpq_series>t).prec)
            if cap > 0:
                fmpq_poly_add((<fmpq_series>u).val, (<fmpq_series>s).val, (<fmpq_series>t).val)
                fmpq_poly_truncate((<fmpq_series>u).val, cap)
            (<fmpq_series>u).prec = cap
            return u
        s, t = fmpq_series_coerce_operands(s, t)
        if s is NotImplemented:
            return s
        return s + t

    def __sub__(s, t):
        cdef long cap
        if type(s) is type(t):
            u = fmpq_series.__new__(fmpq_series)
            cap = getcap()
            cap = min(cap, (<fmpq_series>s).prec)
            cap = min(cap, (<fmpq_series>t).prec)
            if cap > 0:
                fmpq_poly_sub((<fmpq_series>u).val, (<fmpq_series>s).val, (<fmpq_series>t).val)
                fmpq_poly_truncate((<fmpq_series>u).val, cap)
            (<fmpq_series>u).prec = cap
            return u
        s, t = fmpq_series_coerce_operands(s, t)
        if s is NotImplemented:
            return s
        return s - t

    def __mul__(s, t):
        cdef long cap
        if type(s) is type(t):
            u = fmpq_series.__new__(fmpq_series)
            cap = getcap()
            cap = min(cap, (<fmpq_series>s).prec)
            cap = min(cap, (<fmpq_series>t).prec)
            if cap > 0:
                fmpq_poly_mullow((<fmpq_series>u).val, (<fmpq_series>s).val, (<fmpq_series>t).val, cap)
            (<fmpq_series>u).prec = cap
            return u
        s, t = fmpq_series_coerce_operands(s, t)
        if s is NotImplemented:
            return s
        return s * t

    cpdef valuation(self):
        cdef long i
        if fmpq_poly_is_zero(self.val):
            return -1
        i = 0
        while fmpz_is_zero(&(self.val.coeffs[i])):
            i += 1
        return i

    @staticmethod
    def _div_(s, t):
        cdef long cap, sval, tval
        cdef fmpq_poly_t stmp, ttmp
        if type(s) is type(t):
            cap = getcap()
            cap = min(cap, (<fmpq_series>s).prec)
            cap = min(cap, (<fmpq_series>t).prec)

            if fmpq_poly_is_zero((<fmpq_series>t).val):
                raise ZeroDivisionError("power series division")

            u = fmpq_series.__new__(fmpq_series)

            if fmpq_poly_is_zero((<fmpq_series>s).val):
                u.cap = cap
                return u

            sval = (<fmpq_series>s).valuation()
            tval = (<fmpq_series>t).valuation()

            if sval < tval:
                raise ValueError("quotient would not be a power series")

            if fmpz_is_zero(&((<fmpq_series>t).val.coeffs[tval])):
                raise ValueError("leading term in denominator is not a unit")

            if tval == 0:
                fmpq_poly_div_series((<fmpq_series>u).val, (<fmpq_series>s).val, (<fmpq_series>t).val, cap)
            else:
                fmpq_poly_init(stmp)
                fmpq_poly_init(ttmp)
                fmpq_poly_shift_right(stmp, (<fmpq_series>s).val, tval)
                fmpq_poly_shift_right(ttmp, (<fmpq_series>t).val, tval)
                cap -= tval
                fmpq_poly_div_series((<fmpq_series>u).val, stmp, ttmp, cap)
                fmpq_poly_clear(stmp)
                fmpq_poly_clear(ttmp)

            (<fmpq_series>u).prec = cap
            return u

        s, t = fmpq_series_coerce_operands(s, t)
        if s is NotImplemented:
            return s
        return s / t

    def __truediv__(s, t):
        return fmpq_series._div_(s, t)

    def __div__(s, t):
        return fmpq_series._div_(s, t)

    # generic exponentiation (fallback code)
    def __pow__(s, ulong exp, mod):
        if mod is not None:
            raise NotImplementedError("modular exponentiation")
        cdef long i
        if exp == 0:
            return (s * 0) + 1
        if exp == 1:
            return s
        if exp == 2:
            return s * s
        y = s
        bits = exp.bit_length()
        for i in range(bits-2, -1, -1):
            y = y * y
            if (exp >> i) & 1:
                y = y * s
        return y

    def __call__(s, t):
        cdef long cap
        if typecheck(t, fmpq_series):
            u = fmpq_series.__new__(fmpq_series)
            if (<fmpq_series>t).valuation() < 1:
                raise ValueError("power series composition with nonzero constant term")
            cap = getcap()
            cap = min(cap, (<fmpq_series>s).prec)
            cap = min(cap, (<fmpq_series>t).prec)
            fmpq_poly_compose_series((<fmpq_series>u).val, (<fmpq_series>s).val, (<fmpq_series>t).val, cap)
            (<fmpq_series>u).prec = cap
            return u
        raise TypeError("cannot call fmpq_series with input of type %s", type(t))

    def reversion(s):
        """
        Returns the power series reversion (compositional inverse) of *s*.

            >>> x = fmpq_series([0,1]); print((x/2-x**2).reversion())
            2*x + 8*x^2 + 64*x^3 + 640*x^4 + 7168*x^5 + 86016*x^6 + 1081344*x^7 + 14057472*x^8 + 187432960*x^9 + O(x^10)
        """
        cdef long cap
        cap = getcap()
        cap = min(cap, (<fmpq_series>s).prec)
        if (<fmpq_series>s).valuation() != 1:
            raise ValueError("power series reversion must have valuation 1")
        if fmpz_is_zero(&((<fmpq_series>s).val.coeffs[1])):
            raise ValueError("leading term is not a unit")
        u = fmpq_series.__new__(fmpq_series)
        fmpq_poly_revert_series((<fmpq_series>u).val, (<fmpq_series>s).val, cap)
        (<fmpq_series>u).prec = cap
        return u

    def inv(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<fmpq_series>s).prec)
        if (<fmpq_series>s).valuation() != 0:
            raise ValueError("can only invert series with nonzero constant term")
        u = fmpq_series.__new__(fmpq_series)
        fmpq_poly_inv_series((<fmpq_series>u).val, (<fmpq_series>s).val, cap)
        (<fmpq_series>u).prec = cap
        return u

    cdef bint zero_constant_term(s):
        if fmpq_poly_is_zero((<fmpq_series>s).val):
            return True
        if fmpz_is_zero(&((<fmpq_series>s).val.coeffs[0])):
            return True
        return False

    cdef bint one_constant_term(s):
        if fmpq_poly_is_zero((<fmpq_series>s).val):
            return False
        if fmpz_is_one(&((<fmpq_series>s).val.coeffs[0])):
            return True
        return False

    def derivative(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<fmpq_series>s).prec - 1)
        u = fmpq_series.__new__(fmpq_series)
        fmpq_poly_derivative((<fmpq_series>u).val, (<fmpq_series>s).val)
        fmpq_poly_truncate((<fmpq_series>u).val, max(0, cap))
        (<fmpq_series>u).prec = cap
        return u

    def integral(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<fmpq_series>s).prec + 1)
        u = fmpq_series.__new__(fmpq_series)
        fmpq_poly_integral((<fmpq_series>u).val, (<fmpq_series>s).val)
        fmpq_poly_truncate((<fmpq_series>u).val, max(0, cap))
        (<fmpq_series>u).prec = cap
        return u

    def sqrt(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<fmpq_series>s).prec)
        if not s.one_constant_term():
            raise ValueError("sqrt() of power series: constant term != 1")
        u = fmpq_series.__new__(fmpq_series)
        fmpq_poly_sqrt_series((<fmpq_series>u).val, (<fmpq_series>s).val, cap)
        (<fmpq_series>u).prec = cap
        return u

    def rsqrt(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<fmpq_series>s).prec)
        if not s.one_constant_term():
            raise ValueError("rsqrt() of power series: constant term != 1")
        u = fmpq_series.__new__(fmpq_series)
        fmpq_poly_invsqrt_series((<fmpq_series>u).val, (<fmpq_series>s).val, cap)
        (<fmpq_series>u).prec = cap
        return u

    def exp(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<fmpq_series>s).prec)
        if not s.zero_constant_term():
            raise ValueError("exp() of power series: constant term must be zero")
        u = fmpq_series.__new__(fmpq_series)
        fmpq_poly_exp_series((<fmpq_series>u).val, (<fmpq_series>s).val, cap)
        (<fmpq_series>u).prec = cap
        return u

    def log(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<fmpq_series>s).prec)
        if not s.one_constant_term():
            raise ValueError("log() of power series: constant term must be one")
        u = fmpq_series.__new__(fmpq_series)
        fmpq_poly_log_series((<fmpq_series>u).val, (<fmpq_series>s).val, cap)
        (<fmpq_series>u).prec = cap
        return u

    def atan(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<fmpq_series>s).prec)
        if not s.zero_constant_term():
            raise ValueError("atan() of power series: constant term must be zero")
        u = fmpq_series.__new__(fmpq_series)
        fmpq_poly_atan_series((<fmpq_series>u).val, (<fmpq_series>s).val, cap)
        (<fmpq_series>u).prec = cap
        return u

    def atanh(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<fmpq_series>s).prec)
        if not s.zero_constant_term():
            raise ValueError("atanh() of power series: constant term must be zero")
        u = fmpq_series.__new__(fmpq_series)
        fmpq_poly_atanh_series((<fmpq_series>u).val, (<fmpq_series>s).val, cap)
        (<fmpq_series>u).prec = cap
        return u

    def asin(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<fmpq_series>s).prec)
        if not s.zero_constant_term():
            raise ValueError("asin() of power series: constant term must be zero")
        u = fmpq_series.__new__(fmpq_series)
        fmpq_poly_asin_series((<fmpq_series>u).val, (<fmpq_series>s).val, cap)
        (<fmpq_series>u).prec = cap
        return u

    def asinh(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<fmpq_series>s).prec)
        if not s.zero_constant_term():
            raise ValueError("asinh() of power series: constant term must be zero")
        u = fmpq_series.__new__(fmpq_series)
        fmpq_poly_asinh_series((<fmpq_series>u).val, (<fmpq_series>s).val, cap)
        (<fmpq_series>u).prec = cap
        return u

    def sin(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<fmpq_series>s).prec)
        if not s.zero_constant_term():
            raise ValueError("sin() of power series: constant term must be zero")
        u = fmpq_series.__new__(fmpq_series)
        fmpq_poly_sin_series((<fmpq_series>u).val, (<fmpq_series>s).val, cap)
        (<fmpq_series>u).prec = cap
        return u

    def cos(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<fmpq_series>s).prec)
        if not s.zero_constant_term():
            raise ValueError("cos() of power series: constant term must be zero")
        u = fmpq_series.__new__(fmpq_series)
        fmpq_poly_cos_series((<fmpq_series>u).val, (<fmpq_series>s).val, cap)
        (<fmpq_series>u).prec = cap
        return u

    def tan(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<fmpq_series>s).prec)
        if not s.zero_constant_term():
            raise ValueError("tan() of power series: constant term must be zero")
        u = fmpq_series.__new__(fmpq_series)
        fmpq_poly_tan_series((<fmpq_series>u).val, (<fmpq_series>s).val, cap)
        (<fmpq_series>u).prec = cap
        return u

    def sinh(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<fmpq_series>s).prec)
        if not s.zero_constant_term():
            raise ValueError("sinh() of power series: constant term must be zero")
        u = fmpq_series.__new__(fmpq_series)
        fmpq_poly_sinh_series((<fmpq_series>u).val, (<fmpq_series>s).val, cap)
        (<fmpq_series>u).prec = cap
        return u

    def cosh(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<fmpq_series>s).prec)
        if not s.zero_constant_term():
            raise ValueError("cosh() of power series: constant term must be zero")
        u = fmpq_series.__new__(fmpq_series)
        fmpq_poly_cosh_series((<fmpq_series>u).val, (<fmpq_series>s).val, cap)
        (<fmpq_series>u).prec = cap
        return u

    def tanh(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<fmpq_series>s).prec)
        if not s.zero_constant_term():
            raise ValueError("tanh() of power series: constant term must be zero")
        u = fmpq_series.__new__(fmpq_series)
        fmpq_poly_tanh_series((<fmpq_series>u).val, (<fmpq_series>s).val, cap)
        (<fmpq_series>u).prec = cap
        return u

