cdef fmpz_series_coerce_operands(x, y):
    if typecheck(x, fmpz_series):
        if isinstance(y, (int, long, fmpz, fmpz_poly)):
            return x, fmpz_series(y)
        if isinstance(y, (fmpq, fmpq_poly, fmpq_series)):
            return fmpq_series(x), fmpq_series(y)
        #if isinstance(y, (nmod, nmod_poly, nmod_series)):
        #    return nmod_series(x), nmod_series(y)
        if isinstance(y, (float, arb, arb_poly, arb_series)):
            return arb_series(x), arb_series(y)
        if isinstance(y, (complex, acb, acb_poly, acb_series)):
            return acb_series(x), acb_series(y)
    else:
        if isinstance(x, (int, long, fmpz, fmpz_poly)):
            return fmpz_series(x), y
        if isinstance(x, (fmpq, fmpq_poly, fmpq_series)):
            return fmpq_series(x), fmpq_series(y)
        #if isinstance(x, (nmod, nmod_poly, nmod_series)):
        #    return nmod_series(x), nmod_series(y)
        if isinstance(x, (float, arb, arb_poly, arb_series)):
            return arb_series(x), arb_series(y)
        if isinstance(x, (complex, acb, acb_poly, acb_series)):
            return acb_series(x), acb_series(y)
    return NotImplemented, NotImplemented

cdef class fmpz_series(flint_series):

    cdef fmpz_poly_t val
    cdef long prec

    def __cinit__(self):
        fmpz_poly_init(self.val)
        self.prec = 0

    def __dealloc__(self):
        fmpz_poly_clear(self.val)

    def __init__(self, val=None, prec=None):
        if prec is None:
            self.prec = getcap()
        else:
            self.prec = prec
        if self.prec < 0:
            self.prec = -1
        if val is not None:
            if typecheck(val, fmpz_series):
                fmpz_poly_set(self.val, (<fmpz_series>val).val)
                self.prec = min((<fmpz_series>val).prec, getcap())
            elif typecheck(val, fmpz_poly):
                fmpz_poly_set(self.val, (<fmpz_poly>val).val)
            elif typecheck(val, list):
                fmpz_poly_set_list(self.val, val)
            else:
                fmpz_poly_set_list(self.val, [val])
        fmpz_poly_truncate(self.val, max(0, self.prec))

    def __len__(self):
        return fmpz_poly_length(self.val)

    cpdef long length(self):
        return fmpz_poly_length(self.val)

    def __getitem__(self, long i):
        cdef fmpz x
        x = fmpz()
        if i >= 0:
            fmpz_poly_get_coeff_fmpz(x.val, self.val, i)
        return x

    def __setitem__(self, long i, x):
        if i < 0:
            raise ValueError("cannot assign to index < 0 of series")
        if not typecheck(x, fmpz):
            x = fmpz(x)
        fmpz_poly_set_coeff_fmpz(self.val, i, (<fmpz>x).val)

    def repr(self, **kwargs):
        return "fmpz_series([%s], prec=%s)" % (", ".join(map(str, self)), self.prec)

    def str(self, **kwargs):
        if self.prec > 0:
            s = fmpz_poly(list(self)).str(ascending=True)
            return s + (" + O(x^%s)" % self.prec)
        elif self.prec == 0:
            return "O(x^0)"
        else:
            return "(invalid power series)"

    def __pos__(self):
        return self

    def __neg__(s):
        cdef long cap
        u = fmpz_series.__new__(fmpz_series)
        cap = getcap()
        cap = min(cap, (<fmpz_series>s).prec)
        if cap > 0:
            fmpz_poly_neg((<fmpz_series>u).val, (<fmpz_series>s).val)
            fmpz_poly_truncate((<fmpz_series>u).val, cap)
        (<fmpz_series>u).prec = cap
        return u

    def __add__(s, t):
        cdef long cap
        if type(s) is type(t):
            u = fmpz_series.__new__(fmpz_series)
            cap = getcap()
            cap = min(cap, (<fmpz_series>s).prec)
            cap = min(cap, (<fmpz_series>t).prec)
            if cap > 0:
                fmpz_poly_add((<fmpz_series>u).val, (<fmpz_series>s).val, (<fmpz_series>t).val)
                fmpz_poly_truncate((<fmpz_series>u).val, cap)
            (<fmpz_series>u).prec = cap
            return u
        s, t = fmpz_series_coerce_operands(s, t)
        if s is NotImplemented:
            return s
        return s + t

    def __sub__(s, t):
        cdef long cap
        if type(s) is type(t):
            u = fmpz_series.__new__(fmpz_series)
            cap = getcap()
            cap = min(cap, (<fmpz_series>s).prec)
            cap = min(cap, (<fmpz_series>t).prec)
            if cap > 0:
                fmpz_poly_sub((<fmpz_series>u).val, (<fmpz_series>s).val, (<fmpz_series>t).val)
                fmpz_poly_truncate((<fmpz_series>u).val, cap)
            (<fmpz_series>u).prec = cap
            return u
        s, t = fmpz_series_coerce_operands(s, t)
        if s is NotImplemented:
            return s
        return s - t

    def __mul__(s, t):
        cdef long cap
        if type(s) is type(t):
            u = fmpz_series.__new__(fmpz_series)
            cap = getcap()
            cap = min(cap, (<fmpz_series>s).prec)
            cap = min(cap, (<fmpz_series>t).prec)
            if cap > 0:
                fmpz_poly_mullow((<fmpz_series>u).val, (<fmpz_series>s).val, (<fmpz_series>t).val, cap)
            (<fmpz_series>u).prec = cap
            return u
        s, t = fmpz_series_coerce_operands(s, t)
        if s is NotImplemented:
            return s
        return s * t

    cpdef valuation(self):
        cdef long i
        if fmpz_poly_is_zero(self.val):
            return -1
        i = 0
        while fmpz_is_zero(fmpz_poly_get_coeff_ptr(self.val, i)):
            i += 1
        return i

    def __div__(s, t):
        cdef long cap, sval, tval
        cdef fmpz_poly_t stmp, ttmp
        if type(s) is type(t):
            cap = getcap()
            cap = min(cap, (<fmpz_series>s).prec)
            cap = min(cap, (<fmpz_series>t).prec)

            if fmpz_poly_is_zero((<fmpz_series>t).val):
                raise ZeroDivisionError("power series division")

            u = fmpz_series.__new__(fmpz_series)

            if fmpz_poly_is_zero((<fmpz_series>s).val):
                u.cap = cap
                return u

            sval = (<fmpz_series>s).valuation()
            tval = (<fmpz_series>t).valuation()

            if sval < tval:
                raise ValueError("quotient would not be a power series")

            if not fmpz_is_pm1(fmpz_poly_get_coeff_ptr((<fmpz_series>t).val, tval)):
                raise ValueError("leading term in denominator is not a unit")

            if tval == 0:
                fmpz_poly_div_series((<fmpz_series>u).val, (<fmpz_series>s).val, (<fmpz_series>t).val, cap)
            else:
                fmpz_poly_init(stmp)
                fmpz_poly_init(ttmp)
                fmpz_poly_shift_right(stmp, (<fmpz_series>s).val, tval)
                fmpz_poly_shift_right(ttmp, (<fmpz_series>t).val, tval)
                cap -= tval
                fmpz_poly_div_series((<fmpz_series>u).val, stmp, ttmp, cap)
                fmpz_poly_clear(stmp)
                fmpz_poly_clear(ttmp)

            (<fmpz_series>u).prec = cap
            return u

        s, t = fmpz_series_coerce_operands(s, t)
        if s is NotImplemented:
            return s
        return s / t

    def __truediv__(s, t):
        return fmpz_series.__div__(s, t)

    def __pow__(fmpz_series s, ulong exp, mod):
        cdef long cap
        if mod is not None:
            raise NotImplementedError("fmpz_series modular exponentiation")
        cap = getcap()
        cap = min(cap, (<fmpz_series>s).prec)
        u = fmpz_series.__new__(fmpz_series)
        fmpz_poly_pow_trunc((<fmpz_series>u).val, (<fmpz_series>s).val, exp, cap)
        (<fmpz_series>u).prec = cap
        return u

    def __call__(s, t):
        cdef long cap
        if typecheck(t, fmpz_series):
            u = fmpz_series.__new__(fmpz_series)
            if (<fmpz_series>t).valuation() < 1:
                raise ValueError("power series composition with nonzero constant term")
            cap = getcap()
            cap = min(cap, (<fmpz_series>s).prec)
            cap = min(cap, (<fmpz_series>t).prec)
            fmpz_poly_compose_series((<fmpz_series>u).val, (<fmpz_series>s).val, (<fmpz_series>t).val, cap)
            (<fmpz_series>u).prec = cap
            return u
        raise TypeError("cannot call fmpz_series with input of type %s", type(t))

    def reversion(s):
        cdef long cap
        cap = getcap()
        cap = min(cap, (<fmpz_series>s).prec)
        if (<fmpz_series>s).valuation() != 1:
            raise ValueError("power series reversion must have valuation 1")
        if not fmpz_is_pm1(fmpz_poly_get_coeff_ptr((<fmpz_series>s).val, 1)):
            raise ValueError("leading term is not a unit")
        u = fmpz_series.__new__(fmpz_series)
        fmpz_poly_revert_series((<fmpz_series>u).val, (<fmpz_series>s).val, cap)
        (<fmpz_series>u).prec = cap
        return u

