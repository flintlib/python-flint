cdef class arf:

    cdef arf_t val

    def __cinit__(self):
        arf_init(self.val)

    def __dealloc__(self):
        arf_clear(self.val)

    def __init__(self, val=None):
        """
        Create a new arf from an integer, a Python float, an existing arf,
        or a tuple containing (mantissa, exponent)::

            >>> arf(-100)
            -100.000000000000
            >>> arf(15.125)
            15.1250000000000
            >>> arf(arf(3))
            3.00000000000000
            >>> arf((10,-2))
            2.50000000000000

        """
        if val is not None:
            if typecheck(val, fmpz):
                arf_set_fmpz(self.val, (<fmpz>val).val)
            elif typecheck(val, arf):
                arf_set(self.val, (<arf>val).val)
            elif typecheck(val, float):
                arf_set_d(self.val, PyFloat_AS_DOUBLE(<PyObject*>val))
            elif typecheck(val, tuple):
                man = any_as_fmpz(val[0])
                exp = any_as_fmpz(val[1])
                if man is NotImplemented or exp is NotImplemented:
                    raise TypeError("cannot create arf from tuple (%s, %s)" % (type(val[0]), type(val[1])))
                arf_set_fmpz_2exp(self.val, (<fmpz>man).val, (<fmpz>exp).val)
            elif typecheck(val, str):
                if val == "+inf" or val == "inf":
                    arf_pos_inf(self.val)
                elif val == "-inf":
                    arf_neg_inf(self.val)
                elif val == "nan":
                    arf_nan(self.val)
                else:
                    raise TypeError("cannot create arf from type %s" % type(val))
            else:
                v = any_as_fmpz(val)
                if v is not NotImplemented:
                    arf_set_fmpz(self.val, (<fmpz>v).val)
                else:
                    raise TypeError("cannot create arf from type %s" % type(val))

    cpdef bint is_finite(self):
        return arf_is_finite(self.val)

    cpdef bint is_pos_inf(self):
        return arf_is_pos_inf(self.val)

    cpdef bint is_neg_inf(self):
        return arf_is_neg_inf(self.val)

    cpdef bint is_nan(self):
        return arf_is_nan(self.val)

    cpdef bint is_zero(self):
        return arf_is_zero(self.val)

    # todo: exception when not finite!
    def man_exp(self):
        cdef fmpz man, exp
        man = fmpz()
        exp = fmpz()
        arf_get_fmpz_2exp(man.val, exp.val, self.val)
        return man, exp

    def _repr_str(self):
        cdef fmpz man, exp
        if arf_is_zero(self.val):
            return "0.0"
        elif arf_is_finite(self.val):
            return "(%s, %s)" % self.man_exp()
        elif arf_is_pos_inf(self.val):
            return "'+inf'"
        elif arf_is_neg_inf(self.val):
            return "'-inf'"
        else:
            return "'nan'"

    def _dec_str(self, num_digits=None):
        if arf_is_normal(self.val):
            if num_digits is None:
                num_digits = prec_to_dps(getprec())
            return arb(self).str(num_digits, radius=False)
        elif arf_is_zero(self.val):
            return "0.0"
        elif arf_is_pos_inf(self.val):
            return "inf"
        elif arf_is_neg_inf(self.val):
            return "-inf"
        else:
            return "nan"

    def __repr__(self):
        if ctx.pretty:
            return str(self)
        return "arf(%s)" % self._repr_str()

    def __str__(self):
        return self._dec_str()

    def __richcmp__(s, t, int op):
        cdef bint res = 0
        if not typecheck(t, arf):
            t = arf(t)
        if   op == 2: res = arf_equal((<arf>s).val, (<arf>t).val)
        elif op == 3: res = not arf_equal((<arf>s).val, (<arf>t).val)
        elif op == 0: res = arf_cmp((<arf>s).val, (<arf>t).val) < 0
        elif op == 1: res = arf_cmp((<arf>s).val, (<arf>t).val) <= 0
        elif op == 4: res = arf_cmp((<arf>s).val, (<arf>t).val) > 0
        elif op == 5: res = arf_cmp((<arf>s).val, (<arf>t).val) >= 0
        return res

    def __pos__(self):
        res = arf.__new__(arf)
        arf_set_round((<arf>res).val, (<arf>self).val, getprec(), thectx.rnd)
        return res

    def __neg__(self):
        res = arf.__new__(arf)
        arf_neg_round((<arf>res).val, (<arf>self).val, getprec(), thectx.rnd)
        return res

    def __abs__(self):
        res = arf.__new__(arf)
        arf_abs((<arf>res).val, (<arf>self).val)
        arf_set_round((<arf>res).val, (<arf>res).val, getprec(), thectx.rnd)
        return res

    def __add__(s, t):
        if typecheck(s, arf):
            if typecheck(t, arf):
                u = arf.__new__(arf)
                arf_add((<arf>u).val, (<arf>s).val, (<arf>t).val, getprec(), thectx.rnd)
                return u
        return NotImplemented

    def __sub__(s, t):
        if typecheck(s, arf):
            if typecheck(t, arf):
                u = arf.__new__(arf)
                arf_sub((<arf>u).val, (<arf>s).val, (<arf>t).val, getprec(), thectx.rnd)
                return u
        return NotImplemented

    def __mul__(s, t):
        if typecheck(s, arf):
            if typecheck(t, arf):
                u = arf.__new__(arf)
                arf_mul((<arf>u).val, (<arf>s).val, (<arf>t).val, getprec(), thectx.rnd)
                return u
        return NotImplemented

    def __div__(s, t):
        if typecheck(s, arf):
            if typecheck(t, arf):
                u = arf.__new__(arf)
                arf_div((<arf>u).val, (<arf>s).val, (<arf>t).val, getprec(), thectx.rnd)
                return u
        return NotImplemented

