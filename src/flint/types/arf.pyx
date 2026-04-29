from cpython.float cimport PyFloat_AsDouble
from flint.flint_base.flint_context cimport getprec
from flint.flint_base.flint_context cimport thectx
from flint.utils.typecheck cimport typecheck
from flint.utils.conversion cimport prec_to_dps
from flint.types.fmpz cimport fmpz
from flint.types.fmpz cimport any_as_fmpz
from flint.types.fmpq cimport fmpq
from flint.types.arb cimport arb

from flint.flintlib.functions.arf cimport *
from flint.flintlib.functions.fmpq cimport fmpq_numref, fmpq_denref
from flint.flintlib.types.arf cimport ARF_RND_DOWN

ctx = thectx

cdef class arf:

    # cdef arf_t val

    def __cinit__(self):
        arf_init(self.val)

    def __dealloc__(self):
        arf_clear(self.val)

    def __init__(self, val=None):
        """
        Create a new arf from an integer, a Python float, an existing arf,
        or a tuple containing (mantissa, exponent)::

            >>> from flint import arf, ctx
            >>> ctx.prec = 53
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
            elif typecheck(val, fmpq):
                arf_fmpz_div_fmpz(self.val, fmpq_numref((<fmpq>val).val), fmpq_denref((<fmpq>val).val), getprec(), thectx.rnd)
            elif typecheck(val, float):
                arf_set_d(self.val, PyFloat_AsDouble(val))
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
        if arf_is_zero(self.val):
            return "0.0"
        elif arf_is_finite(self.val):
            man, exp = self.man_exp()
            return "(%s, %s)" % (hex(man), hex(exp))
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

    def __bool__(self):
        return not arf_is_zero(self.val)

    def __float__(self):
        return arf_get_d(self.val, thectx.rnd)

    def __int__(self):
        cdef fmpz z
        if arf_is_nan(self.val):
            raise ValueError("cannot convert NaN to integer")
        elif arf_is_inf(self.val):
            raise OverflowError("cannot convert infinity to integer")
        z = fmpz.__new__(fmpz)
        arf_get_fmpz(z.val, self.val, ARF_RND_DOWN)
        return int(z)

    def as_integer_ratio(self):
        if arf_is_nan(self.val):
            raise ValueError("cannot convert NaN to integer ratio")
        elif arf_is_inf(self.val):
            raise OverflowError("cannot convert infinity to integer ratio")
        man, exp = self.man_exp()
        if exp >= 0:
            return int(man << exp), 1
        return int(man), int(fmpz(1) << (-exp))

    def __richcmp__(s, t, int op):
        cdef bint res = 0
        cdef int cmp
        cdef object z
        cdef arf_t v, num_arf

        if typecheck(t, arf):
            if op == 2:
                return bool(arf_equal((<arf>s).val, (<arf>t).val))
            if op == 3:
                return not arf_equal((<arf>s).val, (<arf>t).val)
            cmp = arf_cmp((<arf>s).val, (<arf>t).val)
        elif typecheck(t, float):
            if op == 2:
                return bool(arf_equal_d((<arf>s).val, PyFloat_AsDouble(t)))
            if op == 3:
                return not arf_equal_d((<arf>s).val, PyFloat_AsDouble(t))
            cmp = arf_cmp_d((<arf>s).val, PyFloat_AsDouble(t))
        elif typecheck(t, fmpq):
            if arf_is_nan((<arf>s).val):
                if op == 3:
                    return True
                return False
            if arf_is_pos_inf((<arf>s).val):
                cmp = 1
            elif arf_is_neg_inf((<arf>s).val):
                cmp = -1
            else:
                # Compare exactly by cross-multiplying with positive denominator:
                # s ? (n/d)  <=>  s*d ? n
                arf_init(v)
                arf_init(num_arf)
                arf_mul_fmpz(v, (<arf>s).val, fmpq_denref((<fmpq>t).val), ARF_PREC_EXACT, ARF_RND_DOWN)
                arf_set_fmpz(num_arf, fmpq_numref((<fmpq>t).val))
                cmp = arf_cmp(v, num_arf)
                arf_clear(v)
                arf_clear(num_arf)
            if op == 2:
                return cmp == 0
            if op == 3:
                return cmp != 0
        else:
            z = any_as_fmpz(t)
            if z is NotImplemented:
                return NotImplemented
            arf_init(v)
            arf_set_fmpz(v, (<fmpz>z).val)
            if op == 2:
                res = arf_equal((<arf>s).val, v)
                arf_clear(v)
                return res
            if op == 3:
                res = not arf_equal((<arf>s).val, v)
                arf_clear(v)
                return res
            cmp = arf_cmp((<arf>s).val, v)
            arf_clear(v)

        if op == 0:
            res = cmp < 0
        elif op == 1:
            res = cmp <= 0
        elif op == 4:
            res = cmp > 0
        elif op == 5:
            res = cmp >= 0
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
        cdef slong t_si
        cdef arf_t v
        cdef arf u
        if typecheck(t, arf):
            u = arf.__new__(arf)
            arf_add((<arf>u).val, (<arf>s).val, (<arf>t).val, getprec(), thectx.rnd)
            return u
        elif typecheck(t, int):
            try:
                t_si = <slong>t
            except OverflowError:
                pass
            else:
                u = arf.__new__(arf)
                arf_add_si((<arf>u).val, (<arf>s).val, t_si, getprec(), thectx.rnd)
                return u
        elif typecheck(t, float):
            u = arf.__new__(arf)
            arf_init(v)
            arf_set_d(v, PyFloat_AsDouble(t))
            arf_add((<arf>u).val, (<arf>s).val, v, getprec(), thectx.rnd)
            arf_clear(v)
            return u
        z = any_as_fmpz(t)
        if z is not NotImplemented:
            u = arf.__new__(arf)
            arf_add_fmpz((<arf>u).val, (<arf>s).val, (<fmpz>z).val, getprec(), thectx.rnd)
            return u
        return NotImplemented

    def __radd__(s, t):
        cdef slong t_si
        cdef arf_t v
        cdef arf u
        if typecheck(t, arf):
            u = arf.__new__(arf)
            arf_add((<arf>u).val, (<arf>t).val, (<arf>s).val, getprec(), thectx.rnd)
            return u
        elif typecheck(t, int):
            try:
                t_si = <slong>t
            except OverflowError:
                pass
            else:
                u = arf.__new__(arf)
                arf_add_si((<arf>u).val, (<arf>s).val, t_si, getprec(), thectx.rnd)
                return u
        elif typecheck(t, float):
            u = arf.__new__(arf)
            arf_init(v)
            arf_set_d(v, PyFloat_AsDouble(t))
            arf_add((<arf>u).val, v, (<arf>s).val, getprec(), thectx.rnd)
            arf_clear(v)
            return u
        z = any_as_fmpz(t)
        if z is not NotImplemented:
            u = arf.__new__(arf)
            arf_add_fmpz((<arf>u).val, (<arf>s).val, (<fmpz>z).val, getprec(), thectx.rnd)
            return u
        return NotImplemented

    def __sub__(s, t):
        cdef slong t_si
        cdef arf_t v
        cdef arf u
        if typecheck(t, arf):
            u = arf.__new__(arf)
            arf_sub((<arf>u).val, (<arf>s).val, (<arf>t).val, getprec(), thectx.rnd)
            return u
        elif typecheck(t, int):
            try:
                t_si = <slong>t
            except OverflowError:
                pass
            else:
                u = arf.__new__(arf)
                arf_sub_si((<arf>u).val, (<arf>s).val, t_si, getprec(), thectx.rnd)
                return u
        elif typecheck(t, float):
            u = arf.__new__(arf)
            arf_init(v)
            arf_set_d(v, PyFloat_AsDouble(t))
            arf_sub((<arf>u).val, (<arf>s).val, v, getprec(), thectx.rnd)
            arf_clear(v)
            return u
        z = any_as_fmpz(t)
        if z is not NotImplemented:
            u = arf.__new__(arf)
            arf_sub_fmpz((<arf>u).val, (<arf>s).val, (<fmpz>z).val, getprec(), thectx.rnd)
            return u
        return NotImplemented

    def __rsub__(s, t):
        cdef slong t_si
        cdef arf_t v
        cdef arf u
        if typecheck(t, arf):
            u = arf.__new__(arf)
            arf_sub((<arf>u).val, (<arf>t).val, (<arf>s).val, getprec(), thectx.rnd)
            return u
        elif typecheck(t, int):
            try:
                t_si = <slong>t
            except OverflowError:
                pass
            else:
                u = arf.__new__(arf)
                arf_init(v)
                arf_set_si(v, t_si)
                arf_sub((<arf>u).val, v, (<arf>s).val, getprec(), thectx.rnd)
                arf_clear(v)
                return u
        elif typecheck(t, float):
            u = arf.__new__(arf)
            arf_init(v)
            arf_set_d(v, PyFloat_AsDouble(t))
            arf_sub((<arf>u).val, v, (<arf>s).val, getprec(), thectx.rnd)
            arf_clear(v)
            return u
        z = any_as_fmpz(t)
        if z is not NotImplemented:
            u = arf.__new__(arf)
            arf_init(v)
            arf_set_fmpz(v, (<fmpz>z).val)
            arf_sub((<arf>u).val, v, (<arf>s).val, getprec(), thectx.rnd)
            arf_clear(v)
            return u
        return NotImplemented

    def __mul__(s, t):
        cdef slong t_si
        cdef arf_t v
        cdef arf u
        if typecheck(t, arf):
            u = arf.__new__(arf)
            arf_mul((<arf>u).val, (<arf>s).val, (<arf>t).val, getprec(), thectx.rnd)
            return u
        elif typecheck(t, int):
            try:
                t_si = <slong>t
            except OverflowError:
                pass
            else:
                u = arf.__new__(arf)
                arf_mul_si((<arf>u).val, (<arf>s).val, t_si, getprec(), thectx.rnd)
                return u
        elif typecheck(t, float):
            u = arf.__new__(arf)
            arf_init(v)
            arf_set_d(v, PyFloat_AsDouble(t))
            arf_mul((<arf>u).val, (<arf>s).val, v, getprec(), thectx.rnd)
            arf_clear(v)
            return u
        z = any_as_fmpz(t)
        if z is not NotImplemented:
            u = arf.__new__(arf)
            arf_mul_fmpz((<arf>u).val, (<arf>s).val, (<fmpz>z).val, getprec(), thectx.rnd)
            return u
        return NotImplemented

    def __rmul__(s, t):
        cdef slong t_si
        cdef arf_t v
        cdef arf u
        if typecheck(t, arf):
            u = arf.__new__(arf)
            arf_mul((<arf>u).val, (<arf>t).val, (<arf>s).val, getprec(), thectx.rnd)
            return u
        elif typecheck(t, int):
            try:
                t_si = <slong>t
            except OverflowError:
                pass
            else:
                u = arf.__new__(arf)
                arf_mul_si((<arf>u).val, (<arf>s).val, t_si, getprec(), thectx.rnd)
                return u
        elif typecheck(t, float):
            u = arf.__new__(arf)
            arf_init(v)
            arf_set_d(v, PyFloat_AsDouble(t))
            arf_mul((<arf>u).val, v, (<arf>s).val, getprec(), thectx.rnd)
            arf_clear(v)
            return u
        z = any_as_fmpz(t)
        if z is not NotImplemented:
            u = arf.__new__(arf)
            arf_mul_fmpz((<arf>u).val, (<arf>s).val, (<fmpz>z).val, getprec(), thectx.rnd)
            return u
        return NotImplemented

    def __truediv__(s, t):
        cdef slong t_si
        cdef arf_t v
        cdef arf u
        if typecheck(t, arf):
            u = arf.__new__(arf)
            arf_div((<arf>u).val, (<arf>s).val, (<arf>t).val, getprec(), thectx.rnd)
            return u
        elif typecheck(t, int):
            try:
                t_si = <slong>t
            except OverflowError:
                pass
            else:
                u = arf.__new__(arf)
                arf_div_si((<arf>u).val, (<arf>s).val, t_si, getprec(), thectx.rnd)
                return u
        elif typecheck(t, float):
            u = arf.__new__(arf)
            arf_init(v)
            arf_set_d(v, PyFloat_AsDouble(t))
            arf_div((<arf>u).val, (<arf>s).val, v, getprec(), thectx.rnd)
            arf_clear(v)
            return u
        z = any_as_fmpz(t)
        if z is not NotImplemented:
            u = arf.__new__(arf)
            arf_div_fmpz((<arf>u).val, (<arf>s).val, (<fmpz>z).val, getprec(), thectx.rnd)
            return u
        return NotImplemented

    def __rtruediv__(s, t):
        cdef slong t_si
        cdef arf_t v
        cdef arf u
        if typecheck(t, arf):
            u = arf.__new__(arf)
            arf_div((<arf>u).val, (<arf>t).val, (<arf>s).val, getprec(), thectx.rnd)
            return u
        elif typecheck(t, int):
            try:
                t_si = <slong>t
            except OverflowError:
                pass
            else:
                u = arf.__new__(arf)
                arf_si_div((<arf>u).val, t_si, (<arf>s).val, getprec(), thectx.rnd)
                return u
        elif typecheck(t, float):
            u = arf.__new__(arf)
            arf_init(v)
            arf_set_d(v, PyFloat_AsDouble(t))
            arf_div((<arf>u).val, v, (<arf>s).val, getprec(), thectx.rnd)
            arf_clear(v)
            return u
        z = any_as_fmpz(t)
        if z is not NotImplemented:
            u = arf.__new__(arf)
            arf_fmpz_div((<arf>u).val, (<fmpz>z).val, (<arf>s).val, getprec(), thectx.rnd)
            return u
        return NotImplemented
