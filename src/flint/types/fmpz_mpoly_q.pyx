from flint.flint_base.flint_base cimport flint_rational_function
from flint.utils.typecheck cimport typecheck
from flint.flintlib.fmpz_mpoly cimport fmpz_mpoly_set, fmpz_mpoly_get_str_pretty
from flint.flintlib.fmpz_mpoly_q cimport *
from flint.types.fmpz_mpoly cimport fmpz_mpoly, fmpz_mpoly_ctx

cdef inline init_fmpz_mpoly_q(fmpz_mpoly_q var, fmpz_mpoly_ctx ctx):
    var.ctx = ctx
    fmpz_mpoly_q_init(var.fraction, ctx.val)
    var._init = True

cdef inline create_fmpz_mpoly_q(fmpz_mpoly_ctx ctx):
    cdef fmpz_mpoly_q var
    var = fmpz_mpoly_q.__new__(fmpz_mpoly_q)
    var.ctx = ctx
    fmpz_mpoly_q_init(var.fraction, ctx.val)
    var._init = True
    return var


cdef class fmpz_mpoly_q(flint_rational_function):
    """
    The `fmpz_mpoly_q` represents multivariate rational functions
    over the integers
    """
    def __cinit__(self):
        self._init = False

    def __dealloc__(self):
        if self._init:
            fmpz_mpoly_q_clear(self.fraction, self.ctx.val)
            self._init = False

    def __init__(self, num, den):
        if typecheck(num, fmpz_mpoly) and typecheck(den, fmpz_mpoly):
            if (<fmpz_mpoly>num).ctx == (<fmpz_mpoly>den).ctx:
                self.ctx = (<fmpz_mpoly>num).ctx
                fmpz_mpoly_q_init(self.fraction, self.ctx.val)
                fmpz_mpoly_set(fmpz_mpoly_q_numref(self.fraction),
                               (<fmpz_mpoly>num).val, self.ctx.val)
                fmpz_mpoly_set(fmpz_mpoly_q_denref(self.fraction),
                               (<fmpz_mpoly>den).val, self.ctx.val)
                self._init = True
            else:
                raise ValueError("numerator and denominator must have identical contexts")
        else:
            raise TypeError("fmpz_mpoly_q is a fraction of two fmpz_mpolys fs")


    def __nonzero__(self):
        return not fmpz_mpoly_q_is_zero(self.fraction, self.ctx.val)

    def __bool__(self):
        return not fmpz_mpoly_q_is_zero(self.fraction, self.ctx.val)

    def is_one(self):
        return fmpz_mpoly_q_is_one(self.fraction, self.ctx.val)

    def __richcmp__(self, other, int op):
        if op != 2 and op != 3:
            return NotImplemented
        if typecheck(self, fmpz_mpoly_q) and typecheck(other, fmpz_mpoly_q):
            if (<fmpz_mpoly_q> self).ctx is (<fmpz_mpoly_q> other).ctx:
                if op == 2:
                    return bool(fmpz_mpoly_q_equal((<fmpz_mpoly_q>self).fraction, (<fmpz_mpoly_q>other).fraction, self.ctx.val))
                else:
                    return not bool(fmpz_mpoly_q_equal((<fmpz_mpoly_q>self).fraction, (<fmpz_mpoly_q>other).fraction, self.ctx.val))
            else:
                if op == 2:
                    return False
                else:
                    return True
        if op == 2:
            return not bool(self - other)
        else:
            return bool(self - other)

    def repr(self):
        return self.str() + "  (nvars=%s, ordering=%s names=%s)" % (self.ctx.nvars(), self.ctx.ordering(), self.ctx.py_names)

    def str(self):
        cdef bytes numerator = fmpz_mpoly_get_str_pretty(&(self.fraction.num), self.ctx.c_names, self.ctx.val)
        cdef bytes denominator = fmpz_mpoly_get_str_pretty(&(self.fraction.den), self.ctx.c_names, self.ctx.val)
        res = str(b"(" + numerator + b")/(" + denominator + b")", encoding='utf-8')
        res = res.replace("+", " + ")
        res = res.replace("-", " - ")
        return res

    def __neg__(self):
        cdef fmpz_mpoly_q res
        res = create_fmpz_mpoly_q(self.ctx)
        fmpz_mpoly_q_neg(res.fraction, self.fraction, res.ctx.val)
