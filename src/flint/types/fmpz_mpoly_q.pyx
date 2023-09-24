from flint.flint_base.flint_base cimport flint_rational_function
from flint.utils.typecheck cimport typecheck
from flint.flintlib.fmpz_mpoly cimport fmpz_mpoly_set
from flint.flintlib.fmpz_mpoly_q cimport *
from flint.types.fmpz_mpoly cimport fmpz_mpoly, fmpz_mpoly_ctx



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
