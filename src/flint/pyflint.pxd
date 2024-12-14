from flint.flintlib.functions.arf cimport arf_rnd_t
from flint.flint_base.flint_base cimport flint_mat
from flint.flint_base.flint_base cimport flint_mpoly
from flint.flint_base.flint_base cimport flint_series
from flint.flint_base.flint_base cimport flint_scalar
from flint.flint_base.flint_base cimport flint_poly

from flint.types.fmpz cimport fmpz
from flint.flintlib.types.flint cimport *

cdef flint_rand_t global_random_state

cdef class Context:
    cdef public bint pretty
    cdef public long prec
    cdef arf_rnd_t rnd
