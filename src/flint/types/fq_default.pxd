from flint.flintlib.fq_default cimport *

from flint.types.fmpz cimport fmpz
from flint.flint_base.flint_base cimport flint_scalar

cpdef enum fq_default_type:
    DEFAULT = 0
    FQ_ZECH = 1
    FQ_NMOD = 2
    FQ      = 3

cdef class fq_default_ctx:
    cdef fq_default_ctx_t val
    cdef readonly char *var
    cdef bint initialized

    cdef set_any_as_fq_default(self, fq_default_t val, obj)
    cdef any_as_fmpz_mod(self, obj)
    
    @staticmethod
    cdef fq_default_ctx c_from_order(fmpz p, int d, char *var, fq_default_type type=*)
    @staticmethod
    cdef fq_default_ctx c_from_modulus(modulus, char *var, fq_default_type type=*)


cdef class fq_default(flint_scalar):
    cdef fq_default_ctx ctx
    cdef fq_default_t val
