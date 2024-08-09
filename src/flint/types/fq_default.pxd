from flint.flintlib.fq_default cimport *
from flint.flintlib.fq_zech cimport fq_zech_is_primitive, fq_zech_multiplicative_order
from flint.flintlib.fq_nmod cimport fq_nmod_is_primitive, fq_nmod_multiplicative_order
from flint.flintlib.fq cimport fq_is_primitive, fq_multiplicative_order
from flint.types.fmpz cimport fmpz
from flint.flint_base.flint_base cimport flint_scalar

cpdef enum fq_default_type:
    DEFAULT  = 0
    FQ_ZECH  = 1
    FQ_NMOD  = 2
    FQ       = 3
    NMOD     = 4
    FMPZ_MOD = 5

cdef class fq_default_ctx:
    cdef fq_default_ctx_t val
    cdef readonly char *var
    cdef bint _initialized

    cdef new_ctype_fq_default(self)
    cdef set_list_as_fq_default(self, fq_default_t val, obj)
    cdef set_any_scalar_as_fq_default(self, fq_default_t fq_ele, obj)
    cdef set_any_as_fq_default(self, fq_default_t val, obj)
    cdef any_as_fq_default(self, obj)

    cdef _c_set_from_order(self, fmpz p, int d, char *var, fq_default_type fq_type=*)
    cdef _c_set_from_modulus(self, modulus, char *var, fq_default_type fq_type=*)

cdef class fq_default(flint_scalar):
    cdef fq_default_ctx ctx
    cdef fq_default_t val

    # Arithmetic for flint_scalar base class
    cpdef fq_default _neg_(fq_default self)
    cpdef fq_default _add_(fq_default self, fq_default other)
    cpdef fq_default _mul_(fq_default self, fq_default other)
    cpdef fq_default _sub_(fq_default self, fq_default other)
    cpdef fq_default _rsub_(fq_default self, fq_default other)
    cpdef fq_default _div_(fq_default self, fq_default other)
    cpdef fq_default _rdiv_(fq_default self, fq_default other)
    cpdef fq_default _invert_(fq_default self)
