from flint.flintlib.functions.fq_default_poly cimport *
from flint.flintlib.functions.fq_default_poly_factor cimport *
from flint.flintlib.functions.fq_default cimport fq_default_neg

from flint.flint_base.flint_base cimport flint_poly
from flint.types.fq_default cimport fq_default_ctx


cdef class fq_default_poly_ctx:
    cdef fq_default_ctx field
    cdef any_as_fq_default_poly(self, obj)
    cdef set_any_as_fq_default_poly(self, fq_default_poly_t poly, obj)
    cdef set_list_as_fq_default_poly(self, fq_default_poly_t poly, val)
    cdef new_ctype_poly(self)

cdef class fq_default_poly(flint_poly):
    cdef fq_default_poly_t val
    cdef fq_default_poly_ctx ctx
    cpdef long length(self)
    cpdef long degree(self)
