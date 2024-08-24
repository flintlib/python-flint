from flint.flintlib.nmod cimport nmod_t
from flint.flintlib.nmod_poly cimport nmod_poly_t
from flint.flintlib.flint cimport mp_limb_t

from flint.flint_base.flint_base cimport flint_poly

from flint.types.nmod cimport nmod_ctx


cdef nmod_poly_ctx any_as_nmod_poly_ctx(obj)
cdef nmod_poly nmod_poly_new_init(nmod_poly_ctx ctx)


cdef class nmod_poly_ctx:
    cdef nmod_ctx ctx
    cdef nmod_t mod
    cdef bint _is_prime

    cdef nmod_poly_set_list(self, nmod_poly_t poly, list val)
    cdef int any_as_nmod(self, mp_limb_t * val, obj) except -1
    cdef any_as_nmod_poly(self, obj)


cdef class nmod_poly(flint_poly):
    cdef nmod_poly_t val
    cdef nmod_poly_ctx ctx

    cpdef long length(self)
    cpdef long degree(self)
    cpdef mp_limb_t modulus(self)
