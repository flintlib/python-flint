from flint.flint_base.flint_base cimport flint_poly

from flint.flintlib.nmod_poly cimport nmod_poly_t
from flint.flintlib.flint cimport mp_limb_t

cdef class nmod_poly(flint_poly):
    cdef nmod_poly_t val
    cpdef long length(self)
    cpdef long degree(self)
    cpdef mp_limb_t modulus(self)
