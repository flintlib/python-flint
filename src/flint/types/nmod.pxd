from flint.flint_base.flint_base cimport flint_scalar
from flint.flintlib.types.flint cimport mp_limb_t
from flint.flintlib.functions.nmod cimport nmod_t

cdef int any_as_nmod(mp_limb_t * val, obj, nmod_t mod) except -1

cdef class nmod(flint_scalar):
    cdef mp_limb_t val
    cdef nmod_t mod
