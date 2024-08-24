from flint.flint_base.flint_base cimport flint_scalar
from flint.flintlib.flint cimport mp_limb_t
from flint.flintlib.nmod cimport nmod_t

cdef int any_as_nmod(mp_limb_t * val, obj, nmod_t mod) except -1
cdef nmod_ctx any_as_nmod_ctx(obj)


cdef class nmod_ctx:
    cdef nmod_t mod
    cdef bint _is_prime

    cdef int any_as_nmod(self, mp_limb_t * val, obj) except -1
    cdef nmod _new(self, mp_limb_t * val)


cdef class nmod(flint_scalar):
    cdef mp_limb_t val
    cdef nmod_ctx ctx
