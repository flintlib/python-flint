from flint.flint_base.flint_base cimport flint_scalar
from flint.flintlib.flint cimport mp_limb_t, ulong
from flint.flintlib.nmod cimport nmod_t


cdef class nmod_ctx:
    cdef nmod_t mod
    cdef bint _is_prime

    @staticmethod
    cdef nmod_ctx any_as_nmod_ctx(obj)
    @staticmethod
    cdef _get_ctx(int mod)
    @staticmethod
    cdef _new_ctx(ulong mod)

    cdef int any_as_nmod(self, mp_limb_t * val, obj) except -1
    cdef nmod new_nmod(self)


cdef class nmod(flint_scalar):
    cdef mp_limb_t val
    cdef nmod_ctx ctx
