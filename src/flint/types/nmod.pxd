cimport cython

from flint.flint_base.flint_base cimport flint_scalar
from flint.flintlib.flint cimport mp_limb_t, ulong
from flint.flintlib.nmod cimport nmod_t


@cython.no_gc
cdef class nmod_ctx:
    cdef nmod_t mod
    cdef bint _is_prime

    @staticmethod
    cdef any_as_nmod_ctx(obj)
    @staticmethod
    cdef _get_ctx(int mod)
    @staticmethod
    cdef _new_ctx(ulong mod)

    @cython.final
    cdef int any_as_nmod(self, mp_limb_t * val, obj) except -1
    @cython.final
    cdef nmod new_nmod(self)


cdef class nmod(flint_scalar):
    cdef mp_limb_t val
    cdef nmod_ctx ctx
