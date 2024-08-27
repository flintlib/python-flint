cimport cython

from flint.flintlib.nmod cimport nmod_t
from flint.flintlib.nmod_poly cimport nmod_poly_t
from flint.flintlib.flint cimport mp_limb_t, ulong

from flint.flint_base.flint_base cimport flint_poly

from flint.types.nmod cimport nmod_ctx, nmod


@cython.no_gc
cdef class nmod_poly_ctx:
    cdef nmod_t mod
    cdef bint _is_prime
    cdef nmod_ctx scalar_ctx

    @staticmethod
    cdef any_as_nmod_poly_ctx(obj)
    @staticmethod
    cdef nmod_poly_ctx _get_ctx(int mod)
    @staticmethod
    cdef nmod_poly_ctx _new_ctx(ulong mod)

    @cython.final
    cdef nmod_poly_set_list(self, nmod_poly_t poly, list val)
    @cython.final
    cdef int any_as_nmod(self, mp_limb_t * val, obj) except -1
    @cython.final
    cdef any_as_nmod_poly(self, obj)
    @cython.final
    cdef nmod new_nmod(self)
    @cython.final
    cdef nmod_poly new_nmod_poly(self)


@cython.no_gc
cdef class nmod_poly(flint_poly):
    cdef nmod_poly_t val
    cdef nmod_poly_ctx ctx

    cpdef long length(self)
    cpdef long degree(self)
    cpdef mp_limb_t modulus(self)
