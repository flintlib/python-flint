cimport cython

from flint.flint_base.flint_base cimport flint_mat

from flint.flintlib.nmod cimport nmod_t
from flint.flintlib.nmod_mat cimport nmod_mat_t
from flint.flintlib.flint cimport mp_limb_t, ulong

from flint.types.nmod cimport nmod_ctx, nmod
from flint.types.nmod_poly cimport nmod_poly_ctx, nmod_poly


@cython.no_gc
cdef class nmod_mat_ctx:
    cdef nmod_t mod
    cdef bint _is_prime
    cdef nmod_ctx scalar_ctx
    cdef nmod_poly_ctx poly_ctx

    @staticmethod
    cdef nmod_mat_ctx any_as_nmod_mat_ctx(obj)

    @staticmethod
    cdef nmod_mat_ctx _get_ctx(int mod)

    @staticmethod
    cdef nmod_mat_ctx _new_ctx(ulong mod)

    @cython.final
    cdef int any_as_nmod(self, mp_limb_t * val, obj) except -1

    @cython.final
    cdef any_as_nmod_mat(self, obj)

    @cython.final
    cdef nmod new_nmod(self)

    @cython.final
    cdef nmod_poly new_nmod_poly(self)

    @cython.final
    cdef nmod_mat new_nmod_mat(self, ulong m, ulong n)

    @cython.final
    cdef nmod_mat new_nmod_mat_copy(self, nmod_mat other)


@cython.no_gc
cdef class nmod_mat(flint_mat):
    cdef nmod_mat_t val
    cdef nmod_mat_ctx ctx

    cpdef long nrows(self)
    cpdef long ncols(self)
    cpdef mp_limb_t modulus(self)
    cdef __mul_nmod(self, mp_limb_t c)
