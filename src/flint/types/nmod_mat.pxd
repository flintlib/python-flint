from flint.flint_base.flint_base cimport flint_mat

from flint.flintlib.functions.nmod_mat cimport nmod_mat_t
from flint.flintlib.types.flint cimport mp_limb_t

cdef class nmod_mat(flint_mat):
    cdef nmod_mat_t val
    cpdef long nrows(self)
    cpdef long ncols(self)
    cpdef mp_limb_t modulus(self)
    cdef __mul_nmod(self, mp_limb_t c)
