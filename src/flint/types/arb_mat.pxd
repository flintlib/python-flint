from flint.flint_base.flint_base cimport flint_mat
from flint.flintlib.functions.arb_mat cimport arb_mat_t

cdef class arb_mat(flint_mat):
    cdef arb_mat_t val
    cpdef long nrows(self)
    cpdef long ncols(self)
