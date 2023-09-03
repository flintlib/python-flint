from flint.flint_base.flint_base cimport flint_mat

from flint._flint cimport acb_mat_t

cdef class acb_mat(flint_mat):
    cdef acb_mat_t val
    cpdef long nrows(self)
    cpdef long ncols(self)
