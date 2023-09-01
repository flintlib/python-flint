from flint.flint_base.flint_base cimport flint_mat
from flint._flint cimport fmpz_mat_t
from flint.fmpz cimport fmpz
cdef any_as_fmpz_mat(obj)

cdef class fmpz_mat(flint_mat):

    cdef fmpz_mat_t val
    cpdef long nrows(self)
    cpdef long ncols(self)
    cdef __mul_fmpz(self, fmpz c)
