from flint.flint_base.flint_base cimport flint_mat

from flint.flintlib.functions.fmpq_mat cimport fmpq_mat_t
from flint.types.fmpz cimport fmpz
from flint.types.fmpq cimport fmpq
from flint.types.fmpz_mat cimport fmpz_mat

cdef class fmpq_mat(flint_mat):
    cdef fmpq_mat_t val

    cpdef long nrows(self)
    cpdef long ncols(self)
    cdef __mul_fmpz(self, fmpz c)
    cdef __mul_fmpq(self, fmpq c)
    cdef __mul_fmpq_mat(self, fmpq_mat other)
    cdef __mul_fmpz_mat(self, fmpz_mat other)
    cdef __mul_r_fmpz_mat(self, fmpz_mat other)
