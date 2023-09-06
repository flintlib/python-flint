from flint.flint_base.flint_base cimport flint_scalar

from flint.flintlib.fmpq cimport fmpq_t

cdef any_as_fmpq(obj)
cdef class fmpq(flint_scalar):
    cdef fmpq_t val


