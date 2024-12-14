from flint.flint_base.flint_base cimport flint_scalar
from flint.flintlib.functions.fmpq cimport fmpq_t

cdef int fmpq_set_any_ref(fmpq_t x, obj)
cdef any_as_fmpq(obj)

cdef class fmpq(flint_scalar):
    cdef fmpq_t val
