from flint.flintlib.fmpq cimport fmpq_struct
from flint.flintlib.flint cimport slong

cdef class fmpq_vec:
    cdef fmpq_struct *val
    cdef fmpq_struct **double_indirect
    cdef slong length
