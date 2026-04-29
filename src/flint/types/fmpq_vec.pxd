from flint.flintlib.types.flint cimport slong
from flint.flintlib.types.fmpq cimport fmpq_struct

cdef class fmpq_vec:
    cdef fmpq_struct *val
    cdef fmpq_struct **double_indirect
    cdef slong length
