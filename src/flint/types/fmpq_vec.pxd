from flint.flintlib.functions.fmpq cimport fmpq_struct
from flint.flintlib.types.flint cimport slong

cdef class fmpq_vec:
    cdef fmpq_struct *val
    cdef fmpq_struct **double_indirect
    cdef slong length
