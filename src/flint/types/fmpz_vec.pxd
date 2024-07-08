from flint.flintlib.fmpz cimport fmpz_struct
from flint.flintlib.flint cimport slong

cdef class fmpz_vec:
    cdef fmpz_struct *val
    cdef fmpz_struct **double_indirect
    cdef slong length
