from flint.flintlib.functions.fmpz cimport fmpz_struct
from flint.flintlib.types.flint cimport slong

cdef class fmpz_vec:
    cdef fmpz_struct *val
    cdef fmpz_struct **double_indirect
    cdef slong length
