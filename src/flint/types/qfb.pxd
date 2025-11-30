from flint.flintlib.functions.qfb cimport *
from flint.types.fmpz cimport fmpz

cdef class qfb:
    cdef qfb_t val
    # the discriminant, stored for convenience
    cdef fmpz D
