"""Experimental Cython API for python-flint consumers."""

cdef extern class flint.types.fmpz.fmpz [object PyFlintFMPZObject, check_size ignore]:
    pass

cdef extern from "flint/fmpz.h":
    ctypedef long flint_fmpz "fmpz"

cdef api const flint_fmpz *fmpz_get_value(fmpz value) except NULL
cdef api fmpz fmpz_from_value(const flint_fmpz *value)
