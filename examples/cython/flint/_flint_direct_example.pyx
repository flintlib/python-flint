from flint cimport fmpz, fmpz_get_value

cdef extern from "flint/fmpz.h":
    ctypedef long flint_fmpz "fmpz"
    unsigned long fmpz_bits(const flint_fmpz *)


def bit_length(fmpz value):
    """Call a read-only FLINT operation on borrowed fmpz storage."""
    return fmpz_bits(<const flint_fmpz *>fmpz_get_value(value))
