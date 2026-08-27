from flint._capi.cython cimport (
    flint_fmpz,
    fmpz,
    fmpz_from_value,
    fmpz_get_value,
)

cdef extern from "flint/fmpz.h":
    unsigned long fmpz_bits(const flint_fmpz *)


def bit_length(fmpz value):
    """Call a read-only FLINT operation on borrowed fmpz storage."""
    return fmpz_bits(fmpz_get_value(value))


def copy(fmpz value):
    """Create a new Python fmpz by copying a borrowed FLINT value."""
    return fmpz_from_value(fmpz_get_value(value))
