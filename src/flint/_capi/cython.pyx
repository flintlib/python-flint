"""Implementation of python-flint's experimental Cython API."""

from flint.types.fmpz cimport fmpz as fmpz_internal
from flint.flintlib.functions.fmpz cimport fmpz_set


cdef api const flint_fmpz *fmpz_get_value(fmpz value) except NULL:
    """Return a borrowed, read-only pointer to an fmpz value."""
    return (<fmpz_internal>value).val


cdef api fmpz fmpz_from_value(const flint_fmpz *value):
    """Return a new fmpz containing a copy of a FLINT fmpz value."""
    cdef fmpz_internal result

    if value == NULL:
        raise ValueError("fmpz_from_value does not accept NULL")

    result = fmpz_internal.__new__(fmpz_internal)
    fmpz_set(result.val, value)
    return result
