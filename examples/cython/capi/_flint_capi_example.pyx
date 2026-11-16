from flint cimport fmpz, fmpz_add


def add(fmpz a, fmpz b):
    """Add through python-flint's capsule; no FLINT linkage is needed."""
    return fmpz_add(a, b)
