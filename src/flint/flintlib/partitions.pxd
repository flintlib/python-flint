from flint.flintlib.fmpz cimport fmpz_t

cdef extern from "flint/partitions.h":
    void partitions_fmpz_fmpz(fmpz_t, const fmpz_t, int)
