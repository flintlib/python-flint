from flint.flintlib.flint cimport (
    mp_ptr,
    fmpz_struct,
    slong,
    flint_bitcnt_t,
)


cdef extern from "flint/fmpz.h":
    ctypedef fmpz_struct fmpz_t[1]

    ctypedef struct fmpz_preinvn_struct:
        mp_ptr dinv
        slong n
        flint_bitcnt_t norm
    ctypedef fmpz_preinvn_struct fmpz_preinvn_t[1]
