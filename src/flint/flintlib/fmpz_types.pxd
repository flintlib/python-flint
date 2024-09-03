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


cdef extern from "flint/fmpz_poly.h":
    ctypedef struct fmpz_poly_struct:
        fmpz_struct * coeffs
        slong alloc
        slong length
    ctypedef fmpz_poly_struct fmpz_poly_t[1]

    ctypedef struct fmpz_poly_factor_struct:
        fmpz_struct c
        fmpz_poly_struct *p
        slong *exp
        slong num
        slong alloc
    ctypedef fmpz_poly_factor_struct fmpz_poly_factor_t[1]
