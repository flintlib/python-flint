from flint.flintlib.types.flint cimport (
    mp_ptr,
    fmpz_struct,
    ulong,
    slong,
    flint_bitcnt_t,
)


cdef extern from "flint/fmpz_types.h":

    ctypedef fmpz_struct fmpz_t[1]

    ctypedef struct fmpz_preinvn_struct:
        mp_ptr dinv
        slong n
        flint_bitcnt_t norm
    ctypedef fmpz_preinvn_struct fmpz_preinvn_t[1]

    ctypedef struct fmpz_factor_struct:
        int sign
        fmpz_struct * p
        ulong * exp
        slong alloc
        slong num
    ctypedef fmpz_factor_struct fmpz_factor_t[1]

    ctypedef struct fmpz_mat_struct:
        fmpz_struct * entries
        slong r
        slong c
        fmpz_struct ** rows
    ctypedef fmpz_mat_struct fmpz_mat_t[1]
    long fmpz_mat_nrows(fmpz_mat_t mat)
    long fmpz_mat_ncols(fmpz_mat_t mat)

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


cdef extern from "flint/fmpz_factor.h":
    # XXX: Missing from docs...
    void fmpz_factor_expand(fmpz_t n, const fmpz_factor_t factor)


cdef extern from "flint/fmpz_lll.h":

    cdef enum _rep_type:
        GRAM
        Z_BASIS
    ctypedef _rep_type rep_type

    cdef enum _gram_type:
        APPROX
        EXACT
    ctypedef _gram_type gram_type

    ctypedef struct fmpz_lll_struct:
        double delta
        double eta
        rep_type rt
        gram_type gt
    ctypedef fmpz_lll_struct fmpz_lll_t[1]