from flint.flintlib.types.flint cimport slong, fmpz_struct, fmpz_t
from flint.flintlib.types.fmpz cimport fmpz_mat_struct
from flint.flintlib.types.fmpz_mod cimport fmpz_mod_ctx_struct


cdef extern from "flint/fmpz_mod_poly.h":

    #  Type definitions *********************************************************/
    ctypedef struct fmpz_mod_poly_struct:
        fmpz_struct * coeffs
        slong alloc
        slong length

    ctypedef fmpz_mod_poly_struct fmpz_mod_poly_t[1]

    ctypedef struct fmpz_mod_poly_res_struct:
        fmpz_t res
        fmpz_t lc
        slong len0
        slong len1
        slong off

    ctypedef fmpz_mod_poly_res_struct fmpz_mod_poly_res_t[1]

    ctypedef struct fmpz_mod_poly_frobenius_powers_2exp_struct:
        fmpz_mod_poly_struct * pow
        slong len

    ctypedef fmpz_mod_poly_frobenius_powers_2exp_struct fmpz_mod_poly_frobenius_powers_2exp_t[1]

    ctypedef struct fmpz_mod_poly_frobenius_powers_struct:
        fmpz_mod_poly_struct * pow
        slong len

    ctypedef fmpz_mod_poly_frobenius_powers_struct fmpz_mod_poly_frobenius_powers_t[1]

    ctypedef struct fmpz_mod_poly_matrix_precompute_arg_t:
        fmpz_mat_struct * A
        fmpz_mod_poly_struct * poly1
        fmpz_mod_poly_struct * poly2
        fmpz_mod_poly_struct * poly2inv
        const fmpz_mod_ctx_struct * ctx

    ctypedef struct fmpz_mod_poly_compose_mod_precomp_preinv_arg_t:
        fmpz_mat_struct * A
        fmpz_mod_poly_struct * res
        fmpz_mod_poly_struct * poly1
        fmpz_mod_poly_struct * poly3
        fmpz_mod_poly_struct * poly3inv
        const fmpz_mod_ctx_struct * ctx

    # Radix conversion *********************************************************/
    ctypedef struct fmpz_mod_poly_radix_struct:
        fmpz_struct *V
        fmpz_struct *W
        fmpz_struct **Rpow
        fmpz_struct **Rinv
        slong degR
        slong k
        fmpz_struct invL

    ctypedef fmpz_mod_poly_radix_struct fmpz_mod_poly_radix_t[1]

    # Berlekamp-Massey Algorithm - see fmpz_mod_poly/berlekamp_massey.c for more info ********/
    ctypedef struct fmpz_mod_berlekamp_massey_struct:
        slong npoints
        fmpz_mod_poly_t R0, R1
        fmpz_mod_poly_t V0, V1
        fmpz_mod_poly_t qt, rt
        fmpz_mod_poly_t points

    ctypedef fmpz_mod_berlekamp_massey_struct fmpz_mod_berlekamp_massey_t[1]
