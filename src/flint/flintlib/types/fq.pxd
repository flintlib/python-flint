from flint.flintlib.types.flint cimport fmpz_struct, slong
from flint.flintlib.types.fmpz cimport fmpz_poly_t, fmpz_poly_struct
from flint.flintlib.functions.fmpz_mod cimport fmpz_mod_ctx_t
from flint.flintlib.functions.fmpz_mod_poly cimport fmpz_mod_poly_t


cdef extern from "flint/fq.h":
    # Type definitions **********************************************/
    ctypedef fmpz_poly_t fq_t
    ctypedef fmpz_poly_struct fq_struct

    ctypedef struct fq_ctx_struct:
        fmpz_mod_ctx_t ctxp

        int sparse_modulus
        int is_conway  # whether field was initialized with the Flint Conway tables  (assures primitivity)

        fmpz_struct * a
        slong * j
        slong len

        fmpz_mod_poly_t modulus
        fmpz_mod_poly_t inv

        char * var
    ctypedef fq_ctx_struct fq_ctx_t[1]

    ctypedef struct fq_mat_struct:
        fq_struct * entries
        slong r
        slong s
        fq_struct ** rows
    ctypedef fq_mat_struct fq_mat_t[1]

    ctypedef struct fq_poly_struct:
        fq_struct * coeffs
        slong alloc
        slong length
    ctypedef fq_poly_struct fq_poly_t[1]

    ctypedef struct fq_poly_factor_struct:
        fq_poly_struct * poly
        slong * exp
        slong num
        slong alloc

    ctypedef fq_poly_factor_struct fq_poly_factor_t[1]
