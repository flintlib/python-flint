from flint.flintlib.types.flint cimport mp_limb_t, slong
from flint.flintlib.functions.fmpz cimport fmpz_struct
from flint.flintlib.functions.nmod cimport nmod_t
from flint.flintlib.functions.nmod_poly cimport nmod_poly_struct, nmod_poly_t
from flint.flintlib.types.mpoly cimport mpoly_ctx_t


cdef extern from "flint/fq_nmod.h":
    # Type definitions **********************************************/
    ctypedef nmod_poly_t fq_nmod_t
    ctypedef nmod_poly_struct fq_nmod_struct

    ctypedef struct fq_nmod_ctx_struct:
        fmpz_struct p
        nmod_t mod

        int sparse_modulus
        int is_conway  # whether field was generated using Flint Conway table (assures primitivity

        mp_limb_t *a
        slong *j
        slong len

        nmod_poly_t modulus
        nmod_poly_t inv

        char *var
    ctypedef fq_nmod_ctx_struct fq_nmod_ctx_t[1]

    ctypedef struct fq_nmod_poly_struct:
        fq_nmod_struct * coeffs
        slong alloc
        slong length
    ctypedef fq_nmod_poly_struct fq_nmod_poly_t[1]

    ctypedef struct fq_nmod_poly_factor_struct:
        fq_nmod_poly_struct * poly
        slong * exp
        slong num
        slong alloc
    ctypedef fq_nmod_poly_factor_struct fq_nmod_poly_factor_t[1]

    ctypedef struct fq_nmod_mat_struct:
        fq_nmod_struct * entries
        slong r
        slong s
        fq_nmod_struct ** rows
    ctypedef fq_nmod_mat_struct fq_nmod_mat_t[1]

    ctypedef struct fq_nmod_mpoly_ctx_struct:
        mpoly_ctx_t minfo
        fq_nmod_ctx_t fqctx

    ctypedef fq_nmod_mpoly_ctx_struct fq_nmod_mpoly_ctx_t[1]
