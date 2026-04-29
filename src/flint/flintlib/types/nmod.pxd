from flint.flintlib.types.flint cimport ulong, slong, nmod_t, fmpz_struct, nn_ptr, flint_bitcnt_t
from flint.flintlib.types.mpoly cimport mpoly_ctx_t


cdef extern from "flint/nmod_types.h":

    ctypedef struct nmod_mat_struct:
        ulong * entries
        slong r
        slong c
        ulong ** rows
        nmod_t mod

    ctypedef nmod_mat_struct nmod_mat_t[1]

    # Macros:
    ulong nmod_mat_entry(nmod_mat_t mat, slong i, slong j)

    ctypedef struct nmod_poly_struct:
        nn_ptr coeffs
        slong alloc
        slong length
        nmod_t mod

    ctypedef nmod_poly_struct nmod_poly_t[1]

    ctypedef struct nmod_poly_factor_struct:
        nmod_poly_struct * p
        slong *exp
        slong num
        slong alloc

    ctypedef nmod_poly_factor_struct nmod_poly_factor_t[1]

    ctypedef struct nmod_poly_mat_struct:
        nmod_poly_struct * entries
        slong r
        slong c
        nmod_poly_struct ** rows
        ulong modulus

    ctypedef nmod_poly_mat_struct nmod_poly_mat_t[1]

    ctypedef struct nmod_mpoly_ctx_struct:
        mpoly_ctx_t minfo
        nmod_t mod

    ctypedef nmod_mpoly_ctx_struct nmod_mpoly_ctx_t[1]

    ctypedef struct nmod_mpoly_struct:
        ulong * coeffs
        ulong * exps
        slong length
        flint_bitcnt_t bits
        slong coeffs_alloc
        slong exps_alloc

    ctypedef nmod_mpoly_struct nmod_mpoly_t[1]

    ctypedef struct nmod_mpoly_univar_struct:
        nmod_mpoly_struct * coeffs
        fmpz_struct * exps
        slong alloc
        slong length

    ctypedef nmod_mpoly_univar_struct nmod_mpoly_univar_t[1]

    ctypedef struct nmod_mpoly_factor_struct:
        ulong constant
        nmod_mpoly_struct * poly
        fmpz_struct * exp
        slong num
        slong alloc

    ctypedef nmod_mpoly_factor_struct nmod_mpoly_factor_t[1]
