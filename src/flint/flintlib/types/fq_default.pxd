from flint.flintlib.types.flint cimport ulong, fmpz_t
from flint.flintlib.types.nmod cimport nmod_t, nmod_mat_t, nmod_poly_t, nmod_poly_factor_t
from flint.flintlib.types.fmpz_mod cimport fmpz_mod_ctx_t, fmpz_mod_mat_t, fmpz_mod_poly_t, fmpz_mod_poly_factor_t
from flint.flintlib.types.fq_nmod cimport fq_nmod_t, fq_nmod_ctx_t, fq_nmod_mat_t, fq_nmod_poly_t, fq_nmod_poly_factor_t
from flint.flintlib.types.fq_zech cimport fq_zech_t, fq_zech_ctx_t, fq_zech_mat_t, fq_zech_poly_t, fq_zech_poly_factor_t
from flint.flintlib.types.fq cimport fq_t, fq_ctx_t, fq_mat_t, fq_poly_t, fq_poly_factor_t


cdef extern from "flint/fq_default.h":

    ctypedef union fq_default_struct:
        fq_t fq
        fq_nmod_t fq_nmod
        fq_zech_t fq_zech
        ulong nmod
        fmpz_t fmpz_mod
    ctypedef fq_default_struct fq_default_t[1]

    ctypedef struct fq_default_fmpz_mod_ctx_struct:
        fmpz_mod_ctx_t mod
        fmpz_t a       # minpoly is x - a

    ctypedef struct fq_default_nmod_ctx_struct:
        nmod_t mod
        ulong a    # minpoly is x - a

    # This is how it is actually defined now:
    # ctypedef gr_ctx_struct fq_default_ctx_struct;
    # ctypedef fq_default_ctx_struct fq_default_ctx_t[1];

    ctypedef union fq_default_ctx_struct:
        fq_ctx_t fq
        fq_nmod_ctx_t fq_nmod
        fq_zech_ctx_t fq_zech
        fq_default_nmod_ctx_struct nmod
        fq_default_fmpz_mod_ctx_struct fmpz_mod

    ctypedef fq_default_ctx_struct fq_default_ctx_t[1]

    ctypedef union fq_default_mat_struct:
        fq_mat_t fq
        fq_nmod_mat_t fq_nmod
        fq_zech_mat_t fq_zech
        nmod_mat_t nmod
        fmpz_mod_mat_t fmpz_mod
    ctypedef fq_default_mat_struct fq_default_mat_t[1]

    ctypedef union fq_default_poly_struct:
        fq_poly_t fq
        fq_nmod_poly_t fq_nmod
        fq_zech_poly_t fq_zech
        nmod_poly_t nmod
        fmpz_mod_poly_t fmpz_mod

    ctypedef fq_default_poly_struct fq_default_poly_t[1]

    ctypedef union fq_default_poly_factor_struct:
        fq_poly_factor_t fq
        fq_nmod_poly_factor_t fq_nmod
        fq_zech_poly_factor_t fq_zech
        nmod_poly_factor_t nmod
        fmpz_mod_poly_factor_t fmpz_mod

    ctypedef fq_default_poly_factor_struct fq_default_poly_factor_t[1]
