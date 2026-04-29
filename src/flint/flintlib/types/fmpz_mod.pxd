from flint.flintlib.types.flint cimport ulong, slong, fmpz_struct, fmpz_t, nmod_t, flint_bitcnt_t
from flint.flintlib.types.fmpz cimport fmpz_struct, fmpz_preinvn_struct, fmpz_mat_struct
from flint.flintlib.types.mpoly cimport mpoly_ctx_t


cdef extern from "flint/fmpz_mod_types.h":

    ctypedef struct fmpz_mod_ctx_struct:
        fmpz_t n
        void (* add_fxn)(fmpz_t, const fmpz_t, const fmpz_t, const fmpz_mod_ctx_struct*)
        void (* sub_fxn)(fmpz_t, const fmpz_t, const fmpz_t, const fmpz_mod_ctx_struct*)
        void (* mul_fxn)(fmpz_t, const fmpz_t, const fmpz_t, const fmpz_mod_ctx_struct*)
        nmod_t mod
        ulong n_limbs[3]
        ulong ninv_limbs[3]
        fmpz_preinvn_struct * ninv_huge

    ctypedef fmpz_mod_ctx_struct fmpz_mod_ctx_t[1]

    ctypedef fmpz_mat_struct fmpz_mod_mat_struct

    ctypedef fmpz_mod_mat_struct fmpz_mod_mat_t[1]

    ctypedef struct fmpz_mod_poly_struct:
        fmpz_struct * coeffs
        slong alloc
        slong length

    ctypedef fmpz_mod_poly_struct fmpz_mod_poly_t[1]

    ctypedef struct fmpz_mod_poly_factor_struct:
        fmpz_mod_poly_struct * poly
        slong *exp
        slong num
        slong alloc

    ctypedef fmpz_mod_poly_factor_struct fmpz_mod_poly_factor_t[1]

    ctypedef struct fmpz_mod_mpoly_ctx_struct:
        mpoly_ctx_t minfo
        fmpz_mod_ctx_t ffinfo

    ctypedef fmpz_mod_mpoly_ctx_struct fmpz_mod_mpoly_ctx_t[1]

    ctypedef struct fmpz_mod_mpoly_struct:
        fmpz_struct * coeffs
        ulong * exps
        slong length
        flint_bitcnt_t bits
        slong coeffs_alloc
        slong exps_alloc

    ctypedef fmpz_mod_mpoly_struct fmpz_mod_mpoly_t[1]

    ctypedef struct fmpz_mod_mpoly_factor_struct:
        fmpz_t constant
        fmpz_mod_mpoly_struct * poly
        fmpz_struct * exp
        slong num
        slong alloc

    ctypedef fmpz_mod_mpoly_factor_struct fmpz_mod_mpoly_factor_t[1]


cdef extern from "flint/fmpz_mod.h":

    ctypedef struct fmpz_mod_discrete_log_pohlig_hellman_table_entry_struct:
        fmpz_t gammapow
        ulong cm

    ctypedef struct fmpz_mod_discrete_log_pohlig_hellman_entry_struct:
        slong exp
        ulong prime
        fmpz_t gamma
        fmpz_t gammainv
        fmpz_t startingbeta
        fmpz_t co
        fmpz_t startinge
        fmpz_t idem
        ulong cbound
        ulong dbound
        fmpz_mod_discrete_log_pohlig_hellman_table_entry_struct * table  # length cbound */

    ctypedef struct fmpz_mod_discrete_log_pohlig_hellman_struct:
        fmpz_mod_ctx_t fpctx
        fmpz_t pm1         # p - 1 */
        fmpz_t alpha       # p.r. of p */
        fmpz_t alphainv
        slong num_factors  # factors of p - 1
        fmpz_mod_discrete_log_pohlig_hellman_entry_struct * entries

    ctypedef fmpz_mod_discrete_log_pohlig_hellman_struct fmpz_mod_discrete_log_pohlig_hellman_t[1]
