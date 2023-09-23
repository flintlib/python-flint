from flint.flintlib.flint cimport ulong, slong
from flint.flintlib.fmpz cimport fmpz_t, fmpz_struct, fmpz_preinvn_struct
from flint.flintlib.nmod cimport nmod_t

cdef extern from "flint/fmpz_mod.h":
    #
    # fmpz_mod structs, a la Pohlig - Hellman
    #
    ctypedef struct fmpz_mod_ctx_struct:
        fmpz_t n
        nmod_t mod
        ulong n_limbs[3]
        ulong ninv_limbs[3]
        fmpz_preinvn_struct * ninv_huge
    ctypedef fmpz_mod_ctx_struct fmpz_mod_ctx_t[1]

    #
    # discrete logs structs, a la Pohlig - Hellman
    #

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
        fmpz_mod_discrete_log_pohlig_hellman_table_entry_struct * table # length cbound */

    ctypedef struct fmpz_mod_discrete_log_pohlig_hellman_struct:
        fmpz_mod_ctx_t fpctx
        fmpz_t pm1         # p - 1 */
        fmpz_t alpha       # p.r. of p */
        fmpz_t alphainv
        slong num_factors  # factors of p - 1
        fmpz_mod_discrete_log_pohlig_hellman_entry_struct * entries
    ctypedef fmpz_mod_discrete_log_pohlig_hellman_struct fmpz_mod_discrete_log_pohlig_hellman_t[1]

    # Parsed from here
    void fmpz_mod_ctx_init(fmpz_mod_ctx_t ctx, const fmpz_t n)
    void fmpz_mod_ctx_clear(fmpz_mod_ctx_t ctx)
    void fmpz_mod_ctx_set_modulus(fmpz_mod_ctx_t ctx, const fmpz_t n)
    void fmpz_mod_set_fmpz(fmpz_t a, const fmpz_t b, const fmpz_mod_ctx_t ctx)
    int fmpz_mod_is_canonical(const fmpz_t a, const fmpz_mod_ctx_t ctx)
    int fmpz_mod_is_one(const fmpz_t a, const fmpz_mod_ctx_t ctx)
    void fmpz_mod_add(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx)
    void fmpz_mod_add_fmpz(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx)
    void fmpz_mod_add_ui(fmpz_t a, const fmpz_t b, ulong c, const fmpz_mod_ctx_t ctx)
    void fmpz_mod_add_si(fmpz_t a, const fmpz_t b, slong c, const fmpz_mod_ctx_t ctx)
    void fmpz_mod_sub(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx)
    void fmpz_mod_sub_fmpz(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx)
    void fmpz_mod_sub_ui(fmpz_t a, const fmpz_t b, ulong c, const fmpz_mod_ctx_t ctx)
    void fmpz_mod_sub_si(fmpz_t a, const fmpz_t b, slong c, const fmpz_mod_ctx_t ctx)
    void fmpz_mod_fmpz_sub(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx)
    void fmpz_mod_ui_sub(fmpz_t a, ulong b, const fmpz_t c, const fmpz_mod_ctx_t ctx)
    void fmpz_mod_si_sub(fmpz_t a, slong b, const fmpz_t c, const fmpz_mod_ctx_t ctx)
    void fmpz_mod_neg(fmpz_t a, const fmpz_t b, const fmpz_mod_ctx_t ctx)
    void fmpz_mod_mul(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx)
    void fmpz_mod_inv(fmpz_t a, const fmpz_t b, const fmpz_mod_ctx_t ctx)
    int fmpz_mod_divides(fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_mod_ctx_t ctx)
    void fmpz_mod_pow_ui(fmpz_t a, const fmpz_t b, ulong e, const fmpz_mod_ctx_t ctx)
    int fmpz_mod_pow_fmpz(fmpz_t a, const fmpz_t b, const fmpz_t e, const fmpz_mod_ctx_t ctx)
    void fmpz_mod_discrete_log_pohlig_hellman_init(fmpz_mod_discrete_log_pohlig_hellman_t L)
    void fmpz_mod_discrete_log_pohlig_hellman_clear(fmpz_mod_discrete_log_pohlig_hellman_t L)
    double fmpz_mod_discrete_log_pohlig_hellman_precompute_prime(fmpz_mod_discrete_log_pohlig_hellman_t L, const fmpz_t p)
    const fmpz_struct * fmpz_mod_discrete_log_pohlig_hellman_primitive_root(const fmpz_mod_discrete_log_pohlig_hellman_t L)
    void fmpz_mod_discrete_log_pohlig_hellman_run(fmpz_t x, const fmpz_mod_discrete_log_pohlig_hellman_t L, const fmpz_t y)
    int fmpz_next_smooth_prime(fmpz_t a, const fmpz_t b)
