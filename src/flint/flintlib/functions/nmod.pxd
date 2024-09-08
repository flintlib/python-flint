from flint.flintlib.types.flint cimport fmpz_t, mp_limb_t, nmod_t, ulong

# unknown type nmod_discrete_log_pohlig_hellman_t

# .. macro:: NMOD_BITS(mod)
# .. macro:: NMOD_CAN_USE_SHOUP(mod)
# .. macro:: NMOD_RED2(r, a_hi, a_lo, mod)
# .. macro:: NMOD_RED(r, a, mod)
# .. macro:: NMOD2_RED2(r, a_hi, a_lo, mod)
# .. macro:: NMOD_RED3(r, a_hi, a_me, a_lo, mod)
# .. macro:: NMOD_MUL_PRENORM(res, a, b, mod)
# .. macro:: NMOD_MUL_FULLWORD(res, a, b, mod)
# .. macro:: NMOD_ADDMUL(r, a, b, mod)

cdef extern from "flint/nmod.h":
    void nmod_init(nmod_t * mod, mp_limb_t n)
    mp_limb_t _nmod_add(mp_limb_t a, mp_limb_t b, nmod_t mod)
    mp_limb_t nmod_add(mp_limb_t a, mp_limb_t b, nmod_t mod)
    mp_limb_t _nmod_sub(mp_limb_t a, mp_limb_t b, nmod_t mod)
    mp_limb_t nmod_sub(mp_limb_t a, mp_limb_t b, nmod_t mod)
    mp_limb_t nmod_neg(mp_limb_t a, nmod_t mod)
    mp_limb_t nmod_mul(mp_limb_t a, mp_limb_t b, nmod_t mod)
    mp_limb_t _nmod_mul_fullword(mp_limb_t a, mp_limb_t b, nmod_t mod)
    mp_limb_t nmod_inv(mp_limb_t a, nmod_t mod)
    mp_limb_t nmod_div(mp_limb_t a, mp_limb_t b, nmod_t mod)
    int nmod_divides(mp_limb_t * a, mp_limb_t b, mp_limb_t c, nmod_t mod)
    mp_limb_t nmod_pow_ui(mp_limb_t a, ulong e, nmod_t mod)
    mp_limb_t nmod_pow_fmpz(mp_limb_t a, const fmpz_t e, nmod_t mod)
    # void nmod_discrete_log_pohlig_hellman_init(nmod_discrete_log_pohlig_hellman_t L)
    # void nmod_discrete_log_pohlig_hellman_clear(nmod_discrete_log_pohlig_hellman_t L)
    # double nmod_discrete_log_pohlig_hellman_precompute_prime(nmod_discrete_log_pohlig_hellman_t L, mp_limb_t p)
    # mp_limb_t nmod_discrete_log_pohlig_hellman_primitive_root(const nmod_discrete_log_pohlig_hellman_t L)
    # ulong nmod_discrete_log_pohlig_hellman_run(const nmod_discrete_log_pohlig_hellman_t L, mp_limb_t y)
