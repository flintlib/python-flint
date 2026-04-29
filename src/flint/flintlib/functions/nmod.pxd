from flint.flintlib.types.flint cimport fmpz_t, nmod_t, ulong

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
    void nmod_init(nmod_t * mod, ulong n)
    ulong _nmod_add(ulong a, ulong b, nmod_t mod)
    ulong nmod_add(ulong a, ulong b, nmod_t mod)
    ulong _nmod_sub(ulong a, ulong b, nmod_t mod)
    ulong nmod_sub(ulong a, ulong b, nmod_t mod)
    ulong nmod_neg(ulong a, nmod_t mod)
    ulong nmod_mul(ulong a, ulong b, nmod_t mod)
    ulong _nmod_mul_fullword(ulong a, ulong b, nmod_t mod)
    ulong nmod_inv(ulong a, nmod_t mod)
    ulong nmod_div(ulong a, ulong b, nmod_t mod)
    int nmod_divides(ulong * a, ulong b, ulong c, nmod_t mod)
    ulong nmod_pow_ui(ulong a, ulong e, nmod_t mod)
    ulong nmod_pow_fmpz(ulong a, const fmpz_t e, nmod_t mod)
    # void nmod_discrete_log_pohlig_hellman_init(nmod_discrete_log_pohlig_hellman_t L)
    # void nmod_discrete_log_pohlig_hellman_clear(nmod_discrete_log_pohlig_hellman_t L)
    # double nmod_discrete_log_pohlig_hellman_precompute_prime(nmod_discrete_log_pohlig_hellman_t L, ulong p)
    # ulong nmod_discrete_log_pohlig_hellman_primitive_root(const nmod_discrete_log_pohlig_hellman_t L)
    # ulong nmod_discrete_log_pohlig_hellman_run(const nmod_discrete_log_pohlig_hellman_t L, ulong y)
