from flint.flintlib.types.flint cimport flint_rand_t, fmpz_t, mp_limb_t, mp_ptr, slong, ulong
from flint.flintlib.types.fmpz cimport fmpz_factor_t

# unknown type FILE
# unknown type ecm_t


cdef extern from "flint/fmpz_factor.h":
    void fmpz_factor_init(fmpz_factor_t factor)
    void fmpz_factor_clear(fmpz_factor_t factor)
    void _fmpz_factor_append_ui(fmpz_factor_t factor, mp_limb_t p, ulong exp)
    void _fmpz_factor_append(fmpz_factor_t factor, const fmpz_t p, ulong exp)
    void fmpz_factor(fmpz_factor_t factor, const fmpz_t n)
    int fmpz_factor_smooth(fmpz_factor_t factor, const fmpz_t n, slong bits, int proved)
    void fmpz_factor_si(fmpz_factor_t factor, slong n)
    int fmpz_factor_trial_range(fmpz_factor_t factor, const fmpz_t n, ulong start, ulong num_primes)
    int fmpz_factor_trial(fmpz_factor_t factor, const fmpz_t n, slong num_primes)
    void fmpz_factor_refine(fmpz_factor_t res, const fmpz_factor_t f)
    void fmpz_factor_expand_iterative(fmpz_t n, const fmpz_factor_t factor)
    int fmpz_factor_pp1(fmpz_t factor, const fmpz_t n, ulong B1, ulong B2_sqrt, ulong c)
    int fmpz_factor_pollard_brent_single(fmpz_t p_factor, fmpz_t n_in, fmpz_t yi, fmpz_t ai, mp_limb_t max_iters)
    int fmpz_factor_pollard_brent(fmpz_t factor, flint_rand_t state, fmpz_t n, mp_limb_t max_tries, mp_limb_t max_iters)
    # int fmpz_factor_fprint(FILE * fs, const fmpz_factor_t factor)
    int fmpz_factor_print(const fmpz_factor_t factor)
    # void fmpz_factor_ecm_init(ecm_t ecm_inf, mp_limb_t sz)
    # void fmpz_factor_ecm_clear(ecm_t ecm_inf)
    void fmpz_factor_ecm_addmod(mp_ptr a, mp_ptr b, mp_ptr c, mp_ptr n, mp_limb_t n_size)
    void fmpz_factor_ecm_submod(mp_ptr x, mp_ptr a, mp_ptr b, mp_ptr n, mp_limb_t n_size)
    # void fmpz_factor_ecm_double(mp_ptr x, mp_ptr z, mp_ptr x0, mp_ptr z0, mp_ptr n, ecm_t ecm_inf)
    # void fmpz_factor_ecm_add(mp_ptr x, mp_ptr z, mp_ptr x1, mp_ptr z1, mp_ptr x2, mp_ptr z2, mp_ptr x0, mp_ptr z0, mp_ptr n, ecm_t ecm_inf)
    # void fmpz_factor_ecm_mul_montgomery_ladder(mp_ptr x, mp_ptr z, mp_ptr x0, mp_ptr z0, mp_limb_t k, mp_ptr n, ecm_t ecm_inf)
    # int fmpz_factor_ecm_select_curve(mp_ptr f, mp_ptr sigma, mp_ptr n, ecm_t ecm_inf)
    # int fmpz_factor_ecm_stage_I(mp_ptr f, const mp_limb_t * prime_array, mp_limb_t num, mp_limb_t B1, mp_ptr n, ecm_t ecm_inf)
    # int fmpz_factor_ecm_stage_II(mp_ptr f, mp_limb_t B1, mp_limb_t B2, mp_limb_t P, mp_ptr n, ecm_t ecm_inf)
    int fmpz_factor_ecm(fmpz_t f, mp_limb_t curves, mp_limb_t B1, mp_limb_t B2, flint_rand_t state, const fmpz_t n_in)
