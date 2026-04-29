from flint.flintlib.types.flint cimport flint_rand_t, fmpz_t, nn_ptr, slong, ulong
from flint.flintlib.types.fmpz cimport fmpz_factor_t

# unknown type FILE
# unknown type ecm_t


cdef extern from "flint/fmpz_factor.h":
    void fmpz_factor_init(fmpz_factor_t factor)
    void fmpz_factor_clear(fmpz_factor_t factor)
    void _fmpz_factor_append_ui(fmpz_factor_t factor, ulong p, ulong exp)
    void _fmpz_factor_append(fmpz_factor_t factor, const fmpz_t p, ulong exp)
    void fmpz_factor(fmpz_factor_t factor, const fmpz_t n)
    int fmpz_factor_smooth(fmpz_factor_t factor, const fmpz_t n, slong bits, int proved)
    void fmpz_factor_si(fmpz_factor_t factor, slong n)
    int fmpz_factor_trial_range(fmpz_factor_t factor, const fmpz_t n, ulong start, ulong num_primes)
    int fmpz_factor_trial(fmpz_factor_t factor, const fmpz_t n, slong num_primes)
    void fmpz_factor_refine(fmpz_factor_t res, const fmpz_factor_t f)
    void fmpz_factor_expand_iterative(fmpz_t n, const fmpz_factor_t factor)
    int fmpz_factor_pp1(fmpz_t factor, const fmpz_t n, ulong B1, ulong B2_sqrt, ulong c)
    int fmpz_factor_pollard_brent_single(fmpz_t p_factor, fmpz_t n_in, fmpz_t yi, fmpz_t ai, ulong max_iters)
    int fmpz_factor_pollard_brent(fmpz_t factor, flint_rand_t state, fmpz_t n, ulong max_tries, ulong max_iters)
    # int fmpz_factor_fprint(FILE * fs, const fmpz_factor_t factor)
    int fmpz_factor_print(const fmpz_factor_t factor)
    # void fmpz_factor_ecm_init(ecm_t ecm_inf, ulong sz)
    # void fmpz_factor_ecm_clear(ecm_t ecm_inf)
    # void fmpz_factor_ecm_double(nn_ptr x, nn_ptr z, nn_ptr x0, nn_ptr z0, nn_ptr n, ecm_t ecm_inf)
    # void fmpz_factor_ecm_add(nn_ptr x, nn_ptr z, nn_ptr x1, nn_ptr z1, nn_ptr x2, nn_ptr z2, nn_ptr x0, nn_ptr z0, nn_ptr n, ecm_t ecm_inf)
    # void fmpz_factor_ecm_mul_montgomery_ladder(nn_ptr x, nn_ptr z, nn_ptr x0, nn_ptr z0, ulong k, nn_ptr n, ecm_t ecm_inf)
    # int fmpz_factor_ecm_select_curve(nn_ptr f, nn_ptr sigma, nn_ptr n, ecm_t ecm_inf)
    # int fmpz_factor_ecm_stage_I(nn_ptr f, const ulong * prime_array, ulong num, ulong B1, nn_ptr n, ecm_t ecm_inf)
    # int fmpz_factor_ecm_stage_II(nn_ptr f, ulong B1, ulong B2, ulong P, nn_ptr n, ecm_t ecm_inf)
    int fmpz_factor_ecm(fmpz_t f, ulong curves, ulong B1, ulong B2, flint_rand_t state, const fmpz_t n_in)
