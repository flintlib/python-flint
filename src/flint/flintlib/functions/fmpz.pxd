from flint.flintlib.types.flint cimport flint_bitcnt_t, flint_rand_t, fmpz_struct, fmpz_t, nmod_t, nn_ptr, nn_srcptr, slong, ulong
from flint.flintlib.types.fmpz cimport fmpz_factor_t, fmpz_preinvn_t

# unknown type FILE
# unknown type fmpz_comb_t
# unknown type fmpz_comb_temp_t
# unknown type fmpz_multi_CRT_t
# unknown type mpf_t
# unknown type mpfr_rnd_t
# unknown type mpfr_t
# unknown type mpz_ptr
# unknown type mpz_t
# unknown type size_t

# .. macro:: COEFF_MAX
# .. macro:: COEFF_MIN
# .. macro:: COEFF_IS_MPZ(f)
# .. macro:: MPZ_MIN_ALLOC

cdef extern from "flint/fmpz.h":
    # fmpz_struct PTR_TO_COEFF(mpz_ptr ptr)
    # mpz_ptr COEFF_TO_PTR(fmpz_struct f)
    # mpz_ptr _fmpz_new_mpz(void)
    void _fmpz_clear_mpz(fmpz_struct f)
    void _fmpz_cleanup_mpz_content()
    void _fmpz_cleanup()
    # mpz_ptr _fmpz_promote(fmpz_t f)
    # mpz_ptr _fmpz_promote_val(fmpz_t f)
    void _fmpz_demote(fmpz_t f)
    void _fmpz_demote_val(fmpz_t f)
    int _fmpz_is_canonical(const fmpz_t f)
    void fmpz_init(fmpz_t f)
    void fmpz_init2(fmpz_t f, ulong limbs)
    void fmpz_clear(fmpz_t f)
    void fmpz_init_set(fmpz_t f, const fmpz_t g)
    void fmpz_init_set_ui(fmpz_t f, ulong g)
    void fmpz_init_set_si(fmpz_t f, slong g)
    void fmpz_randbits_unsigned(fmpz_t f, flint_rand_t state, flint_bitcnt_t bits)
    void fmpz_randbits(fmpz_t f, flint_rand_t state, flint_bitcnt_t bits)
    void fmpz_randtest_unsigned(fmpz_t f, flint_rand_t state, flint_bitcnt_t bits)
    void fmpz_randtest(fmpz_t f, flint_rand_t state, flint_bitcnt_t bits)
    void fmpz_randtest_not_zero(fmpz_t f, flint_rand_t state, flint_bitcnt_t bits)
    void fmpz_randm(fmpz_t f, flint_rand_t state, const fmpz_t m)
    void fmpz_randtest_mod(fmpz_t f, flint_rand_t state, const fmpz_t m)
    void fmpz_randtest_mod_signed(fmpz_t f, flint_rand_t state, const fmpz_t m)
    void fmpz_randprime(fmpz_t f, flint_rand_t state, flint_bitcnt_t bits, int proved)
    slong fmpz_get_si(const fmpz_t f)
    ulong fmpz_get_ui(const fmpz_t f)
    void fmpz_get_uiui(ulong * hi, ulong * low, const fmpz_t f)
    ulong fmpz_get_nmod(const fmpz_t f, nmod_t mod)
    double fmpz_get_d(const fmpz_t f)
    # void fmpz_set_mpf(fmpz_t f, const mpf_t x)
    # void fmpz_get_mpf(mpf_t x, const fmpz_t f)
    # void fmpz_get_mpfr(mpfr_t x, const fmpz_t f, mpfr_rnd_t rnd)
    double fmpz_get_d_2exp(slong * exp, const fmpz_t f)
    # void fmpz_get_mpz(mpz_t x, const fmpz_t f)
    int fmpz_get_mpn(nn_ptr * n, fmpz_t n_in)
    char * fmpz_get_str(char * str, int b, const fmpz_t f)
    void fmpz_set_si(fmpz_t f, slong val)
    void fmpz_set_ui(fmpz_t f, ulong val)
    void fmpz_set_d(fmpz_t f, double c)
    void fmpz_set_d_2exp(fmpz_t f, double d, slong exp)
    void fmpz_neg_ui(fmpz_t f, ulong val)
    void fmpz_set_uiui(fmpz_t f, ulong hi, ulong lo)
    void fmpz_neg_uiui(fmpz_t f, ulong hi, ulong lo)
    void fmpz_set_signed_uiui(fmpz_t f, ulong hi, ulong lo)
    void fmpz_set_signed_uiuiui(fmpz_t f, ulong hi, ulong mid, ulong lo)
    void fmpz_set_ui_array(fmpz_t out, const ulong * in_, slong n)
    void fmpz_set_signed_ui_array(fmpz_t out, const ulong * in_, slong n)
    void fmpz_get_ui_array(ulong * out, slong n, const fmpz_t in_)
    void fmpz_get_signed_ui_array(ulong * out, slong n, const fmpz_t in_)
    void fmpz_set_mpn_large(fmpz_t z, nn_srcptr src, slong n, int negative)
    void fmpz_get_signed_uiui(ulong * hi, ulong * lo, const fmpz_t in_)
    # void fmpz_set_mpz(fmpz_t f, const mpz_t x)
    int fmpz_set_str(fmpz_t f, const char * str, int b)
    void fmpz_set_ui_smod(fmpz_t f, ulong x, ulong m)
    # void flint_mpz_init_set_readonly(mpz_t z, const fmpz_t f)
    # void flint_mpz_clear_readonly(mpz_t z)
    # void fmpz_init_set_readonly(fmpz_t f, const mpz_t z)
    void fmpz_clear_readonly(fmpz_t f)
    int fmpz_read(fmpz_t f)
    # int fmpz_fread(FILE * file, fmpz_t f)
    # size_t fmpz_inp_raw(fmpz_t x, FILE * fin)
    # int fmpz_fprint(FILE * fs, const fmpz_t x)
    int fmpz_print(const fmpz_t x)
    # size_t fmpz_out_raw(FILE * fout, const fmpz_t x )
    # size_t fmpz_sizeinbase(const fmpz_t f, int b)
    flint_bitcnt_t fmpz_bits(const fmpz_t f)
    slong fmpz_size(const fmpz_t f)
    int fmpz_sgn(const fmpz_t f)
    flint_bitcnt_t fmpz_val2(const fmpz_t f)
    void fmpz_swap(fmpz_t f, fmpz_t g)
    void fmpz_set(fmpz_t f, const fmpz_t g)
    void fmpz_zero(fmpz_t f)
    void fmpz_one(fmpz_t f)
    int fmpz_abs_fits_ui(const fmpz_t f)
    int fmpz_fits_si(const fmpz_t f)
    void fmpz_setbit(fmpz_t f, ulong i)
    int fmpz_tstbit(const fmpz_t f, ulong i)
    ulong fmpz_abs_lbound_ui_2exp(slong * exp, const fmpz_t x, int bits)
    ulong fmpz_abs_ubound_ui_2exp(slong * exp, const fmpz_t x, int bits)
    int fmpz_cmp(const fmpz_t f, const fmpz_t g)
    int fmpz_cmp_ui(const fmpz_t f, ulong g)
    int fmpz_cmp_si(const fmpz_t f, slong g)
    int fmpz_cmpabs(const fmpz_t f, const fmpz_t g)
    int fmpz_cmp2abs(const fmpz_t f, const fmpz_t g)
    int fmpz_equal(const fmpz_t f, const fmpz_t g)
    int fmpz_equal_ui(const fmpz_t f, ulong g)
    int fmpz_equal_si(const fmpz_t f, slong g)
    int fmpz_is_zero(const fmpz_t f)
    int fmpz_is_one(const fmpz_t f)
    int fmpz_is_pm1(const fmpz_t f)
    int fmpz_is_even(const fmpz_t f)
    int fmpz_is_odd(const fmpz_t f)
    void fmpz_neg(fmpz_t f1, const fmpz_t f2)
    void fmpz_abs(fmpz_t f1, const fmpz_t f2)
    void fmpz_add(fmpz_t f, const fmpz_t g, const fmpz_t h)
    void fmpz_add_ui(fmpz_t f, const fmpz_t g, ulong h)
    void fmpz_add_si(fmpz_t f, const fmpz_t g, slong h)
    void fmpz_sub(fmpz_t f, const fmpz_t g, const fmpz_t h)
    void fmpz_sub_ui(fmpz_t f, const fmpz_t g, ulong h)
    void fmpz_sub_si(fmpz_t f, const fmpz_t g, slong h)
    void fmpz_mul(fmpz_t f, const fmpz_t g, const fmpz_t h)
    void fmpz_mul_ui(fmpz_t f, const fmpz_t g, ulong h)
    void fmpz_mul_si(fmpz_t f, const fmpz_t g, slong h)
    void fmpz_mul2_uiui(fmpz_t f, const fmpz_t g, ulong x, ulong y)
    void fmpz_mul_2exp(fmpz_t f, const fmpz_t g, ulong e)
    void fmpz_one_2exp(fmpz_t f, ulong e)
    void fmpz_addmul(fmpz_t f, const fmpz_t g, const fmpz_t h)
    void fmpz_addmul_ui(fmpz_t f, const fmpz_t g, ulong h)
    void fmpz_addmul_si(fmpz_t f, const fmpz_t g, slong h)
    void fmpz_submul(fmpz_t f, const fmpz_t g, const fmpz_t h)
    void fmpz_submul_ui(fmpz_t f, const fmpz_t g, ulong h)
    void fmpz_submul_si(fmpz_t f, const fmpz_t g, slong h)
    void fmpz_fmma(fmpz_t f, const fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_t d)
    void fmpz_fmms(fmpz_t f, const fmpz_t a, const fmpz_t b, const fmpz_t c, const fmpz_t d)
    void fmpz_cdiv_qr(fmpz_t f, fmpz_t s, const fmpz_t g, const fmpz_t h)
    void fmpz_fdiv_qr(fmpz_t f, fmpz_t s, const fmpz_t g, const fmpz_t h)
    void fmpz_tdiv_qr(fmpz_t f, fmpz_t s, const fmpz_t g, const fmpz_t h)
    void fmpz_ndiv_qr(fmpz_t f, fmpz_t s, const fmpz_t g, const fmpz_t h)
    void fmpz_cdiv_q(fmpz_t f, const fmpz_t g, const fmpz_t h)
    void fmpz_fdiv_q(fmpz_t f, const fmpz_t g, const fmpz_t h)
    void fmpz_tdiv_q(fmpz_t f, const fmpz_t g, const fmpz_t h)
    void fmpz_cdiv_q_si(fmpz_t f, const fmpz_t g, slong h)
    void fmpz_fdiv_q_si(fmpz_t f, const fmpz_t g, slong h)
    void fmpz_tdiv_q_si(fmpz_t f, const fmpz_t g, slong h)
    void fmpz_cdiv_q_ui(fmpz_t f, const fmpz_t g, ulong h)
    void fmpz_fdiv_q_ui(fmpz_t f, const fmpz_t g, ulong h)
    void fmpz_tdiv_q_ui(fmpz_t f, const fmpz_t g, ulong h)
    void fmpz_cdiv_q_2exp(fmpz_t f, const fmpz_t g, ulong exp)
    void fmpz_fdiv_q_2exp(fmpz_t f, const fmpz_t g, ulong exp)
    void fmpz_tdiv_q_2exp(fmpz_t f, const fmpz_t g, ulong exp)
    void fmpz_fdiv_r(fmpz_t s, const fmpz_t g, const fmpz_t h)
    void fmpz_cdiv_r_2exp(fmpz_t s, const fmpz_t g, ulong exp)
    void fmpz_fdiv_r_2exp(fmpz_t s, const fmpz_t g, ulong exp)
    void fmpz_tdiv_r_2exp(fmpz_t s, const fmpz_t g, ulong exp)
    ulong fmpz_cdiv_ui(const fmpz_t g, ulong h)
    ulong fmpz_fdiv_ui(const fmpz_t g, ulong h)
    ulong fmpz_tdiv_ui(const fmpz_t g, ulong h)
    void fmpz_divexact(fmpz_t f, const fmpz_t g, const fmpz_t h)
    void fmpz_divexact_si(fmpz_t f, const fmpz_t g, slong h)
    void fmpz_divexact_ui(fmpz_t f, const fmpz_t g, ulong h)
    void fmpz_divexact2_uiui(fmpz_t f, const fmpz_t g, ulong x, ulong y)
    int fmpz_divisible(const fmpz_t f, const fmpz_t g)
    int fmpz_divisible_si(const fmpz_t f, slong g)
    int fmpz_divides(fmpz_t q, const fmpz_t f, const fmpz_t g)
    void fmpz_mod(fmpz_t f, const fmpz_t g, const fmpz_t h)
    ulong fmpz_mod_ui(fmpz_t f, const fmpz_t g, ulong h)
    void fmpz_smod(fmpz_t f, const fmpz_t g, const fmpz_t h)
    void fmpz_preinvn_init(fmpz_preinvn_t inv, const fmpz_t f)
    void fmpz_preinvn_clear(fmpz_preinvn_t inv)
    void fmpz_fdiv_qr_preinvn(fmpz_t f, fmpz_t s, const fmpz_t g, const fmpz_t h, const fmpz_preinvn_t hinv)
    void fmpz_pow_ui(fmpz_t f, const fmpz_t g, ulong x)
    void fmpz_ui_pow_ui(fmpz_t f, ulong g, ulong x)
    int fmpz_pow_fmpz(fmpz_t f, const fmpz_t g, const fmpz_t x)
    void fmpz_powm_ui(fmpz_t f, const fmpz_t g, ulong e, const fmpz_t m)
    void fmpz_powm(fmpz_t f, const fmpz_t g, const fmpz_t e, const fmpz_t m)
    slong fmpz_clog(const fmpz_t x, const fmpz_t b)
    slong fmpz_clog_ui(const fmpz_t x, ulong b)
    slong fmpz_flog(const fmpz_t x, const fmpz_t b)
    slong fmpz_flog_ui(const fmpz_t x, ulong b)
    double fmpz_dlog(const fmpz_t x)
    int fmpz_sqrtmod(fmpz_t b, const fmpz_t a, const fmpz_t p)
    void fmpz_sqrt(fmpz_t f, const fmpz_t g)
    void fmpz_sqrtrem(fmpz_t f, fmpz_t r, const fmpz_t g)
    int fmpz_is_square(const fmpz_t f)
    int fmpz_root(fmpz_t r, const fmpz_t f, slong n)
    int fmpz_is_perfect_power(fmpz_t root, const fmpz_t f)
    void fmpz_fac_ui(fmpz_t f, ulong n)
    void fmpz_fib_ui(fmpz_t f, ulong n)
    void fmpz_bin_uiui(fmpz_t f, ulong n, ulong k)
    void _fmpz_rfac_ui(fmpz_t r, const fmpz_t x, ulong a, ulong b)
    void fmpz_rfac_ui(fmpz_t r, const fmpz_t x, ulong k)
    void fmpz_rfac_uiui(fmpz_t r, ulong x, ulong k)
    void fmpz_mul_tdiv_q_2exp(fmpz_t f, const fmpz_t g, const fmpz_t h, ulong exp)
    void fmpz_mul_si_tdiv_q_2exp(fmpz_t f, const fmpz_t g, slong x, ulong exp)
    void fmpz_gcd_ui(fmpz_t f, const fmpz_t g, ulong h)
    void fmpz_gcd(fmpz_t f, const fmpz_t g, const fmpz_t h)
    void fmpz_gcd3(fmpz_t f, const fmpz_t a, const fmpz_t b, const fmpz_t c)
    void fmpz_lcm(fmpz_t f, const fmpz_t g, const fmpz_t h)
    void fmpz_gcdinv(fmpz_t d, fmpz_t a, const fmpz_t f, const fmpz_t g)
    void fmpz_xgcd(fmpz_t d, fmpz_t a, fmpz_t b, const fmpz_t f, const fmpz_t g)
    void fmpz_xgcd_canonical_bezout(fmpz_t d, fmpz_t a, fmpz_t b, const fmpz_t f, const fmpz_t g)
    void fmpz_xgcd_partial(fmpz_t co2, fmpz_t co1, fmpz_t r2, fmpz_t r1, const fmpz_t L)
    slong _fmpz_remove(fmpz_t x, const fmpz_t f, double finv)
    slong fmpz_remove(fmpz_t rop, const fmpz_t op, const fmpz_t f)
    int fmpz_invmod(fmpz_t f, const fmpz_t g, const fmpz_t h)
    void fmpz_negmod(fmpz_t f, const fmpz_t g, const fmpz_t h)
    int fmpz_jacobi(const fmpz_t a, const fmpz_t n)
    int fmpz_kronecker(const fmpz_t a, const fmpz_t n)
    void fmpz_divides_mod_list(fmpz_t xstart, fmpz_t xstride, fmpz_t xlength, const fmpz_t a, const fmpz_t b, const fmpz_t n)
    int fmpz_bit_pack(ulong * arr, flint_bitcnt_t shift, flint_bitcnt_t bits, const fmpz_t coeff, int negate, int borrow)
    int fmpz_bit_unpack(fmpz_t coeff, ulong * arr, flint_bitcnt_t shift, flint_bitcnt_t bits, int negate, int borrow)
    void fmpz_bit_unpack_unsigned(fmpz_t coeff, const ulong * arr, flint_bitcnt_t shift, flint_bitcnt_t bits)
    void fmpz_complement(fmpz_t r, const fmpz_t f)
    void fmpz_clrbit(fmpz_t f, ulong i)
    void fmpz_combit(fmpz_t f, ulong i)
    void fmpz_and(fmpz_t r, const fmpz_t a, const fmpz_t b)
    void fmpz_or(fmpz_t r, const fmpz_t a, const fmpz_t b)
    void fmpz_xor(fmpz_t r, const fmpz_t a, const fmpz_t b)
    ulong fmpz_popcnt(const fmpz_t a)
    void fmpz_CRT_ui(fmpz_t out, const fmpz_t r1, const fmpz_t m1, ulong r2, ulong m2, int sign)
    void fmpz_CRT(fmpz_t out, const fmpz_t r1, const fmpz_t m1, const fmpz_t r2, const fmpz_t m2, int sign)
    # void fmpz_multi_mod_ui(ulong * out, const fmpz_t in_, const fmpz_comb_t comb, fmpz_comb_temp_t temp)
    # void fmpz_multi_CRT_ui(fmpz_t output, nn_srcptr residues, const fmpz_comb_t comb, fmpz_comb_temp_t ctemp, int sign)
    # void fmpz_comb_init(fmpz_comb_t comb, nn_srcptr primes, slong num_primes)
    # void fmpz_comb_temp_init(fmpz_comb_temp_t temp, const fmpz_comb_t comb)
    # void fmpz_comb_clear(fmpz_comb_t comb)
    # void fmpz_comb_temp_clear(fmpz_comb_temp_t temp)
    # void fmpz_multi_CRT_init(fmpz_multi_CRT_t CRT)
    # int fmpz_multi_CRT_precompute(fmpz_multi_CRT_t CRT, const fmpz_struct * moduli, slong len)
    # void fmpz_multi_CRT_precomp(fmpz_t output, const fmpz_multi_CRT_t P, const fmpz_struct * inputs, int sign)
    int fmpz_multi_CRT(fmpz_t output, const fmpz_struct * moduli, const fmpz_struct * values, slong len, int sign)
    # void fmpz_multi_CRT_clear(fmpz_multi_CRT_t P)
    int fmpz_is_strong_probabprime(const fmpz_t n, const fmpz_t a)
    int fmpz_is_probabprime_lucas(const fmpz_t n)
    int fmpz_is_probabprime_BPSW(const fmpz_t n)
    int fmpz_is_probabprime(const fmpz_t p)
    int fmpz_is_prime_pseudosquare(const fmpz_t n)
    int fmpz_is_prime_pocklington(fmpz_t F, fmpz_t R, const fmpz_t n, nn_ptr pm1, slong num_pm1)
    void _fmpz_nm1_trial_factors(const fmpz_t n, nn_ptr pm1, slong * num_pm1, ulong limit)
    int fmpz_is_prime_morrison(fmpz_t F, fmpz_t R, const fmpz_t n, nn_ptr pp1, slong num_pp1)
    void _fmpz_np1_trial_factors(const fmpz_t n, nn_ptr pp1, slong * num_pp1, ulong limit)
    int fmpz_is_prime(const fmpz_t n)
    void fmpz_lucas_chain(fmpz_t Vm, fmpz_t Vm1, const fmpz_t A, const fmpz_t m, const fmpz_t n)
    void fmpz_lucas_chain_full(fmpz_t Vm, fmpz_t Vm1, const fmpz_t A, const fmpz_t B, const fmpz_t m, const fmpz_t n)
    void fmpz_lucas_chain_double(fmpz_t U2m, fmpz_t U2m1, const fmpz_t Um, const fmpz_t Um1, const fmpz_t A, const fmpz_t B, const fmpz_t n)
    void fmpz_lucas_chain_add(fmpz_t Umn, fmpz_t Umn1, const fmpz_t Um, const fmpz_t Um1, const fmpz_t Un, const fmpz_t Un1, const fmpz_t A, const fmpz_t B, const fmpz_t n)
    void fmpz_lucas_chain_mul(fmpz_t Ukm, fmpz_t Ukm1, const fmpz_t Um, const fmpz_t Um1, const fmpz_t A, const fmpz_t B, const fmpz_t k, const fmpz_t n)
    void fmpz_lucas_chain_VtoU(fmpz_t Um, fmpz_t Um1, const fmpz_t Vm, const fmpz_t Vm1, const fmpz_t A, const fmpz_t B, const fmpz_t Dinv, const fmpz_t n)
    int fmpz_divisor_in_residue_class_lenstra(fmpz_t fac, const fmpz_t n, const fmpz_t r, const fmpz_t s)
    void fmpz_nextprime(fmpz_t res, const fmpz_t n, int proved)
    void fmpz_primorial(fmpz_t res, ulong n)
    void fmpz_factor_euler_phi(fmpz_t res, const fmpz_factor_t fac)
    void fmpz_euler_phi(fmpz_t res, const fmpz_t n)
    int fmpz_factor_moebius_mu(const fmpz_factor_t fac)
    int fmpz_moebius_mu(const fmpz_t n)
    void fmpz_factor_divisor_sigma(fmpz_t res, ulong k, const fmpz_factor_t fac)
    void fmpz_divisor_sigma(fmpz_t res, ulong k, const fmpz_t n)
