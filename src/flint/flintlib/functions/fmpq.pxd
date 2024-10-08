from flint.flintlib.types.flint cimport flint_bitcnt_t, flint_rand_t, fmpz_struct, fmpz_t, slong, ulong
from flint.flintlib.types.fmpq cimport fmpq_t

# unknown type FILE
# unknown type mpfr_rnd_t
# unknown type mpfr_t
# unknown type mpq_t
# unknown type mpz_t


cdef extern from "flint/fmpq.h":
    fmpz_struct * fmpq_numref(const fmpq_t x)
    fmpz_struct * fmpq_denref(const fmpq_t x)
    void fmpq_init(fmpq_t x)
    void fmpq_clear(fmpq_t x)
    void fmpq_canonicalise(fmpq_t res)
    void _fmpq_canonicalise(fmpz_t num, fmpz_t den)
    int fmpq_is_canonical(const fmpq_t x)
    int _fmpq_is_canonical(const fmpz_t num, const fmpz_t den)
    void fmpq_set(fmpq_t dest, const fmpq_t src)
    void fmpq_swap(fmpq_t op1, fmpq_t op2)
    void fmpq_neg(fmpq_t dest, const fmpq_t src)
    void fmpq_abs(fmpq_t dest, const fmpq_t src)
    void fmpq_zero(fmpq_t res)
    void fmpq_one(fmpq_t res)
    int fmpq_is_zero(const fmpq_t res)
    int fmpq_is_one(const fmpq_t res)
    int fmpq_is_pm1(const fmpq_t res)
    int fmpq_equal(const fmpq_t x, const fmpq_t y)
    int fmpq_equal_fmpz(const fmpq_t x, const fmpz_t y)
    int fmpq_equal_si(fmpq_t x, slong y)
    int fmpq_equal_ui(fmpq_t x, ulong y)
    int fmpq_sgn(const fmpq_t x)
    int fmpq_cmp(const fmpq_t x, const fmpq_t y)
    int fmpq_cmp_fmpz(const fmpq_t x, const fmpz_t y)
    int fmpq_cmp_si(const fmpq_t x, slong y)
    int fmpq_cmp_ui(const fmpq_t x, ulong y)
    void fmpq_height(fmpz_t height, const fmpq_t x)
    flint_bitcnt_t fmpq_height_bits(const fmpq_t x)
    void fmpq_set_fmpz_frac(fmpq_t res, const fmpz_t p, const fmpz_t q)
    # void fmpq_get_mpz_frac(mpz_t a, mpz_t b, fmpq_t c)
    void fmpq_set_si(fmpq_t res, slong p, ulong q)
    void _fmpq_set_si(fmpz_t rnum, fmpz_t rden, slong p, ulong q)
    void fmpq_set_ui(fmpq_t res, ulong p, ulong q)
    void _fmpq_set_ui(fmpz_t rnum, fmpz_t rden, ulong p, ulong q)
    # void fmpq_set_mpq(fmpq_t dest, const mpq_t src)
    int fmpq_set_str(fmpq_t dest, const char * s, int base)
    double fmpq_get_d(const fmpq_t f)
    # void fmpq_get_mpq(mpq_t dest, const fmpq_t src)
    # int fmpq_get_mpfr(mpfr_t dest, const fmpq_t src, mpfr_rnd_t rnd)
    char * _fmpq_get_str(char * str, int b, const fmpz_t num, const fmpz_t den)
    char * fmpq_get_str(char * str, int b, const fmpq_t x)
    # void flint_mpq_init_set_readonly(mpq_t z, const fmpq_t f)
    # void flint_mpq_clear_readonly(mpq_t z)
    # void fmpq_init_set_readonly(fmpq_t f, const mpq_t z)
    void fmpq_clear_readonly(fmpq_t f)
    # int fmpq_fprint(FILE * file, const fmpq_t x)
    # int _fmpq_fprint(FILE * file, const fmpz_t num, const fmpz_t den)
    int fmpq_print(const fmpq_t x)
    int _fmpq_print(const fmpz_t num, const fmpz_t den)
    void fmpq_randtest(fmpq_t res, flint_rand_t state, flint_bitcnt_t bits)
    void _fmpq_randtest(fmpz_t num, fmpz_t den, flint_rand_t state, flint_bitcnt_t bits)
    void fmpq_randtest_not_zero(fmpq_t res, flint_rand_t state, flint_bitcnt_t bits)
    void fmpq_randbits(fmpq_t res, flint_rand_t state, flint_bitcnt_t bits)
    void _fmpq_randbits(fmpz_t num, fmpz_t den, flint_rand_t state, flint_bitcnt_t bits)
    void fmpq_add(fmpq_t res, const fmpq_t op1, const fmpq_t op2)
    void fmpq_sub(fmpq_t res, const fmpq_t op1, const fmpq_t op2)
    void fmpq_mul(fmpq_t res, const fmpq_t op1, const fmpq_t op2)
    void fmpq_div(fmpq_t res, const fmpq_t op1, const fmpq_t op2)
    void _fmpq_add(fmpz_t rnum, fmpz_t rden, const fmpz_t op1num, const fmpz_t op1den, const fmpz_t op2num, const fmpz_t op2den)
    void _fmpq_sub(fmpz_t rnum, fmpz_t rden, const fmpz_t op1num, const fmpz_t op1den, const fmpz_t op2num, const fmpz_t op2den)
    void _fmpq_mul(fmpz_t rnum, fmpz_t rden, const fmpz_t op1num, const fmpz_t op1den, const fmpz_t op2num, const fmpz_t op2den)
    void _fmpq_div(fmpz_t rnum, fmpz_t rden, const fmpz_t op1num, const fmpz_t op1den, const fmpz_t op2num, const fmpz_t op2den)
    void _fmpq_add_si(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q, slong r)
    void _fmpq_sub_si(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q, slong r)
    void _fmpq_add_ui(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q, ulong r)
    void _fmpq_sub_ui(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q, ulong r)
    void _fmpq_add_fmpz(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q, const fmpz_t r)
    void _fmpq_sub_fmpz(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q, const fmpz_t r)
    void fmpq_add_si(fmpq_t res, const fmpq_t op1, slong c)
    void fmpq_sub_si(fmpq_t res, const fmpq_t op1, slong c)
    void fmpq_add_ui(fmpq_t res, const fmpq_t op1, ulong c)
    void fmpq_sub_ui(fmpq_t res, const fmpq_t op1, ulong c)
    void fmpq_add_fmpz(fmpq_t res, const fmpq_t op1, const fmpz_t c)
    void fmpq_sub_fmpz(fmpq_t res, const fmpq_t op1, const fmpz_t c)
    void _fmpq_mul_si(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q, slong r)
    void fmpq_mul_si(fmpq_t res, const fmpq_t op1, slong c)
    void _fmpq_mul_ui(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q, ulong r)
    void fmpq_mul_ui(fmpq_t res, const fmpq_t op1, ulong c)
    void fmpq_addmul(fmpq_t res, const fmpq_t op1, const fmpq_t op2)
    void fmpq_submul(fmpq_t res, const fmpq_t op1, const fmpq_t op2)
    void _fmpq_addmul(fmpz_t rnum, fmpz_t rden, const fmpz_t op1num, const fmpz_t op1den, const fmpz_t op2num, const fmpz_t op2den)
    void _fmpq_submul(fmpz_t rnum, fmpz_t rden, const fmpz_t op1num, const fmpz_t op1den, const fmpz_t op2num, const fmpz_t op2den)
    void fmpq_inv(fmpq_t dest, const fmpq_t src)
    void _fmpq_pow_si(fmpz_t rnum, fmpz_t rden, const fmpz_t opnum, const fmpz_t opden, slong e)
    void fmpq_pow_si(fmpq_t res, const fmpq_t op, slong e)
    int fmpq_pow_fmpz(fmpq_t a, const fmpq_t b, const fmpz_t e)
    void fmpq_mul_fmpz(fmpq_t res, const fmpq_t op, const fmpz_t x)
    void fmpq_div_fmpz(fmpq_t res, const fmpq_t op, const fmpz_t x)
    void fmpq_mul_2exp(fmpq_t res, const fmpq_t x, flint_bitcnt_t exp)
    void fmpq_div_2exp(fmpq_t res, const fmpq_t x, flint_bitcnt_t exp)
    void _fmpq_gcd(fmpz_t rnum, fmpz_t rden, const fmpz_t p, const fmpz_t q, const fmpz_t r, const fmpz_t s)
    void fmpq_gcd(fmpq_t res, const fmpq_t op1, const fmpq_t op2)
    void _fmpq_gcd_cofactors(fmpz_t gnum, fmpz_t gden, fmpz_t abar, fmpz_t bbar, const fmpz_t anum, const fmpz_t aden, const fmpz_t bnum, const fmpz_t bden)
    void fmpq_gcd_cofactors(fmpq_t g, fmpz_t abar, fmpz_t bbar, const fmpq_t a, const fmpq_t b)
    void _fmpq_add_small(fmpz_t rnum, fmpz_t rden, slong p1, ulong q1, slong p2, ulong q2)
    void _fmpq_mul_small(fmpz_t rnum, fmpz_t rden, slong p1, ulong q1, slong p2, ulong q2)
    int _fmpq_mod_fmpz(fmpz_t res, const fmpz_t num, const fmpz_t den, const fmpz_t mod)
    int fmpq_mod_fmpz(fmpz_t res, const fmpq_t x, const fmpz_t mod)
    int _fmpq_reconstruct_fmpz_2_naive(fmpz_t n, fmpz_t d, const fmpz_t a, const fmpz_t m, const fmpz_t N, const fmpz_t D)
    int _fmpq_reconstruct_fmpz_2(fmpz_t n, fmpz_t d, const fmpz_t a, const fmpz_t m, const fmpz_t N, const fmpz_t D)
    int fmpq_reconstruct_fmpz_2(fmpq_t res, const fmpz_t a, const fmpz_t m, const fmpz_t N, const fmpz_t D)
    int _fmpq_reconstruct_fmpz(fmpz_t n, fmpz_t d, const fmpz_t a, const fmpz_t m)
    int fmpq_reconstruct_fmpz(fmpq_t res, const fmpz_t a, const fmpz_t m)
    void _fmpq_next_minimal(fmpz_t rnum, fmpz_t rden, const fmpz_t num, const fmpz_t den)
    void fmpq_next_minimal(fmpq_t res, const fmpq_t x)
    void _fmpq_next_signed_minimal(fmpz_t rnum, fmpz_t rden, const fmpz_t num, const fmpz_t den)
    void fmpq_next_signed_minimal(fmpq_t res, const fmpq_t x)
    void _fmpq_next_calkin_wilf(fmpz_t rnum, fmpz_t rden, const fmpz_t num, const fmpz_t den)
    void fmpq_next_calkin_wilf(fmpq_t res, const fmpq_t x)
    void _fmpq_next_signed_calkin_wilf(fmpz_t rnum, fmpz_t rden, const fmpz_t num, const fmpz_t den)
    void fmpq_next_signed_calkin_wilf(fmpq_t res, const fmpq_t x)
    void fmpq_farey_neighbors(fmpq_t l, fmpq_t r, const fmpq_t x, const fmpz_t Q)
    void _fmpq_simplest_between(fmpz_t x_num, fmpz_t x_den, const fmpz_t l_num, const fmpz_t l_den, const fmpz_t r_num, const fmpz_t r_den)
    void fmpq_simplest_between(fmpq_t x, const fmpq_t l, const fmpq_t r)
    slong fmpq_get_cfrac(fmpz_struct * c, fmpq_t rem, const fmpq_t x, slong n)
    slong fmpq_get_cfrac_naive(fmpz_struct * c, fmpq_t rem, const fmpq_t x, slong n)
    void fmpq_set_cfrac(fmpq_t x, const fmpz_struct * c, slong n)
    slong fmpq_cfrac_bound(const fmpq_t x)
    void _fmpq_harmonic_ui(fmpz_t num, fmpz_t den, ulong n)
    void fmpq_harmonic_ui(fmpq_t x, ulong n)
    void fmpq_dedekind_sum(fmpq_t s, const fmpz_t h, const fmpz_t k)
    void fmpq_dedekind_sum_naive(fmpq_t s, const fmpz_t h, const fmpz_t k)
