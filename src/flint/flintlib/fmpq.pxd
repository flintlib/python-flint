from flint.flintlib.flint cimport ulong, flint_rand_t, mp_bitcnt_t
from flint.flintlib.fmpz cimport fmpz_struct, fmpz_t

cdef extern from "flint/fmpq.h":
    ctypedef struct fmpq_struct:
        fmpz_struct num
        fmpz_struct den
    ctypedef fmpq_struct fmpq_t[1]
    fmpz_struct * fmpq_numref(fmpq_t x)
    fmpz_struct * fmpq_denref(fmpq_t x)
    void fmpq_init(fmpq_t x)
    void fmpq_clear(fmpq_t x)
    void fmpq_zero(fmpq_t res)
    void fmpq_one(fmpq_t res)
    int fmpq_equal(fmpq_t x, fmpq_t y)
    int fmpq_sgn(fmpq_t x)
    int fmpq_is_zero(fmpq_t x)
    int fmpq_is_one(fmpq_t x)
    void fmpq_set(fmpq_t dest, fmpq_t src)
    void fmpq_neg(fmpq_t dest, fmpq_t src)
    void fmpq_canonicalise(fmpq_t res)
    int fmpq_is_canonical(fmpq_t x)
    void fmpq_set_si(fmpq_t res, long p, ulong q)
    #void fmpq_set_mpq(fmpq_t dest, mpq_t src)
    #void fmpq_get_mpq(mpq_t dest, fmpq_t src)
    void fmpq_print(fmpq_t x)
    void fmpq_randtest(fmpq_t res, flint_rand_t state, mp_bitcnt_t bits)
    void fmpq_randbits(fmpq_t res, flint_rand_t state, mp_bitcnt_t bits)
    void fmpq_add(fmpq_t res, fmpq_t op1, fmpq_t op2)
    void fmpq_sub(fmpq_t res, fmpq_t op1, fmpq_t op2)
    void fmpq_mul(fmpq_t res, fmpq_t op1, fmpq_t op2)
    void fmpq_mul_fmpz(fmpq_t res, fmpq_t op, fmpz_t x)
    void fmpq_addmul(fmpq_t res, fmpq_t op1, fmpq_t op2)
    void fmpq_submul(fmpq_t res, fmpq_t op1, fmpq_t op2)
    void fmpq_inv(fmpq_t dest, fmpq_t src)
    void fmpq_div(fmpq_t res, fmpq_t op1, fmpq_t op2)
    void fmpq_div_fmpz(fmpq_t res, fmpq_t op, fmpz_t x)
    int fmpq_mod_fmpz(fmpz_t res, fmpq_t x, fmpz_t mod)
    int fmpq_pow_fmpz(fmpq_t a, const fmpq_t b, const fmpz_t e)
    int fmpq_reconstruct_fmpz(fmpq_t res, fmpz_t a, fmpz_t m)
    int fmpq_reconstruct_fmpz_2(fmpq_t res, fmpz_t a, fmpz_t m, fmpz_t N, fmpz_t D)
    mp_bitcnt_t fmpq_height_bits(fmpq_t x)
    void fmpq_height(fmpz_t height, fmpq_t x)
    void fmpq_next_calkin_wilf(fmpq_t res, fmpq_t x)
    void fmpq_next_signed_calkin_wilf(fmpq_t res, fmpq_t x)
    void fmpq_next_minimal(fmpq_t res, fmpq_t x)
    void fmpq_next_signed_minimal(fmpq_t res, fmpq_t x)
    void fmpq_harmonic_ui(fmpq_t res, ulong n)
    void fmpq_dedekind_sum(fmpq_t s, const fmpz_t h, const fmpz_t k)
