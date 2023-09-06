from flint.flintlib.fmpz cimport fmpz_t, fmpz_struct
from flint._flint cimport mp_limb_t, ulong, flint_rand_t
from flint.flintlib.fmpq cimport fmpq_t
from flint.flintlib.mag cimport mag_t

cdef extern from "flint/arf.h":
    ctypedef struct arf_struct:
        fmpz_struct exp
        long size
        mp_limb_t d0
        mp_limb_t d1

    ctypedef arf_struct arf_t[1]
    ctypedef arf_struct * arf_ptr
    ctypedef const arf_struct * arf_srcptr

    ctypedef int arf_rnd_t
    cdef arf_rnd_t ARF_RND_DOWN
    cdef arf_rnd_t ARF_RND_NEAR
    cdef arf_rnd_t ARF_RND_FLOOR
    cdef arf_rnd_t ARF_RND_CEIL
    cdef arf_rnd_t ARF_RND_UP

    void arf_init(arf_t x)
    void arf_clear(arf_t x)
    void arf_zero(arf_t x)
    void arf_pos_inf(arf_t x)
    void arf_neg_inf(arf_t x)
    void arf_nan(arf_t x)
    int arf_is_special(const arf_t x)
    int arf_is_zero(const arf_t x)
    int arf_is_pos_inf(const arf_t x)
    int arf_is_neg_inf(const arf_t x)
    int arf_is_nan(const arf_t x)
    int arf_is_normal(const arf_t x)
    int arf_is_finite(const arf_t x)
    int arf_is_inf(const arf_t x)
    void arf_one(arf_t x)
    int arf_is_one(const arf_t x)
    int arf_sgn(const arf_t x)
    int arf_cmp(const arf_t x, const arf_t y)
    int arf_cmpabs(const arf_t x, const arf_t y)
    void arf_swap(arf_t y, arf_t x)
    void arf_set(arf_t y, const arf_t x)
    void arf_neg(arf_t y, const arf_t x)
    void arf_init_set_ui(arf_t x, ulong v)
    void arf_init_set_si(arf_t x, long v)
    void arf_set_ui(arf_t x, ulong v)
    void arf_set_si(arf_t x, long v)
    int arf_cmpabs_ui(const arf_t x, ulong y)
    void arf_init_set_shallow(arf_t z, const arf_t x)
    void arf_init_neg_shallow(arf_t z, const arf_t x)
    void arf_init_set_mag_shallow(arf_t y, const mag_t x)
    void arf_init_neg_mag_shallow(arf_t z, const mag_t x)
    int arf_cmpabs_mag(const arf_t x, const mag_t y)
    int arf_mag_cmpabs(const mag_t x, const arf_t y)
    void arf_set_fmpz(arf_t y, const fmpz_t x)
    int arf_set_round_ui(arf_t x, ulong v, long prec, arf_rnd_t rnd)
    int arf_set_round_si(arf_t x, long v, long prec, arf_rnd_t rnd)
    int arf_set_round_fmpz(arf_t y, const fmpz_t x, long prec, arf_rnd_t rnd)
    int arf_set_round(arf_t y, const arf_t x, long prec, arf_rnd_t rnd)
    int arf_neg_round(arf_t y, const arf_t x, long prec, arf_rnd_t rnd)
    int arf_equal(const arf_t x, const arf_t y)
    void arf_min(arf_t z, const arf_t a, const arf_t b)
    void arf_max(arf_t z, const arf_t a, const arf_t b)
    void arf_abs(arf_t y, const arf_t x)
    long arf_bits(const arf_t x)
    void arf_bot(fmpz_t e, const arf_t x)
    int arf_is_int(const arf_t x)
    int arf_is_int_2exp_si(const arf_t x, long e)
    int arf_cmp_2exp_si(const arf_t x, long e)
    int arf_cmpabs_2exp_si(const arf_t x, long e)
    void arf_set_si_2exp_si(arf_t x, long man, long exp)
    void arf_set_ui_2exp_si(arf_t x, ulong man, long exp)
    void arf_mul_2exp_si(arf_t y, const arf_t x, long e)
    void arf_mul_2exp_fmpz(arf_t y, const arf_t x, const fmpz_t e)
    int arf_set_round_fmpz_2exp(arf_t y, const fmpz_t x, const fmpz_t exp, long prec, arf_rnd_t rnd)
    void arf_abs_bound_lt_2exp_fmpz(fmpz_t b, const arf_t x)
    void arf_abs_bound_le_2exp_fmpz(fmpz_t b, const arf_t x)
    long arf_abs_bound_lt_2exp_si(const arf_t x)
    void arf_get_fmpz_2exp(fmpz_t man, fmpz_t exp, const arf_t x)
    void arf_get_fmpz(fmpz_t z, const arf_t x, arf_rnd_t rnd)
    long arf_get_si(const arf_t x, arf_rnd_t rnd)
    int arf_get_fmpz_fixed_fmpz(fmpz_t y, const arf_t x, const fmpz_t e)
    int arf_get_fmpz_fixed_si(fmpz_t y, const arf_t x, long e)
    void arf_set_fmpz_2exp(arf_t x, const fmpz_t man, const fmpz_t exp)
    void arf_floor(arf_t z, const arf_t x)
    void arf_ceil(arf_t z, const arf_t x)
    void arf_debug(const arf_t x)
    void arf_print(const arf_t x)
    void arf_printd(const arf_t y, long d)
    void arf_randtest(arf_t x, flint_rand_t state, long bits, long mag_bits)
    void arf_randtest_not_zero(arf_t x, flint_rand_t state, long bits, long mag_bits)
    void arf_randtest_special(arf_t x, flint_rand_t state, long bits, long mag_bits)
    int arf_mul(arf_t z, const arf_t x, const arf_t y, long prec, arf_rnd_t rnd)
    int arf_neg_mul(arf_t z, const arf_t x, const arf_t y, long prec, arf_rnd_t rnd)
    int arf_mul_ui(arf_ptr z, arf_srcptr x, ulong y, long prec, arf_rnd_t rnd)
    int arf_mul_si(arf_ptr z, arf_srcptr x, long y, long prec, arf_rnd_t rnd)
    int arf_mul_fmpz(arf_ptr z, arf_srcptr x, const fmpz_t y, long prec, arf_rnd_t rnd)
    int arf_add(arf_ptr z, arf_srcptr x, arf_srcptr y, long prec, arf_rnd_t rnd)
    int arf_add_si(arf_ptr z, arf_srcptr x, long y, long prec, arf_rnd_t rnd)
    int arf_add_ui(arf_ptr z, arf_srcptr x, ulong y, long prec, arf_rnd_t rnd)
    int arf_add_fmpz(arf_ptr z, arf_srcptr x, const fmpz_t y, long prec, arf_rnd_t rnd)
    int arf_add_fmpz_2exp(arf_ptr z, arf_srcptr x, const fmpz_t y, const fmpz_t exp, long prec, arf_rnd_t rnd)
    int arf_sub(arf_ptr z, arf_srcptr x, arf_srcptr y, long prec, arf_rnd_t rnd)
    int arf_sub_si(arf_ptr z, arf_srcptr x, long y, long prec, arf_rnd_t rnd)
    int arf_sub_ui(arf_ptr z, arf_srcptr x, ulong y, long prec, arf_rnd_t rnd)
    int arf_sub_fmpz(arf_ptr z, arf_srcptr x, const fmpz_t y, long prec, arf_rnd_t rnd)
    int arf_addmul(arf_ptr z, arf_srcptr x, arf_srcptr y, long prec, arf_rnd_t rnd)
    int arf_addmul_ui(arf_ptr z, arf_srcptr x, ulong y, long prec, arf_rnd_t rnd)
    int arf_addmul_si(arf_ptr z, arf_srcptr x, long y, long prec, arf_rnd_t rnd)
    int arf_addmul_fmpz(arf_ptr z, arf_srcptr x, const fmpz_t y, long prec, arf_rnd_t rnd)
    int arf_submul(arf_ptr z, arf_srcptr x, arf_srcptr y, long prec, arf_rnd_t rnd)
    int arf_submul_ui(arf_ptr z, arf_srcptr x, ulong y, long prec, arf_rnd_t rnd)
    int arf_submul_si(arf_ptr z, arf_srcptr x, long y, long prec, arf_rnd_t rnd)
    int arf_submul_fmpz(arf_ptr z, arf_srcptr x, const fmpz_t y, long prec, arf_rnd_t rnd)
    int arf_div(arf_ptr z, arf_srcptr x, arf_srcptr y, long prec, arf_rnd_t rnd)
    int arf_div_ui(arf_ptr z, arf_srcptr x, ulong y, long prec, arf_rnd_t rnd)
    int arf_ui_div(arf_ptr z, ulong x, arf_srcptr y, long prec, arf_rnd_t rnd)
    int arf_div_si(arf_ptr z, arf_srcptr x, long y, long prec, arf_rnd_t rnd)
    int arf_si_div(arf_ptr z, long x, arf_srcptr y, long prec, arf_rnd_t rnd)
    int arf_div_fmpz(arf_ptr z, arf_srcptr x, const fmpz_t y, long prec, arf_rnd_t rnd)
    int arf_fmpz_div(arf_ptr z, const fmpz_t x, arf_srcptr y, long prec, arf_rnd_t rnd)
    int arf_fmpz_div_fmpz(arf_ptr z, const fmpz_t x, const fmpz_t y, long prec, arf_rnd_t rnd)
    int arf_sqrt(arf_ptr z, arf_srcptr x, long prec, arf_rnd_t rnd)
    int arf_sqrt_ui(arf_t z, ulong x, long prec, arf_rnd_t rnd)
    int arf_sqrt_fmpz(arf_t z, const fmpz_t x, long prec, arf_rnd_t rnd)
    int arf_rsqrt(arf_ptr z, arf_srcptr x, long prec, arf_rnd_t rnd)
    void arf_get_mag(mag_t y, const arf_t x)
    void arf_get_mag_lower(mag_t y, const arf_t x)
    void arf_set_mag(arf_t y, const mag_t x)
    void mag_init_set_arf(mag_t y, const arf_t x)
    void arf_mag_add_ulp(mag_t z, const mag_t x, const arf_t y, long prec)
    void arf_mag_set_ulp(mag_t z, const arf_t y, long prec)
    void arf_get_fmpq(fmpq_t y, const arf_t x)
    int arf_set_fmpq(arf_t y, const fmpq_t x, long prec, arf_rnd_t rnd)
    double arf_get_d(const arf_t x, arf_rnd_t rnd)
    void arf_set_d(arf_t x, double v)
