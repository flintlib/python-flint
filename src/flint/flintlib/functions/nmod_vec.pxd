from flint.flintlib.types.flint cimport flint_bitcnt_t, flint_rand_t, mp_limb_t, mp_ptr, mp_srcptr, nmod_t, slong

# unknown type FILE

# .. macro:: NMOD_VEC_DOT(res, i, len, expr1, expr2, mod, nlimbs)

cdef extern from "flint/nmod_vec.h":
    mp_ptr _nmod_vec_init(slong len)
    void _nmod_vec_clear(mp_ptr vec)
    void _nmod_vec_randtest(mp_ptr vec, flint_rand_t state, slong len, nmod_t mod)
    void _nmod_vec_set(mp_ptr res, mp_srcptr vec, slong len)
    void _nmod_vec_zero(mp_ptr vec, slong len)
    void _nmod_vec_swap(mp_ptr a, mp_ptr b, slong length)
    void _nmod_vec_reduce(mp_ptr res, mp_srcptr vec, slong len, nmod_t mod)
    flint_bitcnt_t _nmod_vec_max_bits(mp_srcptr vec, slong len)
    int _nmod_vec_equal(mp_srcptr vec, mp_srcptr vec2, slong len)
    void _nmod_vec_print_pretty(mp_srcptr vec, slong len, nmod_t mod)
    # int _nmod_vec_fprint_pretty(FILE * file, mp_srcptr vec, slong len, nmod_t mod)
    int _nmod_vec_print(mp_srcptr vec, slong len, nmod_t mod)
    # int _nmod_vec_fprint(FILE * f, mp_srcptr vec, slong len, nmod_t mod)
    void _nmod_vec_add(mp_ptr res, mp_srcptr vec1, mp_srcptr vec2, slong len, nmod_t mod)
    void _nmod_vec_sub(mp_ptr res, mp_srcptr vec1, mp_srcptr vec2, slong len, nmod_t mod)
    void _nmod_vec_neg(mp_ptr res, mp_srcptr vec, slong len, nmod_t mod)
    void _nmod_vec_scalar_mul_nmod(mp_ptr res, mp_srcptr vec, slong len, mp_limb_t c, nmod_t mod)
    void _nmod_vec_scalar_mul_nmod_shoup(mp_ptr res, mp_srcptr vec, slong len, mp_limb_t c, nmod_t mod)
    void _nmod_vec_scalar_addmul_nmod(mp_ptr res, mp_srcptr vec, slong len, mp_limb_t c, nmod_t mod)
    int _nmod_vec_dot_bound_limbs(slong len, nmod_t mod)
    mp_limb_t _nmod_vec_dot(mp_srcptr vec1, mp_srcptr vec2, slong len, nmod_t mod, int nlimbs)
    mp_limb_t _nmod_vec_dot_rev(mp_srcptr vec1, mp_srcptr vec2, slong len, nmod_t mod, int nlimbs)
    mp_limb_t _nmod_vec_dot_ptr(mp_srcptr vec1, const mp_ptr * vec2, slong offset, slong len, nmod_t mod, int nlimbs)
