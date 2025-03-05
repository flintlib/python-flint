from flint.flintlib.types.flint cimport flint_bitcnt_t, flint_rand_t, nmod_t, nn_ptr, nn_srcptr, slong, ulong

# unknown type FILE
# unknown type dot_params_t

# .. macro:: NMOD_VEC_DOT(res, i, len, expr1, expr2, mod, params)

cdef extern from "flint/nmod_vec.h":
    nn_ptr _nmod_vec_init(slong len)
    void _nmod_vec_clear(nn_ptr vec)
    void _nmod_vec_randtest(nn_ptr vec, flint_rand_t state, slong len, nmod_t mod)
    void _nmod_vec_set(nn_ptr res, nn_srcptr vec, slong len)
    void _nmod_vec_zero(nn_ptr vec, slong len)
    void _nmod_vec_swap(nn_ptr a, nn_ptr b, slong length)
    void _nmod_vec_reduce(nn_ptr res, nn_srcptr vec, slong len, nmod_t mod)
    flint_bitcnt_t _nmod_vec_max_bits(nn_srcptr vec, slong len)
    int _nmod_vec_equal(nn_srcptr vec, nn_srcptr vec2, slong len)
    void _nmod_vec_print_pretty(nn_srcptr vec, slong len, nmod_t mod)
    # int _nmod_vec_fprint_pretty(FILE * file, nn_srcptr vec, slong len, nmod_t mod)
    int _nmod_vec_print(nn_srcptr vec, slong len, nmod_t mod)
    # int _nmod_vec_fprint(FILE * f, nn_srcptr vec, slong len, nmod_t mod)
    void _nmod_vec_add(nn_ptr res, nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod)
    void _nmod_vec_sub(nn_ptr res, nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod)
    void _nmod_vec_neg(nn_ptr res, nn_srcptr vec, slong len, nmod_t mod)
    void _nmod_vec_scalar_mul_nmod(nn_ptr res, nn_srcptr vec, slong len, ulong c, nmod_t mod)
    void _nmod_vec_scalar_mul_nmod_shoup(nn_ptr res, nn_srcptr vec, slong len, ulong c, nmod_t mod)
    void _nmod_vec_scalar_addmul_nmod(nn_ptr res, nn_srcptr vec, slong len, ulong c, nmod_t mod)
    # dot_params_t _nmod_vec_dot_params(slong len, nmod_t mod)
    # ulong _nmod_vec_dot(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod, dot_params_t params)
    # ulong _nmod_vec_dot_rev(nn_srcptr vec1, nn_srcptr vec2, slong len, nmod_t mod, dot_params_t params)
    # ulong _nmod_vec_dot_ptr(nn_srcptr vec1, const nn_ptr * vec2, slong offset, slong len, nmod_t mod, dot_params_t params)
    int _nmod_vec_dot_bound_limbs(slong len, nmod_t mod)
    # int _nmod_vec_dot_bound_limbs_from_params(slong len, nmod_t mod, dot_params_t params)
