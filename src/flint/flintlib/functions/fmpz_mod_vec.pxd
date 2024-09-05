from flint.flintlib.types.flint cimport fmpz_struct, fmpz_t, slong
from flint.flintlib.types.fmpz_mod cimport fmpz_mod_ctx_t



cdef extern from "flint/fmpz_mod_vec.h":
    void _fmpz_mod_vec_set_fmpz_vec(fmpz_struct * A, const fmpz_struct * B, slong len, const fmpz_mod_ctx_t ctx)
    void _fmpz_mod_vec_neg(fmpz_struct * A, const fmpz_struct * B, slong len, const fmpz_mod_ctx_t ctx)
    void _fmpz_mod_vec_add(fmpz_struct * a, const fmpz_struct * b, const fmpz_struct * c, slong n, const fmpz_mod_ctx_t ctx)
    void _fmpz_mod_vec_sub(fmpz_struct * a, const fmpz_struct * b, const fmpz_struct * c, slong n, const fmpz_mod_ctx_t ctx)
    void _fmpz_mod_vec_scalar_mul_fmpz_mod(fmpz_struct * A, const fmpz_struct * B, slong len, const fmpz_t c, const fmpz_mod_ctx_t ctx)
    void _fmpz_mod_vec_scalar_addmul_fmpz_mod(fmpz_struct * A, const fmpz_struct * B, slong len, const fmpz_t c, const fmpz_mod_ctx_t ctx)
    void _fmpz_mod_vec_scalar_div_fmpz_mod(fmpz_struct * A, const fmpz_struct * B, slong len, const fmpz_t c, const fmpz_mod_ctx_t ctx)
    void _fmpz_mod_vec_dot(fmpz_t d, const fmpz_struct * A, const fmpz_struct * B, slong len, const fmpz_mod_ctx_t ctx)
    void _fmpz_mod_vec_dot_rev(fmpz_t d, const fmpz_struct * A, const fmpz_struct * B, slong len, const fmpz_mod_ctx_t ctx)
    void _fmpz_mod_vec_mul(fmpz_struct * A, const fmpz_struct * B, const fmpz_struct * C, slong len, const fmpz_mod_ctx_t ctx)
