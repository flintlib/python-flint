from flint.flintlib.types.flint cimport flint_bitcnt_t, flint_rand_t, fmpz_struct, fmpz_t, slong
from flint.flintlib.types.fmpq cimport fmpq_struct, fmpq_t

# unknown type FILE


cdef extern from "flint/fmpq_vec.h":
    fmpq_struct * _fmpq_vec_init(slong n)
    void _fmpq_vec_clear(fmpq_struct * vec, slong n)
    void _fmpq_vec_randtest(fmpq_struct * f, flint_rand_t state, slong len, flint_bitcnt_t bits)
    void _fmpq_vec_randtest_uniq_sorted(fmpq_struct * vec, flint_rand_t state, slong len, flint_bitcnt_t bits)
    void _fmpq_vec_sort(fmpq_struct * vec, slong len)
    void _fmpq_vec_set_fmpz_vec(fmpq_struct * res, const fmpz_struct * vec, slong len)
    void _fmpq_vec_get_fmpz_vec_fmpz(fmpz_struct * num, fmpz_t den, const fmpq_struct * a, slong len)
    void _fmpq_vec_dot(fmpq_t res, const fmpq_struct * vec1, const fmpq_struct * vec2, slong len)
    # int _fmpq_vec_fprint(FILE * file, const fmpq_struct * vec, slong len)
    int _fmpq_vec_print(const fmpq_struct * vec, slong len)
