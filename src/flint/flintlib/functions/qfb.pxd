from flint.flintlib.types.flint cimport fmpz_t, slong, ulong
from flint.flintlib.types.qfb cimport qfb_t

# unknown type qfb
# unknown type qfb_hash_t


cdef extern from "flint/qfb.h":
    void qfb_init(qfb_t q)
    void qfb_clear(qfb_t q)
    # void qfb_array_clear(qfb ** forms, slong num)
    # qfb_hash_t * qfb_hash_init(slong depth)
    # void qfb_hash_clear(qfb_hash_t * qhash, slong depth)
    # void qfb_hash_insert(qfb_hash_t * qhash, qfb_t q, qfb_t q2, slong iter, slong depth)
    # slong qfb_hash_find(qfb_hash_t * qhash, qfb_t q, slong depth)
    void qfb_set(qfb_t f, qfb_t g)
    int qfb_equal(qfb_t f, qfb_t g)
    void qfb_print(qfb_t q)
    void qfb_discriminant(fmpz_t D, qfb_t f)
    void qfb_reduce(qfb_t r, qfb_t f, fmpz_t D)
    int qfb_is_reduced(qfb_t r)
    # slong qfb_reduced_forms(qfb ** forms, slong d)
    # slong qfb_reduced_forms_large(qfb ** forms, slong d)
    void qfb_nucomp(qfb_t r, const qfb_t f, const qfb_t g, fmpz_t D, fmpz_t L)
    void qfb_nudupl(qfb_t r, const qfb_t f, fmpz_t D, fmpz_t L)
    void qfb_pow_ui(qfb_t r, qfb_t f, fmpz_t D, ulong exp)
    void qfb_pow(qfb_t r, qfb_t f, fmpz_t D, fmpz_t exp)
    void qfb_inverse(qfb_t r, qfb_t f)
    int qfb_is_principal_form(qfb_t f, fmpz_t D)
    void qfb_principal_form(qfb_t f, fmpz_t D)
    int qfb_is_primitive(qfb_t f)
    void qfb_prime_form(qfb_t r, fmpz_t D, fmpz_t p)
    int qfb_exponent_element(fmpz_t exponent, qfb_t f, fmpz_t n, ulong B1, ulong B2_sqrt)
    int qfb_exponent(fmpz_t exponent, fmpz_t n, ulong B1, ulong B2_sqrt, slong c)
    int qfb_exponent_grh(fmpz_t exponent, fmpz_t n, ulong B1, ulong B2_sqrt)
