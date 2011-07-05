from flint cimport fmpz_t, fmpz_poly_t, fmpz_mat_t

cdef class fmpz:
    cdef fmpz_t val

cdef class fmpz_poly:
    cdef fmpz_poly_t val

cdef class fmpz_mat:
    cdef fmpz_mat_t val
