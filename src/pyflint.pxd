from flint cimport *

cdef class Context:
    cpdef public bint pretty
    cpdef public long prec
    cpdef arf_rnd_t rnd

cdef class fmpz:
    cdef fmpz_t val

cdef class fmpz_poly:
    cdef fmpz_poly_t val

cdef class fmpz_mat:
    cdef fmpz_mat_t val

cdef class fmpq:
    cdef fmpq_t val

cdef class fmpq_poly:
    cdef fmpq_poly_t val

cdef class fmpq_mat:
    cdef fmpq_mat_t val

cdef class nmod:
    cdef mp_limb_t val
    cdef nmod_t mod

cdef class nmod_poly:
    cdef nmod_poly_t val

cdef class nmod_mat:
    cdef nmod_mat_t val

cdef class arf:
    cdef arf_t val

cdef class arb:
    cdef arb_t val

cdef class acb:
    cdef acb_t val

