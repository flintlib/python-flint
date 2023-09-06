from flint._flint cimport *

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

cdef class fmpz_series:
    cdef fmpz_poly_t val
    cdef long prec

cdef class fmpq:
    cdef fmpq_t val

cdef class fmpq_poly:
    cdef fmpq_poly_t val

cdef class fmpq_mat:
    cdef fmpq_mat_t val

cdef class fmpq_series:
    cdef fmpq_poly_t val
    cdef long prec

cdef class nmod:
    cdef mp_limb_t val
    cdef nmod_t mod

cdef class nmod_poly:
    cdef nmod_poly_t val

cdef class nmod_mat:
    cdef nmod_mat_t val

cdef class nmod_series:
    cdef nmod_poly_t val
    cdef long prec

cdef class arf:
    cdef arf_t val

cdef class arb:
    cdef arb_t val

cdef class acb:
    cdef acb_t val

cdef class arb_poly:
    cdef arb_poly_t val

cdef class acb_poly:
    cdef acb_poly_t val

cdef class arb_mat:
    cdef arb_mat_t val

cdef class acb_mat:
    cdef acb_mat_t val

cdef class arb_series:
    cdef arb_poly_t val
    cdef long prec

cdef class acb_series:
    cdef acb_poly_t val
    cdef long prec

cdef class fmpz_mpoly_ctx:
    cdef fmpz_mpoly_ctx_t val

cdef class fmpz_mpoly:
    cdef fmpz_mpoly_t val
    cdef bint initialized

