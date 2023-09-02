from flint._flint cimport *
from flint.flint_base.flint_base cimport flint_mat
from flint.flint_base.flint_base cimport flint_mpoly
from flint.flint_base.flint_base cimport flint_series
from flint.flint_base.flint_base cimport flint_scalar
from flint.flint_base.flint_base cimport flint_poly

from flint.fmpz cimport fmpz
from flint._flint cimport *

cdef FMPZ_UNKNOWN = 0
cdef FMPZ_REF = 1
cdef FMPZ_TMP = 2


cdef flint_rand_t global_random_state

cdef class Context:
    cdef public bint pretty
    cdef public long prec
    cdef arf_rnd_t rnd




#cdef class fmpz:
#    cdef fmpz_t val

#cdef class fmpz_poly:
#    cdef fmpz_poly_t val

# cdef class fmpz_mat(flint_mat):
#     cdef fmpz_mat_t val
#     cpdef long nrows(self)
#     cpdef long ncols(self)
#     cdef __mul_fmpz(self, fmpz c)

# cdef class fmpz_series(flint_series):
#     cdef fmpz_poly_t val
#     cdef long prec
#     cpdef long length(self)
#     cpdef valuation(self)

# cdef any_as_fmpq(obj)
# cdef class fmpq(flint_scalar):
#     cdef fmpq_t val
from flint.fmpq cimport fmpq

# cdef any_as_fmpq_poly(obj)
# cdef class fmpq_poly(flint_poly):
#     cdef fmpq_poly_t val
#     cpdef long length(self)
#     cpdef long degree(self)

# from flint.fmpz_mat cimport fmpz_mat
# cdef class fmpq_mat(flint_mat):
#     cdef fmpq_mat_t val
#     cpdef long nrows(self)
#     cpdef long ncols(self)
#     cdef __mul_fmpz(self, fmpz c)
#     cdef __mul_fmpq(self, fmpq c)
#     cdef __mul_fmpq_mat(self, fmpq_mat other)
#     cdef __mul_fmpz_mat(self, fmpz_mat other)
#     cdef __mul_r_fmpz_mat(self, fmpz_mat other)

# cdef class fmpq_series(flint_series):
#     cdef fmpq_poly_t val
#     cdef long prec
#     cpdef long length(self)
#     cpdef valuation(self)
#     cdef bint zero_constant_term(s)
#     cdef bint one_constant_term(s)

# cdef class nmod(flint_scalar):
#     cdef mp_limb_t val
#     cdef nmod_t mod

# cdef class nmod_poly(flint_poly):
#     cdef nmod_poly_t val
#     cpdef long length(self)
#     cpdef long degree(self)
#     cpdef mp_limb_t modulus(self)

# cdef class nmod_mat:
#     cdef nmod_mat_t val
#     cpdef long nrows(self)
#     cpdef long ncols(self)
#     cpdef mp_limb_t modulus(self)
#     cdef __mul_nmod(self, mp_limb_t c)

# cdef class nmod_series(flint_series):
#     pass
    # cdef nmod_poly_t val
    # cdef long prec

# cdef class arf:
#     cdef arf_t val
#     cpdef bint is_finite(self)
#     cpdef bint is_pos_inf(self)
#     cpdef bint is_neg_inf(self)
#     cpdef bint is_nan(self)
#     cpdef bint is_zero(self)

# cdef any_as_arb_or_notimplemented(x)
# cdef class arb(flint_scalar):
#     cdef arb_t val

#     cpdef bint is_zero(self)
#     cpdef bint is_finite(self)
#     cpdef bint is_nan(self)
#     cpdef bint is_exact(self)
#     cpdef bint is_integer(self)


cdef any_as_acb_or_notimplemented(x)
cdef class acb(flint_scalar):
    cdef acb_t val
    cpdef bint is_zero(self)
    cpdef bint is_finite(self)
    cpdef bint is_exact(self)

cdef class arb_poly(flint_poly):
    cdef arb_poly_t val
    cpdef long length(self)
    cpdef long degree(self)

cdef class acb_poly(flint_poly):
    cdef acb_poly_t val
    cpdef long length(self)
    cpdef long degree(self)


cdef class arb_mat(flint_mat):
    cdef arb_mat_t val
    cpdef long nrows(self)
    cpdef long ncols(self)

cdef class acb_mat(flint_mat):
    cdef acb_mat_t val
    cpdef long nrows(self)
    cpdef long ncols(self)

cdef class arb_series(flint_series):
    cdef arb_poly_t val
    cdef long prec
    cpdef long length(self)
    cpdef valuation(self)

cdef class acb_series(flint_series):
    cdef acb_poly_t val
    cdef long prec
    cpdef long length(self)
    cpdef valuation(self)

cdef class fmpz_mpoly_ctx:
    cdef fmpz_mpoly_ctx_t val
    cpdef slong nvars(self)

    cpdef ordering(self)

cdef class fmpz_mpoly(flint_mpoly):
    cdef fmpz_mpoly_t val
    cdef fmpz_mpoly_ctx ctx
    cdef bint _init

