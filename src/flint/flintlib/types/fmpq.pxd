from flint.flintlib.types.flint cimport slong, fmpz_struct, fmpz_t
from flint.flintlib.types.flint cimport fmpq as fmpq_struct
from flint.flintlib.types.fmpz cimport fmpz_mpoly_ctx_t, fmpz_mpoly_t

cdef extern from "flint/fmpq_types.h":

    ctypedef fmpq_struct fmpq_t[1]

    ctypedef struct fmpq_mat_struct:
        fmpq_struct * entries
        slong r
        slong c
        fmpq_struct ** rows

    ctypedef fmpq_mat_struct fmpq_mat_t[1]

    ctypedef struct fmpq_poly_struct:
        fmpz_struct * coeffs
        slong alloc
        slong length
        fmpz_t den

    ctypedef fmpq_poly_struct fmpq_poly_t[1]

    ctypedef struct fmpq_mpoly_ctx_struct:
        fmpz_mpoly_ctx_t zctx

    ctypedef fmpq_mpoly_ctx_struct fmpq_mpoly_ctx_t[1]

    ctypedef struct fmpq_mpoly_struct:
        fmpq_t content
        fmpz_mpoly_t zpoly

    ctypedef fmpq_mpoly_struct fmpq_mpoly_t[1]

    ctypedef struct fmpq_mpoly_factor_struct:
        fmpq_t constant
        fmpq_mpoly_struct * poly
        fmpz_struct * exp
        slong num
        slong alloc

    ctypedef fmpq_mpoly_factor_struct fmpq_mpoly_factor_t[1]
