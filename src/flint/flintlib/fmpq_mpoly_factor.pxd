from flint.flintlib.fmpq cimport fmpq_t
from flint.flintlib.fmpz cimport fmpz_t
from flint.flintlib.fmpq_mpoly cimport fmpq_mpoly_ctx_t, fmpq_mpoly_t, fmpq_mpoly_struct
from flint.flintlib.types.flint cimport slong, fmpz_struct

cdef extern from "flint/fmpq_mpoly_factor.h":
    ctypedef struct fmpq_mpoly_factor_struct:
        fmpq_t constant
        fmpq_mpoly_struct * poly
        fmpz_struct * exp
        slong num
        slong alloc
    ctypedef fmpq_mpoly_factor_struct fmpq_mpoly_factor_t[1]

    void fmpq_mpoly_factor_init(fmpq_mpoly_factor_t f, const fmpq_mpoly_ctx_t ctx)
    void fmpq_mpoly_factor_clear(fmpq_mpoly_factor_t f, const fmpq_mpoly_ctx_t ctx)
    void fmpq_mpoly_factor_swap(fmpq_mpoly_factor_t f, fmpq_mpoly_factor_t g, const fmpq_mpoly_ctx_t ctx)
    slong fmpq_mpoly_factor_length(const fmpq_mpoly_factor_t f, const fmpq_mpoly_ctx_t ctx)
    void fmpq_mpoly_factor_get_constant_fmpq(fmpq_t c, const fmpq_mpoly_factor_t f, const fmpq_mpoly_ctx_t ctx)
    void fmpq_mpoly_factor_get_base(fmpq_mpoly_t B, const fmpq_mpoly_factor_t f, slong i, const fmpq_mpoly_ctx_t ctx)
    void fmpq_mpoly_factor_swap_base(fmpq_mpoly_t B, fmpq_mpoly_factor_t f, slong i, const fmpq_mpoly_ctx_t ctx)
    slong fmpq_mpoly_factor_get_exp_si(fmpq_mpoly_factor_t f, slong i, const fmpq_mpoly_ctx_t ctx)
    void fmpq_mpoly_factor_sort(fmpq_mpoly_factor_t f, const fmpq_mpoly_ctx_t ctx)
    int fmpq_mpoly_factor_make_monic(fmpq_mpoly_factor_t f, const fmpq_mpoly_ctx_t ctx)
    int fmpq_mpoly_factor_make_integral(fmpq_mpoly_factor_t f, const fmpq_mpoly_ctx_t ctx)
    int fmpq_mpoly_factor_squarefree(fmpq_mpoly_factor_t f, const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)
    int fmpq_mpoly_factor(fmpq_mpoly_factor_t f, const fmpq_mpoly_t A, const fmpq_mpoly_ctx_t ctx)
