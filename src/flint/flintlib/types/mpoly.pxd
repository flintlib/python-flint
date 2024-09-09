from flint.flintlib.types.flint cimport slong, FLINT_BITS

cdef extern from "flint/mpoly.h":

    ctypedef enum ordering_t:
        ORD_LEX
        ORD_DEGLEX
        ORD_DEGREVLEX

    ctypedef struct mpoly_ctx_struct:
        slong nvars
        slong nfields
        ordering_t ord
        int deg
        int rev
        slong lut_words_per_exp[FLINT_BITS]
        unsigned char lut_fix_bits[FLINT_BITS]

    ctypedef mpoly_ctx_struct mpoly_ctx_t[1]

    void mpoly_ctx_init(mpoly_ctx_t ctx, slong nvars, const ordering_t ord)
    void mpoly_ctx_clear(mpoly_ctx_t mctx)
