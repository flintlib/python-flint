from flint.flintlib.types.nmod cimport nmod_mat_t

from flint.flintlib.types.flint cimport fmpz_struct
from flint.flintlib.types.nmod cimport nmod_mpoly_ctx_t, nmod_mpoly_t

cdef extern from "flint/nmod_types.h":
    int nmod_mat_is_square(const nmod_mat_t mat)

cdef extern from "flint/nmod_mpoly.h":
    void nmod_mpoly_deflation(fmpz_struct * shift, fmpz_struct * stride, const nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_deflate(nmod_mpoly_t A, const nmod_mpoly_t B, const fmpz_struct * shift, const fmpz_struct * stride, const nmod_mpoly_ctx_t ctx)
    void nmod_mpoly_inflate(nmod_mpoly_t A, const nmod_mpoly_t B, const fmpz_struct * shift, const fmpz_struct * stride, const nmod_mpoly_ctx_t ctx)
