from flint.flintlib.types.flint cimport slong
from flint.flintlib.functions.fmpz cimport fmpz_t
from flint.flintlib.functions.fmpz_mod_mat cimport fmpz_mod_mat_t

from flint.flint_base.flint_base cimport flint_mat
from flint.types.fmpz_mod cimport fmpz_mod_ctx, fmpz_mod


cdef class fmpz_mod_mat(flint_mat):
    cdef fmpz_mod_mat_t val
    cdef bint _initialized
    cdef fmpz_mod_ctx ctx

    cdef void _init_empty(self, slong m, slong n, list args)
    cdef void _init_empty_ctx(self, slong m, slong n, fmpz_mod_ctx ctx)
    cdef void _init_from_list(self, slong m, slong n, list entries, list args)
    # cdef void _init_from_matrix(self, flint_mat M, list args)
    cdef void _init_from_matrix(self, M, list args)
    cdef fmpz_mod_ctx _parse_args(self, list args)
    cdef fmpz_mod_mat _new(self, slong m, slong n, fmpz_mod_ctx ctx)
    cdef fmpz_mod_mat _newlike(self)
    cdef fmpz_mod_mat _copy(self)

    cpdef slong nrows(self)
    cpdef slong ncols(self)
    cdef fmpz_mod _getitem(self, slong i, slong j)
    cdef void _setitem(self, slong i, slong j, fmpz_t e)
