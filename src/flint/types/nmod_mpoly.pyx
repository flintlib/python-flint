from flint.flint_base.flint_base cimport (
    flint_mpoly,
    flint_mpoly_context,
    ordering_py_to_c,
    ordering_c_to_py,
)

from flint.utils.typecheck cimport typecheck
from flint.utils.flint_exceptions import DomainError, IncompatibleContextError

from flint.types.fmpz cimport any_as_fmpz, fmpz
from flint.types.fmpz_vec cimport fmpz_vec

from cpython.object cimport Py_EQ, Py_NE
cimport libc.stdlib

cdef dict _nmod_mpoly_ctx_cache = {}

cdef class nmod_mpoly_ctx(flint_mpoly_context):
    pass

cdef class nmod_mpoly(flint_mpoly):
    pass

cdef class nmod_mpoly_vec:
    pass
