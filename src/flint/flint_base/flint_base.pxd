from flint.flintlib.types.mpoly cimport ordering_t
from flint.flintlib.types.flint cimport slong

cdef class flint_ctx:
    pass

cdef class flint_elem:
    pass

cdef class flint_scalar(flint_elem):
    pass

cdef class flint_poly(flint_elem):
    pass

cdef class flint_mpoly_context(flint_elem):
    cdef public object py_names
    cdef const char ** c_names

cdef class flint_mod_mpoly_context(flint_mpoly_context):
    cdef readonly bint __prime_modulus

cdef class flint_mpoly(flint_elem):
    cdef _add_scalar_(self, other)
    cdef _sub_scalar_(self, other)
    cdef _mul_scalar_(self, other)

    cdef _add_mpoly_(self, other)
    cdef _sub_mpoly_(self, other)
    cdef _mul_mpoly_(self, other)

    cdef _divmod_mpoly_(self, other)
    cdef _truediv_scalar_(self, other)
    cdef _divexact_scalar_(self, other)
    cdef _floordiv_mpoly_(self, other)
    cdef _truediv_mpoly_(self, other)
    cdef _mod_mpoly_(self, other)

    cdef _rsub_scalar_(self, other)
    cdef _rsub_mpoly_(self, other)

    cdef _rdivmod_mpoly_(self, other)
    cdef _rfloordiv_mpoly_(self, other)
    cdef _rtruediv_mpoly_(self, other)
    cdef _rmod_mpoly_(self, other)

    cdef _pow_(self, other)

    cdef _iadd_scalar_(self, other)
    cdef _isub_scalar_(self, other)
    cdef _imul_scalar_(self, other)

    cdef _iadd_mpoly_(self, other)
    cdef _isub_mpoly_(self, other)
    cdef _imul_mpoly_(self, other)

    cdef _compose_gens_(self, ctx, slong *mapping)

cdef class flint_mat(flint_elem):
    pass

cdef class flint_series(flint_elem):
    pass

cdef ordering_t ordering_py_to_c(ordering)
cdef ordering_c_to_py(ordering_t ordering)
