from flint.flintlib.mpoly cimport ordering_t

cdef class flint_elem:
    pass

cdef class flint_scalar(flint_elem):
    pass

cdef class flint_poly(flint_elem):
    pass

cdef class flint_mpoly_context(flint_elem):
    cdef public object py_names
    cdef const char ** c_names

cdef class flint_mpoly(flint_elem):
    pass

cdef class flint_mat(flint_elem):
    pass

cdef class flint_series(flint_elem):
    pass

cdef class flint_rational_function(flint_elem):
    pass

cpdef enum Ordering:
    lex, deglex, degrevlex

cdef ordering_t ordering_py_to_c(ordering: Ordering)
cdef ordering_c_to_py(ordering_t ordering)
