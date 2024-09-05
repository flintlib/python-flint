from flint.flintlib.functions.dirichlet cimport dirichlet_group_t
from flint.flintlib.functions.dirichlet cimport dirichlet_char_t

cdef class dirichlet_group(object):
    cdef dirichlet_group_t val
    cdef int _init

    cpdef long size(self)

cdef class dirichlet_char(object):

    cdef dirichlet_char_t val
    cdef dirichlet_group G
