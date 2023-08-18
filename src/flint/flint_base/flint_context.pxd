from flint._flint cimport (
    arf_rnd_t,
)

cdef class FlintContext:
    cdef public bint pretty
    cdef public long _prec
    cdef public long _dps
    cdef arf_rnd_t rnd
    cdef public bint unicode
    cdef public long _cap

    pass
