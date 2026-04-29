from flint.flintlib.types.arf cimport (
    arf_rnd_t,
)

cdef class FlintContext:
    cdef public bint pretty
    cdef public long _prec
    cdef public long _dps
    cdef arf_rnd_t rnd
    cdef public bint unicode
    cdef public long _cap

cdef FlintContext thectx

cdef inline long getprec(long prec=0):
    if prec > 1:
        return prec
    else:
        return thectx._prec

cdef inline long getcap(long cap=-1):
    if cap >= 0:
        return cap
    else:
        return thectx._cap
