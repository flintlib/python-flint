from flint.flint_base.flint_context cimport FlintContext

cdef FlintContext thectx = FlintContext()

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
