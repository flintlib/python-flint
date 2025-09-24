from flint.flintlib.types.flint cimport fmpz_t

cdef extern from "flint/qfb.h":

    ctypedef struct qfb_struct:
        fmpz_t a
        fmpz_t b
        fmpz_t c

    ctypedef qfb_struct qfb_t[1]
