from flint.flintlib.acb cimport acb_ptr, acb_srcptr

cdef extern from "flint/acb_dft.h":
    void acb_dft(acb_ptr w, acb_srcptr v, long n, long prec)
    void acb_dft_inverse(acb_ptr w, acb_srcptr v, long n, long prec)
