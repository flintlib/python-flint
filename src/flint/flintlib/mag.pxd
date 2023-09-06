from flint.flintlib.flint cimport ulong, mp_limb_t
from flint.flintlib.fmpz cimport fmpz_struct

cdef extern from "mag.h":
    ctypedef struct mag_struct:
        fmpz_struct exp
        mp_limb_t man
    ctypedef mag_struct mag_t[1]
    ctypedef mag_struct * mag_ptr
    ctypedef const mag_struct * mag_srcptr

    void mag_init(mag_t x)
    void mag_clear(mag_t x)
    void mag_zero(mag_t x)
    void mag_set(mag_t x, const mag_t y)
    void mag_set_ui_2exp_si(mag_t x, ulong v, long e)
    void mag_hypot(mag_t x, const mag_t y, const mag_t z)
