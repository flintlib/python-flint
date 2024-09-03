from flint.flintlib.flint cimport ulong, slong, fmpz_struct
from flint.flintlib.arf_types cimport arf_struct, arf_ptr

cdef extern from "flint/arb_types.h":
    ctypedef struct mag_struct:
        fmpz_struct exp
        ulong man

    ctypedef mag_struct mag_t[1]
    ctypedef mag_struct * mag_ptr
    ctypedef const mag_struct * mag_srcptr

    ctypedef struct arb_struct:
        arf_struct mid
        mag_struct rad

    ctypedef arb_struct arb_t[1]
    ctypedef arb_struct * arb_ptr
    ctypedef const arb_struct * arb_srcptr

    # Macros
    cdef arf_ptr arb_midref(const arb_t x)
    cdef mag_ptr arb_radref(const arb_t x)

    cdef const ulong ARB_STR_MORE
    cdef const ulong ARB_STR_NO_RADIUS
    cdef const ulong ARB_STR_CONDENSE

    ctypedef struct arb_mat_struct:
        arb_ptr entries
        slong r
        slong c
        arb_ptr * rows

    ctypedef arb_mat_struct arb_mat_t[1]

    ctypedef struct arb_poly_struct:
        arb_ptr coeffs
        slong alloc
        slong length

    ctypedef arb_poly_struct arb_poly_t[1]
