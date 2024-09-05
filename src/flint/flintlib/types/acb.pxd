from flint.flintlib.types.flint cimport slong, ulong
from flint.flintlib.types.arb cimport arb_struct, arb_ptr, arb_poly_t, mag_struct


cdef extern from "flint/acb_types.h":

    ctypedef struct acb_struct:
        arb_struct real
        arb_struct imag

    ctypedef acb_struct acb_t[1]
    ctypedef acb_struct * acb_ptr
    ctypedef const acb_struct * acb_srcptr

    ctypedef struct acb_mat_struct:
        acb_ptr entries
        slong r
        slong c
        acb_ptr * rows

    ctypedef acb_mat_struct acb_mat_t[1]

    ctypedef struct acb_poly_struct:
        acb_ptr coeffs
        slong alloc
        slong length

    ctypedef acb_poly_struct acb_poly_t[1]


cdef extern from "flint/acb.h":

    arb_ptr acb_realref(const acb_t x)
    arb_ptr acb_imagref(const acb_t x)


cdef extern from "flint/acb_mat.h":

    acb_struct * acb_mat_entry(acb_mat_t mat, slong i, slong j)
    slong acb_mat_nrows(const acb_mat_t x)
    slong acb_mat_ncols(const acb_mat_t x)


cdef extern from "flint/acb_poly.h":

    acb_ptr acb_poly_get_coeff_ptr(arb_poly_t poly, slong n)
