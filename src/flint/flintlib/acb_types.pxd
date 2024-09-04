from flint.flintlib.flint cimport slong, ulong
from flint.flintlib.arb_types cimport arb_struct, arb_ptr, arb_poly_t, mag_struct


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


cdef extern from "flint/acb_calc.h":

    ctypedef int (*acb_calc_func_t)(acb_ptr out, const acb_t inp, void * param, slong order, slong prec)

    ctypedef struct acb_calc_integrate_opt_struct:
        slong deg_limit
        slong eval_limit
        slong depth_limit
        int use_heap
        int verbose

    ctypedef acb_calc_integrate_opt_struct acb_calc_integrate_opt_t[1]


cdef extern from "flint/acb_dirichlet.h":

    ctypedef struct acb_dirichlet_roots_struct:
        ulong order
        ulong reduced_order
        acb_t z
        slong size
        slong depth
        acb_ptr * Z
        int use_pow

    ctypedef acb_dirichlet_roots_struct acb_dirichlet_roots_t[1]

    ctypedef struct acb_dirichlet_hurwitz_precomp_struct:
        acb_struct s
        mag_struct err
        acb_ptr coeffs
        int deflate
        slong A
        slong N
        slong K

    ctypedef acb_dirichlet_hurwitz_precomp_struct acb_dirichlet_hurwitz_precomp_t[1]


cdef extern from "flint/acb_theta.h":

    ctypedef struct acb_theta_eld_struct:
        slong dim
        slong ambient_dim
        slong * last_coords
        slong min
        slong mid
        slong max
        slong nr
        slong nl
        acb_theta_eld_struct * rchildren
        acb_theta_eld_struct * lchildren
        slong nb_pts, nb_border
        slong * box

    ctypedef acb_theta_eld_struct acb_theta_eld_t[1]
    ctypedef void (*acb_theta_naive_worker_t)(acb_ptr, acb_srcptr, acb_srcptr, const slong *, slong, const acb_t, const slong *, slong, slong, slong, slong)
    ctypedef int (*acb_theta_ql_worker_t)(acb_ptr, acb_srcptr, acb_srcptr, arb_srcptr, arb_srcptr, const acb_mat_t, slong, slong)
