from flint.flintlib.types.flint cimport slong
from flint.flintlib.types.acb cimport acb_ptr, acb_srcptr, acb_t, acb_mat_t

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
