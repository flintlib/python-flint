from flint.flintlib.types.flint cimport slong
from flint.flintlib.types.acb cimport acb_ptr, acb_t

cdef extern from "flint/acb_calc.h":

    ctypedef int (*acb_calc_func_t)(acb_ptr out, const acb_t inp, void * param, slong order, slong prec)

    ctypedef struct acb_calc_integrate_opt_struct:
        slong deg_limit
        slong eval_limit
        slong depth_limit
        int use_heap
        int verbose

    ctypedef acb_calc_integrate_opt_struct acb_calc_integrate_opt_t[1]
