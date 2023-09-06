from flint.flintlib.acb cimport acb_t, acb_ptr
from flint.flintlib.mag cimport mag_t

cdef extern from "flint/acb_calc.h":
    ctypedef int (*acb_calc_func_t)(acb_ptr out, const acb_t inp, void * param, long order, long prec)

    ctypedef struct acb_calc_integrate_opt_struct:
        long deg_limit
        long eval_limit
        long depth_limit
        int use_heap
        int verbose

    ctypedef acb_calc_integrate_opt_struct acb_calc_integrate_opt_t[1]

    void acb_calc_integrate_opt_init(acb_calc_integrate_opt_t options)

    int acb_calc_integrate(acb_t res, acb_calc_func_t f, void * param,
        const acb_t a, const acb_t b,
        long goal, const mag_t tol,
        const acb_calc_integrate_opt_t options,
        long prec)
