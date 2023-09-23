from flint.flintlib.acb cimport acb_t, acb_ptr
from flint.flintlib.mag cimport mag_t
from flint.flintlib.flint cimport slong
from flint.flintlib.arb cimport arb_t
from flint.flintlib.arf cimport arf_t

cdef extern from "flint/acb_calc.h":
    ctypedef int (*acb_calc_func_t)(acb_ptr out, const acb_t inp, void * param, long order, long prec)

    ctypedef struct acb_calc_integrate_opt_struct:
        long deg_limit
        long eval_limit
        long depth_limit
        int use_heap
        int verbose

    ctypedef acb_calc_integrate_opt_struct acb_calc_integrate_opt_t[1]

# from /Users/davideinstein/projects/arb/doc/source/acb_calc.rst
    int acb_calc_integrate(acb_t res, acb_calc_func_t func, void * param, const acb_t a, const acb_t b, slong rel_goal, const mag_t abs_tol, const acb_calc_integrate_opt_t options, slong prec)
    void acb_calc_integrate_opt_init(acb_calc_integrate_opt_t options)
    int acb_calc_integrate_gl_auto_deg(acb_t res, slong * num_eval, acb_calc_func_t func, void * param, const acb_t a, const acb_t b, const mag_t tol, slong deg_limit, int flags, slong prec)
    void acb_calc_cauchy_bound(arb_t bound, acb_calc_func_t func, void * param, const acb_t x, const arb_t radius, slong maxdepth, slong prec)
    int acb_calc_integrate_taylor(acb_t res, acb_calc_func_t func, void * param, const acb_t a, const acb_t b, const arf_t inner_radius, const arf_t outer_radius, slong accuracy_goal, slong prec)

