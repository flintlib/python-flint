from flint.flintlib.types.acb cimport acb_t
from flint.flintlib.types.acb_calc cimport acb_calc_func_t, acb_calc_integrate_opt_t
from flint.flintlib.types.arb cimport arb_t, mag_t
from flint.flintlib.types.arf cimport arf_t
from flint.flintlib.types.flint cimport slong



cdef extern from "flint/acb_calc.h":
    int acb_calc_integrate(acb_t res, acb_calc_func_t func, void * param, const acb_t a, const acb_t b, slong rel_goal, const mag_t abs_tol, const acb_calc_integrate_opt_t options, slong prec)
    void acb_calc_integrate_opt_init(acb_calc_integrate_opt_t options)
    int acb_calc_integrate_gl_auto_deg(acb_t res, slong * num_eval, acb_calc_func_t func, void * param, const acb_t a, const acb_t b, const mag_t tol, slong deg_limit, int flags, slong prec)
    void acb_calc_cauchy_bound(arb_t bound, acb_calc_func_t func, void * param, const acb_t x, const arb_t radius, slong maxdepth, slong prec)
    int acb_calc_integrate_taylor(acb_t res, acb_calc_func_t func, void * param, const acb_t a, const acb_t b, const arf_t inner_radius, const arf_t outer_radius, slong accuracy_goal, slong prec)
