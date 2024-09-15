from flint.flintlib.types.flint cimport slong
from flint.flintlib.types.gr cimport gr_ctx_t, gr_funcptr, truth_t

# unknown type gr_method_tab_input


cdef extern from "flint/gr_implementing.h":
    # void gr_method_tab_init(gr_funcptr * methods, gr_method_tab_input * tab)
    int gr_not_implemented(void)
    int gr_not_in_domain(void)
    truth_t gr_generic_ctx_predicate(gr_ctx_t ctx)
    truth_t gr_generic_ctx_predicate_true(gr_ctx_t ctx)
    truth_t gr_generic_ctx_predicate_false(gr_ctx_t ctx)
    void gr_test_ring(gr_ctx_t R, slong iters, int test_flags)
