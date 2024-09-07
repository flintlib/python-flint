cimport cython

from flint.flintlib.types.gr cimport (
    truth_t,
    T_TRUE,
    T_FALSE,
    T_UNKNOWN,

    GR_SUCCESS,
    GR_DOMAIN,
    GR_UNABLE,

    gr_ctx_t,
    gr_ptr,
)
from flint.flintlib.functions.gr cimport (
    gr_heap_init,
    gr_set_str,
    gr_get_str,

    gr_neg,
    gr_add,
    gr_add_si,
    gr_sub,
    gr_sub_si,
    gr_mul,
    gr_mul_si,
    gr_inv,
    gr_div,
    gr_div_si,
    gr_pow_si,
)

from flint.flint_base.flint_base cimport flint_ctx, flint_scalar

from flint.utils.flint_exceptions import DomainError, UnableError


cdef inline truth_to_py(truth_t t):
    if t == T_TRUE:
        return True
    elif t == T_FALSE:
        return False
    else:
        return None


@cython.no_gc
cdef class gr_ctx(flint_ctx):
    cdef gr_ctx_t ctx_t
    cdef bint _init

    @cython.final
    cdef inline gr new_gr(self):
        cdef gr py_val
        cdef gr_ptr pval
        pval = gr_heap_init(self.ctx_t)
        if pval == NULL:
            raise MemoryError("Failed to allocate memory for gr object")
        py_val = gr.__new__(gr)
        py_val.pval = pval
        py_val.ctx = self
        py_val._init = True
        return py_val

    @cython.final
    cdef inline gr from_str(self, s: str):
        cdef gr py_val
        cdef bytes b
        b = s.encode('utf-8')
        py_val = self.new_gr()
        err = gr_set_str(py_val.pval, b, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Failed to parse string")
        return py_val

    @cython.final
    cdef inline str to_str(self, val: gr):
        cdef str py_str
        cdef char *s
        err = gr_get_str(&s, val.pval, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Failed to convert to string")
        py_str = (<bytes>s).decode("utf-8")
        # flint_free(s)
        return py_str

    @cython.final
    cdef inline _error(self, int err, str msg):
        if err & GR_DOMAIN:
            return DomainError(msg)
        elif err & GR_UNABLE:
            return UnableError(msg)
        else:
            return AssertionError("Bad error code")


@cython.no_gc
cdef class _gr_fmpz_ctx(gr_ctx):

    @staticmethod
    cdef _gr_fmpz_ctx _new()


@cython.no_gc
cdef class _gr_fmpq_ctx(gr_ctx):

    @staticmethod
    cdef _gr_fmpq_ctx _new()


# The global contexts for use in cython code:
cdef _gr_fmpz_ctx gr_fmpz_ctx_c
cdef _gr_fmpq_ctx gr_fmpq_ctx_c


@cython.no_gc
cdef class gr(flint_scalar):
    cdef gr_ptr pval
    cdef gr_ctx ctx
    cdef bint _init

    @cython.final
    cdef inline _error(self, int err, str msg):
        return self.ctx._error(err, msg)
