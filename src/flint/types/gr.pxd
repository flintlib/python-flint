cimport cython

from flint.flintlib.types.flint cimport (
    slong,
)
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
    gr_vec_t,
)
from flint.flintlib.functions.gr cimport (
    gr_heap_init,
    gr_set_str,
    gr_get_str,
    gr_set,
    gr_get_si,

    gr_zero,
    gr_one,
    gr_i,

    gr_is_zero,
    gr_is_one,
    gr_is_neg_one,
    # gr_is_integer,
    # gr_is_rational,

    gr_equal,

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

    gr_is_square,
    gr_sqrt,
    gr_rsqrt,

    gr_gcd,
    gr_lcm,
    gr_factor,

    gr_floor,
    gr_ceil,
    gr_trunc,
    gr_nint,

    gr_abs,
    gr_conj,
    gr_re,
    gr_im,
    gr_sgn,
    gr_csgn,
    gr_arg,

    # gr_le,
    # gr_lt,
    # gr_ge,
    # gr_gt,

    gr_numerator,
    gr_denominator,
)
from flint.flintlib.functions.gr_vec cimport (
    gr_vec_init,
    gr_vec_clear,
    gr_vec_length,
    gr_vec_entry_ptr,
)

from flint.flint_base.flint_base cimport flint_ctx, flint_scalar

# XXX: Can't import from Python modules in a .pxd file
# from flint.utils.flint_exceptions import DomainError, UnableError, UknownError


cdef inline truth_to_py(truth_t t):
    if t == T_TRUE:
        return True
    elif t == T_FALSE:
        return False
    else:
        return None


cdef inline truth_to_bool(truth_t t):
    if t == T_TRUE:
        return True
    elif t == T_FALSE:
        return False
    else:
        # raise UnknownError("Unknown truth value")
        raise AssertionError("Unknown truth value")


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
        # XXX: Is this a memory leak?
        # flint_free(s)
        return py_str

    @cython.final
    cdef inline _error(self, int err, str msg):
        if err & GR_DOMAIN:
            # return DomainError(msg)
            return AssertionError(msg)
        elif err & GR_UNABLE:
            # return UnableError(msg)
            return AssertionError(msg)
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


@cython.no_gc
cdef class _gr_fmpzi_ctx(gr_ctx):

    @staticmethod
    cdef _gr_fmpzi_ctx _new()


# The global contexts for use in cython code:
cdef _gr_fmpz_ctx gr_fmpz_ctx_c
cdef _gr_fmpq_ctx gr_fmpq_ctx_c
cdef _gr_fmpzi_ctx gr_fmpzi_ctx_c


@cython.no_gc
cdef class gr(flint_scalar):
    cdef gr_ptr pval
    cdef gr_ctx ctx
    cdef bint _init

    @cython.final
    cdef inline _error(self, int err, str msg):
        return self.ctx._error(err, msg)

    @cython.final
    cdef inline truth_t _equal(self, other: gr):
        return gr_equal(self.pval, other.pval, self.ctx.ctx_t)

    @cython.final
    cdef inline truth_t _is_zero(self):
        return gr_is_zero(self.pval, self.ctx.ctx_t)

    @cython.final
    cdef inline truth_t _is_one(self):
        return gr_is_one(self.pval, self.ctx.ctx_t)

    @cython.final
    cdef inline truth_t _is_neg_one(self):
        return gr_is_neg_one(self.pval, self.ctx.ctx_t)

    # @cython.final
    # cdef inline truth_t _is_integer(self):
    #     return gr_is_integer(self.pval, self.ctx.ctx_t)

    # @cython.final
    # cdef inline truth_t _is_rational(self):
    #     return gr_is_rational(self.pval, self.ctx.ctx_t)

    @cython.final
    cdef inline gr _neg(self):
        cdef int err
        cdef gr res = self.ctx.new_gr()
        err = gr_neg(res.pval, self.pval, self.ctx.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Failed to negate gr object")
        return res

    @cython.final
    cdef inline gr _add(self, other: gr):
        cdef int err
        cdef gr res = self.ctx.new_gr()
        err = gr_add(res.pval, self.pval, other.pval, self.ctx.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Failed to add gr objects")
        return res

    @cython.final
    cdef inline gr _add_si(self, other: slong):
        cdef int err
        cdef gr res = self.ctx.new_gr()
        err = gr_add_si(res.pval, self.pval, other, self.ctx.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Failed to add gr object and int")
        return res

    @cython.final
    cdef inline gr _sub(self, other: gr):
        cdef int err
        cdef gr res = self.ctx.new_gr()
        err = gr_sub(res.pval, self.pval, other.pval, self.ctx.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Failed to subtract gr objects")
        return res

    @cython.final
    cdef inline gr _sub_si(self, other: int):
        cdef int err
        cdef gr res = self.ctx.new_gr()
        err = gr_sub_si(res.pval, self.pval, other, self.ctx.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Failed to subtract gr object and int")
        return res

    @cython.final
    cdef inline gr _mul(self, other: gr):
        cdef int err
        cdef gr res = self.ctx.new_gr()
        err = gr_mul(res.pval, self.pval, other.pval, self.ctx.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Failed to multiply gr objects")
        return res

    @cython.final
    cdef inline gr _mul_si(self, other: slong):
        cdef int err
        cdef gr res = self.ctx.new_gr()
        err = gr_mul_si(res.pval, self.pval, other, self.ctx.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Failed to multiply gr object and int")
        return res

    @cython.final
    cdef inline gr _inv(self):
        cdef int err
        cdef gr res = self.ctx.new_gr()
        err = gr_inv(res.pval, self.pval, self.ctx.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Failed to invert gr object")
        return res

    @cython.final
    cdef inline gr _div(self, other: gr):
        cdef int err
        cdef gr res = self.ctx.new_gr()
        err = gr_div(res.pval, self.pval, other.pval, self.ctx.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Failed to divide gr objects")
        return res

    @cython.final
    cdef inline gr _div_si(self, other: slong):
        cdef int err
        cdef gr res = self.ctx.new_gr()
        err = gr_div_si(res.pval, self.pval, other, self.ctx.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Failed to divide gr object and int")
        return res

    @cython.final
    cdef inline gr _pow_si(self, other: slong):
        cdef int err
        cdef gr res = self.ctx.new_gr()
        err = gr_pow_si(res.pval, self.pval, other, self.ctx.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Failed to raise gr object to an integer power")
        return res

    @cython.final
    cdef inline truth_t _is_square(self):
        return gr_is_square(self.pval, self.ctx.ctx_t)

    @cython.final
    cdef inline gr _sqrt(self):
        cdef int err
        cdef gr res = self.ctx.new_gr()
        err = gr_sqrt(res.pval, self.pval, self.ctx.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Failed to compute square root of gr object")
        return res

    @cython.final
    cdef inline gr _rsqrt(self):
        cdef int err
        cdef gr res = self.ctx.new_gr()
        err = gr_rsqrt(res.pval, self.pval, self.ctx.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Failed to compute reciprocal square root of gr object")
        return res

    @cython.final
    cdef inline gr _gcd(self, other: gr):
        cdef int err
        cdef gr res = self.ctx.new_gr()
        err = gr_gcd(res.pval, self.pval, other.pval, self.ctx.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Failed to compute gcd of gr objects")
        return res

    @cython.final
    cdef inline gr _lcm(self, other: gr):
        cdef int err
        cdef gr res = self.ctx.new_gr()
        err = gr_lcm(res.pval, self.pval, other.pval, self.ctx.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Failed to compute lcm of gr objects")
        return res

    # XXX: Needs gr_vec_t
    cdef inline _factor(self):
        cdef int err, i
        cdef slong length, exp_s
        cdef gr c, f
        cdef gr_ptr fac_i, exp_i
        cdef gr_vec_t factors, exponents
        cdef int flags = 0  # XXX: What is flags?
        c = self.ctx.new_gr()
        gr_vec_init(factors, 0, self.ctx.ctx_t)
        gr_vec_init(exponents, 0, gr_fmpz_ctx_c.ctx_t)
        err = gr_factor(c.pval, factors, exponents, self.pval, flags, self.ctx.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Failed to factor gr object")
        length = gr_vec_length(factors, self.ctx.ctx_t)
        py_factors = [None] * length
        for 0 <= i < length:
            f = self.ctx.new_gr()
            fac_i = gr_vec_entry_ptr(factors, i, self.ctx.ctx_t)
            err = gr_set(f.pval, fac_i, self.ctx.ctx_t)
            if err != GR_SUCCESS:
                raise self._error(err, "Failed to copy factor.")
            exp_i = gr_vec_entry_ptr(exponents, i, gr_fmpz_ctx_c.ctx_t)
            err = gr_get_si(&exp_s, exp_i, gr_fmpz_ctx_c.ctx_t)
            if err != GR_SUCCESS:
                raise self._error(err, "Failed to get integer value of exponent.")
            exp = exp_s
            py_factors[i] = (f, exp)
        gr_vec_clear(factors, self.ctx.ctx_t)
        gr_vec_clear(exponents, gr_fmpz_ctx_c.ctx_t)
        return c, py_factors

    cdef inline gr _numerator(self):
        cdef int err
        cdef gr res = self.ctx.new_gr()
        err = gr_numerator(res.pval, self.pval, self.ctx.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Failed to compute numerator of gr object")
        return res

    cdef inline gr _denominator(self):
        cdef int err
        cdef gr res = self.ctx.new_gr()
        err = gr_denominator(res.pval, self.pval, self.ctx.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Failed to compute denominator of gr object")
        return res

    cdef inline gr _floor(self):
        cdef int err
        cdef gr res = self.ctx.new_gr()
        err = gr_floor(res.pval, self.pval, self.ctx.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Failed to compute floor of gr object")
        return res

    cdef inline gr _ceil(self):
        cdef int err
        cdef gr res = self.ctx.new_gr()
        err = gr_ceil(res.pval, self.pval, self.ctx.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Failed to compute ceil of gr object")
        return res

    cdef inline gr _trunc(self):
        cdef int err
        cdef gr res = self.ctx.new_gr()
        err = gr_trunc(res.pval, self.pval, self.ctx.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Failed to compute trunc of gr object")
        return res

    cdef inline gr _nint(self):
        cdef int err
        cdef gr res = self.ctx.new_gr()
        err = gr_nint(res.pval, self.pval, self.ctx.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Failed to compute nint of gr object")
        return res

    cdef inline gr _abs(self):
        cdef int err
        cdef gr res = self.ctx.new_gr()
        err = gr_abs(res.pval, self.pval, self.ctx.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Failed to compute absolute value of gr object")
        return res

    cdef inline gr _conj(self):
        cdef int err
        cdef gr res = self.ctx.new_gr()
        err = gr_conj(res.pval, self.pval, self.ctx.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Failed to compute conjugate of gr object")
        return res

    cdef inline gr _re(self):
        cdef int err
        cdef gr res = self.ctx.new_gr()
        err = gr_re(res.pval, self.pval, self.ctx.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Failed to compute real part of gr object")
        return res

    cdef inline gr _im(self):
        cdef int err
        cdef gr res = self.ctx.new_gr()
        err = gr_im(res.pval, self.pval, self.ctx.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Failed to compute imaginary part of gr object")
        return res

    cdef inline gr _sgn(self):
        cdef int err
        cdef gr res = self.ctx.new_gr()
        err = gr_sgn(res.pval, self.pval, self.ctx.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Failed to compute sign of gr object")
        return res

    cdef inline gr _csgn(self):
        cdef int err
        cdef gr res = self.ctx.new_gr()
        err = gr_csgn(res.pval, self.pval, self.ctx.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Failed to compute complex sign of gr object")
        return res

    cdef inline gr _arg(self):
        cdef int err
        cdef gr res = self.ctx.new_gr()
        err = gr_arg(res.pval, self.pval, self.ctx.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Failed to compute argument of gr object")
        return res

    # cdef inline truth_t _le(self, other: gr):
    #     return gr_le(self.pval, other.pval, self.ctx.ctx_t)

    # cdef inline truth_t _lt(self, other: gr):
    #     return gr_lt(self.pval, other.pval, self.ctx.ctx_t)

    # cdef inline truth_t _ge(self, other: gr):
    #     return gr_ge(self.pval, other.pval, self.ctx.ctx_t)

    # cdef inline truth_t _gt(self, other: gr):
    #     return gr_gt(self.pval, other.pval, self.ctx.ctx_t)

    # cdef inline truth_t _abs_le(self, other: gr):
    #    return gr_abs_le(self.pval, other.pval, self.ctx.ctx_t)

    # cdef inline truth_t _abs_lt(self, other: gr):
    #     return gr_abs_lt(self.pval, other.pval, self.ctx.ctx_t)

    # cdef inline truth_t _abs_ge(self, other: gr):
    #     return gr_abs_ge(self.pval, other.pval, self.ctx.ctx_t)

    # cdef inline truth_t _abs_gt(self, other: gr):
    #     return gr_abs_gt(self.pval, other.pval, self.ctx.ctx_t)

    # cdef inline gr _min(self, other: gr):
    #     cdef int err
    #     cdef gr res = self.ctx.new_gr()
    #     err = gr_min(res.pval, self.pval, other.pval, self.ctx.ctx_t)
    #     if err != GR_SUCCESS:
    #         raise self._error(err, "Failed to compute minimum of gr objects")
    #     return res

    # cdef inline gr _max(self, other: gr):
    #     cdef int err
    #     cdef gr res = self.ctx.new_gr()
    #     err = gr_max(res.pval, self.pval, other.pval, self.ctx.ctx_t)
    #     if err != GR_SUCCESS:
    #         raise self._error(err, "Failed to compute maximum of gr objects")
    #     return res
