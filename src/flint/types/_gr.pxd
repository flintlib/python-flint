cimport cython
cimport libc.stdlib

from flint.flintlib.types.flint cimport (
    slong,
    ulong,
    fmpz_t,
)
from flint.flintlib.types.fmpz cimport (
    fmpz_poly_t,
)
from flint.flintlib.types.fmpq cimport (
    fmpq_poly_t,
)
from flint.flintlib.types.mpoly cimport (
    ordering_t,
)
from flint.flintlib.functions.ulong_extras cimport (
    n_is_prime,
)
from flint.flintlib.functions.fmpz cimport (
    fmpz_init_set,
)
from flint.flintlib.functions.fmpz_poly cimport (
    fmpz_poly_init,
    fmpz_poly_set,
)
from flint.flintlib.functions.fmpq_poly cimport (
    fmpq_poly_init,
    fmpq_poly_set,
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
from flint.flintlib.functions.gr_domains cimport (
    gr_ctx_init_fmpz,
    gr_ctx_init_fmpq,
    gr_ctx_init_fmpzi,
    gr_ctx_init_fexpr,

    gr_ctx_init_nmod,
    gr_ctx_init_fmpz_mod,
    # gr_ctx_set_is_field,

    gr_ctx_init_fq,
    gr_ctx_init_fq_nmod,
    gr_ctx_init_fq_zech,

    gr_ctx_init_nf,
    gr_ctx_init_nf_fmpz_poly,

    gr_ctx_init_real_qqbar,
    gr_ctx_init_complex_qqbar,
    # _gr_ctx_qqbar_set_limits,

    gr_ctx_init_real_ca,
    gr_ctx_init_complex_ca,
    gr_ctx_init_real_algebraic_ca,
    gr_ctx_init_complex_algebraic_ca,
    gr_ctx_init_complex_extended_ca,

    gr_ctx_init_real_float_arf,
    gr_ctx_init_complex_float_acf,

    gr_ctx_init_real_arb,
    gr_ctx_init_complex_acb,

    gr_ctx_init_gr_poly,
    gr_ctx_init_gr_mpoly,
    # gr_ctx_init_fmpz_mpoly_q,
    # gr_ctx_init_series_mod_gr_poly,
    gr_ctx_init_gr_series,
)
from flint.flintlib.functions.gr cimport (
    gr_heap_init,
    gr_set_d,
    gr_set_other,
    gr_set_str,
    gr_get_str,
    gr_set,
    gr_set_si,
    gr_get_si,

    gr_zero,
    gr_one,
    gr_gen,
    gr_gens,
    # gr_gens_recursive,
    gr_ctx_set_gen_names,

    gr_i,
    gr_pos_inf,
    gr_neg_inf,
    gr_uinf,
    gr_undefined,
    gr_unknown,

    gr_is_zero,
    gr_is_one,
    gr_is_neg_one,
    gr_is_integer,
    gr_is_rational,

    gr_equal,

    gr_neg,
    gr_add,
    gr_add_si,
    gr_add_other,
    gr_other_add,
    gr_sub,
    gr_sub_si,
    gr_sub_other,
    gr_other_sub,
    gr_mul,
    gr_mul_si,
    gr_mul_other,
    gr_other_mul,
    gr_is_invertible,
    gr_inv,
    gr_div,
    gr_div_si,
    gr_div_other,
    gr_div_nonunique,
    gr_divides,
    gr_other_div,
    gr_divexact,
    gr_divexact_si,
    gr_divexact_other,
    gr_other_divexact,
    gr_euclidean_div,
    gr_euclidean_rem,
    gr_euclidean_divrem,
    gr_pow,
    gr_pow_si,
    gr_pow_other,
    gr_other_pow,

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

    gr_cmp,
    gr_cmp_other,
    gr_cmpabs,
    gr_cmpabs_other,
    gr_le,
    gr_lt,
    gr_ge,
    gr_gt,
    gr_abs_le,
    gr_abs_lt,
    gr_abs_ge,
    gr_abs_gt,
    gr_min,
    gr_max,

    gr_numerator,
    gr_denominator,
)
from flint.flintlib.functions.gr_vec cimport (
    gr_vec_init,
    gr_vec_clear,
    gr_vec_length,
    gr_vec_entry_ptr,
)

from flint.flint_base.flint_base cimport (
    flint_ctx,
    flint_scalar,
    ordering_py_to_c,
)

from flint.types.fmpz cimport fmpz
from flint.types.fmpz_poly cimport fmpz_poly
from flint.types.fmpq_poly cimport fmpq_poly


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
    cdef inline gr from_d(self, double d):
        cdef gr py_val
        py_val = self.new_gr()
        err = gr_set_d(py_val.pval, d, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Incorrect conversion from a double")
        return py_val

    @cython.final
    cdef inline gr from_other(self, gr x):
        cdef gr py_val
        py_val = self.new_gr()
        err = gr_set_other(py_val.pval, x.pval, x.ctx.ctx_t, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot convert x to the current context")
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
    cdef inline gr from_si(self, slong n):
        cdef gr py_val
        py_val = self.new_gr()
        err = gr_set_si(py_val.pval, n, self.ctx_t)
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

    @cython.final
    cdef inline gr _zero(self):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_zero(res.pval, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot create zero")
        return res

    @cython.final
    cdef inline gr _one(self):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_one(res.pval, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot create one")
        return res

    @cython.final
    cdef inline gr _i(self):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_i(res.pval, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot create sqrt(-1)")
        return res

    @cython.final
    cdef inline gr _pos_inf(self):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_pos_inf(res.pval, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot create +inf")
        return res

    @cython.final
    cdef inline gr _neg_inf(self):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_neg_inf(res.pval, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot create -inf")
        return res

    @cython.final
    cdef inline gr _uinf(self):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_uinf(res.pval, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot create uinf")
        return res

    @cython.final
    cdef inline gr _undefined(self):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_undefined(res.pval, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot create undefined")
        return res

    @cython.final
    cdef inline gr _unknown(self):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_unknown(res.pval, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot create unknown")
        return res

    @cython.final
    cdef inline gr _gen(self):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_gen(res.pval, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot create generator")
        return res

    @cython.final
    cdef inline list _gens(self):
        cdef int err
        cdef gr g
        cdef gr_vec_t gens
        gr_vec_init(gens, 0, self.ctx_t)
        err = gr_gens(gens, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot get generators")
        length = gr_vec_length(gens, self.ctx_t)
        py_gens = [None] * length
        for 0 <= i < length:
            g = self.new_gr()
            err = gr_set(g.pval, gr_vec_entry_ptr(gens, i, self.ctx_t), self.ctx_t)
            if err != GR_SUCCESS:
                raise self._error(err, "Failed to copy generator.")
            py_gens[i] = g
        gr_vec_clear(gens, self.ctx_t)
        return py_gens

    @cython.final
    cdef inline truth_t _is_zero(self, gr x):
        return gr_is_zero(x.pval, self.ctx_t)

    @cython.final
    cdef inline truth_t _is_one(self, gr x):
        return gr_is_one(x.pval, self.ctx_t)

    @cython.final
    cdef inline truth_t _is_neg_one(self, gr x):
        return gr_is_neg_one(x.pval, self.ctx_t)

    @cython.final
    cdef inline truth_t _equal(self, gr x, gr y):
        return gr_equal(x.pval, y.pval, self.ctx_t)

    # @cython.final
    # cdef inline truth_t _is_integer(self, gr x):
    #     return gr_is_integer(x.pval, self.ctx_t)
    #
    # @cython.final
    # cdef inline truth_t _is_rational(self, gr x):
    #     return gr_is_rational(x.pval, self.ctx_t)

    @cython.final
    cdef inline gr _neg(self, gr x):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_neg(res.pval, x.pval, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot negate x in this context")
        return res

    @cython.final
    cdef inline gr _add(self, gr x, gr y):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_add(res.pval, x.pval, y.pval, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot add x and y in this context")
        return res

    @cython.final
    cdef inline gr _add_other(self, gr x, gr y):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_add_other(res.pval, x.pval, y.pval, y.ctx.ctx_t, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot add x and y in this context")
        return res

    @cython.final
    cdef inline gr _other_add(self, gr x, gr y):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_other_add(res.pval, x.pval, x.ctx.ctx_t, y.pval, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot add x and y in this context")
        return res

    @cython.final
    cdef inline gr _add_si(self, gr x, slong y):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_add_si(res.pval, x.pval, y, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot add x and y in this context")
        return res

    @cython.final
    cdef inline gr _sub(self, gr x, gr y):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_sub(res.pval, x.pval, y.pval, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot sub x and y in this context")
        return res

    @cython.final
    cdef inline gr _sub_other(self, gr x, gr y):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_sub_other(res.pval, x.pval, y.pval, y.ctx.ctx_t, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot sub x and y in this context")
        return res

    @cython.final
    cdef inline gr _other_sub(self, gr x, gr y):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_other_sub(res.pval, x.pval, x.ctx.ctx_t, y.pval, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot sub x and y in this context")
        return res

    @cython.final
    cdef inline gr _sub_si(self, gr x, slong y):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_sub_si(res.pval, x.pval, y, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot sub x and y in this context")
        return res

    @cython.final
    cdef inline gr _mul(self, gr x, gr y):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_mul(res.pval, x.pval, y.pval, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot mul x and y in this context")
        return res

    @cython.final
    cdef inline gr _mul_other(self, gr x, gr y):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_mul_other(res.pval, x.pval, y.pval, y.ctx.ctx_t, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot mul x and y in this context")
        return res

    @cython.final
    cdef inline gr _other_mul(self, gr x, gr y):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_other_mul(res.pval, x.pval, x.ctx.ctx_t, y.pval, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot mul x and y in this context")
        return res

    @cython.final
    cdef inline gr _mul_si(self, gr x, slong y):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_mul_si(res.pval, x.pval, y, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot mul x and y in this context")
        return res

    ###
    # Division

    @cython.final
    cdef inline truth_t _is_invertible(self, gr x):
        return gr_is_invertible(x.pval, self.ctx_t)

    @cython.final
    cdef inline gr _inv(self, gr x):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_inv(res.pval, x.pval, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "x in not invertible in this context")
        return res

    @cython.final
    cdef inline gr _div(self, gr x, gr y):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_div(res.pval, x.pval, y.pval, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot div x and y in this context")
        return res

    @cython.final
    cdef inline gr _div_other(self, gr x, gr y):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_div_other(res.pval, x.pval, y.pval, y.ctx.ctx_t, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot div x and y in this context")
        return res

    @cython.final
    cdef inline gr _other_div(self, gr x, gr y):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_other_div(res.pval, x.pval, x.ctx.ctx_t, y.pval, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot div x and y in this context")
        return res

    @cython.final
    cdef inline gr _div_si(self, gr x, slong y):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_div_si(res.pval, x.pval, y, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot div x and y in this context")
        return res

    @cython.final
    cdef inline gr _divexact(self, gr x, gr y):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_divexact(res.pval, x.pval, y.pval, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot divexact x and y in this context")
        return res

    @cython.final
    cdef inline gr _divexact_other(self, gr x, gr y):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_divexact_other(res.pval, x.pval, y.pval, y.ctx.ctx_t, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot divexact x and y in this context")
        return res

    @cython.final
    cdef inline gr _other_divexact(self, gr x, gr y):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_other_divexact(res.pval, x.pval, x.ctx.ctx_t, y.pval, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot divexact x and y in this context")
        return res

    @cython.final
    cdef inline gr _divexact_si(self, gr x, slong y):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_divexact_si(res.pval, x.pval, y, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot divexact x and y in this context")
        return res

    @cython.final
    cdef inline gr _div_nonunique(self, gr x, gr y):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_div_nonunique(res.pval, x.pval, y.pval, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "No solution q to x = qy has been found")
        return res

    @cython.final
    cdef inline truth_t _divides(self, gr x, gr d):
        return gr_divides(x.pval, d.pval, self.ctx_t)

    @cython.final
    cdef inline gr _euclidean_div(self, gr x, gr y):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_euclidean_div(res.pval, x.pval, y.pval, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "No solution q to x = qy has been found")
        return res

    @cython.final
    cdef inline gr _euclidean_rem(self, gr x, gr y):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_euclidean_rem(res.pval, x.pval, y.pval, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "No solution q to x = qy has been found")
        return res

    @cython.final
    cdef inline tuple[gr, gr] _euclidean_divrem(self, gr x, gr y):
        cdef int err
        cdef gr rem = self.new_gr()
        cdef gr div = self.new_gr()
        err = gr_euclidean_divrem(div.pval, rem.pval, x.pval, y.pval, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "No solution q to x = qy has been found")
        return (div, rem)

    ###
    # Powering

    @cython.final
    cdef inline gr _pow(self, gr x, gr y):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_pow(res.pval, x.pval, y.pval, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot pow x and y in this context")
        return res

    @cython.final
    cdef inline gr _pow_other(self, gr x, gr y):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_pow_other(res.pval, x.pval, y.pval, y.ctx.ctx_t, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot pow x and y in this context")
        return res

    @cython.final
    cdef inline gr _other_pow(self, gr x, gr y):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_other_pow(res.pval, x.pval, x.ctx.ctx_t, y.pval, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot pow x and y in this context")
        return res

    @cython.final
    cdef inline gr _pow_si(self, gr x, slong y):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_pow_si(res.pval, x.pval, y, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot pow x and y in this context")
        return res

    ###
    # Square roots

    @cython.final
    cdef inline truth_t _is_square(self, gr x):
        return gr_is_square(x.pval, self.ctx_t)

    @cython.final
    cdef inline gr _sqrt(self, gr x):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_sqrt(res.pval, x.pval, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot pow x and y in this context")
        return res

    @cython.final
    cdef inline gr _rsqrt(self, gr x):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_rsqrt(res.pval, x.pval, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot pow x and y in this context")
        return res

    ###
    # Greates Common Divisors

    @cython.final
    cdef inline gr _gcd(self, gr x, gr y):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_gcd(res.pval, x.pval, y.pval, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot pow x and y in this context")
        return res

    @cython.final
    cdef inline gr _lcm(self, gr x, gr y):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_lcm(res.pval, x.pval, y.pval, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot pow x and y in this context")
        return res

    ###
    # Factorization

    @cython.final
    cdef inline tuple[gr, list[tuple[gr, int]]] _factor(self, gr x):
        cdef int err, i
        cdef slong length, exp_s
        cdef gr c, f
        cdef gr_ptr fac_i, exp_i
        cdef gr_vec_t factors, exponents
        cdef int flags = 0  # XXX: What is flags?
        c = self.new_gr()
        gr_vec_init(factors, 0, self.ctx_t)
        gr_vec_init(exponents, 0, gr_fmpz_ctx_c.ctx_t)
        err = gr_factor(c.pval, factors, exponents, x.pval, flags, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Failed to factor gr object")
        length = gr_vec_length(factors, self.ctx_t)
        py_factors = [None] * length
        for 0 <= i < length:
            f = self.new_gr()
            fac_i = gr_vec_entry_ptr(factors, i, self.ctx_t)
            err = gr_set(f.pval, fac_i, self.ctx_t)
            if err != GR_SUCCESS:
                raise self._error(err, "Failed to copy factor.")
            exp_i = gr_vec_entry_ptr(exponents, i, gr_fmpz_ctx_c.ctx_t)
            err = gr_get_si(&exp_s, exp_i, gr_fmpz_ctx_c.ctx_t)
            if err != GR_SUCCESS:
                raise self._error(err, "Failed to get integer value of exponent.")
            exp = exp_s
            py_factors[i] = (f, exp)
        gr_vec_clear(factors, self.ctx_t)
        gr_vec_clear(exponents, gr_fmpz_ctx_c.ctx_t)
        return c, py_factors

    ###
    # Fractions

    @cython.final
    cdef inline gr _numerator(self, gr x):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_numerator(res.pval, x.pval, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot pow x and y in this context")
        return res

    @cython.final
    cdef inline gr _denominator(self, gr x):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_denominator(res.pval, x.pval, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot pow x and y in this context")
        return res

    ###
    # Integer and Complex parts

    @cython.final
    cdef inline gr _floor(self, gr x):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_floor(res.pval, x.pval, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot compute floor(x) in this context")
        return res

    @cython.final
    cdef inline gr _ceil(self, gr x):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_ceil(res.pval, x.pval, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot compute ceil(x) in this context")
        return res

    @cython.final
    cdef inline gr _trunc(self, gr x):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_trunc(res.pval, x.pval, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot compute trunc(x) in this context")
        return res

    @cython.final
    cdef inline gr _nint(self, gr x):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_nint(res.pval, x.pval, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot compute nint(x) in this context")
        return res

    @cython.final
    cdef inline gr _abs(self, gr x):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_abs(res.pval, x.pval, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot compute abs(x) in this context")
        return res

    @cython.final
    cdef inline gr _conj(self, gr x):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_conj(res.pval, x.pval, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot compute conj(x) in this context")
        return res

    @cython.final
    cdef inline gr _re(self, gr x):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_re(res.pval, x.pval, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot compute re(x) in this context")
        return res

    @cython.final
    cdef inline gr _im(self, gr x):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_im(res.pval, x.pval, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot compute im(x) in this context")
        return res

    @cython.final
    cdef inline gr _sgn(self, gr x):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_sgn(res.pval, x.pval, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot compute sgn(x) in this context")
        return res

    @cython.final
    cdef inline gr _csgn(self, gr x):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_csgn(res.pval, x.pval, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot compute csgn(x) in this context")
        return res

    @cython.final
    cdef inline gr _arg(self, gr x):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_arg(res.pval, x.pval, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot compute arg(x) in this context")
        return res

    ###
    # Ordering methods

    @cython.final
    cdef inline int _cmp(self, gr x, gr y):
        cdef int err
        cdef int res
        err = gr_cmp(&res, x.pval, y.pval, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot compare x and y")
        return res

    @cython.final
    cdef inline int _cmp_other(self, gr x, gr y):
        cdef int err
        cdef int res
        err = gr_cmp_other(&res, x.pval, y.pval, y.ctx.ctx_t, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot compare x and y")
        return res

    @cython.final
    cdef inline int _cmpabs(self, gr x, gr y):
        cdef int err
        cdef int res
        err = gr_cmpabs(&res, x.pval, y.pval, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot compare x and y")
        return res

    @cython.final
    cdef inline int _cmpabs_other(self, gr x, gr y):
        cdef int err
        cdef int res
        err = gr_cmpabs_other(&res, x.pval, y.pval, y.ctx.ctx_t, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot compare x and y")
        return res

    @cython.final
    cdef inline truth_t _le(self, gr x, gr y):
        return gr_le(x.pval, y.pval, self.ctx_t)

    @cython.final
    cdef inline truth_t _abs_le(self, gr x, gr y):
        return gr_abs_le(x.pval, y.pval, self.ctx_t)

    @cython.final
    cdef inline truth_t _lt(self, gr x, gr y):
        return gr_lt(x.pval, y.pval, self.ctx_t)

    @cython.final
    cdef inline truth_t _abs_lt(self, gr x, gr y):
        return gr_abs_lt(x.pval, y.pval, self.ctx_t)

    @cython.final
    cdef inline truth_t _ge(self, gr x, gr y):
        return gr_ge(x.pval, y.pval, self.ctx_t)

    @cython.final
    cdef inline truth_t _abs_ge(self, gr x, gr y):
        return gr_abs_ge(x.pval, y.pval, self.ctx_t)

    @cython.final
    cdef inline truth_t _gt(self, gr x, gr y):
        return gr_gt(x.pval, y.pval, self.ctx_t)

    @cython.final
    cdef inline truth_t _abs_gt(self, gr x, gr y):
        return gr_abs_gt(x.pval, y.pval, self.ctx_t)

    @cython.final
    cdef inline gr _min(self, gr x, gr y):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_min(res.pval, x.pval, y.pval, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot compute min(x) in this context")
        return res

    @cython.final
    cdef inline gr _max(self, gr x, gr y):
        cdef int err
        cdef gr res = self.new_gr()
        err = gr_max(res.pval, x.pval, y.pval, self.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Cannot compute min(x) in this context")
        return res

    # @cython.final
    # cdef inline list _gens_recursive(self):
    #     cdef int err
    #     cdef gr g
    #     cdef gr_vec_t gens
    #     gr_vec_init(gens, 0, self.ctx_t)
    #     err = gr_gens_recursive(gens, self.ctx_t)
    #     if err != GR_SUCCESS:
    #         raise self._error(err, "Cannot get recursive generators")
    #     length = gr_vec_length(gens, self.ctx_t)
    #     py_gens = [None] * length
    #     for 0 <= i < length:
    #         g = self.new_gr()
    #         err = gr_set(g.pval, gr_vec_entry_ptr(gens, i, self.ctx_t), self.ctx_t)
    #         if err != GR_SUCCESS:
    #             raise self._error(err, "Failed to copy generator.")
    #         py_gens[i] = g
    #     gr_vec_clear(gens, self.ctx_t)
    #     return py_gens


cdef class gr_scalar_ctx(gr_ctx):
    pass


cdef class gr_poly_ctx(gr_ctx):
    pass


cdef class gr_mpoly_ctx(gr_ctx):
    pass


# cdef class gr_matrix_domain_ctx(gr_ctx):
#     pass


# cdef class gr_matrix_space_ctx(gr_ctx):
#     pass


# cdef class gr_matrix_ring_ctx(gr_ctx):
#     pass


@cython.no_gc
cdef class _gr_fmpz_ctx(gr_scalar_ctx):

    @staticmethod
    cdef inline _gr_fmpz_ctx _new():
        cdef _gr_fmpz_ctx ctx
        ctx = _gr_fmpz_ctx.__new__(_gr_fmpz_ctx)
        gr_ctx_init_fmpz(ctx.ctx_t)
        ctx._init = True
        return ctx


@cython.no_gc
cdef class _gr_fmpq_ctx(gr_scalar_ctx):

    @staticmethod
    cdef inline _gr_fmpq_ctx _new():
        cdef _gr_fmpq_ctx ctx
        ctx = _gr_fmpq_ctx.__new__(_gr_fmpq_ctx)
        gr_ctx_init_fmpq(ctx.ctx_t)
        ctx._init = True
        return ctx


@cython.no_gc
cdef class _gr_fmpzi_ctx(gr_scalar_ctx):

    @staticmethod
    cdef inline _gr_fmpzi_ctx _new():
        cdef _gr_fmpzi_ctx ctx
        ctx = _gr_fmpzi_ctx.__new__(_gr_fmpzi_ctx)
        gr_ctx_init_fmpzi(ctx.ctx_t)
        ctx._init = True
        return ctx


@cython.no_gc
cdef class _gr_fexpr_ctx(gr_scalar_ctx):

    @staticmethod
    cdef inline _gr_fexpr_ctx _new():
        cdef _gr_fexpr_ctx ctx
        ctx = _gr_fexpr_ctx.__new__(_gr_fexpr_ctx)
        gr_ctx_init_fexpr(ctx.ctx_t)
        ctx._init = True
        return ctx


# The global contexts for use in cython code:
cdef _gr_fmpz_ctx gr_fmpz_ctx_c
cdef _gr_fmpq_ctx gr_fmpq_ctx_c
cdef _gr_fmpzi_ctx gr_fmpzi_ctx_c
cdef _gr_fexpr_ctx gr_fexpr_ctx_c


@cython.no_gc
cdef class gr_nmod_ctx(gr_scalar_ctx):
    cdef ulong n
    cdef bint is_field

    @staticmethod
    cdef inline gr_nmod_ctx _new(ulong n):
        cdef gr_nmod_ctx ctx
        cdef bint is_prime = n_is_prime(n)
        ctx = gr_nmod_ctx.__new__(gr_nmod_ctx)
        ctx.n = n
        ctx.is_field = is_prime
        gr_ctx_init_nmod(ctx.ctx_t, n)
        ctx._init = True
        # if is_prime:
        #     gr_ctx_set_is_field(ctx.ctx_t, T_TRUE)
        # else:
        #     gr_ctx_set_is_field(ctx.ctx_t, T_FALSE)
        return ctx


@cython.no_gc
cdef class gr_fmpz_mod_ctx(gr_scalar_ctx):
    cdef fmpz_t n
    cdef bint is_field

    @staticmethod
    cdef inline gr_fmpz_mod_ctx _new(fmpz n):
        cdef gr_fmpz_mod_ctx ctx
        cdef bint is_prime = n.is_prime()
        ctx = gr_fmpz_mod_ctx.__new__(gr_fmpz_mod_ctx)
        ctx.is_field = is_prime
        fmpz_init_set(ctx.n, n.val)
        gr_ctx_init_fmpz_mod(ctx.ctx_t, ctx.n)
        ctx._init = True
        # if is_prime:
        #     gr_ctx_set_is_field(ctx.ctx_t, T_TRUE)
        # else:
        #     gr_ctx_set_is_field(ctx.ctx_t, T_FALSE)
        return ctx


@cython.no_gc
cdef class gr_fq_ctx(gr_scalar_ctx):
    cdef fmpz_t p
    cdef slong d

    @staticmethod
    cdef inline gr_fq_ctx _new(fmpz p, slong d, char* name):
        cdef gr_fq_ctx ctx
        ctx = gr_fq_ctx.__new__(gr_fq_ctx)
        ctx.d = d
        fmpz_init_set(ctx.p, p.val)
        gr_ctx_init_fq(ctx.ctx_t, ctx.p, d, name)
        ctx._init = True
        # XXX: free name_c?
        return ctx


@cython.no_gc
cdef class gr_fq_nmod_ctx(gr_scalar_ctx):
    cdef ulong p
    cdef slong d

    @staticmethod
    cdef inline gr_fq_nmod_ctx _new(ulong p, slong d, char* name):
        cdef gr_fq_nmod_ctx ctx
        ctx = gr_fq_nmod_ctx.__new__(gr_fq_nmod_ctx)
        ctx.p = p
        ctx.d = d
        gr_ctx_init_fq_nmod(ctx.ctx_t, p, d, name)
        ctx._init = True
        return ctx


@cython.no_gc
cdef class gr_fq_zech_ctx(gr_scalar_ctx):
    cdef ulong p
    cdef slong d

    @staticmethod
    cdef inline gr_fq_zech_ctx _new(ulong p, slong d, char* name):
        cdef gr_fq_zech_ctx ctx
        ctx = gr_fq_zech_ctx.__new__(gr_fq_zech_ctx)
        ctx.p = p
        ctx.d = d
        gr_ctx_init_fq_zech(ctx.ctx_t, p, d, name)
        ctx._init = True
        return ctx


@cython.no_gc
cdef class gr_nf_ctx(gr_scalar_ctx):
    cdef fmpq_poly_t poly

    @staticmethod
    cdef inline gr_nf_ctx _new(fmpq_poly poly):
        cdef gr_nf_ctx ctx
        ctx = gr_nf_ctx.__new__(gr_nf_ctx)
        fmpq_poly_init(ctx.poly)
        fmpq_poly_set(ctx.poly, poly.val)
        gr_ctx_init_nf(ctx.ctx_t, ctx.poly)
        ctx._init = True
        return ctx


@cython.no_gc
cdef class gr_nf_fmpz_poly_ctx(gr_scalar_ctx):
    cdef fmpz_poly_t poly

    @staticmethod
    cdef inline gr_nf_fmpz_poly_ctx _new(fmpz_poly poly):
        cdef gr_nf_fmpz_poly_ctx ctx
        ctx = gr_nf_fmpz_poly_ctx.__new__(gr_nf_fmpz_poly_ctx)
        fmpz_poly_init(ctx.poly)
        fmpz_poly_set(ctx.poly, poly.val)
        gr_ctx_init_nf_fmpz_poly(ctx.ctx_t, ctx.poly)
        ctx._init = True
        return ctx


@cython.no_gc
cdef class gr_real_qqbar_ctx(gr_scalar_ctx):
    cdef slong deg_limit
    cdef slong bits_limit

    @staticmethod
    cdef inline gr_real_qqbar_ctx _new(slong deg_limit, slong bits_limit):
        cdef gr_real_qqbar_ctx ctx
        ctx = gr_real_qqbar_ctx.__new__(gr_real_qqbar_ctx)
        gr_ctx_init_real_qqbar(ctx.ctx_t)
        ctx._init = True
        ctx.deg_limit = deg_limit
        ctx.bits_limit = bits_limit
        if deg_limit != -1 or bits_limit != -1:
            # Maybe add setters for these?
            # XXX: Not available in FLINT 3.0.1
            # _gr_ctx_qqbar_set_limits(ctx.ctx_t, deg_limit, bits_limit)
            pass
        return ctx


@cython.no_gc
cdef class gr_complex_qqbar_ctx(gr_scalar_ctx):
    cdef slong deg_limit
    cdef slong bits_limit

    @staticmethod
    cdef inline gr_complex_qqbar_ctx _new(slong deg_limit, slong bits_limit):
        cdef gr_complex_qqbar_ctx ctx
        ctx = gr_complex_qqbar_ctx.__new__(gr_complex_qqbar_ctx)
        gr_ctx_init_complex_qqbar(ctx.ctx_t)
        ctx._init = True
        ctx.deg_limit = deg_limit
        ctx.bits_limit = bits_limit
        if deg_limit != -1 or bits_limit != -1:
            # Maybe add setters for these?
            # XXX: Not available in FLINT 3.0.1
            # _gr_ctx_qqbar_set_limits(ctx.ctx_t, deg_limit, bits_limit)
            pass
        return ctx


@cython.no_gc
cdef class gr_real_ca_ctx(gr_scalar_ctx):
    cdef dict options

    @staticmethod
    cdef inline gr_real_ca_ctx _new(dict options):
        cdef gr_real_ca_ctx ctx
        if options:
            raise NotImplementedError("Options not implemented yet")
        ctx = gr_real_ca_ctx.__new__(gr_real_ca_ctx)
        gr_ctx_init_real_ca(ctx.ctx_t)
        ctx._init = True
        ctx.options = options.copy()
        return ctx


@cython.no_gc
cdef class gr_complex_ca_ctx(gr_scalar_ctx):
    cdef dict options

    @staticmethod
    cdef inline gr_complex_ca_ctx _new(dict options):
        cdef gr_complex_ca_ctx ctx
        if options:
            raise NotImplementedError("Options not implemented yet")
        ctx = gr_complex_ca_ctx.__new__(gr_complex_ca_ctx)
        gr_ctx_init_complex_ca(ctx.ctx_t)
        ctx._init = True
        ctx.options = options.copy()
        return ctx


@cython.no_gc
cdef class gr_real_algebraic_ca_ctx(gr_scalar_ctx):
    cdef dict options

    @staticmethod
    cdef inline gr_real_algebraic_ca_ctx _new(dict options):
        cdef gr_real_algebraic_ca_ctx ctx
        if options:
            raise NotImplementedError("Options not implemented yet")
        ctx = gr_real_algebraic_ca_ctx.__new__(gr_real_algebraic_ca_ctx)
        gr_ctx_init_real_algebraic_ca(ctx.ctx_t)
        ctx._init = True
        ctx.options = options.copy()
        return ctx


@cython.no_gc
cdef class gr_complex_algebraic_ca_ctx(gr_scalar_ctx):
    cdef dict options

    @staticmethod
    cdef inline gr_complex_algebraic_ca_ctx _new(dict options):
        cdef gr_complex_algebraic_ca_ctx ctx
        if options:
            raise NotImplementedError("Options not implemented yet")
        ctx = gr_complex_algebraic_ca_ctx.__new__(gr_complex_algebraic_ca_ctx)
        gr_ctx_init_complex_algebraic_ca(ctx.ctx_t)
        ctx._init = True
        ctx.options = options.copy()
        return ctx


@cython.no_gc
cdef class gr_complex_extended_ca_ctx(gr_scalar_ctx):
    cdef dict options

    @staticmethod
    cdef inline gr_complex_extended_ca_ctx _new(dict options):
        cdef gr_complex_extended_ca_ctx ctx
        if options:
            raise NotImplementedError("Options not implemented yet")
        ctx = gr_complex_extended_ca_ctx.__new__(gr_complex_extended_ca_ctx)
        gr_ctx_init_complex_extended_ca(ctx.ctx_t)
        ctx._init = True
        ctx.options = options.copy()
        return ctx


@cython.no_gc
cdef class gr_real_float_arf_ctx(gr_scalar_ctx):
    cdef slong prec

    @staticmethod
    cdef inline gr_real_float_arf_ctx _new(slong prec):
        cdef gr_real_float_arf_ctx ctx
        ctx = gr_real_float_arf_ctx.__new__(gr_real_float_arf_ctx)
        gr_ctx_init_real_float_arf(ctx.ctx_t, prec)
        ctx._init = True
        ctx.prec = prec
        return ctx


@cython.no_gc
cdef class gr_complex_float_acf_ctx(gr_scalar_ctx):
    cdef slong prec

    @staticmethod
    cdef inline gr_complex_float_acf_ctx _new(slong prec):
        cdef gr_complex_float_acf_ctx ctx
        ctx = gr_complex_float_acf_ctx.__new__(gr_complex_float_acf_ctx)
        gr_ctx_init_complex_float_acf(ctx.ctx_t, prec)
        ctx._init = True
        ctx.prec = prec
        return ctx


@cython.no_gc
cdef class gr_real_arb_ctx(gr_scalar_ctx):
    cdef slong prec

    @staticmethod
    cdef inline gr_real_arb_ctx _new(slong prec):
        cdef gr_real_arb_ctx ctx
        ctx = gr_real_arb_ctx.__new__(gr_real_arb_ctx)
        gr_ctx_init_real_arb(ctx.ctx_t, prec)
        ctx._init = True
        ctx.prec = prec
        return ctx


@cython.no_gc
cdef class gr_complex_acb_ctx(gr_scalar_ctx):
    cdef slong prec

    @staticmethod
    cdef inline gr_complex_acb_ctx _new(slong prec):
        cdef gr_complex_acb_ctx ctx
        ctx = gr_complex_acb_ctx.__new__(gr_complex_acb_ctx)
        gr_ctx_init_complex_acb(ctx.ctx_t, prec)
        ctx._init = True
        ctx.prec = prec
        return ctx


# @cython.no_gc
# cdef class _gr_fmpz_poly_ctx(gr_poly_ctx):
#
#     @staticmethod
#     cdef _gr_fmpz_poly_ctx _new()


# @cython.no_gc
# cdef class _gr_fmpq_poly_ctx(gr_poly_ctx):
#
#     @staticmethod
#     cdef _gr_fmpq_poly_ctx _new()


# @cython.no_gc
# cdef class _gr_gr_poly_ctx(gr_poly_ctx):
#
#     @staticmethod
#     cdef _gr_gr_poly_ctx _new()


@cython.no_gc
cdef class gr_gr_poly_ctx(gr_poly_ctx):
    cdef gr_ctx base_ctx

    @staticmethod
    cdef inline gr_gr_poly_ctx _new(gr_ctx base_ctx):
        cdef gr_gr_poly_ctx ctx
        ctx = gr_gr_poly_ctx.__new__(gr_gr_poly_ctx)
        gr_ctx_init_gr_poly(ctx.ctx_t, base_ctx.ctx_t)
        ctx._init = True
        ctx.base_ctx = base_ctx
        return ctx


@cython.no_gc
cdef class gr_gr_mpoly_ctx(gr_mpoly_ctx):
    cdef gr_ctx base_ctx
    cdef ordering_t _order
    cdef tuple _names
    cdef slong _nvars

    @staticmethod
    cdef inline gr_gr_mpoly_ctx _new(gr_ctx base_ctx, tuple names, order):
        cdef gr_gr_mpoly_ctx ctx
        cdef ordering_t ord_c
        cdef slong nvars
        cdef const char **names_c
        cdef int status

        nvars = len(names)
        names_b = [name.encode('utf-8') for name in names]

        ord_c = ordering_py_to_c(order)
        ctx = gr_gr_mpoly_ctx.__new__(gr_gr_mpoly_ctx)
        gr_ctx_init_gr_mpoly(ctx.ctx_t, base_ctx.ctx_t, nvars, ord_c)
        ctx._init = True
        ctx.base_ctx = base_ctx
        ctx._order = ord_c
        ctx._nvars = nvars
        ctx._names = names

        names_c = <const char **>libc.stdlib.malloc(nvars * sizeof(const char *))
        if names_c == NULL:
            raise MemoryError("Failed to allocate memory for generator names")
        try:
            for i in range(nvars):
                names_c[i] = names_b[i]
            status = gr_ctx_set_gen_names(ctx.ctx_t, names_c)
        finally:
            libc.stdlib.free(names_c)
        if status != GR_SUCCESS:
            raise MemoryError("Failed to set generator names")

        return ctx


# @cython.no_gc
# cdef class gr_fmpz_mpoly_q_ctx(gr_mpoly_ctx):
#     cdef ordering_t _order
#     cdef slong _nvars
#
#     @staticmethod
#     cdef inline gr_fmpz_mpoly_q_ctx _new(slong nvars, order):
#         cdef gr_fmpz_mpoly_q_ctx ctx
#         cdef ordering_t ord_c
#         ord_c = ordering_py_to_c(order)
#         ctx = gr_fmpz_mpoly_q_ctx.__new__(gr_fmpz_mpoly_q_ctx)
#         gr_ctx_init_fmpz_mpoly_q(ctx.ctx_t, nvars, ord_c)
#         ctx._init = True
#         ctx._order = ord_c
#         ctx._nvars = nvars
#         return ctx


# @cython.no_gc
# cdef class gr_series_mod_gr_poly_ctx(gr_ctx):
#     cdef gr_ctx base_ctx
#     cdef slong _n
#
#     @staticmethod
#     cdef inline gr_series_mod_gr_poly_ctx _new(gr_ctx base_ctx, slong n):
#         cdef gr_series_mod_gr_poly_ctx ctx
#         ctx = gr_series_mod_gr_poly_ctx.__new__(gr_series_mod_gr_poly_ctx)
#         gr_ctx_init_series_mod_gr_poly(ctx.ctx_t, base_ctx.ctx_t, n)
#         ctx._init = True
#         ctx.base_ctx = base_ctx
#         ctx._n = n
#         return ctx


@cython.no_gc
cdef class gr_series_ctx(gr_ctx):
    cdef gr_ctx base_ctx
    cdef slong _prec

    @staticmethod
    cdef inline gr_series_ctx _new(gr_ctx base_ctx, slong prec):
        cdef gr_series_ctx ctx
        ctx = gr_series_ctx.__new__(gr_series_ctx)
        gr_ctx_init_gr_series(ctx.ctx_t, base_ctx.ctx_t, prec)
        ctx._init = True
        ctx.base_ctx = base_ctx
        ctx._prec = prec
        return ctx


@cython.no_gc
cdef class gr(flint_scalar):
    cdef gr_ptr pval
    cdef public gr_ctx ctx
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
