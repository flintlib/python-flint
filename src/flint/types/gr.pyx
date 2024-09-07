cimport cython

from flint.flintlib.types.flint cimport slong
from flint.flintlib.types.gr cimport GR_SUCCESS
from flint.flintlib.functions.gr cimport (
    gr_ctx_clear,
    gr_heap_clear,

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
from flint.flintlib.functions.gr_domains cimport (
    gr_ctx_init_fmpz,
    gr_ctx_init_fmpq,

    gr_ctx_is_finite,
    gr_ctx_is_multiplicative_group,
    gr_ctx_is_ring,
    gr_ctx_is_commutative_ring,
    gr_ctx_is_integral_domain,
    gr_ctx_is_unique_factorization_domain,
    gr_ctx_is_field,
    gr_ctx_is_algebraically_closed,
    gr_ctx_is_finite_characteristic,
    gr_ctx_is_ordered_ring,
    gr_ctx_is_zero_ring,
    gr_ctx_is_exact,
    gr_ctx_is_canonical,
    gr_ctx_has_real_prec,
    gr_ctx_set_real_prec,
    gr_ctx_get_real_prec,
)


@cython.no_gc
cdef class gr_ctx(flint_ctx):
    """Base class for all gr contexts.

    This class should not be instantiated directly. Instead, use one of the
    derived classes, e.g. gr_nmod_ctx or gr_fmpz_ctx.
    """
    def __init__(self, *args, **kwargs):
        raise TypeError("Cannot create gr_ctx object directly."
                        " Use e.g. gr_nmod_ctx.new(n) instead.")

    def new(self):
        raise NotImplementedError("Cannot create gr_ctx object directly."
                                  " Use e.g. gr_nmod_ctx.new(n) instead.")

    def __dealloc__(self):
        if self._init:
            self._init = False
            gr_ctx_clear(self.ctx_t)

    @property
    def is_finite(self) -> bool | None:
        return truth_to_py(gr_ctx_is_finite(self.ctx_t))

    @property
    def is_multiplicative_group(self) -> bool | None:
        return truth_to_py(gr_ctx_is_multiplicative_group(self.ctx_t))

    @property
    def is_ring(self) -> bool | None:
        return truth_to_py(gr_ctx_is_ring(self.ctx_t))

    @property
    def is_commutative_ring(self) -> bool | None:
        return truth_to_py(gr_ctx_is_commutative_ring(self.ctx_t))

    @property
    def is_integral_domain(self) -> bool | None:
        return truth_to_py(gr_ctx_is_integral_domain(self.ctx_t))

    @property
    def is_unique_factorization_domain(self) -> bool | None:
        return truth_to_py(gr_ctx_is_unique_factorization_domain(self.ctx_t))

    @property
    def is_field(self) -> bool | None:
        return truth_to_py(gr_ctx_is_field(self.ctx_t))

    @property
    def is_algebraically_closed(self) -> bool | None:
        return truth_to_py(gr_ctx_is_algebraically_closed(self.ctx_t))

    @property
    def is_finite_characteristic(self) -> bool | None:
        return truth_to_py(gr_ctx_is_finite_characteristic(self.ctx_t))

    @property
    def is_ordered_ring(self) -> bool | None:
        return truth_to_py(gr_ctx_is_ordered_ring(self.ctx_t))

    @property
    def is_zero_ring(self) -> bool | None:
        return truth_to_py(gr_ctx_is_zero_ring(self.ctx_t))

    @property
    def is_exact(self) -> bool | None:
        return truth_to_py(gr_ctx_is_exact(self.ctx_t))

    @property
    def is_canonical(self) -> bool | None:
        return truth_to_py(gr_ctx_is_canonical(self.ctx_t))

    @property
    def has_real_prec(self) -> bool | None:
        return truth_to_py(gr_ctx_has_real_prec(self.ctx_t))

    @property
    def real_prec(self) -> int:
        cdef slong prec
        if gr_ctx_has_real_prec(self.ctx_t) != truth_t.T_TRUE:
            raise ValueError("Real precision is not available")
        if gr_ctx_get_real_prec(&prec, self.ctx_t) != GR_SUCCESS:
            raise ValueError("Failed to get real precision")
        return prec

    @real_prec.setter
    def real_prec(self, prec: slong) -> None:
        if gr_ctx_has_real_prec(self.ctx_t) != truth_t.T_TRUE:
            raise ValueError("Real precision is not available")
        if gr_ctx_set_real_prec(self.ctx_t, prec) != GR_SUCCESS:
            raise ValueError("Failed to set real precision")

    def __call__(self, arg) -> gr:
        return self.from_str(str(arg))


cdef _gr_fmpz_ctx gr_fmpz_ctx_c
cdef _gr_fmpq_ctx gr_fmpq_ctx_c


@cython.no_gc
cdef class _gr_fmpz_ctx(gr_ctx):

    @staticmethod
    def new() -> _gr_fmpz_ctx:
        return gr_fmpz_ctx_c

    @staticmethod
    cdef _gr_fmpz_ctx _new():
        cdef _gr_fmpz_ctx ctx
        ctx = _gr_fmpz_ctx.__new__(_gr_fmpz_ctx)
        gr_ctx_init_fmpz(ctx.ctx_t)
        ctx._init = True
        return ctx

    def __repr__(self):
        return "gr_fmpz_ctx"


@cython.no_gc
cdef class _gr_fmpq_ctx(gr_ctx):

    @staticmethod
    def new() -> _gr_fmpq_ctx:
        return gr_fmpq_ctx_c

    @staticmethod
    cdef _gr_fmpq_ctx _new():
        cdef _gr_fmpq_ctx ctx
        ctx = _gr_fmpq_ctx.__new__(_gr_fmpq_ctx)
        gr_ctx_init_fmpq(ctx.ctx_t)
        ctx._init = True
        return ctx

    def __repr__(self):
        return "gr_fmpq_ctx"


# Global contexts for Cython code
gr_fmpz_ctx_c = _gr_fmpz_ctx._new()
gr_fmpq_ctx_c = _gr_fmpq_ctx._new()

# Global contexts for Python code
gr_fmpz_ctx = gr_fmpz_ctx_c
gr_fmpq_ctx = gr_fmpq_ctx_c


@cython.no_gc
cdef class gr(flint_scalar):

    def __dealloc__(self):
        if self._init:
            self._init = False
            gr_heap_clear(self.pval, self.ctx.ctx_t)

    def __repr__(self):
        return self.ctx.to_str(self)

    def _neg(self) -> gr:
        cdef int err
        cdef gr res = self.ctx.new_gr()
        err = gr_neg(res.pval, self.pval, self.ctx.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Failed to negate gr object")
        return res

    def _add(self, other: gr) -> gr:
        cdef int err
        cdef gr res = self.ctx.new_gr()
        err = gr_add(res.pval, self.pval, other.pval, self.ctx.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Failed to add gr objects")
        return res

    def _add_si(self, other: int) -> gr:
        cdef int err
        cdef gr res = self.ctx.new_gr()
        err = gr_add_si(res.pval, self.pval, other, self.ctx.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Failed to add gr object and int")
        return res

    def _sub(self, other: gr) -> gr:
        cdef int err
        cdef gr res = self.ctx.new_gr()
        err = gr_sub(res.pval, self.pval, other.pval, self.ctx.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Failed to subtract gr objects")
        return res

    def _sub_si(self, other: int) -> gr:
        cdef int err
        cdef gr res = self.ctx.new_gr()
        err = gr_sub_si(res.pval, self.pval, other, self.ctx.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Failed to subtract gr object and int")
        return res

    def _mul(self, other: gr) -> gr:
        cdef int err
        cdef gr res = self.ctx.new_gr()
        err = gr_mul(res.pval, self.pval, other.pval, self.ctx.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Failed to multiply gr objects")
        return res

    def _mul_si(self, other: int) -> gr:
        cdef int err
        cdef gr res = self.ctx.new_gr()
        err = gr_mul_si(res.pval, self.pval, other, self.ctx.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Failed to multiply gr object and int")
        return res

    def _inv(self) -> gr:
        cdef int err
        cdef gr res = self.ctx.new_gr()
        err = gr_inv(res.pval, self.pval, self.ctx.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Failed to invert gr object")
        return res

    def _div(self, other: gr) -> gr:
        cdef int err
        cdef gr res = self.ctx.new_gr()
        err = gr_div(res.pval, self.pval, other.pval, self.ctx.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Failed to divide gr objects")
        return res

    def _div_si(self, other: int) -> gr:
        cdef int err
        cdef gr res = self.ctx.new_gr()
        err = gr_div_si(res.pval, self.pval, other, self.ctx.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Failed to divide gr object and int")
        return res

    def _pow_si(self, other: int) -> gr:
        cdef int err
        cdef gr res = self.ctx.new_gr()
        err = gr_pow_si(res.pval, self.pval, other, self.ctx.ctx_t)
        if err != GR_SUCCESS:
            raise self._error(err, "Failed to raise gr object to an integer power")
        return res

    def __neg__(self) -> gr:
        return self._neg()

    # XXX: Maybe +a should return a copy or for arb it should round...
    # def __pos__(self) -> gr:
    #     return self

    def __add__(self, other) -> gr:
        cdef gr_ctx ctx
        cdef gr other_gr, res
        if isinstance(other, gr):
            other_gr = other
            if self.ctx != other_gr.ctx:
                raise ValueError("Cannot add elements from different contexts")
            return self._add(other_gr)
        elif isinstance(other, int):
            return self._add_si(other)
        else:
            return NotImplemented

    def __radd__(self, other) -> gr:
        cdef gr res
        if isinstance(other, int):
            return self._add_si(other)
        else:
            return NotImplemented

    def __sub__(self, other) -> gr:
        cdef gr other_gr, res
        if isinstance(other, gr):
            other_gr = other
            if self.ctx != other_gr.ctx:
                raise ValueError("Cannot subtract elements from different contexts")
            return self._sub(other_gr)
        elif isinstance(other, int):
            return self._sub_si(other)
        else:
            return NotImplemented

    def __rsub__(self, other) -> gr:
        cdef gr res
        if isinstance(other, int):
            return self._neg()._add_si(other)
        else:
            return NotImplemented

    def __mul__(self, other) -> gr:
        cdef gr other_gr, res
        if isinstance(other, gr):
            other_gr = other
            if self.ctx != other._grctx:
                raise ValueError("Cannot multiply elements from different contexts")
            return self._mul(other_gr)
        elif isinstance(other, int):
            return self._mul_si(other)
        else:
            return NotImplemented

    def __rmul__(self, other) -> gr:
        cdef gr res
        if isinstance(other, int):
            return self._mul_si(other)
        else:
            return NotImplemented

    def __truediv__(self, other) -> gr:
        cdef gr other_gr, res
        if isinstance(other, gr):
            other_gr = other
            if self.ctx != other_gr.ctx:
                raise ValueError("Cannot divide elements from different contexts")
            return self._div(other_gr)
        elif isinstance(other, int):
            return self._div_si(other)
        else:
            return NotImplemented

    def __rtruediv__(self, other) -> gr:
        cdef gr res
        if isinstance(other, int):
            return self._div_si(other)._inv()
        else:
            return NotImplemented

    def __pow__(self, other) -> gr:
        if isinstance(other, int):
            return self._pow_si(other)
        else:
            return NotImplemented
