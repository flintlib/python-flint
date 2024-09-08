cimport cython

from flint.flintlib.types.flint cimport slong
from flint.flintlib.types.gr cimport GR_SUCCESS
from flint.flintlib.functions.gr cimport (
    gr_ctx_clear,
    gr_heap_clear,

    gr_set_si,
)
from flint.flintlib.functions.gr_domains cimport (
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
        if gr_ctx_has_real_prec(self.ctx_t) != T_TRUE:
            raise ValueError("Real precision is not available")
        if gr_ctx_get_real_prec(&prec, self.ctx_t) != GR_SUCCESS:
            raise ValueError("Failed to get real precision")
        return prec

    @real_prec.setter
    def real_prec(self, prec: slong) -> None:
        if gr_ctx_has_real_prec(self.ctx_t) != T_TRUE:
            raise ValueError("Real precision is not available")
        if gr_ctx_set_real_prec(self.ctx_t, prec) != GR_SUCCESS:
            raise ValueError("Failed to set real precision")

    def __call__(self, arg) -> gr:
        return self.from_str(str(arg))

    def zero(self) -> gr:
        return self._zero()

    def one(self) -> gr:
        return self._one()

    def i(self) -> gr:
        return self._i()

    def pos_inf(self) -> gr:
        return self._pos_inf()

    def neg_inf(self) -> gr:
        return self._neg_inf()

    def uinf(self) -> gr:
        return self._uinf()

    def undefined(self) -> gr:
        return self._undefined()

    def unknown(self) -> gr:
        return self._unknown()

    def gen(self) -> gr:
        return self._gen()

    def gens(self) -> list[gr]:
        return self._gens()

    def gens_recursive(self) -> list[gr]:
        return self._gens_recursive()


cdef class gr_scalar_ctx(gr_ctx):
    pass


cdef class gr_poly_ctx(gr_ctx):
    pass


cdef class gr_mpoly_ctx(gr_ctx):
    pass


cdef class gr_matrix_domain_ctx(gr_ctx):
    pass


cdef class gr_matrix_space_ctx(gr_ctx):
    pass


cdef class gr_matrix_ring_ctx(gr_ctx):
    pass


# Global contexts for Cython code
cdef _gr_fmpz_ctx gr_fmpz_ctx_c = _gr_fmpz_ctx._new()
cdef _gr_fmpq_ctx gr_fmpq_ctx_c = _gr_fmpq_ctx._new()
cdef _gr_fmpzi_ctx gr_fmpzi_ctx_c = _gr_fmpzi_ctx._new()
cdef _gr_fexpr_ctx gr_fexpr_ctx_c = _gr_fexpr_ctx._new()

# Global contexts for Python code
gr_fmpz_ctx = gr_fmpz_ctx_c
gr_fmpq_ctx = gr_fmpq_ctx_c
gr_fmpzi_ctx = gr_fmpzi_ctx_c
gr_fexpr_ctx = gr_fexpr_ctx_c


@cython.no_gc
cdef class _gr_fmpz_ctx(gr_scalar_ctx):

    @staticmethod
    def new() -> _gr_fmpz_ctx:
        return gr_fmpz_ctx_c

    def __repr__(self):
        return "gr_fmpz_ctx"


@cython.no_gc
cdef class _gr_fmpq_ctx(gr_scalar_ctx):

    @staticmethod
    def new() -> _gr_fmpq_ctx:
        return gr_fmpq_ctx_c

    def __repr__(self):
        return "gr_fmpq_ctx"


@cython.no_gc
cdef class _gr_fmpzi_ctx(gr_scalar_ctx):

    @staticmethod
    def new() -> _gr_fmpzi_ctx:
        return gr_fmpzi_ctx_c

    def __repr__(self):
        return "gr_fmpzi_ctx"


@cython.no_gc
cdef class _gr_fexpr_ctx(gr_scalar_ctx):

    @staticmethod
    def new() -> _gr_fexpr_ctx:
        return gr_fexpr_ctx_c

    def __repr__(self):
        return "gr_fexpr_ctx"


@cython.no_gc
cdef class gr_nmod_ctx(gr_scalar_ctx):

    @staticmethod
    def new(ulong n) -> gr_nmod_ctx:
        return gr_nmod_ctx._new(n)

    def __repr__(self):
        return f"gr_nmod_ctx({self.n})"


@cython.no_gc
cdef class gr_fmpz_mod_ctx(gr_scalar_ctx):

    @staticmethod
    def new(n) -> gr_fmpz_mod_ctx:
        n = fmpz(n)
        return gr_fmpz_mod_ctx._new(n)

    def modulus(self):
        cdef fmpz n = fmpz.__new__(fmpz)
        fmpz_init_set(n.val, self.n)
        return n

    def __repr__(self):
        return f"gr_fmpz_mod_ctx({self.modulus()})"


@cython.no_gc
cdef class gr_fq_ctx(gr_scalar_ctx):

    @staticmethod
    def new(p, d, name=None) -> gr_fq_ctx:
        cdef bytes name_b
        cdef char *name_c
        if name is not None:
            name_b = name.encode('utf-8')
            name_c = name_b
        else:
            name_c = NULL
        p = fmpz(p)
        return gr_fq_ctx._new(p, d, name_c)

    def characteristic(self):
        cdef fmpz p = fmpz.__new__(fmpz)
        fmpz_init_set(p.val, self.p)
        return p

    def degree(self):
        return self.d

    def __repr__(self):
        return f"gr_fq_ctx({self.characteristic()}, {self.degree()})"


@cython.no_gc
cdef class gr_fq_nmod_ctx(gr_scalar_ctx):

    @staticmethod
    def new(p, d, name=None) -> gr_fq_nmod_ctx:
        cdef bytes name_b
        cdef char *name_c
        if name is not None:
            name_b = name.encode('utf-8')
            name_c = name_b
        else:
            name_c = NULL
        return gr_fq_nmod_ctx._new(p, d, name_c)

    def characteristic(self):
        return self.p

    def degree(self):
        return self.d

    def __repr__(self):
        return f"gr_fq_nmod_ctx({self.characteristic()}, {self.degree()})"


@cython.no_gc
cdef class gr_fq_zech_ctx(gr_scalar_ctx):

    @staticmethod
    def new(p, d, name=None) -> gr_fq_zech_ctx:
        cdef bytes name_b
        cdef char *name_c
        if name is not None:
            name_b = name.encode('utf-8')
            name_c = name_b
        else:
            name_c = NULL
        return gr_fq_zech_ctx._new(p, d, name_c)

    def characteristic(self):
        return self.p

    def degree(self):
        return self.d

    def __repr__(self):
        return f"gr_fq_zech_ctx({self.characteristic()}, {self.degree()})"


@cython.no_gc
cdef class gr_nf_ctx(gr_scalar_ctx):

    @staticmethod
    def new(poly) -> gr_nf_ctx:
        poly = fmpq_poly(poly)
        return gr_nf_ctx._new(poly)

    def modulus(self):
        cdef fmpq_poly poly = fmpq_poly.__new__(fmpq_poly)
        fmpq_poly_init(poly.val)
        fmpq_poly_set(poly.val, self.poly)
        return poly

    def __repr__(self):
        return f"gr_nf_ctx({self.modulus()})"


@cython.no_gc
cdef class gr_nf_fmpz_poly_ctx(gr_scalar_ctx):

    @staticmethod
    def new(poly) -> gr_nf_fmpz_poly_ctx:
        poly = fmpz_poly(poly)
        return gr_nf_fmpz_poly_ctx._new(poly)

    def modulus(self):
        cdef fmpz_poly poly = fmpz_poly.__new__(fmpz_poly)
        fmpz_poly_init(poly.val)
        fmpz_poly_set(poly.val, self.poly)
        return poly

    def __repr__(self):
        return f"gr_nf_fmpz_poly_ctx({self.modulus()})"


@cython.no_gc
cdef class gr_real_qqbar_ctx(gr_scalar_ctx):

    @staticmethod
    def new(deg_limit=-1, bits_limit=-1) -> gr_real_qqbar_ctx:
        assert deg_limit == -1 and bits_limit == -1 or deg_limit > 0 and bits_limit > 0
        return gr_real_qqbar_ctx._new(deg_limit, bits_limit)

    def limits(self):
        return {'deg_limit': self.deg_limit, 'bits_limit': self.bits_limit}

    def __repr__(self):
        return f"gr_real_qqbar_ctx({self.deg_limit}, {self.bits_limit})"


@cython.no_gc
cdef class gr_complex_qqbar_ctx(gr_scalar_ctx):

    @staticmethod
    def new(deg_limit=-1, bits_limit=-1) -> gr_complex_qqbar_ctx:
        assert deg_limit == -1 and bits_limit == -1 or deg_limit > 0 and bits_limit > 0
        return gr_complex_qqbar_ctx._new(deg_limit, bits_limit)

    def limits(self):
        return {'deg_limit': self.deg_limit, 'bits_limit': self.bits_limit}

    def __repr__(self):
        return f"gr_complex_qqbar_ctx({self.deg_limit}, {self.bits_limit})"


@cython.no_gc
cdef class gr_real_ca_ctx(gr_scalar_ctx):

    @staticmethod
    def new(**options) -> gr_real_ca_ctx:
        return gr_real_ca_ctx._new(options)

    def __repr__(self):
        return f"gr_real_ca_ctx({self.options})"


@cython.no_gc
cdef class gr_complex_ca_ctx(gr_scalar_ctx):

    @staticmethod
    def new(**options) -> gr_complex_ca_ctx:
        return gr_complex_ca_ctx._new(options)

    def __repr__(self):
        return f"gr_complex_ca_ctx({self.options})"


@cython.no_gc
cdef class gr_real_algebraic_ca_ctx(gr_scalar_ctx):

    @staticmethod
    def new(**options) -> gr_real_algebraic_ca_ctx:
        return gr_real_algebraic_ca_ctx._new(options)

    def __repr__(self):
        return f"gr_real_algebraic_ca_ctx({self.options})"


@cython.no_gc
cdef class gr_complex_algebraic_ca_ctx(gr_scalar_ctx):

    @staticmethod
    def new(**options) -> gr_complex_algebraic_ca_ctx:
        return gr_complex_algebraic_ca_ctx._new(options)

    def __repr__(self):
        return f"gr_complex_algebraic_ca_ctx({self.options})"


@cython.no_gc
cdef class gr_complex_extended_ca_ctx(gr_scalar_ctx):

    @staticmethod
    def new(**options) -> gr_complex_extended_ca_ctx:
        return gr_complex_extended_ca_ctx._new(options)

    def __repr__(self):
        return f"gr_complex_extended_ca_ctx({self.options})"


@cython.no_gc
cdef class gr_real_float_arf_ctx(gr_scalar_ctx):

    @staticmethod
    def new(prec) -> gr_real_float_arf_ctx:
        return gr_real_float_arf_ctx._new(prec)

    def __repr__(self):
        return f"gr_real_float_arf_ctx({self.prec})"


@cython.no_gc
cdef class gr_complex_float_acf_ctx(gr_scalar_ctx):

    @staticmethod
    def new(prec) -> gr_complex_float_acf_ctx:
        return gr_complex_float_acf_ctx._new(prec)

    def __repr__(self):
        return f"gr_complex_float_acf_ctx({self.prec})"


@cython.no_gc
cdef class gr_real_arb_ctx(gr_scalar_ctx):

    @staticmethod
    def new(prec) -> gr_real_arb_ctx:
        return gr_real_arb_ctx._new(prec)

    def __repr__(self):
        return f"gr_real_arb_ctx({self.prec})"


@cython.no_gc
cdef class gr_complex_acb_ctx(gr_scalar_ctx):

    @staticmethod
    def new(prec) -> gr_complex_acb_ctx:
        return gr_complex_acb_ctx._new(prec)

    def __repr__(self):
        return f"gr_complex_acb_ctx({self.prec})"


@cython.no_gc
cdef class gr(flint_scalar):

    def __dealloc__(self):
        if self._init:
            self._init = False
            gr_heap_clear(self.pval, self.ctx.ctx_t)

    def __repr__(self):
        return self.ctx.to_str(self)

    def is_zero(self):
        return truth_to_py(self._is_zero())

    def is_one(self):
        return truth_to_py(self._is_one())

    def is_neg_one(self):
        return truth_to_py(self._is_neg_one())

    def is_integer(self):
        return truth_to_py(self._is_integer())

    def is_rational(self):
        return truth_to_py(self._is_rational())

    def __bool__(self):
        return not truth_to_bool(self._is_zero())

    def __richcmp__(self, other, int op):
        cdef int err
        cdef gr other_gr
        cdef truth_t res

        if isinstance(other, int):
            other_gr = self.ctx.new_gr()
            err = gr_set_si(other_gr.pval, other, self.ctx.ctx_t)
            if err != GR_SUCCESS:
                raise self.ctx._error(err, "Cannot set gr from int")
        elif isinstance(other, gr):
            other_gr = other
        else:
            raise NotImplementedError("Cannot compare gr with non-gr objects")

        if self.ctx != other_gr.ctx:
            raise NotImplementedError("Cannot compare gr with different contexts")

        if op == 0:
            raise NotImplementedError("Cannot compare gr with <")
            # res = self._lt(other_gr)
        elif op == 1:
            raise NotImplementedError("Cannot compare gr with <=")
            # res = self._le(other_gr)
        elif op == 2:
            res = self._equal(other_gr)
        elif op == 3:
            res = self._equal(other_gr)
            if res == T_TRUE:
                res = T_FALSE
            elif res == T_FALSE:
                res = T_TRUE
        elif op == 4:
            raise NotImplementedError("Cannot compare gr with >")
            # res = self._gt(other_gr)
        elif op == 5:
            raise NotImplementedError("Cannot compare gr with >=")
            # res = self._ge(other_gr)
        else:
            assert False, "Invalid rich comparison operator"

        return truth_to_py(res)

    def __neg__(self) -> gr:
        return self._neg()

    # XXX: Maybe +a should return a copy or for arb it should round...
    # def __pos__(self) -> gr:
    #     return self

    def __add__(self, other) -> gr:
        cdef gr other_gr
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
        if isinstance(other, int):
            return self._add_si(other)
        else:
            return NotImplemented

    def __sub__(self, other) -> gr:
        cdef gr other_gr
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
        if isinstance(other, int):
            return self._neg()._add_si(other)
        else:
            return NotImplemented

    def __mul__(self, other) -> gr:
        cdef gr other_gr
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
        if isinstance(other, int):
            return self._mul_si(other)
        else:
            return NotImplemented

    def __truediv__(self, other) -> gr:
        cdef gr other_gr
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
        if isinstance(other, int):
            return self._div_si(other)._inv()
        else:
            return NotImplemented

    def __pow__(self, other) -> gr:
        if isinstance(other, int):
            return self._pow_si(other)
        else:
            return NotImplemented

    def is_square(self):
        return truth_to_py(self._is_square())

    def sqrt(self):
        return self._sqrt()

    def rsqrt(self):
        return self._rqsrt()

    def gcd(self, other):
        cdef gr other_gr
        if not isinstance(other, gr):
            raise TypeError("gcd when other is not gr.")
        other_gr = other
        if not self.ctx == other_gr.ctx:
            raise TypeError("gcd of gr with different contexts.")
        return self._gcd(other_gr)

    def lcm(self, other):
        cdef gr other_gr
        if not isinstance(other, gr):
            raise TypeError("gcd when other is not gr.")
        other_gr = other
        if not self.ctx == other_gr.ctx:
            raise TypeError("gcd of gr with different contexts.")
        return self._gcd(other_gr)

    def factor(self):
        return self._factor()

    def numer(self) -> gr:
        return self._numerator()

    def denom(self) -> gr:
        return self._denominator()

    def __floor__(self) -> gr:
        return self._floor()

    def __ceil__(self) -> gr:
        return self._ceil()

    def __trunc__(self) -> gr:
        return self._trunc()

    def __round__(self, ndigits: int = 0) -> gr:
        if ndigits != 0:
            raise NotImplementedError("Rounding to a specific number of digits is not supported")
        return self._nint()

    # def __int__(self) -> int:
    #     return self._floor().to_int()

    # def __float__(self) -> float:
    #     return ...

    def __abs__(self) -> gr:
        return self._abs()

    def conjugate(self) -> gr:
        return self._conj()

    @property
    def real(self) -> gr:
        return self._re()

    @property
    def imag(self) -> gr:
        return self._im()

    # XXX: Return -1, 0, 1 as int?
    def sgn(self) -> gr:
        return self._sgn()

    def csgn(self) -> gr:
        return self._csgn()

    def arg(self) -> gr:
        return self._arg()
