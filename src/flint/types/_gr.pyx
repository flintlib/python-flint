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
    # gr_ctx_is_zero_ring,
    gr_ctx_is_exact,
    gr_ctx_is_canonical,
    gr_ctx_has_real_prec,
    gr_ctx_set_real_prec,
    gr_ctx_get_real_prec,
)

from flint.flint_base.flint_base cimport ordering_c_to_py
from flint.flint_base.flint_base import Ordering


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
        """Create a context of a given type.

        This method should be called on a derived class:

        >>> from flint.types._gr import gr_real_float_arf_ctx
        >>> ctx10 = gr_real_float_arf_ctx.new(10)
        >>> ctx10
        gr_real_float_arf_ctx(10)
        >>> ctx10.real_prec
        10
        >>> ctx10(2).sqrt()  # 10 bits
        1.414
        >>> ctx100 = gr_real_float_arf_ctx.new(100)
        >>> ctx100(2).sqrt()  # 100 bits
        1.414213562373095048801688724209
        """
        raise NotImplementedError("Cannot create gr_ctx object directly."
                                  " Use e.g. gr_nmod_ctx.new(n) instead.")

    def __dealloc__(self):
        if self._init:
            self._init = False
            gr_ctx_clear(self.ctx_t)

    @property
    def is_finite(self) -> bool | None:
        """True if the domain is finite (can be ``None`` if unknown).

        >>> from flint.types._gr import gr_fmpz_ctx, gr_nmod_ctx
        >>> R1 = gr_fmpz_ctx
        >>> R2 = gr_nmod_ctx.new(5)
        >>> R1.is_finite, R2.is_finite
        (False, True)
        """
        return truth_to_py(gr_ctx_is_finite(self.ctx_t))

    @property
    def is_multiplicative_group(self) -> bool | None:
        """True if the domain is a multiplicative group (can be ``None`` if unknown)."""
        return truth_to_py(gr_ctx_is_multiplicative_group(self.ctx_t))

    @property
    def is_ring(self) -> bool | None:
        """True if the domain is a ring (can be ``None`` if unknown)."""
        return truth_to_py(gr_ctx_is_ring(self.ctx_t))

    @property
    def is_commutative_ring(self) -> bool | None:
        """True if the domain is a commutative ring (can be ``None`` if unknown)."""
        return truth_to_py(gr_ctx_is_commutative_ring(self.ctx_t))

    @property
    def is_integral_domain(self) -> bool | None:
        """True if the domain is an integral domain (can be ``None`` if unknown).

        >>> from flint.types._gr import gr_nmod_ctx
        >>> R1 = gr_nmod_ctx.new(5)
        >>> R2 = gr_nmod_ctx.new(6)
        >>> R1.is_integral_domain, R2.is_integral_domain
        (True, False)
        """
        return truth_to_py(gr_ctx_is_integral_domain(self.ctx_t))

    @property
    def is_unique_factorization_domain(self) -> bool | None:
        """True if the domain is a unique factorization domain (can be ``None`` if unknown).

        >>> from flint.types._gr import gr_fmpz_ctx
        >>> R1 = gr_fmpz_ctx
        >>> R1.is_unique_factorization_domain
        True
        """
        return truth_to_py(gr_ctx_is_unique_factorization_domain(self.ctx_t))

    @property
    def is_field(self) -> bool | None:
        """True if the domain is a field (can be ``None`` if unknown).

        >>> from flint.types._gr import gr_fmpz_ctx, gr_fmpq_ctx
        >>> R1 = gr_fmpz_ctx
        >>> R2 = gr_fmpq_ctx
        >>> R1.is_field, R2.is_field
        (False, True)
        """
        return truth_to_py(gr_ctx_is_field(self.ctx_t))

    @property
    def is_algebraically_closed(self) -> bool | None:
        """True if the domain is algebraically closed (can be ``None`` if unknown).

        >>> from flint.types._gr import gr_complex_qqbar_ctx, gr_real_qqbar_ctx
        >>> R1 = gr_complex_qqbar_ctx.new()
        >>> R2 = gr_real_qqbar_ctx.new()
        >>> R1.is_algebraically_closed, R2.is_algebraically_closed
        (True, False)
        """
        return truth_to_py(gr_ctx_is_algebraically_closed(self.ctx_t))

    @property
    def is_finite_characteristic(self) -> bool | None:
        """True if the domain has finite characteristic (can be ``None`` if unknown).

        >>> from flint.types._gr import gr_nmod_ctx, gr_fmpz_ctx
        >>> R1 = gr_nmod_ctx.new(5)
        >>> R2 = gr_fmpz_ctx
        >>> R1.is_finite_characteristic, R2.is_finite_characteristic
        (True, False)
        """
        return truth_to_py(gr_ctx_is_finite_characteristic(self.ctx_t))

    @property
    def is_ordered_ring(self) -> bool | None:
        """True if the domain is an ordered ring (can be ``None`` if unknown).

        >>> from flint.types._gr import gr_fmpz_ctx, gr_nmod_ctx
        >>> R1 = gr_fmpz_ctx
        >>> R2 = gr_nmod_ctx.new(5)
        >>> R1.is_ordered_ring, R2.is_ordered_ring
        (True, None)
        """
        return truth_to_py(gr_ctx_is_ordered_ring(self.ctx_t))

    # @property
    # def is_zero_ring(self) -> bool | None:
    #     """True if the domain is a zero ring (can be ``None`` if unknown)."""
    #     return truth_to_py(gr_ctx_is_zero_ring(self.ctx_t))

    @property
    def is_exact(self) -> bool | None:
        """True if the domain is exact (can be ``None`` if unknown).

        >>> from flint.types._gr import gr_fmpz_ctx, gr_real_float_arf_ctx
        >>> R1 = gr_fmpz_ctx
        >>> R2 = gr_real_float_arf_ctx.new(10)
        >>> R1.is_exact, R2.is_exact
        (True, False)
        """
        return truth_to_py(gr_ctx_is_exact(self.ctx_t))

    @property
    def is_canonical(self) -> bool | None:
        """True if elements of the domain are always canonical."""
        return truth_to_py(gr_ctx_is_canonical(self.ctx_t))

    @property
    def has_real_prec(self) -> bool | None:
        """True if the domain has real precision (can be ``None`` if unknown).

        See also :attr:`real_prec`.
        """
        return truth_to_py(gr_ctx_has_real_prec(self.ctx_t))

    @property
    def real_prec(self) -> int:
        """Real precision of the domain.

        >>> from flint.types._gr import gr_real_float_arf_ctx
        >>> ctx = gr_real_float_arf_ctx.new(10)
        >>> ctx.real_prec
        10
        >>> ctx.real_prec = 20
        >>> ctx.real_prec
        20
        """
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
        """Create a new element of the domain.

            >>> from flint.types._gr import gr_fmpz_ctx
            >>> ctx = gr_fmpz_ctx
            >>> ctx(2)
            2
            >>> ctx(18446744073709551615)
            18446744073709551615
        """
        if isinstance(arg, gr):
            if arg.ctx == self:
                return arg
            return self.from_other(arg)
        if type(arg) is int:
            try:
                return self.from_si(arg)
            except OverflowError:
                pass
        if type(arg) is float:
            return self.from_d(arg)
        # TODO: fmpz & fmpq ?
        try:
            return self.from_str(str(arg))
        except AssertionError:
            return self.from_si(int(arg))

    def zero(self) -> gr:
        """Return the zero element of the domain.

        >>> from flint.types._gr import gr_fmpz_ctx
        >>> ctx = gr_fmpz_ctx
        >>> ctx.zero()
        0
        """
        return self._zero()

    def one(self) -> gr:
        """Return the one element of the domain.

        >>> from flint.types._gr import gr_fmpz_ctx
        >>> ctx = gr_fmpz_ctx
        >>> ctx.one()
        1
        """
        return self._one()

    def i(self) -> gr:
        """Return the imaginary unit of the domain (if available).

        >>> from flint.types._gr import gr_fmpzi_ctx
        >>> ctx = gr_fmpzi_ctx
        >>> ctx.i()
        I
        """
        return self._i()

    def pos_inf(self) -> gr:
        """Return the positive infinity element of the domain (if available).

        >>> from flint.types._gr import gr_real_float_arf_ctx
        >>> ctx = gr_real_float_arf_ctx.new(10)
        >>> ctx.pos_inf()
        +inf
        """
        return self._pos_inf()

    def neg_inf(self) -> gr:
        """Return the negative infinity element of the domain (if available).

        >>> from flint.types._gr import gr_real_float_arf_ctx
        >>> ctx = gr_real_float_arf_ctx.new(10)
        >>> ctx.neg_inf()
        -inf
        """
        return self._neg_inf()

    def uinf(self) -> gr:
        """Return the unsigned infinity element of the domain (if available).

        >>> from flint.types._gr import gr_complex_extended_ca_ctx
        >>> ctx = gr_complex_extended_ca_ctx.new()
        >>> ctx.uinf()
        UnsignedInfinity
        """
        return self._uinf()

    def undefined(self) -> gr:
        """Return the undefined element of the domain (if available).

        >>> from flint.types._gr import gr_real_float_arf_ctx
        >>> ctx = gr_real_float_arf_ctx.new(10)
        >>> ctx.undefined()
        nan
        """
        return self._undefined()

    def unknown(self) -> gr:
        """Return the unknown element of the domain (if available).

        >>> from flint.types._gr import gr_real_float_arf_ctx
        >>> ctx = gr_real_float_arf_ctx.new(10)
        >>> ctx.unknown()
        nan
        """
        return self._unknown()

    def gen(self) -> gr:
        """Return the generator of the domain (if available).

        >>> from flint.types._gr import gr_fmpzi_ctx, gr_fq_ctx
        >>> ctx = gr_fmpzi_ctx
        >>> ctx.gen()
        I
        >>> ctx = gr_fq_ctx.new(5, 2)
        >>> ctx.gen()
        a
        """
        return self._gen()

    def gens(self) -> list[gr]:
        """Return the top-level generators of the domain

        # XXX: Does not work with FLINT < 3.1

        # >>> from flint.types._gr import gr_fmpzi_ctx, gr_gr_mpoly_ctx
        # >>> ctx = gr_gr_mpoly_ctx.new(gr_fmpzi_ctx, ['x', 'y'])
        # >>> ctx.gens()
        # [x, y]
        # >>> gr_fmpzi_ctx.gens()
        # [I]
        # >>> ctx.gens_recursive()
        # [I, x, y]
        """
        return self._gens()

    def is_zero(self, x) -> bool | None:
        """
        Returns whether x is equal to the ring element 0.

            >>> from flint.types._gr import gr_real_arb_ctx
            >>> ctx = gr_real_arb_ctx.new(10)
            >>> ctx.is_zero(ctx(0))
            True
            >>> ctx.is_zero(ctx(1))
            False
            >>> ctx.is_zero(ctx("[0 +/- 0.1]"))
        """
        return truth_to_py(self._is_zero(x))

    def is_one(self, x) -> bool | None:
        """
        Returns whether x is equal to the ring element 1.

            >>> from flint.types._gr import gr_real_arb_ctx
            >>> ctx = gr_real_arb_ctx.new(10)
            >>> ctx.is_one(ctx(0))
            False
            >>> ctx.is_one(ctx(1))
            True
            >>> ctx.is_one(ctx("[1 +/- 0.1]"))
        """
        return truth_to_py(self._is_one(x))

    def is_neg_one(self, x) -> bool | None:
        """
        Returns whether x is equal to the ring element -1.

            >>> from flint.types._gr import gr_real_arb_ctx
            >>> ctx = gr_real_arb_ctx.new(10)
            >>> ctx.is_neg_one(ctx(1))
            False
            >>> ctx.is_neg_one(ctx(-1))
            True
            >>> ctx.is_neg_one(ctx("[-1 +/- 0.1]"))
        """
        return truth_to_py(self._is_neg_one(x))

    def equal(self, x, y) -> bool | None:
        """
        Returns whether the elements `x` and `y` are equal.
        """
        return self._equal(self(x), self(y))

    # def is_integer(self, x):
    #     """
    #     Returns whether `x` represents an integer.
    #     """
    #     return self._is_integer(self(x))
    #
    # def is_rational(self, x):
    #     """
    #     Returns whether x represents a rational number.
    #     """
    #     return self._is_rational(self(x))

    def neg(self, x) -> gr:
        """
        Returns `-x`.

            >>> from flint.types._gr import gr_complex_acb_ctx, gr_real_arb_ctx
            >>> arb = gr_real_arb_ctx.new(53); acb = gr_complex_acb_ctx.new(106)
            >>> c = acb("2 + I").sqrt(); c
            ([1.4553466902253548081226618397097 +/- 3.48e-32] + [0.3435607497225124641385657439146 +/- 5.23e-32]*I)
            >>> acb.neg(c)
            ([-1.4553466902253548081226618397097 +/- 3.48e-32] + [-0.3435607497225124641385657439146 +/- 5.23e-32]*I)
        """
        return self._neg(self(x))

    ###
    # Arithmetic methods

    def add(self, x, y) -> gr:
        """
        Returns `x + y`

        """
        if isinstance(x, gr) and isinstance(y, gr):
            if x.ctx == self and y.ctx == self:
                return self._add(x, y)
            if x.ctx == self:
                return self._add_other(x, y)
            if y.ctx == self:
                return self._other_add(x, y)

        if isinstance(x, gr) and x.ctx == self and type(y) is int:
            try:
                return self._add_si(x, y)
            except OverflowError:
                pass

        # NOTE: By default, convert everything to the required
        # context before performing the operation
        return self._add(self(x), self(y))

    def sub(self, x, y) -> gr:
        if isinstance(x, gr) and isinstance(y, gr):
            if x.ctx == self and y.ctx == self:
                return self._sub(x, y)
            if x.ctx == self:
                return self._sub_other(x, y)
            if y.ctx == self:
                return self._other_sub(x, y)

        if isinstance(x, gr) and x.ctx == self and type(y) is int:
            try:
                return self._sub_si(x, y)
            except OverflowError:
                pass

        # NOTE: By default, convert everything to the required
        # context before performing the operation
        return self._sub(self(x), self(y))

    def mul(self, x, y) -> gr:
        if isinstance(x, gr) and isinstance(y, gr):
            if x.ctx == self and y.ctx == self:
                return self._mul(x, y)
            if x.ctx == self:
                return self._mul_other(x, y)
            if y.ctx == self:
                return self._other_mul(x, y)

        if isinstance(x, gr) and x.ctx == self and type(y) is int:
            try:
                return self._mul_si(x, y)
            except OverflowError:
                pass

        # NOTE: By default, convert everything to the required
        # context before performing the operation
        return self._mul(self(x), self(y))

    ###
    # Division

    def is_invertible(self, x) -> bool | None:
        """
        Returns whether `x` has a multiplicative inverse in the present ring, i.e. whether `x` is a unit.
        """
        return self._is_invertible(self(x))

    def inv(self, x) -> gr:
        """
        Returns the multiplicative inverse of `x` in the present ring, if such an element exists.
        """
        return self._inv(self(x))

    def div(self, x, y) -> gr:
        """
        Returns the quotient `x/y`
        """
        if isinstance(x, gr) and isinstance(y, gr):
            if x.ctx == self and y.ctx == self:
                return self._div(x, y)
            if x.ctx == self:
                return self._div_other(x, y)
            if y.ctx == self:
                return self._other_div(x, y)

        if isinstance(x, gr) and x.ctx == self and type(y) is int:
            try:
                return self._div_si(x, y)
            except OverflowError:
                pass

        # NOTE: By default, convert everything to the required
        # context before performing the operation
        return self._div(self(x), self(y))

    def divexact(self, x, y) -> gr:
        """
        Returns the quotient `x/y`, assuming that the quotient is exact in the current context.
        """
        if isinstance(x, gr) and isinstance(y, gr):
            if x.ctx == self and y.ctx == self:
                return self._divexact(x, y)
            if x.ctx == self:
                return self._divexact_other(x, y)
            if y.ctx == self:
                return self._other_divexact(x, y)

        if isinstance(x, gr) and x.ctx == self and type(y) is int:
            try:
                return self._divexact_si(x, y)
            except OverflowError:
                pass

        # NOTE: By default, convert everything to the required
        # context before performing the operation
        return self._divexact(self(x), self(y))

    def div_nonunique(self, x, y) -> gr:
        """
        Returns an arbitrary solution `q` of the equation `x = qy`.
        """
        return self._div_nonunique(self(x), self(y))

    def divides(self, d, x) -> bool | None:
        """
        Returns whether `d | x`; that is, whether there is an element `q` such that `x = qd`.
        """
        return self._divides(self(d), self(x))

    def euclidean_div(self, x, y) -> gr:
        return self._euclidean_div(self(x), self(y))

    def euclidean_rem(self, x, y) -> gr:
        return self._euclidean_rem(self(x), self(y))

    def euclidean_divrem(self, x, y) -> tuple[gr, gr]:
        return self._euclidean_divrem(self(x), self(y))

    ###
    # Powering

    def pow(self, x, y) -> gr:
        """
        Returns the power x^y, the interpretation of which depends on the ring when y ∉ ℤ
        """
        if isinstance(x, gr) and isinstance(y, gr):
            if x.ctx == self and y.ctx == self:
                return self._pow(x, y)
            if x.ctx == self:
                return self._pow_other(x, y)
            if y.ctx == self:
                return self._other_pow(x, y)

        if isinstance(x, gr) and x.ctx == self and type(y) is int:
            try:
                return self._pow_si(x, y)
            except OverflowError:
                pass

        # NOTE: By default, convert everything to the required
        # context before performing the operation
        return self._pow(self(x), self(y))

    ###
    # Square Roots

    def is_square(self, x) -> bool | None:
        """
        Returns whether x is a perfect square in the context.
        """
        return self._is_square(self(x))

    def sqrt(self, x) -> gr:
        """
        Returns the square root of `x`.
        """
        return self._sqrt(self(x))

    def rsqrt(self, x) -> gr:
        """
        Returns the reciprocal square root of `x`.
        """
        return self._rsqrt(self(x))

    ###
    # Greatest Common Divisors

    def gcd(self, x, y) -> gr:
        """
        Returns a greatest common divisor (GCD) of `x` and `y`.
        Since the GCD is unique only up to multiplication by a unit,
        an implementation-defined representative is chosen.
        """
        return self._gcd(self(x), self(y))

    def lcm(self, x, y) -> gr:
        """
        Returns a least common multiple (LCM) of x and y.
        Since the LCM is unique only up to multiplication by a unit,
        an implementation-defined representative is chosen.
        """
        return self._lcm(self(x), self(y))

    ###
    # Factorization

    def factor(self, x) -> tuple[gr, list[tuple[gr, int]]]:
        """
        Given an element of the context, this returns a factorization ``(c, f, e)``:
            ``x = c f₁^e₁ ··· fₙ^eₙ``, where ``fₖ`` will be irreducible or prime depending on the ring.
        The prefactor c stores a unit, sign or coefficient.
        Note that c is an element of the same ring as x.
        """
        return self._factor(self(x))

    ###
    # Fractions

    def numerator(self, x) -> gr:
        """
        Return a numerator p such that x = p / q.
        For typical fraction fields, the denominator will be minimal and canonical.
        However, some rings may return an arbitrary denominator as long as the numerator matches.
        The default implementations simply returns p = x.
        """
        return self._numerator(self(x))

    def denominator(self, x) -> gr:
        """
        Return a denominator q such that x = p / q.
        For typical fraction fields, the denominator will be minimal and canonical.
        However, some rings may return an arbitrary denominator as long as the numerator matches.
        The default implementations simply returns q = 1.
        """
        return self._denominator(self(x))

    ###
    # Integer and Complex parts

    def floor(self, x) -> gr:
        return self._floor(self(x))

    def ceil(self, x) -> gr:
        return self._ceil(self(x))

    def trunc(self, x) -> gr:
        return self._trunc(self(x))

    def nint(self, x) -> gr:
        return self._nint(self(x))

    def abs(self, x) -> gr:
        return self._abs(self(x))

    def conj(self, x) -> gr:
        return self._conj(self(x))

    def re(self, x) -> gr:
        return self._re(self(x))

    def im(self, x) -> gr:
        return self._im(self(x))

    def sgn(self, x) -> gr:
        return self._sgn(self(x))

    def csgn(self, x) -> gr:
        return self._csgn(self(x))

    def arg(self, x) -> gr:
        return self._arg(self(x))

    ###
    # Ordering methods

    def cmp(self, x, y) -> int:
        """
        Returns:
        - -1 if x < y
        -  0 if x = y
        -  1 if x > y
        """
        if isinstance(x, gr) and isinstance(y, gr):
            if x.ctx == self and y.ctx != self:
                return self._cmp_other(x, y)

            if y.ctx == self and x.ctx != self:
                return -self._cmp_other(y, x)

            if x.ctx == self == y.ctx:
                return self._cmp(x, y)

        return self._cmp(self(x), self(y))

    def cmpabs(self, x, y) -> int:
        """
        Returns:
        - -1 if `|x| < |y|`
        -  0 if `|x| = |y|`
        -  1 if `|x| > |y|`
        """
        if isinstance(x, gr) and isinstance(y, gr):
            if x.ctx == self and y.ctx != self:
                return self._cmpabs_other(x, y)

            if y.ctx == self and x.ctx != self:
                return -self._cmpabs_other(y, x)

            if x.ctx == self == y.ctx:
                return self._cmpabs(x, y)

        return self._cmpabs(self(x), self(y))

    def le(self, x, y) -> bool | None:
        return truth_to_py(self._le(self(x), self(y)))

    def abs_le(self, x, y) -> bool | None:
        return truth_to_py(self._abs_le(self(x), self(y)))

    def lt(self, x, y) -> bool | None:
        return truth_to_py(self._lt(self(x), self(y)))

    def abs_lt(self, x, y) -> bool | None:
        return truth_to_py(self._abs_lt(self(x), self(y)))

    def ge(self, x, y) -> bool | None:
        return truth_to_py(self._ge(self(x), self(y)))

    def abs_ge(self, x, y) -> bool | None:
        return truth_to_py(self._abs_ge(self(x), self(y)))

    def gt(self, x, y) -> bool | None:
        return truth_to_py(self._gt(self(x), self(y)))

    def abs_gt(self, x, y) -> bool | None:
        return truth_to_py(self._abs_gt(self(x), self(y)))

    def min(self, x, y) -> gr:
        return self._min(self(x), self(y))

    def max(self, x, y) -> gr:
        return self._max(self(x), self(y))

    ###
    # Array-API wrappers

    def divide(self, x, y) -> gr:
        return self.div(x, y)

    def greater(self, x, y):
        return self.gt(x, y)

    def greater_equal(self, x, y):
        return self.ge(x, y)

    def less(self, x, y):
        return self.lt(x, y)

    def less_equal(self, x, y):
        return self.le(x, y)

    def imag(self, x):
        return self.im(x)

    def real(self, x):
        return self.re(x)

    def maximum(self, x, y):
        return self.max(x, y)

    def minimum(self, x, y):
        return self.min(x, y)

    def multiply(self, x, y):
        return self.mul(x, y)

    def negative(self, x):
        return self.neg(x)

    def not_equal(self, x, y):
        eq = self.equal(x, y)
        if eq is None:
            return None
        return not eq

cdef class gr_scalar_ctx(gr_ctx):
    """Base class for all scalar contexts."""
    pass


cdef class gr_poly_ctx(gr_ctx):
    """Base class for dense univariate polynomial contexts."""
    pass


cdef class gr_mpoly_ctx(gr_ctx):
    """Base class for dense multivariate polynomial contexts."""
    pass


# cdef class gr_matrix_domain_ctx(gr_ctx):
#     pass


# cdef class gr_matrix_space_ctx(gr_ctx):
#     pass


# cdef class gr_matrix_ring_ctx(gr_ctx):
#     pass


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
    r"""Context for arbitrary precision integers, `\mathbb{Z}`."""

    @staticmethod
    def new() -> _gr_fmpz_ctx:
        """The global context for arbitrary precision integers.

        >>> from flint.types._gr import gr_fmpz_ctx as Z
        >>> Z
        gr_fmpz_ctx
        >>> Z(10)
        10
        >>> Z(10) + Z(20)
        30
        >>> Z(10) * Z(20)
        200
        >>> Z(10) / Z(20)
        Traceback (most recent call last):
            ...
        AssertionError: Failed to divide gr objects
        >>> Z(10).factor()
        (1, [(2, 1), (5, 1)])
        >>> Z(9).sqrt()
        3
        """
        return gr_fmpz_ctx_c

    def __repr__(self):
        return "gr_fmpz_ctx"


@cython.no_gc
cdef class _gr_fmpq_ctx(gr_scalar_ctx):
    r"""Context for arbitrary precision rationals, `\mathbb{Q}`."""

    @staticmethod
    def new() -> _gr_fmpq_ctx:
        """The global context for arbitrary precision rationals.

        >>> from flint.types._gr import gr_fmpq_ctx as Q
        >>> Q
        gr_fmpq_ctx
        >>> Q(2) / 3
        2/3
        """
        return gr_fmpq_ctx_c

    def __repr__(self):
        return "gr_fmpq_ctx"


@cython.no_gc
cdef class _gr_fmpzi_ctx(gr_scalar_ctx):
    r"""Context for Gaussian integers, `\mathbb{Z}[i]`. """

    @staticmethod
    def new() -> _gr_fmpzi_ctx:
        """The global context for Gaussian integers.

        >>> from flint.types._gr import gr_fmpzi_ctx as ZI
        >>> ZI
        gr_fmpzi_ctx
        >>> ZI.i()
        I
        >>> ZI.gen()
        I
        >>> I = ZI.i()
        >>> (2 + 3*I) ** 2
        (-5+12*I)
        """
        return gr_fmpzi_ctx_c

    def __repr__(self):
        return "gr_fmpzi_ctx"


@cython.no_gc
cdef class _gr_fexpr_ctx(gr_scalar_ctx):
    """Context for symbolic expressions. """

    @staticmethod
    def new() -> _gr_fexpr_ctx:
        """The global context for symbolic expressions.

        >>> from flint.types._gr import gr_fexpr_ctx as EX
        >>> EX
        gr_fexpr_ctx
        >>> (EX(1) + EX(2)) / 3
        Div(Add(1, 2), 3)
        """
        return gr_fexpr_ctx_c

    def __repr__(self):
        return "gr_fexpr_ctx"


@cython.no_gc
cdef class gr_nmod_ctx(gr_scalar_ctx):
    r"""Context for integers modulo `n`, `\mathbb{Z}/n\mathbb{Z}`.

    The modulus `n` must be word-sized. Use :class:`gr_fmpz_mod_ctx` for
    arbitrary precision modulus.
    """
    @staticmethod
    def new(ulong n) -> gr_nmod_ctx:
        """Create a new context for integers modulo `n`.

        >>> from flint.types._gr import gr_nmod_ctx
        >>> Z5 = gr_nmod_ctx.new(5)
        >>> Z5
        gr_nmod_ctx(5)
        >>> Z5(2) + Z5(3)
        0
        >>> Z5.modulus()
        5
        >>> Z5.is_prime
        True
        """
        return gr_nmod_ctx._new(n)

    @property
    def is_prime(self) -> bool:
        return self.is_field

    def modulus(self):
        return self.n

    def __repr__(self):
        return f"gr_nmod_ctx({self.n})"


@cython.no_gc
cdef class gr_fmpz_mod_ctx(gr_scalar_ctx):
    r"""Context for integers modulo `n`, `\mathbb{Z}/n\mathbb{Z}`.

    The modulus `n` can be arbitrary precision. See also :class:`gr_nmod_ctx`
    for word-sized modulus.
    """

    @staticmethod
    def new(n) -> gr_fmpz_mod_ctx:
        """Create a new context for integers modulo `n`.

        >>> from flint.types._gr import gr_fmpz_mod_ctx
        >>> Z64 = gr_fmpz_mod_ctx.new(2**64 + 1)
        >>> Z64
        gr_fmpz_mod_ctx(18446744073709551617)
        >>> Z64(2)**64 + 2
        1
        >>> Z64.is_prime
        False
        >>> Z64.modulus()
        18446744073709551617
        """
        n = fmpz(n)
        return gr_fmpz_mod_ctx._new(n)

    @property
    def is_prime(self) -> bool:
        return self.is_field

    def modulus(self):
        cdef fmpz n = fmpz.__new__(fmpz)
        fmpz_init_set(n.val, self.n)
        return n

    def __repr__(self):
        return f"gr_fmpz_mod_ctx({self.modulus()})"


@cython.no_gc
cdef class gr_fq_ctx(gr_scalar_ctx):
    r"""Context for finite fields, `\mathbb{F}_q` where `q = p^d`.

    The characteristic `p` must be a prime number and the degree `d` must be
    positive. The characteristic `p` can be arbitrary precision. Use
    :class:`gr_fq_nmod_ctx` for word-sized characteristic. For very small
    characteristic, use :class:`gr_fq_zech_ctx`.
    """

    @staticmethod
    def new(p, d, name=None) -> gr_fq_ctx:
        """Create a new context for finite fields.

        >>> from flint.types._gr import gr_fq_ctx
        >>> F9 = gr_fq_ctx.new(3, 2)
        >>> F9
        gr_fq_ctx(3, 2)
        >>> F9(2) + F9(3)
        2
        >>> F9.characteristic()
        3
        >>> F9.degree()
        2
        >>> F9.gen()
        a
        >>> a = F9.gen()
        >>> (1 + a) ** 2 + a
        a+2
        """
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
    r"""Context for finite fields, `\mathbb{F}_q` where `q = p^d`.

    The characteristic `p` must be a prime number and the degree `d` must be
    positive. The characteristic `p` must be word-sized. Use :class:`gr_fq_ctx`
    for arbitrary precision characteristic. For very small characteristic, use
    :class:`gr_fq_zech_ctx`.
    """

    @staticmethod
    def new(p, d, name=None) -> gr_fq_nmod_ctx:
        """Create a new context for finite fields.

        # XXX: Does not work with FLINT < 3.1
        # >>> from flint.types._gr import gr_fq_nmod_ctx
        # >>> F9 = gr_fq_nmod_ctx.new(3, 2)
        # >>> F9
        # gr_fq_nmod_ctx(3, 2)
        # >>> F9(2) + F9(3)
        # 2
        # >>> F9.characteristic()
        # 3
        # >>> F9.degree()
        # 2
        # >>> F9.gen()
        # a
        # >>> a = F9.gen()
        # >>> (1 + a) ** 2 + a
        # a+2
        """
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
    r"""Context for finite fields, `\mathbb{F}_q` where `q = p^d`.

    The characteristic `p` must be a prime number and the degree `d` must be
    positive. The characteristic `p` must be small. Use :class:`gr_nmod_ctx`
    for word-sized characteristic and :class:`gr_fq_ctx` for arbitrary
    precision characteristic.

    .. warning::
        This type is possibly not working correctly. Use with caution.
    """

    @staticmethod
    def new(p, d, name=None) -> gr_fq_zech_ctx:
        """Create a new context for finite fields with small characteristic.

        # XXX: Does not work with FLINT < 3.1
        # >>> from flint.types._gr import gr_fq_zech_ctx
        # >>> F9 = gr_fq_zech_ctx.new(3, 2)
        # >>> F9
        # gr_fq_zech_ctx(3, 2)
        # >>> F9(2) + F9(3) # XXX: Is this correct?
        # a^4
        # >>> F9.characteristic()
        # 3
        # >>> F9.degree()
        # 2
        # >>> F9.gen()
        # a^1
        # >>> a = F9.gen()
        # >>> (1 + a) ** 2 + a  # doctest: +SKIP
        # a+2
        """
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
    r"""Context for number fields, `\mathbb{Q}(\alpha)`.

    The number field is defined by a minimal polynomial `f(x)` in
    `\mathbb{Q}[x]` and represents the field `\mathbb{Q}(\alpha)` where
    `\alpha` is a root of `f(x)`. The minimal polynomial `f(x)` must be
    irreducible.
    """

    @staticmethod
    def new(poly) -> gr_nf_ctx:
        """Create a new context for number fields.

        >>> from flint.types._gr import gr_nf_ctx
        >>> Qa = gr_nf_ctx.new([-2, 0, 1])
        >>> Qa
        gr_nf_ctx(x^2 + (-2))
        >>> Qa.modulus()
        x^2 + (-2)
        >>> a = Qa.gen()
        >>> a
        a
        >>> a**2
        2
        >>> (1 + a) ** 2
        2*a+3
        >>> (1 + a) / 2
        1/2*a+1/2
        """
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
    r"""Context for number fields, `\mathbb{Q}(\alpha)`.

    The number field is defined by a minimal polynomial `f(x)` in
    `\mathbb{Z}[x]` and represents the ring `\mathbb{Z}[\alpha]` where
    `\alpha` is a root of `f(x)`. The minimal polynomial `f(x)` must be
    irreducible.
    """

    @staticmethod
    def new(poly) -> gr_nf_fmpz_poly_ctx:
        """Create a new context for number fields.

        >>> from flint.types._gr import gr_nf_fmpz_poly_ctx
        >>> Qa = gr_nf_fmpz_poly_ctx.new([-2, 0, 1])
        >>> Qa
        gr_nf_fmpz_poly_ctx(x^2 + (-2))
        >>> Qa.modulus()
        x^2 + (-2)
        >>> a = Qa.gen()
        >>> a
        a
        >>> a**2
        2
        >>> (1 + a) ** 2
        2*a+3
        """
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
    r"""Context for real algebraic numbers, `\mathbb{R} \cap \overline{\mathbb{Q}}`.

    The real algebraic numbers are a subset of the real numbers that are roots
    of non-zero polynomials with integer coefficients. See
    :class:`gr_complex_qqbar_ctx` for the full algebraic closure of the
    rationals. See :class:`gr_nf_ctx` for number fields with a single
    generator.
    """
    @staticmethod
    def new(deg_limit=-1, bits_limit=-1) -> gr_real_qqbar_ctx:
        """Create a new context for real algebraic numbers.

        >>> from flint.types._gr import gr_real_qqbar_ctx
        >>> R = gr_real_qqbar_ctx.new()
        >>> R
        gr_real_qqbar_ctx(-1, -1)
        >>> R(2).sqrt()
        Root a = 1.41421 of a^2-2
        """
        assert deg_limit == -1 and bits_limit == -1 or deg_limit > 0 and bits_limit > 0
        return gr_real_qqbar_ctx._new(deg_limit, bits_limit)

    def limits(self):
        return {'deg_limit': self.deg_limit, 'bits_limit': self.bits_limit}

    def __repr__(self):
        return f"gr_real_qqbar_ctx({self.deg_limit}, {self.bits_limit})"


@cython.no_gc
cdef class gr_complex_qqbar_ctx(gr_scalar_ctx):
    r"""Context for algebraic numbers, `\overline{\mathbb{Q}}`.

    The algebraic numbers are a field that contains the rationals and all roots
    of non-zero polynomials with integer coefficients. See
    :class:`gr_real_qqbar_ctx` for the real algebraic numbers,
    :class:`gr_nf_ctx` for number fields with a single generator.
    """
    @staticmethod
    def new(deg_limit=-1, bits_limit=-1) -> gr_complex_qqbar_ctx:
        """Create a new context for algebraic numbers.

        >>> from flint.types._gr import gr_complex_qqbar_ctx
        >>> C = gr_complex_qqbar_ctx.new()
        >>> C
        gr_complex_qqbar_ctx(-1, -1)
        >>> C(-2).sqrt()
        Root a = 1.41421*I of a^2+2
        """
        assert deg_limit == -1 and bits_limit == -1 or deg_limit > 0 and bits_limit > 0
        return gr_complex_qqbar_ctx._new(deg_limit, bits_limit)

    def limits(self):
        return {'deg_limit': self.deg_limit, 'bits_limit': self.bits_limit}

    def __repr__(self):
        return f"gr_complex_qqbar_ctx({self.deg_limit}, {self.bits_limit})"


@cython.no_gc
cdef class gr_real_ca_ctx(gr_scalar_ctx):
    r"""Context for calcium exact real numbers, `\mathbb{R}`."""

    @staticmethod
    def new(**options) -> gr_real_ca_ctx:
        """Create a new context for calcium exact real numbers.

        >>> from flint.types._gr import gr_real_ca_ctx
        >>> R = gr_real_ca_ctx.new()
        >>> R
        gr_real_ca_ctx({})
        >>> R(2).sqrt()
        1.41421 {a where a = 1.41421 [a^2-2=0]}
        """
        return gr_real_ca_ctx._new(options)

    def __repr__(self):
        return f"gr_real_ca_ctx({self.options})"


@cython.no_gc
cdef class gr_complex_ca_ctx(gr_scalar_ctx):
    r"""Context for calcium exact complex numbers, `\mathbb{C}`."""

    @staticmethod
    def new(**options) -> gr_complex_ca_ctx:
        """Create a new context for calcium exact complex numbers.

        >>> from flint.types._gr import gr_complex_ca_ctx
        >>> C = gr_complex_ca_ctx.new()
        >>> C
        gr_complex_ca_ctx({})
        >>> C(2).sqrt()
        1.41421 {a where a = 1.41421 [a^2-2=0]}
        """
        return gr_complex_ca_ctx._new(options)

    def __repr__(self):
        return f"gr_complex_ca_ctx({self.options})"


@cython.no_gc
cdef class gr_real_algebraic_ca_ctx(gr_scalar_ctx):
    r"""Context for calcium exact real algebraic numbers, `\mathbb{R} \cap \overline{\mathbb{Q}}`."""

    @staticmethod
    def new(**options) -> gr_real_algebraic_ca_ctx:
        """Create a new context for calcium exact real algebraic numbers.

        >>> from flint.types._gr import gr_real_algebraic_ca_ctx
        >>> R = gr_real_algebraic_ca_ctx.new()
        >>> R
        gr_real_algebraic_ca_ctx({})
        >>> R(2).sqrt()
        1.41421 {a where a = 1.41421 [a^2-2=0]}
        """
        return gr_real_algebraic_ca_ctx._new(options)

    def __repr__(self):
        return f"gr_real_algebraic_ca_ctx({self.options})"


@cython.no_gc
cdef class gr_complex_algebraic_ca_ctx(gr_scalar_ctx):
    r"""Context for calcium exact complex algebraic numbers, `\overline{\mathbb{Q}}`."""

    @staticmethod
    def new(**options) -> gr_complex_algebraic_ca_ctx:
        """Create a new context for calcium exact complex algebraic numbers.

        >>> from flint.types._gr import gr_complex_algebraic_ca_ctx
        >>> C = gr_complex_algebraic_ca_ctx.new()
        >>> C
        gr_complex_algebraic_ca_ctx({})
        >>> C(2).sqrt()
        1.41421 {a where a = 1.41421 [a^2-2=0]}
        """
        return gr_complex_algebraic_ca_ctx._new(options)

    def __repr__(self):
        return f"gr_complex_algebraic_ca_ctx({self.options})"


@cython.no_gc
cdef class gr_complex_extended_ca_ctx(gr_scalar_ctx):
    r"""Context for calcium exact extended complex numbers, `\mathbb{C} \cup \{\infty\}`."""

    @staticmethod
    def new(**options) -> gr_complex_extended_ca_ctx:
        """Create a new context for calcium exact extended complex numbers.

        >>> from flint.types._gr import gr_complex_extended_ca_ctx
        >>> C = gr_complex_extended_ca_ctx.new()
        >>> C
        gr_complex_extended_ca_ctx({})
        >>> C(1) / 0
        UnsignedInfinity
        """
        return gr_complex_extended_ca_ctx._new(options)

    def __repr__(self):
        return f"gr_complex_extended_ca_ctx({self.options})"


@cython.no_gc
cdef class gr_real_float_arf_ctx(gr_scalar_ctx):
    r"""Context for arbitrary precision approximate real numbers, `\mathbb{R}`."""

    @staticmethod
    def new(prec) -> gr_real_float_arf_ctx:
        """Create a new context for arbitrary precision approximate real numbers.

        >>> from flint.types._gr import gr_real_float_arf_ctx
        >>> R = gr_real_float_arf_ctx.new(10)
        >>> R
        gr_real_float_arf_ctx(10)
        >>> R(2).sqrt()
        1.414
        >>> R.real_prec = 20
        >>> R(2).sqrt()
        1.414213
        """
        return gr_real_float_arf_ctx._new(prec)

    def __repr__(self):
        return f"gr_real_float_arf_ctx({self.prec})"


@cython.no_gc
cdef class gr_complex_float_acf_ctx(gr_scalar_ctx):
    r"""Context for arbitrary precision approximate complex numbers, `\mathbb{C}`."""

    @staticmethod
    def new(prec) -> gr_complex_float_acf_ctx:
        """Create a new context for arbitrary precision approximate complex numbers.

        >>> from flint.types._gr import gr_complex_float_acf_ctx
        >>> C = gr_complex_float_acf_ctx.new(10)
        >>> C
        gr_complex_float_acf_ctx(10)
        >>> C(-2).sqrt()
        1.414*I
        >>> C.real_prec = 20
        >>> C(-2).sqrt()
        1.414213*I
        """
        return gr_complex_float_acf_ctx._new(prec)

    def __repr__(self):
        return f"gr_complex_float_acf_ctx({self.prec})"


@cython.no_gc
cdef class gr_real_arb_ctx(gr_scalar_ctx):
    r"""Context for Arb arbitrary precision real ball arithmetic, `\mathbb{R}`."""

    @staticmethod
    def new(prec) -> gr_real_arb_ctx:
        """Create a new context for arbitrary precision real ball arithmetic.

        >>> from flint.types._gr import gr_real_arb_ctx
        >>> R = gr_real_arb_ctx.new(10)
        >>> R
        gr_real_arb_ctx(10)
        >>> R(2).sqrt()
        [1.41 +/- 6.02e-3]
        >>> R.real_prec = 20
        >>> R(2).sqrt()
        [1.41421 +/- 5.09e-6]
        """
        return gr_real_arb_ctx._new(prec)

    def __repr__(self):
        return f"gr_real_arb_ctx({self.prec})"


@cython.no_gc
cdef class gr_complex_acb_ctx(gr_scalar_ctx):
    r"""Context for Arb arbitrary precision complex ball arithmetic, `\mathbb{C}`."""

    @staticmethod
    def new(prec) -> gr_complex_acb_ctx:
        """Create a new context for arbitrary precision complex ball arithmetic.

        >>> from flint.types._gr import gr_complex_acb_ctx
        >>> C = gr_complex_acb_ctx.new(10)
        >>> C
        gr_complex_acb_ctx(10)
        >>> C(-2).sqrt()
        [1.41 +/- 6.02e-3]*I
        >>> C.real_prec = 20
        >>> C(-2).sqrt()
        [1.41421 +/- 5.09e-6]*I
        """
        return gr_complex_acb_ctx._new(prec)

    def __repr__(self):
        return f"gr_complex_acb_ctx({self.prec})"


@cython.no_gc
cdef class gr_gr_poly_ctx(gr_poly_ctx):
    r"""Context for dense univariate polynomial rings, `\mathbb{D}[x]`."""

    @staticmethod
    def new(base_ring) -> gr_gr_poly_ctx:
        """Create a new context for dense univariate polynomial rings.

        >>> from flint.types._gr import gr_fmpz_ctx, gr_gr_poly_ctx
        >>> Z = gr_fmpz_ctx
        >>> R = gr_gr_poly_ctx.new(Z)
        >>> R
        gr_gr_poly_ctx(gr_fmpz_ctx)
        >>> R.base_ring
        gr_fmpz_ctx
        >>> R.gen()
        x
        >>> x = R.gen()
        >>> (1 + x) ** 2
        1 + 2*x + x^2
        """
        return gr_gr_poly_ctx._new(base_ring)

    def __repr__(self):
        return f"gr_gr_poly_ctx({self.base_ring})"

    @property
    def base_ring(self):
        return self.base_ctx


@cython.no_gc
cdef class gr_gr_mpoly_ctx(gr_mpoly_ctx):
    r"""Context for dense multivariate polynomial rings, `\mathbb{D}[x, y, ...]`."""

    @staticmethod
    def new(base_ring, names, order=None) -> gr_gr_mpoly_ctx:
        """Create a new context for dense multivariate polynomial rings.

        # >>> from flint.types._gr import gr_fmpzi_ctx, gr_gr_mpoly_ctx
        # >>> ZI = gr_fmpzi_ctx
        # >>> R = gr_gr_mpoly_ctx.new(ZI, ['x', 'y'])
        # >>> R
        # gr_gr_mpoly_ctx(gr_fmpzi_ctx, ('x', 'y'), Ordering.lex)
        # >>> R.base_ring
        # gr_fmpzi_ctx
        # >>> R.names
        # ('x', 'y')
        # >>> R.nvars
        # 2
        # >>> R.order
        # <Ordering.lex: 'lex'>
        # >>> R.gens()
        # [x, y]
        # >>> I, [x, y] = ZI.gen(), R.gens()
        # >>> (I + x + y) ** 2
        # x^2 + 2*x*y + (2*I)*x + y^2 + (2*I)*y - 1
        """
        if order is None:
            order = Ordering.lex
        return gr_gr_mpoly_ctx._new(base_ring, tuple(names), order)

    def __repr__(self):
        return f"gr_gr_mpoly_ctx({self.base_ring}, {self.names}, {self.order})"

    @property
    def base_ring(self):
        return self.base_ctx

    @property
    def names(self):
        return self._names

    @property
    def nvars(self):
        return self._nvars

    @property
    def order(self):
        return ordering_c_to_py(self._order)


# @cython.no_gc
# cdef class gr_mpoly_q_ctx(gr_mpoly_ctx):
#
#     @staticmethod
#     def new(nvars, order=None) -> gr_gr_mpoly_q_ctx:
#         if order is None:
#             order = Ordering.lex
#         return gr_mpoly_q_ctx._new(nvars, order)
#
#     def __repr__(self):
#         return f"gr_mpoly_q_ctx({self.nvars}, {self.order})"
#
#     @property
#     def nvars(self):
#         return self._nvars
#
#     @property
#     def order(self):
#         return ordering_c_to_py(self._order)


# @cython.no_gc
# cdef class gr_series_mod_gr_poly_ctx(gr_ctx):
#
#     @staticmethod
#     def new(base_ring, n) -> gr_series_mod_gr_poly_ctx:
#         return gr_series_mod_gr_poly_ctx._new(base_ring, n)
#
#     def __repr__(self):
#         return f"gr_series_mod_gr_poly_ctx({self.base_ring}, {self.nvars})"
#
#     @property
#     def base_ring(self):
#         return self.base_ctx
#
#     @property
#     def n(self):
#         return self._n


@cython.no_gc
cdef class gr_series_ctx(gr_ctx):
    r"""Context for series with precision `n`, `\mathbb{D}[[x]]`."""

    @staticmethod
    def new(base_ring, prec) -> gr_series_ctx:
        """Create a new context for series with precision `n`.

        >>> from flint.types._gr import gr_fmpz_ctx, gr_series_ctx
        >>> Z = gr_fmpz_ctx
        >>> R = gr_series_ctx.new(Z, 10)
        >>> R
        gr_series_ctx(gr_fmpz_ctx, 10)
        >>> R.base_ring
        gr_fmpz_ctx
        >>> R.prec
        10
        >>> x = R.gen()
        >>> 1 / (1 - x)
        1 + x + x^2 + x^3 + x^4 + x^5 + x^6 + x^7 + x^8 + x^9 + O(x^10)
        """
        return gr_series_ctx._new(base_ring, prec)

    def __repr__(self):
        return f"gr_series_ctx({self.base_ring}, {self.prec})"

    @property
    def base_ring(self):
        return self.base_ctx

    @property
    def prec(self):
        return self._prec


@cython.no_gc
cdef class gr(flint_scalar):
    """Type of elements of :class:`gr_ctx` domains."""

    def __dealloc__(self):
        if self._init:
            self._init = False
            gr_heap_clear(self.pval, self.ctx.ctx_t)

    def __repr__(self):
        return self.ctx.to_str(self)

    def parent(self) -> gr_ctx:
        """
        Return the parent context.

            >>> from flint.types._gr import gr_complex_acb_ctx
            >>> acb = gr_complex_acb_ctx.new(53)
            >>> x = acb("2")
            >>> x.parent()
            gr_complex_acb_ctx(53)
        """
        return self.ctx

    def is_zero(self):
        """Return whether the element is zero (may return ``None``).

        >>> from flint.types._gr import gr_fmpz_ctx as Z
        >>> Z(0).is_zero()
        True
        """
        return truth_to_py(self._is_zero())

    def is_one(self):
        """Return whether the element is one (may return ``None``).

        >>> from flint.types._gr import gr_fmpz_ctx as Z
        >>> Z(1).is_one()
        True
        """
        return truth_to_py(self._is_one())

    def is_neg_one(self):
        """Return whether the element is negative one (may return ``None``).

        >>> from flint.types._gr import gr_fmpz_ctx as Z
        >>> Z(-1).is_neg_one()
        True
        """
        return truth_to_py(self._is_neg_one())

#     def is_integer(self):
#         """Return whether the element is an integer (may return ``None``).
#
#         >>> from flint.types._gr import gr_fmpq_ctx as Q
#         >>> Q(2).is_integer()
#         True
#         >>> Q(2, 3).is_integer()
#         False
#         """
#         return truth_to_py(self._is_integer())
#
#     def is_rational(self):
#         """Return whether the element is a rational number (may return ``None``).
#
#         >>> from flint.types._gr import gr_nf_ctx
#         >>> Qa = gr_nf_ctx.new([-2, 0, 1])
#         >>> a = Qa.gen()
#         >>> a.is_rational()
#         False
#         >>> (a**2).is_rational()
#         True
#         """
#         return truth_to_py(self._is_rational())

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
            if self.ctx != other_gr.ctx:
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
        return self.ctx.pow(self, other)

    def is_square(self):
        """Return whether the element is a square (may return ``None``).

        >>> from flint.types._gr import gr_fmpq_ctx as Q
        >>> Q(2).is_square()
        False
        >>> Q(4).is_square()
        True
        >>> Q(4).sqrt()
        2
        """
        return truth_to_py(self.ctx.is_square(self))

    def sqrt(self):
        """Return the square root of the element if it exists.

        >>> from flint.types._gr import gr_fmpz_ctx as Z
        >>> Z(4).sqrt()
        2
        """
        return self.ctx.sqrt(self)

    def rsqrt(self):
        """Return the reciprocal square root of the element if it exists.

        >>> from flint.types._gr import gr_fmpq_ctx as Q
        >>> Q(4).rsqrt()
        1/2
        """
        return self.ctx.rsqrt(self)

    def gcd(self, other):
        """Return the greatest common divisor of two elements.

        >>> from flint.types._gr import gr_fmpz_ctx as Z
        >>> Z(4).gcd(Z(6))
        2
        """
        cdef gr other_gr
        if not isinstance(other, gr):
            raise TypeError("gcd when other is not gr.")
        other_gr = other
        if not self.ctx == other_gr.ctx:
            raise TypeError("gcd of gr with different contexts.")
        return self._gcd(other_gr)

    def lcm(self, other):
        """Return the least common multiple of two elements.

        >>> from flint.types._gr import gr_fmpz_ctx as Z
        >>> Z(4).lcm(Z(6))
        12
        """
        cdef gr other_gr
        if not isinstance(other, gr):
            raise TypeError("gcd when other is not gr.")
        other_gr = other
        if not self.ctx == other_gr.ctx:
            raise TypeError("gcd of gr with different contexts.")
        return self.ctx.lcm(self, other_gr)

    def factor(self):
        """Return the factorization of the element.

        >>> from flint.types._gr import gr_fmpz_ctx as Z
        >>> Z(12).factor()
        (1, [(2, 2), (3, 1)])
        """
        return self.ctx.factor(self)

    def numer(self) -> gr:
        """Return the numerator of the element.

        >>> from flint.types._gr import gr_fmpq_ctx as Q
        >>> q = Q(2) / 3
        >>> q
        2/3
        >>> q.numer()
        2

        See also :meth:`denom`.
        """
        return self.ctx.numerator(self)

    def denom(self) -> gr:
        """Return the denominator of the element.

        >>> from flint.types._gr import gr_fmpq_ctx as Q
        >>> q = Q(2) / 3
        >>> q
        2/3
        >>> q.denom()
        3

        See also :meth:`numer`.
        """
        return self.ctx.denominator(self)

    def __floor__(self) -> gr:
        return self.ctx.floor(self)

    def __ceil__(self) -> gr:
        return self.ctx.ceil(self)

    def __trunc__(self) -> gr:
        return self.ctx.trunc(self)

    def __round__(self, ndigits: int = 0) -> gr:
        if ndigits != 0:
            raise NotImplementedError("Rounding to a specific number of digits is not supported")
        return self.ctx.nint(self)

    # def __int__(self) -> int:
    #     return self._floor().to_int()

    # def __float__(self) -> float:
    #     return ...

    def __abs__(self) -> gr:
        return self.ctx.abs(self)

    def conjugate(self) -> gr:
        """Return complex conjugate of the element.

        >>> from flint.types._gr import gr_fmpzi_ctx as ZI
        >>> I = ZI.gen()
        >>> (1 + I).conjugate()
        (1-I)
        """
        return self.ctx.conj(self)

    @property
    def real(self) -> gr:
        """Return the real part of the element.

        >>> from flint.types._gr import gr_fmpzi_ctx as ZI
        >>> I = ZI.gen()
        >>> (1 + I).real
        1
        """
        return self.ctx.re(self)

    @property
    def imag(self) -> gr:
        """Return the imaginary part of the element.

        >>> from flint.types._gr import gr_fmpzi_ctx as ZI
        >>> I = ZI.gen()
        >>> (1 + I).imag
        1
        """
        return self.ctx.im(self)

    # XXX: Return -1, 0, 1 as int?
    def sgn(self) -> gr:
        """Return the sign of the element.

        >>> from flint.types._gr import gr_fmpq_ctx as Q
        >>> Q(2).sgn()
        1
        >>> Q(-2).sgn()
        -1
        >>> Q(0).sgn()
        0
        """
        return self.ctx.sgn(self)

    def csgn(self) -> gr:
        """Return the complex sign of the element.

        >>> from flint.types._gr import gr_complex_acb_ctx
        >>> C = gr_complex_acb_ctx.new(10)
        >>> (1 + C.i()).csgn()  # doctest: +SKIP
        1
        """
        return self.ctx.csgn(self)

    def arg(self) -> gr:
        """Return the argument of the element.

        >>> from flint.types._gr import gr_complex_acb_ctx
        >>> C = gr_complex_acb_ctx.new(10)
        >>> (1 + C.i()).arg()
        [0.785 +/- 6.45e-4]
        """
        return self.ctx.arg(self)
