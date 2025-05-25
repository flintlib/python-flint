from cpython.list cimport PyList_GET_SIZE
from flint.flint_base.flint_base cimport flint_poly

from flint.types.fmpz cimport fmpz, any_as_fmpz
from flint.types.fmpz_poly cimport fmpz_poly
from flint.types.fmpz_mod_poly cimport fmpz_mod_poly
from flint.types.nmod_poly cimport nmod_poly
from flint.types.fq_default cimport fq_default_ctx, fq_default

from flint.pyflint cimport global_random_state

from flint.utils.typecheck cimport typecheck

from flint.utils.flint_exceptions import DomainError

cdef class fq_default_poly_ctx:
    r"""
    Context object for creating :class:`~.fq_default_poly` initialised
    with a finite field `GF(p^d)`.

        >>> fq_default_poly_ctx(163, 3, fq_type="FQ_NMOD")
        fq_default_poly_ctx(fq_default_ctx(163, 3, 'z', x^3 + 7*x + 161, 'FQ_NMOD'))
    """
    def __cinit__(self):
        pass

    def __dealloc__(self):
        pass

    def __init__(self, *args, **kwargs):
        # Allow context to be made from fq_default_ctx
        if len(args) == 1 and typecheck(args[0], fq_default_ctx):
            self.field = args[0]
        else:  # Otherwise attempt to create context from moduli
            self.field = fq_default_ctx(*args, **kwargs)

    cdef set_list_as_fq_default_poly(self, fq_default_poly_t poly, val):
        cdef long i, n

        # Get size of list and create a polynomial of length n
        n = PyList_GET_SIZE(val)
        fq_default_poly_fit_length(poly, n, self.field.val)

        # Iterate through the list and attempt to coerce every element
        # to an fq_type using the context
        for i in range(n):
            # First coerce the list element as an fq_default type
            v = self.field.any_as_fq_default(val[i])
            if v is NotImplemented:
                raise TypeError(f"unsupported coefficient in list: val[i] = {val[i]}, type(val[i] = {type(val[i])}")

            # Set the coefficient of the polynomial
            fq_default_poly_set_coeff(
                poly, i, (<fq_default>v).val, self.field.val
            )

    cdef set_any_as_fq_default_poly(self, fq_default_poly_t poly, obj):
        # First try and coerce polynomials to polynomials
        if typecheck(obj, fq_default_poly):
            if self != (<fq_default_poly>obj).ctx:
                raise ValueError("fields must match")
            fq_default_poly_set(poly, (<fq_default_poly>obj).val, self.field.val)
            return

        if typecheck(obj, fmpz_poly):
            fq_default_poly_set_fmpz_poly(poly, (<fmpz_poly>obj).val, self.field.val)
            return

        if typecheck(obj, fmpz_mod_poly) and self.prime() == (<fmpz_mod_poly>obj).ctx.modulus():
            fq_default_poly_set_fmpz_mod_poly(poly, (<fmpz_mod_poly>obj).val, self.field.val)
            return

        if typecheck(obj, nmod_poly) and self.prime() == (<nmod_poly>obj).modulus():
            fq_default_poly_set_nmod_poly(poly, (<nmod_poly>obj).val, self.field.val)
            return

        # Otherwise attempt to coerce the object to fq_default and
        # set a constant polynomial
        v = self.field.any_as_fq_default(obj)
        if v is NotImplemented:
            return NotImplemented

        fq_default_poly_set_fq_default(poly, (<fq_default>v).val, self.field.val)

    cdef any_as_fq_default_poly(self, obj):
        # If fq_default_poly, return if fields match
        if typecheck(obj, fq_default_poly):
            if self != (<fq_default_poly>obj).ctx:
                raise ValueError("fields must match")
            return obj

        # For all other types attempt to coerce to a fq_default_poly
        cdef fq_default_poly res
        res = self.new_ctype_poly()
        check = self.set_any_as_fq_default_poly(res.val, obj)
        if check is NotImplemented:
            return NotImplemented
        return res

    def base_field(self):
        """
        Return the base field of the polynomial ring

            >>> R = fq_default_poly_ctx(65537, 3)
            >>> R.base_field()
            fq_default_ctx(65537, 3, 'z', x^3 + 3*x^2 + 30077, 'FQ_NMOD')

        """
        return self.field

    def characteristic(self):
        """
        Return the characteristic of the field from the context
        as an ``fmpz`` type

            >>> R = fq_default_poly_ctx(65537, 3)
            >>> R.characteristic()
            65537

        """
        return self.field.characteristic()

    prime = characteristic

    def zero(self):
        """
        Return the zero element of this polynomial ring

            >>> R = fq_default_poly_ctx(163)
            >>> R.zero()
            0
        """
        cdef fq_default_poly res
        res = self.new_ctype_poly()
        fq_default_poly_zero(res.val, self.field.val)
        return res

    def one(self):
        """
        Return the one element of this polynomial ring

            >>> R = fq_default_poly_ctx(163)
            >>> R.one()
            1
        """
        cdef fq_default_poly res
        res = self.new_ctype_poly()
        fq_default_poly_one(res.val, self.field.val)
        return res

    def gen(self):
        """
        Return the generator of the polynomial: `x`

            >>> R = fq_default_poly_ctx(163)
            >>> R.gen()
            x
        """
        cdef fq_default_poly res
        res = self.new_ctype_poly()
        fq_default_poly_gen(res.val, self.field.val)
        return res

    def random_element(self, degree=3, not_zero=False, monic=False, irreducible=False):
        """
        Return a random element of degree ``degree``. If ``not_zero``
        is ``True``, ensures the output is not zero, if ``monic`` is
        ``True``, ensures the output is monic. If ``irreducible`` is
        ``True``, ensures that the output is irreducible.

            >>> R = fq_default_poly_ctx(163, 3)
            >>> f = R.random_element()
            >>> f.degree() <= 3
            True
            >>> f = R.random_element(degree=123)
            >>> f.degree() <= 123
            True
            >>> f = R.random_element(monic=True)
            >>> f.is_monic()
            True
            >>> f = R.random_element(degree=13, monic=True, irreducible=True)
            >>> f.degree() <= 13
            True
            >>> f.is_monic()
            True
            >>> f.is_irreducible()
            True
        """
        cdef slong length
        if not (isinstance(monic, bool) and isinstance(irreducible, bool) and isinstance(not_zero, bool)):
            raise TypeError("All of 'not_zero', 'monic' and 'irreducible' must be of type bool")

        length = degree + 1
        if length <= 0:
            raise ValueError("The degree argument must be non-negative")

        cdef fq_default_poly res
        res = self.new_ctype_poly()
        # all irreducible elements are returned as monic polynomials
        if irreducible:
            fq_default_poly_randtest_irreducible(
                res.val, global_random_state, length, self.field.val
            )
        elif monic:
            fq_default_poly_randtest_monic(
                res.val, global_random_state, length, self.field.val
            )
        elif not_zero:
            fq_default_poly_randtest_not_zero(
                res.val, global_random_state, length, self.field.val
            )
        else:
            fq_default_poly_randtest(
                res.val, global_random_state, length, self.field.val
            )
        return res

    cdef new_ctype_poly(self):
        return fq_default_poly.__new__(fq_default_poly, None, self)

    def __eq__(self, other):
        # Most often, we expect both `fq_default_poly` to be pointing
        # to the same ctx, so this seems the fastest way to check
        if self is other:
            return True

        # If they're not the same object in memory, they may have the
        # same field, which is good enough
        if typecheck(other, fq_default_poly_ctx):
            return self.field == (<fq_default_poly_ctx>other).field
        return False

    def __hash__(self):
        return hash(("polynomial_ring", self.field))

    def __str__(self):
        return f"Context for fq_default_poly with field: {self.field}"

    def __repr__(self):
        return f"fq_default_poly_ctx({repr(self.field)})"

    def __call__(self, val):
        return fq_default_poly(val, self)


cdef class fq_default_poly(flint_poly):
    """
    The *fq_default_poly* type represents univariate polynomials
    over a finite field.

    An *fq_default_poly* element is constructed from an :class:`~.fq_default_poly_ctx`
    either by passing it as an argument to the type, or
    by directly calling the context:

        >>> fq_default_poly([1,-2,3], fq_default_poly_ctx(2**127 - 1))
        3*x^2 + 170141183460469231731687303715884105725*x + 1
        >>> R = fq_default_poly_ctx(2**127 - 1)
        >>> R([4,5,6])
        6*x^2 + 5*x + 4

    """
    def __cinit__(self, val, ctx):
        if not typecheck(ctx, fq_default_poly_ctx):
            raise TypeError
        self.ctx = ctx
        fq_default_poly_init(self.val, self.ctx.field.val)

    def __dealloc__(self):
        if self.ctx is not None:
            fq_default_poly_clear(self.val, self.ctx.field.val)

    def __init__(self, val, ctx):
        if typecheck(val, list):
            self.ctx.set_list_as_fq_default_poly(self.val, val)
            return

        check = self.ctx.set_any_as_fq_default_poly(self.val, val)
        if check is NotImplemented:
            raise TypeError

    def __getitem__(self, long i):
        cdef fq_default x
        x = self.ctx.field.new_ctype_fq_default()
        if i < 0:
            return x
        fq_default_poly_get_coeff(
            x.val, self.val, i, self.ctx.field.val
        )
        return x

    def __setitem__(self, long i, x):
        if i < 0:
            raise ValueError("cannot assign to index < 0 of polynomial")
        v = self.ctx.field.any_as_fq_default(x)
        if v is NotImplemented:
            raise TypeError
        fq_default_poly_set_coeff(
            self.val, i, (<fq_default>v).val, self.ctx.field.val
        )

    def __richcmp__(self, other, int op):
        cdef bint res
        if op != 2 and op != 3:
            raise TypeError("fq_default_poly cannot be ordered")

        if not typecheck(other, fq_default_poly):
            other = self.ctx.any_as_fq_default_poly(other)

        if typecheck(other, fq_default_poly):
            res = (self.ctx == (<fq_default_poly>other).ctx) and \
                  fq_default_poly_equal(self.val, (<fq_default_poly>other).val, self.ctx.field.val)
            if op == 2:
                return res
            else:
                return not res
        else:
            return NotImplemented

    def __len__(self):
        return self.length()

    def __hash__(self):
        return hash(tuple(self.coeffs()))

    cpdef long length(self):
        """
        Return the length of the polynomial

            >>> R = fq_default_poly_ctx(163)
            >>> f = R([1,2,3])
            >>> f.length()
            3

        """
        return fq_default_poly_length(self.val, self.ctx.field.val)

    cpdef long degree(self):
        """
        Return the degree of the polynomial

            >>> R = fq_default_poly_ctx(163)
            >>> f = R([1,2,3])
            >>> f.degree()
            2

        """
        return fq_default_poly_degree(self.val, self.ctx.field.val)

    def constant_coefficient(self):
        """
        Return the constant coefficient of this polynomial.

            >>> R = fq_default_poly_ctx(163, 3)
            >>> f = R([1,2,3])
            >>> f.constant_coefficient()
            1
        """
        return self[0]

    def leading_coefficient(self):
        """
        Return the leading coefficient of this polynomial.

            >>> R = fq_default_poly_ctx(163)
            >>> f = R([1,2,3])
            >>> f.leading_coefficient()
            3
        """
        d = self.degree()
        if d < 0:
            return self.ctx.field.zero()
        return self[self.degree()]

    def reverse(self, degree=None):
        """
        Return a polynomial with the coefficients of this polynomial
        reversed.

        If ``degree`` is not None, the output polynomial will be zero-padded
        or truncated before being reversed. NOTE: degree must be non-negative.

            >>> R = fq_default_poly_ctx(163, 3)
            >>> f = R([1,2,3,4,5])
            >>> f.reverse()
            x^4 + 2*x^3 + 3*x^2 + 4*x + 5
            >>> f.reverse(degree=1)
            x + 2
            >>> f.reverse(degree=100)
            x^100 + 2*x^99 + 3*x^98 + 4*x^97 + 5*x^96
        """
        cdef fq_default_poly res
        cdef slong d

        if degree is not None:
            d = degree
            if d != degree or d < 0:
                raise ValueError(f"degree argument must be a non-negative integer, got {degree}")
        else:
            d = fq_default_poly_degree(self.val, self.ctx.field.val)

        length = d + 1

        res = self.ctx.new_ctype_poly()
        fq_default_poly_reverse(res.val, self.val, length, self.ctx.field.val)
        return res

    def truncate(self, slong n):
        r"""
        Notionally truncate the polynomial to have length ``n``. If
        ``n`` is larger than the length of the input, then ``self`` is
        returned. If ``n`` is not positive, then the zero polynomial
        is returned.

        Effectively returns this polynomial :math:`\mod x^n`.

            >>> R = fq_default_poly_ctx(163, 5)
            >>> f = R([1,2,3])
            >>> f.truncate(3) == f
            True
            >>> f.truncate(2)
            2*x + 1
            >>> f.truncate(1)
            1
            >>> f.truncate(0)
            0
            >>> f.truncate(-1)
            0
        """
        cdef fq_default_poly res

        length = fq_default_poly_degree(self.val, self.ctx.field.val)
        res = self.ctx.new_ctype_poly()

        if n <= 0:  # return zero
            return res
        elif n > length:  # do nothing
            fq_default_poly_set(
                res.val, self.val, self.ctx.field.val
            )
        else:
            fq_default_poly_set_trunc(
                res.val, self.val, n, self.ctx.field.val
            )
        return res

    def monic(self):
        """
        Return this polynomial divided by its leading coefficient.

            >>> R = fq_default_poly_ctx(163)
            >>> f = R([1,2,3])
            >>> f.monic()
            x^2 + 55*x + 109
        """
        cdef fq_default_poly res
        res = self.ctx.new_ctype_poly()
        fq_default_poly_make_monic(
            res.val, self.val, self.ctx.field.val
        )
        return res

    def is_zero(self):
        """
        Return ``True`` if the polynomial is the zero polynomial
        and ``False`` otherwise

            >>> R = fq_default_poly_ctx(163, 3)
            >>> f = R.zero()
            >>> f.is_zero()
            True
        """
        return 0 != fq_default_poly_is_zero(self.val, self.ctx.field.val)

    def is_one(self):
        """
        Return ``True`` if the polynomial is equal to one
        and ``False`` otherwise

            >>> R = fq_default_poly_ctx(163, 3)
            >>> f = R.one()
            >>> f.is_one()
            True
        """
        return 0 != fq_default_poly_is_one(self.val, self.ctx.field.val)

    def is_gen(self):
        """
        Return ``True`` if the polynomial is the generator
        of the polynomial, `x`, and ``False`` otherwise

            >>> R = fq_default_poly_ctx(163, 3)
            >>> f = R([0,1])
            >>> f.is_gen()
            True
        """
        return 0 != fq_default_poly_is_gen(self.val, self.ctx.field.val)

    def is_constant(self):
        """
        Return ``True`` if this is a constant polynomial.

            >>> R = fq_default_poly_ctx(163, 3)
            >>> x = R.gen()
            >>> x.is_constant()
            False
            >>> R(123).is_constant()
            True
        """
        return self.degree() <= 0

    def is_monic(self):
        """
        Return whether this polynomial is monic.

            >>> R = fq_default_poly_ctx(163)
            >>> x = R.gen()
            >>> f = x**2 + 5*x + 3
            >>> f.is_monic()
            True
            >>> f = 5*x**2 + x + 3
            >>> f.is_monic()
            False
        """
        return self.leading_coefficient().is_one()

    def is_irreducible(self):
        """
        Return whether this polynomial is irreducible.

            >>> R = fq_default_poly_ctx(163)
            >>> x = R.gen()
            >>> f = x**2 + 5*x + 3
            >>> f.is_irreducible()
            True
            >>> f = x**2 + x + 3
            >>> f.is_irreducible()
            False
        """
        return 1 == fq_default_poly_is_irreducible(self.val, self.ctx.field.val)

    def is_squarefree(self):
        """
        Return whether this polynomial is squarefree.

            >>> R = fq_default_poly_ctx(163)
            >>> x = R.gen()
            >>> f = (x + 1)**2 * (x + 3)
            >>> f.is_squarefree()
            False
            >>> f = (x + 1) * (x + 3)
            >>> f.is_squarefree()
            True

        """
        return 1 == fq_default_poly_is_squarefree(self.val, self.ctx.field.val)

    # ====================================
    # Native Arithmetic
    # ====================================

    def __pos__(self):
        return self

    def __neg__(self):
        cdef fq_default_poly res
        res = self.ctx.new_ctype_poly()
        fq_default_poly_neg(
            res.val, self.val, self.ctx.field.val
        )
        return res

    def __add__(self, other):
        cdef fq_default_poly res

        other = self.ctx.any_as_fq_default_poly(other)
        if other is NotImplemented:
            return NotImplemented

        res = self.ctx.new_ctype_poly()
        fq_default_poly_add(
            res.val, self.val, (<fq_default_poly>other).val, self.ctx.field.val
        )
        return res

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        cdef fq_default_poly res

        other = self.ctx.any_as_fq_default_poly(other)
        if other is NotImplemented:
            return NotImplemented

        res = self.ctx.new_ctype_poly()
        fq_default_poly_sub(
            res.val, self.val, (<fq_default_poly>other).val, self.ctx.field.val
        )
        return res

    def __rsub__(self, other):
        cdef fq_default_poly res

        other = self.ctx.any_as_fq_default_poly(other)
        if other is NotImplemented:
            return NotImplemented

        res = self.ctx.new_ctype_poly()
        fq_default_poly_sub(
            res.val, (<fq_default_poly>other).val, self.val, self.ctx.field.val
        )
        return res

    def __mul__(self, other):
        cdef fq_default_poly res
        res = self.ctx.new_ctype_poly()

        # First try scalar multiplication
        if not typecheck(other, fq_default_poly):
            scalar = self.ctx.field.any_as_fq_default(other)
            if scalar is not NotImplemented:
                fq_default_poly_scalar_mul_fq_default(
                    res.val, self.val, (<fq_default>scalar).val, self.ctx.field.val
                )
                return res

        # Otherwise treat other as an fq_default_poly
        other = self.ctx.any_as_fq_default_poly(other)
        if other is NotImplemented:
            return NotImplemented
        fq_default_poly_mul(
            res.val, self.val, (<fq_default_poly>other).val, self.ctx.field.val
        )
        return res

    def __rmul__(self, other):
        return self.__mul__(other)

    def __pow__(self, e, mod=None):
        if mod is not None:
            return self.pow_mod(e, mod)

        cdef fq_default_poly res
        if e < 0:
            self = 1 / self
            e = -e

        if e == 2:
            return self.square()

        res = self.ctx.new_ctype_poly()
        fq_default_poly_pow(
            res.val, self.val, <ulong>e, self.ctx.field.val
        )
        return res

    def pow_mod(self, e, modulus):
        r"""
        Returns ``self`` raised to the power ``e`` modulo ``modulus``:
        :math:`f^e \mod g`/

            >>> R = fq_default_poly_ctx(163)
            >>> x = R.gen()
            >>> f = 30*x**6 + 104*x**5 + 76*x**4 + 33*x**3 + 70*x**2 + 44*x + 65
            >>> mod = x**4 + 93*x**3 + 78*x**2 + 72*x + 149
            >>> f.pow_mod(123, mod)
            3*x^3 + 25*x^2 + 115*x + 161
            >>> f.pow_mod(2**64, mod)
            52*x^3 + 96*x^2 + 136*x + 9
        """
        cdef fq_default_poly res

        if e < 0:
            raise ValueError("Exponent must be non-negative")

        modulus = self.ctx.any_as_fq_default_poly(modulus)
        if modulus is NotImplemented:
            raise TypeError(f"Cannot interpret {modulus} as a polynomial")

        # Output polynomial
        res = self.ctx.new_ctype_poly()

        # For small exponents, use a simple binary exponentiation method
        if e.bit_length() < 32:
            fq_default_poly_powmod_ui_binexp(
                res.val, self.val, <ulong>e, (<fq_default_poly>modulus).val, res.ctx.field.val
            )
            return res

        # For larger exponents we need to cast e to an fmpz first
        e_fmpz = any_as_fmpz(e)
        fq_default_poly_powmod_fmpz_binexp(
                res.val, self.val, (<fmpz>e_fmpz).val, (<fq_default_poly>modulus).val, res.ctx.field.val
            )
        return res

    def divmod(self, other):
        """
        Return `Q`, `R` such that for ``self`` = `F` and ``other`` = `G`,
        `F = Q*G + R`

            >>> R = fq_default_poly_ctx(163)
            >>> f = R([123, 129, 63, 14, 51, 76, 133])
            >>> g = R([106, 134, 32, 41, 158, 115, 115])
            >>> f.divmod(g)
            (21, 106*x^5 + 156*x^4 + 131*x^3 + 43*x^2 + 86*x + 16)
        """
        cdef fq_default_poly Q, R

        other = self.ctx.any_as_fq_default_poly(other)
        if other is NotImplemented:
            raise TypeError(f"Cannot interpret {other} as a polynomial")

        if other.is_zero():
            raise ZeroDivisionError(f"Cannot compute divmod as {other =}")

        Q = self.ctx.new_ctype_poly()
        R = self.ctx.new_ctype_poly()
        fq_default_poly_divrem(
            Q.val, R.val, self.val, (<fq_default_poly>other).val, self.ctx.field.val
        )
        return Q, R

    def __divmod__(self, other):
        return self.divmod(other)

    def __rdivmod__(self, other):
        other = self.ctx.any_as_fq_default_poly(other)
        if other is NotImplemented:
            return other
        return other.divmod(self)

    def exact_division(self, other):
        """
        Attempt to compute the exact quotient of self with other.

        Raises a value error if division without remainder is not
        possible.

            >>> R = fq_default_poly_ctx(163)
            >>> f = R([1,2,1])
            >>> g = R([1,1])
            >>> f.exact_division(g)
            x + 1
        """
        cdef bint check
        cdef fq_default_poly res

        # Case when right is not fq_default_poly, try to convert to fmpz
        other = self.ctx.any_as_fq_default_poly(other)
        if other is NotImplemented:
            raise TypeError(f"Cannot convert {other} to 'fq_default_poly' type.")

        if other.is_zero():
            raise ZeroDivisionError("Cannot divide by zero")

        res = self.ctx.new_ctype_poly()
        check = fq_default_poly_divides(
            res.val, self.val, (<fq_default_poly>other).val, res.ctx.field.val
        )
        if check == 0:
            raise DomainError(
                f"{other} does not divide {self}"
            )

        return res

    def __truediv__(self, other):
        cdef fq_default_poly res
        res = self.ctx.new_ctype_poly()

        # First try scalar division
        if not typecheck(other, fq_default_poly):
            scalar = self.ctx.field.any_as_fq_default(other)
            if scalar is not NotImplemented:
                if scalar.is_zero():
                    raise ZeroDivisionError("Cannot divide by zero")
                fq_default_poly_scalar_div_fq_default(
                    res.val, self.val, (<fq_default>scalar).val, self.ctx.field.val
                )
                return res

        # Otherwise treat other as an fq_default_poly
        return self.exact_division(other)

    def __rtruediv__(self, other):
        if self.is_zero():
            raise ZeroDivisionError("Cannot divide by zero")
        other = self.ctx.any_as_fq_default_poly(other)
        if other is NotImplemented:
            return other
        return other.exact_division(self)

    def __floordiv__(self, other):
        cdef fq_default_poly res
        res = self.ctx.new_ctype_poly()

        # First try scalar division
        if not typecheck(other, fq_default_poly):
            scalar = self.ctx.field.any_as_fq_default(other)
            if scalar is not NotImplemented:
                if scalar.is_zero():
                    raise ZeroDivisionError("Cannot divide by zero")
                fq_default_poly_scalar_div_fq_default(
                    res.val, self.val, (<fq_default>scalar).val, self.ctx.field.val
                )
                return res

        # Coerce right element to fq_default_poly
        other = self.ctx.any_as_fq_default_poly(other)
        if other is NotImplemented:
            return other

        # Do not divide by zero
        if other.is_zero():
            raise ZeroDivisionError("Cannot divide by zero")

        # floor division uses divmod but ignores the rem term
        cdef fq_default_poly rem
        rem = self.ctx.new_ctype_poly()

        fq_default_poly_divrem(
            res.val, rem.val, self.val, (<fq_default_poly>other).val, self.ctx.field.val
        )
        return res

    def __rfloordiv__(self, other):
        cdef fq_default_poly res
        cdef fq_default_poly rem

        # Do not divide by zero
        if self.is_zero():
            raise ZeroDivisionError("Cannot divide by zero")

        # Coerce right element to fq_default_poly
        other = self.ctx.any_as_fq_default_poly(other)
        if other is NotImplemented:
            return other

        # floor division uses divmod but ignores the rem term
        res = self.ctx.new_ctype_poly()
        rem = self.ctx.new_ctype_poly()
        fq_default_poly_divrem(
            res.val, rem.val, (<fq_default_poly>other).val, self.val, self.ctx.field.val
        )
        return res

    def __mod__(self, other):
        cdef fq_default_poly res

        other = self.ctx.any_as_fq_default_poly(other)
        if other is NotImplemented:
            return NotImplemented

        if other.is_zero():
            raise ZeroDivisionError("Cannot compute remainder modulo 0")

        res = self.ctx.new_ctype_poly()
        fq_default_poly_rem(
            res.val, self.val, (<fq_default_poly>other).val, self.ctx.field.val
        )
        return res

    def __rmod__(self, other):
        cdef fq_default_poly res

        if self.is_zero():
            raise ZeroDivisionError("Cannot compute remainder modulo 0")

        other = self.ctx.any_as_fq_default_poly(other)
        if other is NotImplemented:
            return NotImplemented

        res = self.ctx.new_ctype_poly()
        fq_default_poly_rem(
            res.val, (<fq_default_poly>other).val, self.val, self.ctx.field.val
        )
        return res

    # ====================================
    # Additional Arithmetic
    # ====================================

    def square(self):
        """
        Returns the square of ``self``

            >>> R = fq_default_poly_ctx(163)
            >>> f = R([1,2,3])
            >>> f.square()
            9*x^4 + 12*x^3 + 10*x^2 + 4*x + 1
        """
        cdef fq_default_poly res

        res = self.ctx.new_ctype_poly()
        fq_default_poly_sqr(
            res.val, self.val, self.ctx.field.val
        )
        return res

    def left_shift(self, slong n):
        """
        Returns ``self`` shifted left by ``n`` coefficients by inserting
        zero coefficients. This is equivalent to multiplying the polynomial
        by x^n

            >>> R = fq_default_poly_ctx(163)
            >>> f = R([1,2,3])
            >>> f.left_shift(0)
            3*x^2 + 2*x + 1
            >>> f.left_shift(1)
            3*x^3 + 2*x^2 + x
            >>> f.left_shift(4)
            3*x^6 + 2*x^5 + x^4

        """
        cdef fq_default_poly res
        res = self.ctx.new_ctype_poly()

        if n < 0:
            raise ValueError("Value must be shifted by a non-negative integer")

        if n > 0:
            fq_default_poly_shift_left(
                res.val, self.val, n, self.ctx.field.val
            )
        else:  # do nothing, just copy self
            fq_default_poly_set(
                res.val, self.val, self.ctx.field.val
            )

        return res

    def right_shift(self, slong n):
        """
        Returns ``self`` shifted right by ``n`` coefficients.
        This is equivalent to the floor division of the polynomial
        by x^n

            >>> R = fq_default_poly_ctx(163)
            >>> f = R([1,2,3])
            >>> f.right_shift(0)
            3*x^2 + 2*x + 1
            >>> f.right_shift(1)
            3*x + 2
            >>> f.right_shift(4)
            0
        """
        cdef fq_default_poly res

        if n < 0:
            raise ValueError("Value must be shifted by a non-negative integer")

        res = self.ctx.new_ctype_poly()

        if n > 0:
            fq_default_poly_shift_right(
                res.val, self.val, n, self.ctx.field.val
            )
        else:  # do nothing, just copy self
            fq_default_poly_set(
                res.val, self.val, self.ctx.field.val
            )

        return res

    def sqrt(self):
        """
        If ``self`` is a perfect square, compute the square root

            >>> R = fq_default_poly_ctx(163)
            >>> f = R([1,1])
            >>> (f*f).sqrt()
            x + 1

        """
        cdef fq_default_poly res
        cdef int check

        res = self.ctx.new_ctype_poly()
        check = fq_default_poly_sqrt(
            res.val, self.val, res.ctx.field.val
        )
        if check != 1:
            raise DomainError(
                f"Cannot compute square-root {self}"
            )
        return res

    def mul_mod(self, other, modulus):
        """
        Computes the multiplication of ``self`` with ``other``
        modulo the polynomial ``modulus``

            >>> R = fq_default_poly_ctx(163)
            >>> x = R.gen()
            >>> f = 30*x**6 + 104*x**5 + 76*x**4 + 33*x**3 + 70*x**2 + 44*x + 65
            >>> g = 43*x**6 + 91*x**5 + 77*x**4 + 113*x**3 + 71*x**2 + 132*x + 60
            >>> mod = x**4 + 93*x**3 + 78*x**2 + 72*x + 149
            >>>
            >>> f.mul_mod(g, mod)
            106*x^3 + 44*x^2 + 53*x + 77
        """
        cdef fq_default_poly res

        other = self.ctx.any_as_fq_default_poly(other)
        if other is NotImplemented:
            raise TypeError(f"Cannot interpret {other} as a polynomial")

        modulus = self.ctx.any_as_fq_default_poly(modulus)
        if modulus is NotImplemented:
            raise TypeError(f"Cannot interpret {modulus} as a polynomial")

        res = self.ctx.new_ctype_poly()

        fq_default_poly_mulmod(
            res.val, self.val, (<fq_default_poly>other).val, (<fq_default_poly>modulus).val, res.ctx.field.val
        )
        return res

    # ====================================
    # Truncated Arithmetic
    # ====================================

    def equal_trunc(self, other, slong n):
        """
        Returns if two polynomials are equal when truncated to the first ``n`` terms

            >>> R = fq_default_poly_ctx(163)
            >>> f = R([1,2,3,4,5])
            >>> h = R([1,2,3])
            >>> f.equal_trunc(h, 3)
            True
            >>> f.equal_trunc(h, 4)
            False
        """
        # Only allow comparison with other fmpz_mod_poly
        if not typecheck(other, fq_default_poly):
            return False

        # Ensure the contexts match
        other_c = <fq_default_poly>other
        if self.ctx != other_c.ctx:
            return False

        return 1 == fq_default_poly_equal_trunc(self.val, other_c.val, n, self.ctx.field.val)

    def add_trunc(self, other, slong n):
        """
        Truncate ``self`` and ``other`` to polynomials of length ``n`` and return their sum

            >>> R = fq_default_poly_ctx(163)
            >>> f = R([1,2,3,4,5])
            >>> h = R([1,2,3])
            >>> f.add_trunc(h, 2)
            4*x + 2
            >>> f.add_trunc(h, 3)
            6*x^2 + 4*x + 2
        """
        # Only allow addition with other fq_default_poly
        if not typecheck(other, fq_default_poly):
            raise TypeError("other polynomial must be of type fq_default_poly")

        # Ensure the contexts match
        other_c = <fq_default_poly>other
        if self.ctx != other_c.ctx:
            raise ValueError("other polynomial's context does not match")

        cdef fq_default_poly res
        res = self.ctx.new_ctype_poly()
        fq_default_poly_add_series(
            res.val, self.val, other_c.val, n, res.ctx.field.val
        )
        return res

    def sub_trunc(self, other, slong n):
        """
        Truncate ``self`` and ``other`` to polynomials of length ``n`` and return their difference

            >>> R = fq_default_poly_ctx(163)
            >>> f = R([2,3,5,7,11])
            >>> h = R([1,2,4,8,16])
            >>> f.sub_trunc(h, 2)
            x + 1
            >>> f.sub_trunc(h, 4)
            162*x^3 + x^2 + x + 1
        """
        # Only allow subtraction with other fq_default_poly
        if not typecheck(other, fq_default_poly):
            raise TypeError("other polynomial must be of type fq_default_poly")

        # Ensure the contexts match
        other_c = <fq_default_poly>other
        if self.ctx != other_c.ctx:
            raise ValueError("other polynomial's context does not match")

        cdef fq_default_poly res
        res = self.ctx.new_ctype_poly()
        fq_default_poly_sub_series(
            res.val, self.val, other_c.val, n, res.ctx.field.val
        )
        return res

    def mul_low(self, other, slong n):
        r"""
        Returns the lowest ``n`` coefficients of the multiplication of ``self`` with ``other``

        Equivalent to computing `f(x) \cdot g(x) \mod x^n`

            >>> R = fq_default_poly_ctx(163)
            >>> f = R([2,3,5,7,11])
            >>> g = R([1,2,4,8,16])
            >>> f.mul_low(g, 5)
            101*x^4 + 45*x^3 + 19*x^2 + 7*x + 2
            >>> f.mul_low(g, 3)
            19*x^2 + 7*x + 2
            >>> f.mul_low(g, 1)
            2
        """
        # Only allow multiplication with other fq_default_poly
        if not typecheck(other, fq_default_poly):
            raise TypeError("other polynomial must be of type fq_default_poly")

        # Ensure the contexts match
        other_c = <fq_default_poly>other
        if self.ctx != other_c.ctx:
            raise ValueError("other polynomial's context does not match")

        cdef fq_default_poly res
        res = self.ctx.new_ctype_poly()
        fq_default_poly_mullow(
            res.val, self.val, other_c.val, n, res.ctx.field.val
        )
        return res

    def pow_trunc(self, slong e, slong n):
        r"""
        Returns ``self`` raised to the power ``e`` modulo `x^n`:
        :math:`f^e \mod x^n`/

        Note: For exponents larger that 2^31 (which do not fit inside a ulong) use the
        method :meth:`~.pow_mod` with the explicit modulus `x^n`.

            >>> R = fq_default_poly_ctx(163)
            >>> x = R.gen()
            >>> f = 30*x**6 + 104*x**5 + 76*x**4 + 33*x**3 + 70*x**2 + 44*x + 65
            >>> f.pow_trunc(2**20, 30) == pow(f, 2**20, x**30)
            True
            >>> f.pow_trunc(2**20, 5)
            132*x^4 + 113*x^3 + 36*x^2 + 48*x + 6
        """
        if e < 0:
            raise ValueError("Exponent must be non-negative")

        cdef fq_default_poly res
        res = self.ctx.new_ctype_poly()
        fq_default_poly_pow_trunc(res.val, self.val, e, n, res.ctx.field.val)
        return res

    def sqrt_trunc(self, slong n):
        """
        Returns the square root of ``self`` modulo `x^n`.

        Requires that the constant coefficient of the polynomial is one.

            >>> R = fq_default_poly_ctx(163, 3)
            >>> x = R.gen()
            >>> z = R.base_field().gen()
            >>> f = (37*z + 54)*x**3 + (8*z + 94)*x**2 + (52*z + 142)*x + 1
            >>> h = f.sqrt_trunc(4)
            >>> h
            (7*z^2 + 17*z + 148)*x^3 + (151*z^2 + 114*z + 53)*x^2 + (26*z + 71)*x + 1
            >>> h.mul_low(h, 4) == f
            True
        """
        cdef fq_default_poly res

        # FLINT requires the constant term is one, so we need to normalise
        c = self.constant_coefficient()
        if c.is_zero():
            raise ZeroDivisionError("constant coefficient must be invertible (and a square in the field)")
        if not c.is_square():
            raise ValueError("constant coefficient of the polynomial must be a square")
        self = self / c

        res = self.ctx.new_ctype_poly()
        fq_default_poly_sqrt_series(
            res.val, self.val, n, res.ctx.field.val
        )
        return res * c.sqrt()

    def inv_sqrt_trunc(self, slong n):
        """
        Returns the inverse of the square root of ``self`` modulo `x^n`.

        Requires that the constant coefficient of the polynomial is one.

            >>> R = fq_default_poly_ctx(65537, 2)
            >>> x = R.gen()
            >>> z = R.base_field().gen()
            >>> f = 28902*x**3 + (49416*z + 58229)*x**2 + 9441*z*x + (7944*z + 57534)
            >>> h = f.inv_sqrt_trunc(3)
            >>> h
            (23030*z + 8965)*x^2 + (43656*z + 7173)*x + (27935*z + 28199)
            >>> (h*h).mul_low(f, 3).is_one()
            True
        """
        cdef fq_default_poly res

        # FLINT requires the constant term is one, so we need to normalise
        c = self.constant_coefficient()
        if c.is_zero():
            raise ZeroDivisionError("constant coefficient must be invertible (and a square in the field)")
        if not c.is_square():
            raise ValueError("constant coefficient of the polynomial must be a square")
        self = self / c

        res = self.ctx.new_ctype_poly()
        fq_default_poly_invsqrt_series(
            res.val, self.val, n, res.ctx.field.val
        )
        return res / c.sqrt()

    def inverse_series_trunc(self, slong n):
        """
        Returns the inverse of ``self`` modulo `x^n`.

            >>> R = fq_default_poly_ctx(163, 3)
            >>> x = R.gen()
            >>> z = R.base_field().gen()
            >>> f = (37*z + 54)*x**3 + (8*z + 94)*x**2 + (52*z + 142)*x + 1
            >>> f.inverse_series_trunc(2)
            (111*z + 21)*x + 1
            >>> f.inverse_series_trunc(3)
            (96*z^2 + 90*z + 21)*x^2 + (111*z + 21)*x + 1
            >>> f.inverse_series_trunc(4)
            (34*z^2 + z + 2)*x^3 + (96*z^2 + 90*z + 21)*x^2 + (111*z + 21)*x + 1
            >>> h = f.inverse_series_trunc(4)
            >>> f.mul_low(h, 4).is_one()
            True
        """
        cdef fq_default_poly res

        if self.constant_coefficient().is_zero():
            raise ZeroDivisionError("constant coefficient must be invertible")

        res = self.ctx.new_ctype_poly()
        fq_default_poly_inv_series(
            res.val, self.val, n, res.ctx.field.val
        )
        return res

    # ====================================
    # GCD and Extended GCD
    # ====================================

    def gcd(self, other):
        """
        Return the greatest common divisor of self and other.

            >>> R = fq_default_poly_ctx(163)
            >>> x = R.gen()
            >>> f = x*(x + 1)
            >>> f.gcd(x+1)
            x + 1
            >>> f.gcd(x*x)
            x

        """
        cdef fq_default_poly res

        other = self.ctx.any_as_fq_default_poly(other)
        if other is NotImplemented:
            raise TypeError(f"Cannot interpret {other} as a polynomial")

        res = self.ctx.new_ctype_poly()
        fq_default_poly_gcd(
            res.val, self.val, (<fq_default_poly>other).val, self.ctx.field.val
        )
        return res

    def xgcd(self, other):
        r"""
        Computes the extended gcd of self and other: (`G`, `S`, `T`)
        where `G` is the ``gcd(self, other)`` and `S`, `T` are such that:

        :math:`G = \textrm{self}*S +  \textrm{other}*T`

            >>> R = fq_default_poly_ctx(163)
            >>> f = R([143, 19, 37, 138, 102, 127, 95])
            >>> g = R([139, 9, 35, 154, 87, 120, 24])
            >>> f.xgcd(g)
            (x^3 + 128*x^2 + 123*x + 91, 17*x^2 + 49*x + 104, 21*x^2 + 5*x + 25)
            >>> G, S, T = f.xgcd(g)
            >>> assert G == f*S + g*T

        """
        cdef fq_default_poly G, S, T

        other = self.ctx.any_as_fq_default_poly(other)
        if other is NotImplemented:
            raise TypeError(f"Cannot interpret {other} as a polynomial")

        G = self.ctx.new_ctype_poly()
        S = self.ctx.new_ctype_poly()
        T = self.ctx.new_ctype_poly()

        fq_default_poly_xgcd(
            G.val, S.val, T.val, self.val, (<fq_default_poly>other).val, self.ctx.field.val
        )

        return (G, S, T)

    def inverse_mod(self, other):
        """
        Returns the inverse of ``self`` modulo ``other``

            >>> R = fq_default_poly_ctx(163)
            >>> f = R([123, 129, 63, 14, 51, 76, 133])
            >>> h = R([139, 9, 35, 154, 87, 120, 24])
            >>> g = f.inverse_mod(h)
            >>> g
            41*x^5 + 121*x^4 + 47*x^3 + 41*x^2 + 6*x + 5
            >>> assert f.mul_mod(g, h).is_one()
        """
        G, S, _ = self.xgcd(other)
        if not G.is_one():
            raise ValueError(f"polynomial has no inverse modulo other = {other}")
        return S

    # ====================================
    # Derivative
    # ====================================

    def derivative(self):
        """
        The formal derivative of this polynomial

            >>> R = fq_default_poly_ctx(163)
            >>> x = R.gen()
            >>> f = 111*x**4 + 58*x**3 + 98*x**2 + 117*x + 7
            >>> f.derivative()
            118*x^3 + 11*x^2 + 33*x + 117

        """
        cdef fq_default_poly res
        res = self.ctx.new_ctype_poly()
        fq_default_poly_derivative(
            res.val, self.val, self.ctx.field.val
        )
        return res

    def radical(self):
        """
        Return the radical of ``self``, the product of the irreducible
        factors of the polynomial. This is also referred to as the
        square-free part of the polynomial.

            >>> R = fq_default_poly_ctx(163)
            >>> x = R.gen()
            >>> f = (x + 1)**3 * (x + 2)
            >>> f.radical()
            x^2 + 3*x + 2

        """
        return self.exact_division(self.gcd(self.derivative()))

    def integral(self):
        raise NotImplementedError("fq_default_integral is not available from FLINT")

    # ====================================
    # Evaluation and Composition
    # ====================================

    def evaluate(self, input):
        """
        Evaluate ``self`` at a point in the base ring. This is
        the same as calling the polynomial directly.

            >>> R = fq_default_poly_ctx(163)
            >>> f = R([1,2,3,4,5,6])
            >>> f.evaluate(-1)
            160
            >>> f.evaluate(-1) == f(-1)
            True
        """
        val = self.ctx.field.any_as_fq_default(input)
        if val is NotImplemented:
            raise TypeError(f"Cannot evaluate the polynomial with input: {input}")

        # When the fq_default type is FMPZ_MOD this function Segfaults for older
        # FLINT versions
        # See: https://github.com/flintlib/flint/issues/2046
        # TODO:
        # Hack for converting self to an fmpz_mod type, needs to be improved
        # and could be done when we decide exactly how to handle fq_default to
        # fq_type conversions more generally
        if self.ctx.field.fq_type == 5:
            from flint import fmpz_mod_poly_ctx
            ring = fmpz_mod_poly_ctx(self.ctx.characteristic())
            poly = ring([int(i) for i in self.coeffs()])
            return poly(int(input))

        cdef fq_default res
        res = self.ctx.field.new_ctype_fq_default()
        fq_default_poly_evaluate_fq_default(res.val, self.val, (<fq_default>val).val, self.ctx.field.val)
        return res

    def compose(self, other):
        """
        Returns the composition of two polynomials

        To be precise about the order of composition, given ``self``, and ``other``
        by `f(x)`, `g(x)`, returns `f(g(x))`.

            >>> R = fq_default_poly_ctx(163)
            >>> f = R([1,2,3])
            >>> g = R([0,0,1])
            >>> f.compose(g)
            3*x^4 + 2*x^2 + 1
            >>> g.compose(f)
            9*x^4 + 12*x^3 + 10*x^2 + 4*x + 1
        """
        cdef fq_default_poly res
        val = self.ctx.any_as_fq_default_poly(other)
        if val is NotImplemented:
            raise TypeError(f"Cannot compose the polynomial with input: {other}")

        res = self.ctx.new_ctype_poly()
        fq_default_poly_compose(res.val, self.val, (<fq_default_poly>val).val, self.ctx.field.val)
        return res

    def compose_mod(self, other, modulus):
        r"""
        Returns the composition of two polynomials modulo a third.

        To be precise about the order of composition, given ``self``, and ``other``
        and ``modulus`` by `f(x)`, `g(x)` and `h(x)`, returns `f(g(x)) \mod h(x)`.
        We require that `h(x)` is non-zero.

            >>> R = fq_default_poly_ctx(163)
            >>> f = R([1,2,3,4,5])
            >>> g = R([3,2,1])
            >>> h = R([1,0,1,0,1])
            >>> f.compose_mod(g, h)
            63*x^3 + 100*x^2 + 17*x + 63
            >>> g.compose_mod(f, h)
            147*x^3 + 159*x^2 + 4*x + 7
        """
        cdef fq_default_poly res
        val = self.ctx.any_as_fq_default_poly(other)
        if val is NotImplemented:
            raise TypeError(f"cannot compose the polynomial with input: {other}")

        h = self.ctx.any_as_fq_default_poly(modulus)
        if h is NotImplemented:
            raise TypeError(f"cannot reduce the polynomial with input: {modulus}")

        if h.is_zero():
            raise ZeroDivisionError("cannot reduce modulo zero")

        res = self.ctx.new_ctype_poly()
        fq_default_poly_compose_mod(res.val, self.val, (<fq_default_poly>val).val, (<fq_default_poly>h).val, self.ctx.field.val)
        return res

    def __call__(self, input):
        if typecheck(input, fq_default_poly):
            return self.compose(input)
        return self.evaluate(input)

    # ====================================
    # Factoring and Root Finding
    # ====================================

    def factor_squarefree(self):
        """
        Factors self into irreducible, squarefree factors, returning a tuple
        ``(c, factors)`` where `c` is the content of the coefficients and
        factors is a list of ``(poly, exp)`` pairs.

            >>> R = fq_default_poly_ctx(163)
            >>> x = R.gen()
            >>> f = (x + 1) * (x + 2)
            >>> f.factor_squarefree()
            (1, [(x^2 + 3*x + 2, 1)])
            >>> f = (x + 1) * (x + 2)**5
            >>> f.factor_squarefree()
            (1, [(x + 1, 1), (x + 2, 5)])
        """
        cdef fq_default_poly_factor_t fac
        cdef int i

        fq_default_poly_factor_init(fac, self.ctx.field.val)
        fq_default_poly_factor_squarefree(fac, self.val, self.ctx.field.val)

        num = fq_default_poly_factor_length(fac, self.ctx.field.val)
        res = [0] * num

        cdef fq_default_poly u
        for i in range(num):
            u = self.ctx.new_ctype_poly()
            fq_default_poly_factor_get_poly(u.val, fac, i, self.ctx.field.val)
            exp = fq_default_poly_factor_exp(fac, i, self.ctx.field.val)
            res[i] = (u, exp)
        return self.leading_coefficient(), res

    def factor(self):
        """
        Factors self into irreducible factors, returning a tuple
        ``(c, factors)`` where `c` is the content of the coefficients and
        factors is a list of ``(poly, exp)`` pairs.

            >>> R = fq_default_poly_ctx(163)
            >>> x = R.gen()
            >>> f = 6*x**4 + 7*x**3 + 7*x**2 + 8*x + 6
            >>> f.factor()
            (6, [(x^4 + 137*x^3 + 137*x^2 + 110*x + 1, 1)])
            >>> f = (x + 1)**3 * (x + 2)
            >>> f.factor()
            (1, [(x + 2, 1), (x + 1, 3)])
        """
        cdef fq_default_poly_factor_t fac
        cdef int i

        fq_default_poly_factor_init(fac, self.ctx.field.val)

        lead = self.leading_coefficient()
        fq_default_poly_factor(fac, (<fq_default>lead).val, self.val, self.ctx.field.val)

        num = fq_default_poly_factor_length(fac, self.ctx.field.val)
        res = [0] * num

        cdef fq_default_poly u
        for i in range(num):
            u = self.ctx.new_ctype_poly()
            fq_default_poly_factor_get_poly(u.val, fac, i, self.ctx.field.val)
            exp = fq_default_poly_factor_exp(fac, i, self.ctx.field.val)
            res[i] = (u, exp)
        return self.leading_coefficient(), res

    def roots(self, multiplicities=True):
        r"""
        Return the roots of the polynomial in the finite field

            >>> R = fq_default_poly_ctx(163)
            >>> x = R.gen()
            >>> f = (x - 1) * (x - 2)**3 * (x - 3)**5
            >>> f.roots()
            [(1, 1), (2, 3), (3, 5)]
            >>> f.roots(multiplicities=False)
            [1, 2, 3]
        """
        cdef fq_default_poly_factor_t fac
        cdef int i, with_multiplicity

        with_multiplicity = 0
        if multiplicities:
            with_multiplicity = 1

        fq_default_poly_factor_init(fac, self.ctx.field.val)
        fq_default_poly_roots(fac, self.val, with_multiplicity, self.ctx.field.val)

        num = fq_default_poly_factor_length(fac, self.ctx.field.val)
        res = [0] * num

        cdef fq_default_poly linear_factor
        cdef fq_default root

        for i in range(num):
            # Get a factor of the form (x - a)
            linear_factor = self.ctx.new_ctype_poly()
            fq_default_poly_factor_get_poly(linear_factor.val, fac, i, self.ctx.field.val)

            # Compute a from (x - a)
            root = self.ctx.field.new_ctype_fq_default()
            fq_default_poly_get_coeff(
                root.val, linear_factor.val, 0, root.ctx.val
            )
            fq_default_neg(
                root.val, root.val, root.ctx.val
            )
            if multiplicities:
                mul = fq_default_poly_factor_exp(fac, i, self.ctx.field.val)
                res[i] = (root, mul)
            else:
                res[i] = root
        return res

    def real_roots(self):
        """
        This method is not implemented for polynomials in finite fields
        """
        raise DomainError("Cannot compute real roots for polynomials over finite fields")

    def complex_roots(self):
        """
        This method is not implemented for polynomials in finite fields
        """
        raise DomainError("Cannot compute complex roots for polynomials over finite fields")

    # ====================================
    # Inflation and Deflation
    # ====================================

    def inflate(self, ulong n):
        r"""
        Returns the result of the polynomial `f = \textrm{self}` to
        `f(x^n)`

            >>> R = fq_default_poly_ctx(163, 3)
            >>> f = R([1,2,3])
            >>> f.inflate(10)
            3*x^20 + 2*x^10 + 1

        """
        cdef fq_default_poly res
        res = self.ctx.new_ctype_poly()
        fq_default_poly_inflate(
            res.val, self.val, n, res.ctx.field.val
        )
        return res

    def deflate(self, ulong n):
        r"""
        Returns the result of the polynomial `f = \textrm{self}` to
        `f(x^{1/n})`

            >>> R = fq_default_poly_ctx(163, 3)
            >>> f = R([1,0,2,0,3])
            >>> f
            3*x^4 + 2*x^2 + 1
            >>> f.deflate(2)
            3*x^2 + 2*x + 1

        """
        cdef fq_default_poly res

        n_max = fq_default_poly_deflation(
            self.val, self.ctx.field.val
        )
        if n > n_max:
            raise ValueError(f"Cannot deflate with n = {n}, maximum allowed value is n_max = {n_max}")

        res = self.ctx.new_ctype_poly()
        fq_default_poly_deflate(
            res.val, self.val, n, res.ctx.field.val
        )
        return res

    def deflation(self):
        r"""
        Returns the tuple (g, n) where for `f = \textrm{self}` to
        `g = f(x^{1/n})` where n is the largest allowed integer

            >>> R = fq_default_poly_ctx(163, 3)
            >>> f = R([1,0,2,0,3])
            >>> f
            3*x^4 + 2*x^2 + 1
            >>> f.deflate(2)
            3*x^2 + 2*x + 1

        """
        cdef fq_default_poly res
        if self.is_zero():
            return self, 1
        n = fq_default_poly_deflation(
            self.val, self.ctx.field.val
        )
        res = self.ctx.new_ctype_poly()
        fq_default_poly_deflate(
            res.val, self.val, n, res.ctx.field.val
        )
        return res, int(n)
