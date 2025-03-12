from cpython.list cimport PyList_GET_SIZE
from cpython.long cimport PyLong_Check

cimport libc.stdlib

from flint.flint_base.flint_context cimport getprec
from flint.flint_base.flint_base cimport flint_poly
from flint.utils.typecheck cimport typecheck
from flint.types.fmpz cimport fmpz_set_python
from flint.types.fmpz cimport any_as_fmpz
from flint.types.fmpz cimport fmpz
from flint.types.fmpq cimport any_as_fmpq
from flint.types.fmpq_poly cimport fmpq_poly
from flint.types.fmpq_poly cimport any_as_fmpq_poly
from flint.types.acb cimport acb
from flint.types.arb cimport any_as_arb_or_notimplemented
from flint.types.arb cimport arb
from flint.types.acb cimport any_as_acb_or_notimplemented
from flint.flintlib.functions.fmpz cimport fmpz_init, fmpz_clear, fmpz_set
from flint.flintlib.functions.fmpz cimport fmpz_is_zero, fmpz_is_one, fmpz_equal_si, fmpz_equal
from flint.flintlib.functions.acb_modular cimport *
from flint.flintlib.functions.ulong_extras cimport n_is_prime
from flint.flintlib.functions.fmpz_poly cimport *
from flint.flintlib.functions.fmpz_poly_factor cimport *
from flint.flintlib.functions.arith cimport *
from flint.flintlib.types.arith cimport arith_chebyshev_t_polynomial, arith_chebyshev_u_polynomial
from flint.flintlib.functions.acb cimport *
from flint.flintlib.functions.arb_poly cimport *
from flint.flintlib.functions.arb_fmpz_poly cimport *
from flint.flintlib.functions.fmpz_vec cimport _fmpz_vec_content

from flint.utils.flint_exceptions import DomainError


cdef any_as_fmpz_poly(x):
    cdef fmpz_poly res
    if typecheck(x, fmpz_poly):
        return x
    elif typecheck(x, fmpz):
        res = fmpz_poly.__new__(fmpz_poly)
        fmpz_poly_set_fmpz(res.val, (<fmpz>x).val)
        return res
    elif PyLong_Check(x):
        res = fmpz_poly.__new__(fmpz_poly)
        t = fmpz(x)
        fmpz_poly_set_fmpz(res.val, (<fmpz>t).val)
        return res
    return NotImplemented

cdef fmpz_poly_set_list(fmpz_poly_t poly, list val):
    cdef long i, n
    cdef fmpz_t x
    n = PyList_GET_SIZE(val)
    fmpz_poly_fit_length(poly, n)
    fmpz_init(x)
    for i from 0 <= i < n:
        if typecheck(val[i], fmpz):
            fmpz_poly_set_coeff_fmpz(poly, i, (<fmpz>(val[i])).val)
        elif fmpz_set_python(x, val[i]):
            fmpz_poly_set_coeff_fmpz(poly, i, x)
        else:
            fmpz_clear(x)
            raise TypeError("unsupported coefficient in list")
    fmpz_clear(x)

cdef class fmpz_poly(flint_poly):
    """
    The *fmpz_poly* type represents dense univariate polynomials over
    the integers.

        >>> fmpz_poly([1,2,3]) ** 3
        27*x^6 + 54*x^5 + 63*x^4 + 44*x^3 + 21*x^2 + 6*x + 1
        >>> divmod(fmpz_poly([2,0,1,1,6]), fmpz_poly([3,5,7]))
        (0, 6*x^4 + x^3 + x^2 + 2)
    """

    def __cinit__(self):
        fmpz_poly_init(self.val)

    def __dealloc__(self):
        fmpz_poly_clear(self.val)

    def __init__(self, *args):
        if not args:
            return
        elif len(args) == 1:
            val = args[0]
        else:
            raise TypeError("fmpz_poly() takes 0 or 1 arguments")
        if typecheck(val, fmpz_poly):
            fmpz_poly_set(self.val, (<fmpz_poly>val).val)
        elif isinstance(val, list):
            fmpz_poly_set_list(self.val, val)
        elif (v := any_as_fmpz(val)) is not NotImplemented:
            fmpz_poly_set_fmpz(self.val, (<fmpz>v).val)
        else:
            raise TypeError("cannot create fmpz_poly from input of type %s", type(val))

    def __len__(self):
        return fmpz_poly_length(self.val)

    cpdef long length(self):
        return fmpz_poly_length(self.val)

    cpdef long degree(self):
        return fmpz_poly_degree(self.val)

    def __richcmp__(self, other, int op):
        cdef bint r
        if op != 2 and op != 3:
            raise TypeError("polynomials cannot be ordered")
        self = any_as_fmpz_poly(self)
        if self is NotImplemented:
            return self
        other = any_as_fmpz_poly(other)
        if other is NotImplemented:
            return other
        r = fmpz_poly_equal((<fmpz_poly>self).val, (<fmpz_poly>other).val)
        if op == 3:
            r = not r
        return r

    def __getitem__(self, long i):
        cdef fmpz x
        x = fmpz()
        if i < 0:
            return x
        fmpz_poly_get_coeff_fmpz(x.val, self.val, i)
        return x

    def __setitem__(self, long i, x):
        if i < 0:
            raise ValueError("cannot assign to index < 0 of polynomial")
        v = fmpz(x)  # XXX
        fmpz_poly_set_coeff_fmpz(self.val, i, (<fmpz>v).val)

    def repr(self):
        return "fmpz_poly([%s])" % (", ".join(map(str, self.coeffs())))

    def __bool__(self):
        return not fmpz_poly_is_zero(self.val)

    def is_zero(self):
        """
        True if this polynomial is the zero polynomial.

        >>> fmpz_poly([]).is_zero()
        True
        """
        return <bint>fmpz_poly_is_zero(self.val)

    def is_one(self):
        """
        True if this polynomial is equal to one.

        >>> fmpz_poly([2]).is_one()
        False
        """
        return <bint>fmpz_poly_is_one(self.val)

    def is_constant(self):
        """
        True if this is a constant polynomial.

        >>> x = fmpz_poly([0, 1])
        >>> two = fmpz_poly([2])
        >>> x.is_constant()
        False
        >>> two.is_constant()
        True
        """
        return fmpz_poly_degree(self.val) <= 0

    def leading_coefficient(self):
        """
        Returns the leading coefficient of the polynomial.

        >>> f = fmpz_poly([1, 2, 3])
        >>> f
        3*x^2 + 2*x + 1
        >>> f.leading_coefficient()
        3
        """
        cdef fmpz x
        cdef slong d
        d = fmpz_poly_degree(self.val)
        x = fmpz.__new__(fmpz)
        if d >= 0:
            fmpz_poly_get_coeff_fmpz(x.val, self.val, d)
        return x

    def __call__(self, other):
        t = any_as_fmpz(other)
        if t is not NotImplemented:
            v = fmpz.__new__(fmpz)
            fmpz_poly_evaluate_fmpz((<fmpz>v).val, self.val, (<fmpz>t).val)
            return v
        t = any_as_fmpz_poly(other)
        if t is not NotImplemented:
            v = fmpz_poly.__new__(fmpz_poly)
            fmpz_poly_compose((<fmpz_poly>v).val, self.val, (<fmpz_poly>t).val)
            return v
        t = any_as_fmpq(other)
        if t is not NotImplemented:
            return fmpq_poly(self)(t)
        t = any_as_fmpq_poly(other)
        if t is not NotImplemented:
            return fmpq_poly(self)(t)
        t = any_as_arb_or_notimplemented(other)
        if t is not NotImplemented:
            v = arb.__new__(arb)
            arb_fmpz_poly_evaluate_arb((<arb>v).val, self.val, (<arb>t).val, getprec())
            return v
        t = any_as_acb_or_notimplemented(other)
        if t is not NotImplemented:
            v = acb.__new__(acb)
            arb_fmpz_poly_evaluate_acb((<acb>v).val, self.val, (<acb>t).val, getprec())
            return v
        raise TypeError("cannot call fmpz_poly with input of type %s", type(other))

    def derivative(self):
        cdef fmpz_poly res = fmpz_poly.__new__(fmpz_poly)
        fmpz_poly_derivative(res.val, self.val)
        return res

    def __pos__(self):
        return self

    def __neg__(self):
        cdef fmpz_poly res = fmpz_poly.__new__(fmpz_poly)
        fmpz_poly_neg(res.val, self.val)
        return res

    def _add_(self, other):
        cdef fmpz_poly res
        other = any_as_fmpz_poly(other)
        if other is NotImplemented:
            return other
        res = fmpz_poly.__new__(fmpz_poly)
        fmpz_poly_add(res.val, (<fmpz_poly>self).val, (<fmpz_poly>other).val)
        return res

    def __add__(self, other):
        return self._add_(other)

    def __radd__(self, other):
        return self._add_(other)

    def __sub__(self, other):
        cdef fmpz_poly res
        other = any_as_fmpz_poly(other)
        if other is NotImplemented:
            return other
        res = fmpz_poly.__new__(fmpz_poly)
        fmpz_poly_sub(res.val, (<fmpz_poly>self).val, (<fmpz_poly>other).val)
        return res

    def __rsub__(self, other):
        cdef fmpz_poly res
        other = any_as_fmpz_poly(other)
        if other is NotImplemented:
            return other
        res = fmpz_poly.__new__(fmpz_poly)
        fmpz_poly_sub(res.val, (<fmpz_poly>other).val, (<fmpz_poly>self).val)
        return res

    def _mul_(self, other):
        cdef fmpz_poly res
        other = any_as_fmpz_poly(other)
        if other is NotImplemented:
            return other
        res = fmpz_poly.__new__(fmpz_poly)
        fmpz_poly_mul(res.val, (<fmpz_poly>self).val, (<fmpz_poly>other).val)
        return res

    def __mul__(self, other):
        return self._mul_(other)

    def __rmul__(self, other):
        return self._mul_(other)

    def __truediv__(fmpz_poly self, other):
        cdef fmpz_poly res
        o = any_as_fmpz(other)
        if o is NotImplemented:
            o = any_as_fmpz_poly(other)
            if o is NotImplemented:
                return NotImplemented
            if fmpz_poly_is_zero((<fmpz_poly>o).val):
                raise ZeroDivisionError("fmpz_poly division by 0")
            res, r = self._divmod_(o)
            if r:
                raise DomainError("fmpz_poly division is not exact")
        else:
            if fmpz_is_zero((<fmpz>o).val):
                raise ZeroDivisionError("fmpz_poly division by 0")
            res = fmpz_poly.__new__(fmpz_poly)
            fmpz_poly_scalar_divexact_fmpz(res.val, self.val, (<fmpz>o).val)
            # Check division is exact - there should be a better way to do this
            if res * o != self:
                raise DomainError("fmpz_poly division is not exact")
        return res

    def __rtruediv__(fmpz_poly self, other):
        o = any_as_fmpz_poly(other)
        if o is NotImplemented:
            return NotImplemented
        return o / self

    def _floordiv_(self, other):
        cdef fmpz_poly res
        if fmpz_poly_is_zero((<fmpz_poly>other).val):
            raise ZeroDivisionError("fmpz_poly division by 0")
        res = fmpz_poly.__new__(fmpz_poly)
        fmpz_poly_div(res.val, (<fmpz_poly>self).val, (<fmpz_poly>other).val)
        return res

    def __floordiv__(self, other):
        other = any_as_fmpz_poly(other)
        if other is NotImplemented:
            return other
        return self._floordiv_(other)

    def __rfloordiv__(self, other):
        other = any_as_fmpz_poly(other)
        if other is NotImplemented:
            return other
        return other._floordiv_(self)

    def _mod_(self, other):
        cdef fmpz_poly res
        if fmpz_poly_is_zero((<fmpz_poly>other).val):
            raise ZeroDivisionError("fmpz_poly division by 0")
        res = fmpz_poly.__new__(fmpz_poly)
        fmpz_poly_rem(res.val, (<fmpz_poly>self).val, (<fmpz_poly>other).val)
        return res

    def __mod__(self, other):
        other = any_as_fmpz_poly(other)
        if other is NotImplemented:
            return other
        return self._mod_(other)

    def __rmod__(self, other):
        other = any_as_fmpz_poly(other)
        if other is NotImplemented:
            return other
        return other._mod_(self)

    def _divmod_(self, other):
        cdef fmpz_poly P, Q
        if fmpz_poly_is_zero((<fmpz_poly>other).val):
            raise ZeroDivisionError("fmpz_poly divmod by 0")
        P = fmpz_poly.__new__(fmpz_poly)
        Q = fmpz_poly.__new__(fmpz_poly)
        fmpz_poly_divrem(P.val, Q.val, (<fmpz_poly>self).val, (<fmpz_poly>other).val)
        return P, Q

    def __divmod__(self, other):
        other = any_as_fmpz_poly(other)
        if other is NotImplemented:
            return other
        return self._divmod_(other)

    def __rdivmod__(self, other):
        other = any_as_fmpz_poly(other)
        if other is NotImplemented:
            return other
        return other._divmod_(self)

    def __pow__(fmpz_poly self, exp, mod):
        cdef fmpz_poly res
        if mod is not None:
            raise NotImplementedError("fmpz_poly modular exponentiation")
        if exp < 0:
            if not fmpz_poly_is_unit(self.val):
                raise DomainError("fmpz_poly negative exponent, non-unit base")
            exp = -exp
        res = fmpz_poly.__new__(fmpz_poly)
        fmpz_poly_pow(res.val, self.val, <ulong>exp)
        return res

    def gcd(self, other):
        """
        Returns the greatest common divisor of self and other.

            >>> A = fmpz_poly([2,0,1,0,5]); B = fmpz_poly([2,3,4])
            >>> (A*B).gcd(B)
            4*x^2 + 3*x + 2

        """
        cdef fmpz_poly res
        other = any_as_fmpz_poly(other)
        if other is NotImplemented:
            raise TypeError("cannot convert input to fmpz_poly")
        res = fmpz_poly.__new__(fmpz_poly)
        fmpz_poly_gcd(res.val, self.val, (<fmpz_poly>other).val)
        return res

    def factor(self):
        """
        Factors self into irreducible factors, returning a tuple
        (c, factors) where c is the content of the coefficients and
        factors is a list of (poly, exp) pairs.

            >>> (-73 * fmpz_poly([1,2,3]) ** 3 * fmpz_poly([5,6,7,8,9]) ** 8).factor()
            (-73, [(3*x^2 + 2*x + 1, 3), (9*x^4 + 8*x^3 + 7*x^2 + 6*x + 5, 8)])
            >>> fmpz_poly.chebyshev_t(6).factor()
            (1, [(2*x^2 + (-1), 1), (16*x^4 + (-16)*x^2 + 1, 1)])
            >>> (fmpz_poly([-1,1])**100).factor()
            (1, [(x + (-1), 100)])
            >>> fmpz_poly([1,2,3,4,5,6]).factor()
            (1, [(6*x^5 + 5*x^4 + 4*x^3 + 3*x^2 + 2*x + 1, 1)])

        """
        return self._factor('irreducible')

    def factor_squarefree(self):
        """
        Factors self into square-free factors, returning a tuple
        (c, factors) where c is the content of the coefficients and
        factors is a list of (poly, exp) pairs.

            >>> x = fmpz_poly([0, 1])
            >>> p = (-3 * x**2 * (x + 1)**2 * (x - 1)**3)
            >>> p.factor_squarefree()
            (-3, [(x^2 + x, 2), (x + (-1), 3)])
            >>> p.factor()
            (-3, [(x, 2), (x + 1, 2), (x + (-1), 3)])

        """
        return self._factor('squarefree')

    def _factor(self, factor_type):
        cdef fmpz_poly_factor_t fac
        cdef int i
        fmpz_poly_factor_init(fac)

        if factor_type == 'squarefree':
            fmpz_poly_factor_squarefree(fac, self.val)
        elif factor_type == 'irreducible':
            fmpz_poly_factor(fac, self.val)
        else:
            assert False

        res = [0] * fac.num
        for 0 <= i < fac.num:
            u = fmpz_poly.__new__(fmpz_poly)
            fmpz_poly_set((<fmpz_poly>u).val, &fac.p[i])
            exp = fac.exp[i]
            res[i] = (u, exp)
        c = fmpz.__new__(fmpz)
        fmpz_set((<fmpz>c).val, &fac.c)
        fmpz_poly_factor_clear(fac)

        return c, res

    def complex_roots(self, bint verbose=False):
        """
        Computes all the complex roots of this polynomial.
        Returns a list of pairs (*c*, *m*) where *c* is the root
        as an *acb* and *m* is the multiplicity of the root.

            >>> fmpz_poly([]).complex_roots()
            []
            >>> fmpz_poly([1]).complex_roots()
            []
            >>> fmpz_poly([2,0,1]).complex_roots()
            [([1.41421356237310 +/- 4.96e-15]j, 1), ([-1.41421356237310 +/- 4.96e-15]j, 1)]
            >>> for c, m in (fmpz_poly([2,3,4]) * fmpz_poly([5,6,7,11])**3).complex_roots():
            ...     print(f'{float(c.real)} + {float(c.imag)}*j : {m}')
            ...
            -0.375 + 0.5994789404140899*j : 1
            -0.375 + -0.5994789404140899*j : 1
            -0.7352847274048426 + 0.0*j : 3
            0.04946054552060311 + 0.7846931676471847*j : 3
            0.04946054552060311 + -0.7846931676471847*j : 3

        """
        cdef fmpz_poly_factor_t fac
        cdef long deg, i, j
        cdef int exp, flags
        cdef acb_ptr croots
        if not self:
            return []
        flags = 0
        if verbose:
            flags = 1
        roots = []
        fmpz_poly_factor_init(fac)
        fmpz_poly_factor_squarefree(fac, self.val)
        for 0 <= i < fac.num:
            deg = fmpz_poly_degree(&fac.p[i])
            exp = fac.exp[i]
            croots = _acb_vec_init(deg)
            arb_fmpz_poly_complex_roots(croots, &fac.p[i], flags, getprec())
            for 0 <= j < deg:
                v = acb()
                acb_set(v.val, &croots[j])
                roots.append((v, exp))
            _acb_vec_clear(croots, deg)
        fmpz_poly_factor_clear(fac)
        return roots

    @staticmethod
    def cyclotomic(ulong n):
        r"""
        Returns the cyclotomic polynomial `\Phi_n(x)` as an *fmpz_poly*.

            >>> fmpz_poly.cyclotomic(12)
            x^4 + (-1)*x^2 + 1

        """
        u = fmpz_poly.__new__(fmpz_poly)
        fmpz_poly_cyclotomic((<fmpz_poly>u).val, n)
        return u

    @staticmethod
    def cos_minpoly(ulong n):
        r"""
        Returns the monic polynomial of `2 \cos(2 \pi / n)` as an *fmpz_poly*.

            >>> fmpz_poly.cos_minpoly(7)
            x^3 + x^2 + (-2)*x + (-1)

        """
        u = fmpz_poly.__new__(fmpz_poly)
        fmpz_poly_cos_minpoly((<fmpz_poly>u).val, n)
        return u

    @staticmethod
    def chebyshev_t(n):
        r"""
        Returns the Chebyshev polynomial of the first kind `T_n(x)`
        as an *fmpz_poly*.

            >>> fmpz_poly.chebyshev_t(3)
            4*x^3 + (-3)*x

        """
        cdef fmpz_poly v = fmpz_poly()
        arith_chebyshev_t_polynomial(v.val, n)
        return v

    @staticmethod
    def chebyshev_u(n):
        r"""
        Returns the Chebyshev polynomial of the second kind `U_n(x)`
        as an *fmpz_poly*.

            >>> fmpz_poly.chebyshev_u(3)
            8*x^3 + (-4)*x

        """
        cdef fmpz_poly v = fmpz_poly()
        arith_chebyshev_u_polynomial(v.val, n)
        return v

    @staticmethod
    def swinnerton_dyer(ulong n, bint use_arb=True):
        r"""
        Returns the Swinnerton-Dyer polynomial `S_n(x)` as an *fmpz_poly*.
        Warning: the degree is `2^n`.

            >>> fmpz_poly.swinnerton_dyer(0)
            x
            >>> fmpz_poly.swinnerton_dyer(1)
            x^2 + (-2)
            >>> fmpz_poly.swinnerton_dyer(2)
            x^4 + (-10)*x^2 + 1
            >>> fmpz_poly.swinnerton_dyer(3)
            x^8 + (-40)*x^6 + 352*x^4 + (-960)*x^2 + 576

        """
        cdef arb_poly_t t
        if n > 20:
            raise OverflowError("that's way too large...")
        u = fmpz_poly.__new__(fmpz_poly)
        if use_arb:
            arb_poly_init(t)
            arb_poly_swinnerton_dyer_ui(t, n, 0)
            if not arb_poly_get_unique_fmpz_poly((<fmpz_poly>u).val, t):
                arb_poly_clear(t)
                raise ValueError("insufficient precision")
            arb_poly_clear(t)
        else:
            fmpz_poly_swinnerton_dyer((<fmpz_poly>u).val, n)
        return u

    @staticmethod
    def hilbert_class_poly(long D):
        r"""
        Returns the Hilbert class polynomial `H_D(x)` as an *fmpz_poly*.

            >>> fmpz_poly.hilbert_class_poly(-3)
            x
            >>> fmpz_poly.hilbert_class_poly(-4)
            x + (-1728)
            >>> fmpz_poly.hilbert_class_poly(-59)
            x^3 + 30197678080*x^2 + (-140811576541184)*x + 374643194001883136
            >>> fmpz_poly.hilbert_class_poly(-5)
            Traceback (most recent call last):
              ...
            ValueError: D must be an imaginary quadratic discriminant
        """
        cdef fmpz_poly v = fmpz_poly()
        acb_modular_hilbert_class_poly(v.val, D)
        if fmpz_poly_length(v.val) == 0:
            raise ValueError("D must be an imaginary quadratic discriminant")
        return v

    def height_bits(self, bint signed=False):
        if signed:
            return fmpz_poly_max_bits(self.val)
        else:
            return abs(fmpz_poly_max_bits(self.val))

    def sqrt(self):
        cdef fmpz_poly v = fmpz_poly()
        if fmpz_poly_sqrt(v.val, self.val):
            return v
        else:
            raise DomainError(f"Cannot compute square root of {self}")

    def inflate(self, n: int) -> fmpz_poly:
        """
        Compute the inflation of ``self`` for a provided ``n``, that is return ``q``
        such that ``q(x) = p(x^n)``.

            >>> f = fmpz_poly([1, 1])
            >>> f.inflate(2)
            x^2 + 1
        """
        cdef fmpz_poly res = fmpz_poly()
        fmpz_poly_inflate(res.val, self.val, n)
        return res

    def deflate(self, n: int) -> fmpz_poly:
        """
        Compute the deflation of ``self`` for a provided ``n``, that is return ``q``
        such that ``q(x) = p(x^(1/n))``.

            >>> f = fmpz_poly([1, 0, 1])
            >>> f.deflate(2)
            x + 1
        """
        cdef fmpz_poly res = fmpz_poly()
        if n > 0:
            fmpz_poly_deflate(res.val, self.val, n)
            return res
        else:
            raise ValueError("deflate requires n > 0")

    def deflation(self) -> tuple[fmpz_poly, int]:
        """
        Compute the deflation of ``self``, that is ``p(x^(1/n))`` for maximal
        n. returns ``q, n`` such that ``self == q.inflate(n)``.

            >>> f = fmpz_poly([1, 0, 1])
            >>> q, n = f.deflation()
            >>> q, n
            (x + 1, 2)
            >>> q.inflate(n) == f
            True
        """
        cdef ulong n
        if fmpz_poly_is_zero(self.val):
            return self, 1
        n = fmpz_poly_deflation(self.val)
        return self if n <= 1 else self.deflate(n), int(n)

    def deflation_monom(self) -> tuple[fmpz_poly, int, fmpz_poly]:
        """
        Compute the exponent ``n`` and monomial ``m`` such that ``p(x^(1/n)) = m *
        q(x^n)`` for maximal n. The returned monomial allows the undo-ing of the
        deflation.

            >>> f = fmpz_poly([1, 0, 1])
            >>> f.deflation_monom()
            (x^2 + 1, 1, x)
        """
        n, m = self.deflation_index()

        cdef fmpz_poly monom = fmpz_poly.__new__(fmpz_poly)
        cdef fmpz_poly res = fmpz_poly.__new__(fmpz_poly)

        fmpz_poly_set_coeff_ui(monom.val, m, 1)
        fmpz_poly_deflate(res.val, self.val, n)

        return res, n, monom

    def deflation_index(self) -> tuple[int, int]:
        """
        Compute the exponent ``n`` and ``i`` such that ``p(x^(1/n)) = x^i *
        q(x^n)`` for maximal ``n``. Importantly the deflation itself is not computed
        here. The returned exponent ``i`` is the shift that was applied to the
        exponents. It is the exponent of the monomial returned by
        ``deflation_monom``.

            >>> f = fmpz_poly([1, 0, 1])
            >>> f.deflation_index()
            (1, 1)
        """
        cdef fmpz_poly res = fmpz_poly.__new__(fmpz_poly)
        cdef slong length = fmpz_poly_length(self.val)

        if length <= 0:
            return self, 0, fmpz_poly([1])

        # Find the smallest non-zero power, that is the gcd of the monomials
        for i in range(1, length + 1):
            if not fmpz_is_zero(&self.val.coeffs[length - i]):
                break

        fmpz_poly_shift_right(res.val, self.val, i)
        return int(fmpz_poly_deflation(res.val)), int(i)

    def is_cyclotomic(self):
        cdef long * phi
        cdef long i, p, q, d, N1, N2
        cdef double U
        d = self.degree()
        if d < 1:
            return 0
        if d == 1:
            if fmpz_is_one(fmpz_poly_get_coeff_ptr(self.val, 1)) and fmpz_equal_si(fmpz_poly_get_coeff_ptr(self.val, 0), -1):
                return 1
            if fmpz_is_one(fmpz_poly_get_coeff_ptr(self.val, 1)) and fmpz_equal_si(fmpz_poly_get_coeff_ptr(self.val, 0), 1):
                return 2
            return 0
        if d % 2 != 0:
            return 0
        if not fmpz_is_one(fmpz_poly_get_coeff_ptr(self.val, 0)):
            return 0
        for i in range(d//2):
            if not fmpz_equal(fmpz_poly_get_coeff_ptr(self.val, i), fmpz_poly_get_coeff_ptr(self.val, d-i)):
                return 0
        U = d
        for p in range(2, d+1):
            if d % (p-1) == 0:
                if n_is_prime(p):
                    U = (U * p) / (p - 1)
        N1 = d + 1
        N2 = int(U + 3)   # +3 as safety for possible float roundoff
        phi = <long *> libc.stdlib.malloc(N2 * sizeof(long))
        for i in range(N2):
            phi[i] = i
        for p in range(2, N2):
            if phi[p] == p:
                phi[p] = p - 1
                for q in range(2*p, N2, p):
                    phi[q] = (phi[q] // p) * (p-1)
        for i in range(N1, N2):
            if phi[i] == d:
                v = fmpz_poly.cyclotomic(i)
                if self == v:
                    return int(i)
        libc.stdlib.free(phi)
        return 0

    def content(self):
        """
        Return the GCD of the coefficients of ``self``.

            >>> fmpz_poly([3, 6, 0]).content()
            3
        """
        cdef fmpz res = fmpz()
        _fmpz_vec_content(res.val, self.val.coeffs, self.val.length)
        return res
