from cpython.list cimport PyList_GET_SIZE
from flint.flint_base.flint_base cimport flint_poly
from flint.utils.typecheck cimport typecheck
from flint.types.fmpz cimport fmpz, any_as_fmpz
from flint.types.fmpz_poly cimport any_as_fmpz_poly
from flint.types.fmpz_poly cimport fmpz_poly
from flint.types.nmod cimport any_as_nmod
from flint.types.nmod cimport nmod

from flint.flintlib.functions.nmod cimport nmod_init
from flint.flintlib.functions.nmod_poly cimport *
from flint.flintlib.functions.nmod_poly_factor cimport *
from flint.flintlib.functions.fmpz_poly cimport fmpz_poly_get_nmod_poly

from flint.utils.flint_exceptions import DomainError


cdef any_as_nmod_poly(obj, nmod_t mod):
    cdef nmod_poly r
    cdef mp_limb_t v
    # XXX: should check that modulus is the same here, and not all over the place
    if typecheck(obj, nmod_poly):
        return obj
    if any_as_nmod(&v, obj, mod):
        r = nmod_poly.__new__(nmod_poly)
        nmod_poly_init(r.val, mod.n)
        nmod_poly_set_coeff_ui(r.val, 0, v)
        return r
    x = any_as_fmpz_poly(obj)
    if x is not NotImplemented:
        r = nmod_poly.__new__(nmod_poly)
        nmod_poly_init(r.val, mod.n)   # XXX: create flint _nmod_poly_set_modulus for this?
        fmpz_poly_get_nmod_poly(r.val, (<fmpz_poly>x).val)
        return r
    return NotImplemented

cdef nmod_poly_set_list(nmod_poly_t poly, list val):
    cdef long i, n
    cdef nmod_t mod
    cdef mp_limb_t v
    nmod_init(&mod, nmod_poly_modulus(poly))  # XXX
    n = PyList_GET_SIZE(val)
    nmod_poly_fit_length(poly, n)
    for i from 0 <= i < n:
        if any_as_nmod(&v, val[i], mod):
            nmod_poly_set_coeff_ui(poly, i, v)
        else:
            raise TypeError("unsupported coefficient in list")

cdef class nmod_poly(flint_poly):
    """
    The nmod_poly type represents dense univariate polynomials
    over Z/nZ for word-size n.

        >>> from flint import ctx
        >>> a = nmod_poly([5,1,10,14,8], 7)
        >>> a
        x^4 + 3*x^2 + x + 5
        >>> -a
        6*x^4 + 4*x^2 + 6*x + 2
        >>> ctx.pretty = False
        >>> list(nmod_poly(list(range(3)), 2))
        [nmod(0, 2), nmod(1, 2)]
        >>> nmod_poly([1, 2, 3], 23) ** 3
        nmod_poly([1, 6, 21, 21, 17, 8, 4], 23)
        >>> divmod(nmod_poly([2,0,1,1,6],23), nmod_poly([3,5,7],23))
        (nmod_poly([4, 0, 14], 23), nmod_poly([13, 3], 23))
        >>> ctx.pretty = True

    """

    # cdef nmod_poly_t val

    def __cinit__(self):
        nmod_poly_init(self.val, 1)

    def __dealloc__(self):
        nmod_poly_clear(self.val)

    def __init__(self, val=None, ulong mod=0):
        cdef ulong m2
        cdef mp_limb_t v
        if typecheck(val, nmod_poly):
            m2 = nmod_poly_modulus((<nmod_poly>val).val)
            if m2 != mod:
                raise ValueError("different moduli!")
            nmod_poly_init(self.val, m2)
            nmod_poly_set(self.val, (<nmod_poly>val).val)
        else:
            if mod == 0:
                raise ValueError("a nonzero modulus is required")
            nmod_poly_init(self.val, mod)
            if typecheck(val, fmpz_poly):
                fmpz_poly_get_nmod_poly(self.val, (<fmpz_poly>val).val)
            elif typecheck(val, list):
                nmod_poly_set_list(self.val, val)
            elif any_as_nmod(&v, val, self.val.mod):
                nmod_poly_fit_length(self.val, 1)
                nmod_poly_set_coeff_ui(self.val, 0, v)
            else:
                raise TypeError("cannot create nmod_poly from input of type %s", type(val))

    def __len__(self):
        return nmod_poly_length(self.val)

    cpdef long length(self):
        return nmod_poly_length(self.val)

    cpdef long degree(self):
        return nmod_poly_degree(self.val)

    cpdef mp_limb_t modulus(self):
        return nmod_poly_modulus(self.val)

    def __richcmp__(s, t, int op):
        cdef bint res
        if op != 2 and op != 3:
            raise TypeError("nmod_polys cannot be ordered")
        if typecheck(s, nmod_poly) and typecheck(t, nmod_poly):
            if (<nmod_poly>s).val.mod.n != (<nmod_poly>t).val.mod.n:
                res = False
            else:
                res = nmod_poly_equal((<nmod_poly>s).val, (<nmod_poly>t).val)
            if op == 2:
                return res
            if op == 3:
                return not res
        else:
            if not typecheck(s, nmod_poly):
                s, t = t, s
            try:
                t = nmod_poly([t], (<nmod_poly>s).val.mod.n)
            except TypeError:
                pass
            if typecheck(s, nmod_poly) and typecheck(t, nmod_poly):
                if (<nmod_poly>s).val.mod.n != (<nmod_poly>t).val.mod.n:
                    res = False
                else:
                    res = nmod_poly_equal((<nmod_poly>s).val, (<nmod_poly>t).val)
                if op == 2:
                    return res
                if op == 3:
                    return not res
        return NotImplemented

    def __iter__(self):
        cdef long i, n
        n = self.length()
        for i from 0 <= i < n:
            yield self[i]

    def coeffs(self):
        cdef long i, n
        cdef list L
        cdef mp_limb_t m
        n = self.length()
        m = self.modulus()
        L = [nmod(0, m) for i in range(n)]   # XXX: speed up
        for i from 0 <= i < n:
            (<nmod>(L[i])).val = nmod_poly_get_coeff_ui(self.val, i)
        return L

    def repr(self):
        return "nmod_poly(%s, %s)" % ([int(c) for c in self.coeffs()], self.modulus())

    def __getitem__(self, long i):
        cdef nmod x
        x = nmod(0, self.modulus())
        if i < 0:
            return x
        x.val = nmod_poly_get_coeff_ui(self.val, i)
        return x

    def __setitem__(self, long i, x):
        cdef mp_limb_t v
        if i < 0:
            raise ValueError("cannot assign to index < 0 of polynomial")
        if any_as_nmod(&v, x, self.val.mod):
            nmod_poly_set_coeff_ui(self.val, i, v)
        else:
            raise TypeError("cannot set element of type %s" % type(x))

    def __bool__(self):
        return not nmod_poly_is_zero(self.val)

    def is_zero(self):
        """
        Returns True if this is the zero polynomial.
        """
        return <bint>nmod_poly_is_zero(self.val)

    def is_one(self):
        """
        Returns True if this polynomial is equal to 1.
        """
        return <bint>nmod_poly_is_one(self.val)

    def is_constant(self):
        """
        Returns True if this is a constant polynomial.

        >>> nmod_poly([0, 1], 3).is_constant()
        False
        >>> nmod_poly([1], 3).is_constant()
        True
        """
        return nmod_poly_degree(self.val) <= 0

    def is_gen(self):
        """
        Returns True if this polynomial is equal to the generator x.

        >>> x = nmod_poly([0, 1], 3)
        >>> x
        x
        >>> x.is_gen()
        True
        >>> (2*x).is_gen()
        False
        """
        return <bint>nmod_poly_is_gen(self.val)

    def reverse(self, degree=None):
        """
        Return a polynomial with the coefficients of this polynomial
        reversed.

        If ``degree`` is not None, the output polynomial will be zero-padded
        or truncated before being reversed. NOTE: degree must be non-negative.

            >>> f = nmod_poly([1,2,3,4,5], 163)
            >>> f.reverse()
            x^4 + 2*x^3 + 3*x^2 + 4*x + 5
            >>> f.reverse(degree=1)
            x + 2
            >>> f.reverse(degree=100)
            x^100 + 2*x^99 + 3*x^98 + 4*x^97 + 5*x^96
        """
        cdef nmod_poly res
        cdef slong d

        if degree is not None:
            d = degree
            if d != degree or d < 0:
                raise ValueError(f"degree argument must be a non-negative integer, got {degree}")
            length = d + 1
        else:
            length = nmod_poly_length(self.val)

        res = nmod_poly.__new__(nmod_poly)
        nmod_poly_init_preinv(res.val, self.val.mod.n, self.val.mod.ninv)
        nmod_poly_reverse(res.val, self.val, length)
        return res

    def leading_coefficient(self):
        """
        Return the leading coefficient of this polynomial.

            >>> f = nmod_poly([123, 129, 63, 14, 51, 76, 133], 163)
            >>> f.leading_coefficient()
            133
        """
        cdef ulong cu
        cdef slong d
        cdef nmod c

        d = self.degree()
        if d < 0:
            cu = 0
        else:
            cu = nmod_poly_get_coeff_ui(self.val, d)

        c = nmod.__new__(nmod)
        c.mod = self.val.mod
        c.val = cu

        return c

    def inverse_series_trunc(self, slong n):
        """
        Returns the inverse of ``self`` modulo `x^n`. Assumes the leading
        coefficient of the polynomial is invertible.

            >>> f = nmod_poly([123, 129, 63, 14, 51, 76, 133], 163)
            >>> f.inverse_series_trunc(3)
            159*x^2 + 151*x + 110
            >>> f.inverse_series_trunc(4)
            23*x^3 + 159*x^2 + 151*x + 110
            >>> f.inverse_series_trunc(5)
            45*x^4 + 23*x^3 + 159*x^2 + 151*x + 110
        """
        if n <= 0:
            raise ValueError(f"n = {n} must be positive")

        if self.is_zero():
            raise ValueError("cannot invert the zero element")

        cdef nmod_poly res
        res = nmod_poly.__new__(nmod_poly)
        nmod_poly_init_preinv(res.val, self.val.mod.n, self.val.mod.ninv)
        nmod_poly_inv_series(res.val, self.val, n)
        return res

    def compose(self, other):
        """
        Returns the composition of two polynomials

        To be precise about the order of composition, given ``self``, and ``other``
        by `f(x)`, `g(x)`, returns `f(g(x))`.

            >>> f = nmod_poly([1,2,3], 163)
            >>> g = nmod_poly([0,0,1], 163)
            >>> f.compose(g)
            3*x^4 + 2*x^2 + 1
            >>> g.compose(f)
            9*x^4 + 12*x^3 + 10*x^2 + 4*x + 1
        """
        cdef nmod_poly res
        other = any_as_nmod_poly(other, (<nmod_poly>self).val.mod)
        if other is NotImplemented:
            raise TypeError("cannot convert input to nmod_poly")
        res = nmod_poly.__new__(nmod_poly)
        nmod_poly_init_preinv(res.val, self.val.mod.n, self.val.mod.ninv)
        nmod_poly_compose(res.val, self.val, (<nmod_poly>other).val)
        return res

    def compose_mod(self, other, modulus):
        r"""
        Returns the composition of two polynomials modulo a third.

        To be precise about the order of composition, given ``self``, and ``other``
        and ``modulus`` by `f(x)`, `g(x)` and `h(x)`, returns `f(g(x)) \mod h(x)`.
        We require that `h(x)` is non-zero.

            >>> f = nmod_poly([1,2,3,4,5], 163)
            >>> g = nmod_poly([3,2,1], 163)
            >>> h = nmod_poly([1,0,1,0,1], 163)
            >>> f.compose_mod(g, h)
            63*x^3 + 100*x^2 + 17*x + 63
            >>> g.compose_mod(f, h)
            147*x^3 + 159*x^2 + 4*x + 7
        """
        cdef nmod_poly res
        g = any_as_nmod_poly(other, self.val.mod)
        if g is NotImplemented:
            raise TypeError(f"cannot convert other = {other} to nmod_poly")

        h = any_as_nmod_poly(modulus, self.val.mod)
        if h is NotImplemented:
            raise TypeError(f"cannot convert modulus = {modulus} to nmod_poly")

        if modulus.is_zero():
            raise ZeroDivisionError("cannot reduce modulo zero")

        res = nmod_poly.__new__(nmod_poly)
        nmod_poly_init_preinv(res.val, self.val.mod.n, self.val.mod.ninv)
        nmod_poly_compose_mod(res.val, self.val, (<nmod_poly>other).val, (<nmod_poly>modulus).val)
        return res

    def __call__(self, other):
        cdef mp_limb_t c
        if any_as_nmod(&c, other, self.val.mod):
            v = nmod(0, self.modulus())
            (<nmod>v).val = nmod_poly_evaluate_nmod(self.val, c)
            return v
        t = any_as_nmod_poly(other, self.val.mod)
        if t is not NotImplemented:
            r = nmod_poly.__new__(nmod_poly)
            nmod_poly_init_preinv((<nmod_poly>r).val, self.val.mod.n, self.val.mod.ninv)
            nmod_poly_compose((<nmod_poly>r).val, self.val, (<nmod_poly>t).val)
            return r
        raise TypeError("cannot call nmod_poly with input of type %s", type(other))

    def derivative(self):
        cdef nmod_poly res = nmod_poly.__new__(nmod_poly)
        nmod_poly_init_preinv(res.val, self.val.mod.n, self.val.mod.ninv)
        nmod_poly_derivative(res.val, self.val)
        return res

    def integral(self):
        cdef nmod_poly res = nmod_poly.__new__(nmod_poly)
        nmod_poly_init_preinv(res.val, self.val.mod.n, self.val.mod.ninv)
        nmod_poly_integral(res.val, self.val)
        return res

    def __pos__(self):
        return self

    def __neg__(self):
        cdef nmod_poly r = nmod_poly.__new__(nmod_poly)
        nmod_poly_init_preinv(r.val, self.val.mod.n, self.val.mod.ninv)
        nmod_poly_neg(r.val, self.val)
        return r

    def _add_(s, t):
        cdef nmod_poly r
        t = any_as_nmod_poly(t, (<nmod_poly>s).val.mod)
        if t is NotImplemented:
            return t
        if (<nmod_poly>s).val.mod.n != (<nmod_poly>t).val.mod.n:
            raise ValueError("cannot add nmod_polys with different moduli")
        r = nmod_poly.__new__(nmod_poly)
        nmod_poly_init_preinv(r.val, (<nmod_poly>t).val.mod.n, (<nmod_poly>t).val.mod.ninv)
        nmod_poly_add(r.val, (<nmod_poly>s).val, (<nmod_poly>t).val)
        return r

    def __add__(s, t):
        return s._add_(t)

    def __radd__(s, t):
        return s._add_(t)

    def _sub_(s, t):
        cdef nmod_poly r
        if (<nmod_poly>s).val.mod.n != (<nmod_poly>t).val.mod.n:
            raise ValueError("cannot subtract nmod_polys with different moduli")
        r = nmod_poly.__new__(nmod_poly)
        nmod_poly_init_preinv(r.val, (<nmod_poly>t).val.mod.n, (<nmod_poly>t).val.mod.ninv)
        nmod_poly_sub(r.val, (<nmod_poly>s).val, (<nmod_poly>t).val)
        return r

    def __sub__(s, t):
        t = any_as_nmod_poly(t, (<nmod_poly>s).val.mod)
        if t is NotImplemented:
            return t
        return s._sub_(t)

    def __rsub__(s, t):
        t = any_as_nmod_poly(t, (<nmod_poly>s).val.mod)
        if t is NotImplemented:
            return t
        return t._sub_(s)

    def _mul_(s, t):
        cdef nmod_poly r
        t = any_as_nmod_poly(t, (<nmod_poly>s).val.mod)
        if t is NotImplemented:
            return t
        if (<nmod_poly>s).val.mod.n != (<nmod_poly>t).val.mod.n:
            raise ValueError("cannot multiply nmod_polys with different moduli")
        r = nmod_poly.__new__(nmod_poly)
        nmod_poly_init_preinv(r.val, (<nmod_poly>t).val.mod.n, (<nmod_poly>t).val.mod.ninv)
        nmod_poly_mul(r.val, (<nmod_poly>s).val, (<nmod_poly>t).val)
        return r

    def __mul__(s, t):
        return s._mul_(t)

    def __rmul__(s, t):
        return s._mul_(t)

    def __truediv__(s, t):
        t = any_as_nmod_poly(t, (<nmod_poly>s).val.mod)
        if t is NotImplemented:
            return t
        res, r = s._divmod_(t)
        if not nmod_poly_is_zero((<nmod_poly>r).val):
            raise DomainError("nmod_poly inexact division")
        return res

    def __rtruediv__(s, t):
        t = any_as_nmod_poly(t, (<nmod_poly>s).val.mod)
        if t is NotImplemented:
            return t
        res, r = t._divmod_(s)
        if not nmod_poly_is_zero((<nmod_poly>r).val):
            raise DomainError("nmod_poly inexact division")
        return res

    def _floordiv_(s, t):
        cdef nmod_poly r
        if (<nmod_poly>s).val.mod.n != (<nmod_poly>t).val.mod.n:
            raise ValueError("cannot divide nmod_polys with different moduli")
        if nmod_poly_is_zero((<nmod_poly>t).val):
            raise ZeroDivisionError("polynomial division by zero")
        r = nmod_poly.__new__(nmod_poly)
        nmod_poly_init_preinv(r.val, (<nmod_poly>t).val.mod.n, (<nmod_poly>t).val.mod.ninv)
        nmod_poly_div(r.val, (<nmod_poly>s).val, (<nmod_poly>t).val)
        return r

    def __floordiv__(s, t):
        t = any_as_nmod_poly(t, (<nmod_poly>s).val.mod)
        if t is NotImplemented:
            return t
        return s._floordiv_(t)

    def __rfloordiv__(s, t):
        t = any_as_nmod_poly(t, (<nmod_poly>s).val.mod)
        if t is NotImplemented:
            return t
        return t._floordiv_(s)

    def _divmod_(s, t):
        cdef nmod_poly P, Q
        if (<nmod_poly>s).val.mod.n != (<nmod_poly>t).val.mod.n:
            raise ValueError("cannot divide nmod_polys with different moduli")
        if nmod_poly_is_zero((<nmod_poly>t).val):
            raise ZeroDivisionError("polynomial division by zero")
        P = nmod_poly.__new__(nmod_poly)
        Q = nmod_poly.__new__(nmod_poly)
        nmod_poly_init_preinv(P.val, (<nmod_poly>t).val.mod.n, (<nmod_poly>t).val.mod.ninv)
        nmod_poly_init_preinv(Q.val, (<nmod_poly>t).val.mod.n, (<nmod_poly>t).val.mod.ninv)
        nmod_poly_divrem(P.val, Q.val, (<nmod_poly>s).val, (<nmod_poly>t).val)
        return P, Q

    def __divmod__(s, t):
        t = any_as_nmod_poly(t, (<nmod_poly>s).val.mod)
        if t is NotImplemented:
            return t
        return s._divmod_(t)

    def __rdivmod__(s, t):
        t = any_as_nmod_poly(t, (<nmod_poly>s).val.mod)
        if t is NotImplemented:
            return t
        return t._divmod_(s)

    def __mod__(s, t):
        return divmod(s, t)[1]      # XXX

    def __rmod__(s, t):
        return divmod(t, s)[1]      # XXX

    def __pow__(nmod_poly self, exp, mod=None):
        cdef nmod_poly res
        if mod is not None:
            return self.pow_mod(exp, mod)
        if exp < 0:
            self = 1 / self
            exp = -exp
        res = nmod_poly.__new__(nmod_poly)
        nmod_poly_init_preinv(res.val, (<nmod_poly>self).val.mod.n, (<nmod_poly>self).val.mod.ninv)
        nmod_poly_pow(res.val, self.val, <ulong>exp)
        return res

    def pow_mod(self, e, modulus, mod_rev_inv=None):
        r"""
        Returns ``self`` raised to the power ``e`` modulo ``modulus``:
        :math:`f^e \mod g`/

        ``mod_rev_inv`` is the inverse of the reverse of the modulus,
        precomputing it and passing it to ``pow_mod()`` can optimise
        powering of polynomials with large exponents.

            >>> x = nmod_poly([0,1], 163)
            >>> f = 30*x**6 + 104*x**5 + 76*x**4 + 33*x**3 + 70*x**2 + 44*x + 65
            >>> g = 43*x**6 + 91*x**5 + 77*x**4 + 113*x**3 + 71*x**2 + 132*x + 60
            >>> mod = x**4 + 93*x**3 + 78*x**2 + 72*x + 149
            >>>
            >>> f.pow_mod(123, mod)
            3*x^3 + 25*x^2 + 115*x + 161
            >>> f.pow_mod(2**64, mod)
            52*x^3 + 96*x^2 + 136*x + 9
            >>> mod_rev_inv = mod.reverse().inverse_series_trunc(4)
            >>> f.pow_mod(2**64, mod, mod_rev_inv)
            52*x^3 + 96*x^2 + 136*x + 9
        """
        cdef nmod_poly res

        if e < 0:
            raise ValueError("Exponent must be non-negative")

        modulus = any_as_nmod_poly(modulus, (<nmod_poly>self).val.mod)
        if modulus is NotImplemented:
            raise TypeError("cannot convert input to nmod_poly")

        # Output polynomial
        res = nmod_poly.__new__(nmod_poly)
        nmod_poly_init_preinv(res.val, self.val.mod.n, self.val.mod.ninv)

        # For small exponents, use a simple binary exponentiation method
        if e.bit_length() < 32:
            nmod_poly_powmod_ui_binexp(
                res.val, self.val, <ulong>e, (<nmod_poly>modulus).val
            )
            return res

        # For larger exponents we need to cast e to an fmpz first
        e_fmpz = any_as_fmpz(e)
        if e_fmpz is NotImplemented:
            raise TypeError(f"exponent cannot be cast to an fmpz type: {e}")

        # To optimise powering, we precompute the inverse of the reverse of the modulus
        if mod_rev_inv is not None:
            mod_rev_inv = any_as_nmod_poly(mod_rev_inv, (<nmod_poly>self).val.mod)
            if mod_rev_inv is NotImplemented:
                raise TypeError(f"Cannot interpret {mod_rev_inv} as a polynomial")
        else:
            mod_rev_inv = modulus.reverse().inverse_series_trunc(modulus.length())

        # Use windowed exponentiation optimisation when self = x
        if self.is_gen():
            nmod_poly_powmod_x_fmpz_preinv(
                res.val, (<fmpz>e_fmpz).val, (<nmod_poly>modulus).val, (<nmod_poly>mod_rev_inv).val
            )
            return res

        # Otherwise using binary exponentiation for all other inputs
        nmod_poly_powmod_fmpz_binexp_preinv(
                res.val, self.val, (<fmpz>e_fmpz).val, (<nmod_poly>modulus).val, (<nmod_poly>mod_rev_inv).val
            )
        return res

    def gcd(self, other):
        """
        Returns the monic greatest common divisor of self and other.

            >>> A = nmod_poly([1,2,3,4], 7); B = nmod_poly([4,1,5], 7)
            >>> (A * B).gcd(B) * 5
            5*x^2 + x + 4

        """
        cdef nmod_poly res
        other = any_as_nmod_poly(other, (<nmod_poly>self).val.mod)
        if other is NotImplemented:
            raise TypeError("cannot convert input to nmod_poly")
        if self.val.mod.n != (<nmod_poly>other).val.mod.n:
            raise ValueError("moduli must be the same")
        res = nmod_poly.__new__(nmod_poly)
        nmod_poly_init_preinv(res.val, self.val.mod.n, self.val.mod.ninv)
        nmod_poly_gcd(res.val, self.val, (<nmod_poly>other).val)
        return res

    def xgcd(self, other):
        cdef nmod_poly res1, res2, res3
        other = any_as_nmod_poly(other, (<nmod_poly>self).val.mod)
        if other is NotImplemented:
            raise TypeError("cannot convert input to fmpq_poly")
        res1 = nmod_poly.__new__(nmod_poly)
        res2 = nmod_poly.__new__(nmod_poly)
        res3 = nmod_poly.__new__(nmod_poly)
        nmod_poly_init(res1.val, (<nmod_poly>self).val.mod.n)
        nmod_poly_init(res2.val, (<nmod_poly>self).val.mod.n)
        nmod_poly_init(res3.val, (<nmod_poly>self).val.mod.n)
        nmod_poly_xgcd(res1.val, res2.val, res3.val, self.val, (<nmod_poly>other).val)
        return (res1, res2, res3)

    def factor(self, algorithm=None):
        """
        Factors self into irreducible factors, returning a tuple
        (c, factors) where c is the leading coefficient and
        factors is a list of (poly, exp) pairs with all polynomials
        monic.

            >>> nmod_poly(list(range(10)), 3).factor()
            (2, [(x, 1), (x + 2, 7)])
            >>> nmod_poly(list(range(10)), 19).factor()
            (9, [(x, 1), (x^4 + 15*x^3 + 2*x^2 + 7*x + 3, 1), (x^4 + 7*x^3 + 12*x^2 + 15*x + 12, 1)])
            >>> nmod_poly(list(range(10)), 53).factor()
            (9, [(x, 1), (x^8 + 48*x^7 + 42*x^6 + 36*x^5 + 30*x^4 + 24*x^3 + 18*x^2 + 12*x + 6, 1)])

        Algorithm can be None (default), 'berlekamp', or
        'cantor-zassenhaus'.

            >>> nmod_poly([3,2,1,2,3], 7).factor(algorithm='berlekamp')
            (3, [(x + 2, 1), (x + 4, 1), (x^2 + 4*x + 1, 1)])
            >>> nmod_poly([3,2,1,2,3], 7).factor(algorithm='cantor-zassenhaus')
            (3, [(x + 4, 1), (x + 2, 1), (x^2 + 4*x + 1, 1)])

        """
        if algorithm is None:
            algorithm = 'irreducible'
        elif algorithm not in ('berlekamp', 'cantor-zassenhaus'):
            raise ValueError(f"unknown factorization algorithm: {algorithm}")
        return self._factor(algorithm)

    def factor_squarefree(self):
        """
        Factors *self* into square-free polynomials. Returns (*c*, *factors*)
        where *c* is the leading coefficient and *factors* is a list of
        (*poly*, *exp*).

            >>> x = nmod_poly([0, 1], 7)
            >>> p = x**2 * (x/2 - 1)**2 * (x + 1)**3
            >>> p
            2*x^7 + 5*x^6 + 4*x^5 + 2*x^4 + 2*x^3 + x^2
            >>> p.factor_squarefree()
            (2, [(x^2 + 5*x, 2), (x + 1, 3)])
            >>> p.factor()
            (2, [(x, 2), (x + 5, 2), (x + 1, 3)])

        """
        return self._factor('squarefree')

    def _factor(self, factor_type):
        cdef nmod_poly_factor_t fac
        cdef mp_limb_t lead
        cdef int i

        nmod_poly_factor_init(fac)

        if factor_type == 'berlekamp':
            lead = nmod_poly_factor_with_berlekamp(fac, self.val)
        elif factor_type == 'cantor-zassenhaus':
            lead = nmod_poly_factor_with_cantor_zassenhaus(fac, self.val)
        elif factor_type == 'irreducible':
            lead = nmod_poly_factor(fac, self.val)
        elif factor_type == 'squarefree':
            nmod_poly_factor_squarefree(fac, self.val)
            lead = (<nmod>self.leading_coefficient()).val
        else:
            assert False

        res = [None] * fac.num
        for 0 <= i < fac.num:
            u = nmod_poly.__new__(nmod_poly)
            nmod_poly_init_preinv((<nmod_poly>u).val,
                                  (<nmod_poly>self).val.mod.n, (<nmod_poly>self).val.mod.ninv)
            nmod_poly_set((<nmod_poly>u).val, &fac.p[i])
            exp = fac.exp[i]
            res[i] = (u, exp)

        c = nmod.__new__(nmod)
        (<nmod>c).mod = self.val.mod
        (<nmod>c).val = lead

        nmod_poly_factor_clear(fac)

        return c, res

    def sqrt(nmod_poly self):
        """Return exact square root or ``None``. """
        cdef nmod_poly res
        res = nmod_poly.__new__(nmod_poly)
        nmod_poly_init_preinv(res.val, self.val.mod.n, self.val.mod.ninv)
        if nmod_poly_sqrt(res.val, self.val):
            return res
        else:
            raise DomainError(f"Cannot compute square root of {self}")

    def deflation(self):
        cdef nmod_poly v
        cdef ulong n
        if nmod_poly_is_zero(self.val):
            return self, 1
        n = nmod_poly_deflation(self.val)
        if n == 1:
            return self, int(n)
        else:
            v = nmod_poly.__new__(nmod_poly)
            nmod_poly_init(v.val, self.val.mod.n)
            nmod_poly_deflate(v.val, self.val, n)
            return v, int(n)

    def real_roots(self):
        r"""
        This method is not implemented for polynomials in
        :math:`(\mathbb{Z}/N\mathbb{Z})[X]`
        """
        raise DomainError("Cannot compute real roots for polynomials over integers modulo N")

    def complex_roots(self):
        r"""
        This method is not implemented for polynomials in
        :math:`(\mathbb{Z}/N\mathbb{Z})[X]`
        """
        raise DomainError("Cannot compute complex roots for polynomials over integers modulo N")
