from cpython.list cimport PyList_GET_SIZE
from flint.flint_base.flint_base cimport flint_poly
from flint.utils.typecheck cimport typecheck
from flint.types.fmpz cimport fmpz, any_as_fmpz
from flint.types.fmpz_poly cimport any_as_fmpz_poly
from flint.types.fmpz_poly cimport fmpz_poly
from flint.types.nmod cimport nmod, nmod_ctx

from flint.flintlib.nmod_vec cimport *
from flint.flintlib.nmod_poly cimport *
from flint.flintlib.nmod_poly_factor cimport *
from flint.flintlib.fmpz_poly cimport fmpz_poly_get_nmod_poly
from flint.flintlib.ulong_extras cimport n_gcdinv, n_is_prime

from flint.utils.flint_exceptions import DomainError


cdef dict _nmod_poly_ctx_cache = {}


cdef class nmod_poly_ctx:
    """
    Context object for creating :class:`~.nmod_poly` initalised 
    with modulus :math:`N`.

        >>> nmod_poly_ctx.new(17)
        nmod_poly_ctx(17)

    """
    def __init__(self, *args, **kwargs):
        raise TypeError("cannot create nmod_poly_ctx directly: use nmod_poly_ctx.new()")

    @staticmethod
    def new(mod):
        """Get an ``nmod_poly`` context with modulus ``mod``."""
        return nmod_poly_ctx._get_ctx(mod)

    @staticmethod
    cdef any_as_nmod_poly_ctx(obj):
        """Convert an ``nmod_poly_ctx`` or ``int`` to an ``nmod_poly_ctx``."""
        if typecheck(obj, nmod_poly_ctx):
            return obj
        if typecheck(obj, int):
            return nmod_poly_ctx._get_ctx(obj)
        return NotImplemented

    @staticmethod
    cdef nmod_poly_ctx _get_ctx(int mod):
        """Retrieve an nmod_poly context from the cache or create a new one."""
        ctx = _nmod_poly_ctx_cache.get(mod)
        if ctx is None:
            ctx = _nmod_poly_ctx_cache.setdefault(mod, nmod_poly_ctx._new_ctx(mod))
        return ctx

    @staticmethod
    cdef nmod_poly_ctx _new_ctx(ulong mod):
        """Create a new nmod_poly context."""
        cdef nmod_ctx scalar_ctx
        cdef nmod_poly_ctx ctx
        scalar_ctx = nmod_ctx.new(mod)

        ctx = nmod_poly_ctx.__new__(nmod_poly_ctx)
        ctx.mod = scalar_ctx.mod
        ctx._is_prime = scalar_ctx._is_prime
        ctx.scalar_ctx = scalar_ctx

        return ctx

    cdef int any_as_nmod(self, mp_limb_t * val, obj) except -1:
        return self.scalar_ctx.any_as_nmod(val, obj)

    cdef any_as_nmod_poly(self, obj):
        cdef nmod_poly r
        cdef mp_limb_t v
        # XXX: should check that modulus is the same here, and not all over the place
        if typecheck(obj, nmod_poly):
            return obj
        if self.any_as_nmod(&v, obj):
            r = self.new_nmod_poly()
            nmod_poly_set_coeff_ui(r.val, 0, v)
            return r
        x = any_as_fmpz_poly(obj)
        if x is not NotImplemented:
            r = self.new_nmod_poly()
            fmpz_poly_get_nmod_poly(r.val, (<fmpz_poly>x).val)
            return r
        return NotImplemented

    cdef nmod new_nmod(self):
        return self.scalar_ctx.new_nmod()

    cdef nmod_poly new_nmod_poly(self):
        cdef nmod_poly p = nmod_poly.__new__(nmod_poly)
        nmod_poly_init_preinv(p.val, self.mod.n, self.mod.ninv)
        p.ctx = self
        return p

    cdef nmod_poly_set_list(self, nmod_poly_t poly, list val):
        cdef long i, n
        cdef mp_limb_t v
        n = PyList_GET_SIZE(val)
        nmod_poly_fit_length(poly, n)
        for i from 0 <= i < n:
            c = val[i]
            if self.any_as_nmod(&v, val[i]):
                nmod_poly_set_coeff_ui(poly, i, v)
            else:
                raise TypeError("unsupported coefficient in list")

    def modulus(self):
        """Get the modulus of the context.

        >>> ctx = nmod_poly_ctx.new(17)
        >>> ctx.modulus()
        17

        """
        return fmpz(self.mod.n)

    def is_prime(self):
        """Check if the modulus is prime.

        >>> ctx = nmod_poly_ctx.new(17)
        >>> ctx.is_prime()
        True

        """
        return self._is_prime

    def zero(self):
        """Return the zero ``nmod_poly``.

        >>> ctx = nmod_poly_ctx.new(17)
        >>> ctx.zero()
        0

        """
        cdef nmod_poly r = self.new_nmod_poly()
        nmod_poly_zero(r.val)
        return r

    def one(self):
        """Return the one ``nmod_poly``.

        >>> ctx = nmod_poly_ctx.new(17)
        >>> ctx.one()
        1

        """
        cdef nmod_poly r = self.new_nmod_poly()
        nmod_poly_set_coeff_ui(r.val, 0, 1)
        return r

    def __str__(self):
        return f"Context for nmod_poly with modulus: {self.mod.n}"

    def __repr__(self):
        return f"nmod_poly_ctx({self.mod.n})"

    def __call__(self, arg):
        """Create an ``nmod_poly``.

        >>> ctx = nmod_poly_ctx.new(17)
        >>> ctx(10)
        10
        >>> ctx([1,2,3])
        3*x^2 + 2*x + 1

        """
        return nmod_poly(arg, self)


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

    def __init__(self, val=None, mod=0):
        cdef ulong m2
        cdef mp_limb_t v
        cdef nmod_poly_ctx ctx

        if typecheck(val, nmod_poly):
            m2 = nmod_poly_modulus((<nmod_poly>val).val)
            if m2 != mod:
                raise ValueError("different moduli!")
            nmod_poly_init(self.val, m2)
            nmod_poly_set(self.val, (<nmod_poly>val).val)
            self.ctx = (<nmod_poly>val).ctx
        else:
            if mod == 0:
                raise ValueError("a nonzero modulus is required")
            ctx = nmod_poly_ctx.any_as_nmod_poly_ctx(mod)
            if ctx is NotImplemented:
                raise TypeError("cannot create nmod_poly_ctx from input of type %s", type(mod))

            self.ctx = ctx
            nmod_poly_init(self.val, ctx.mod.n)
            if typecheck(val, fmpz_poly):
                fmpz_poly_get_nmod_poly(self.val, (<fmpz_poly>val).val)
            elif typecheck(val, list):
                ctx.nmod_poly_set_list(self.val, val)
            elif ctx.any_as_nmod(&v, val):
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
        if self.ctx.any_as_nmod(&v, x):
            nmod_poly_set_coeff_ui(self.val, i, v)
        else:
            raise TypeError("cannot set element of type %s" % type(x))

    def __bool__(self):
        return not nmod_poly_is_zero(self.val)

    def is_zero(self):
        return <bint>nmod_poly_is_zero(self.val)

    def is_one(self):
        return <bint>nmod_poly_is_one(self.val)

    def is_gen(self):
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

        res = self.ctx.new_nmod_poly()
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

        c = self.ctx.new_nmod()
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

        cdef nmod_poly res = self.ctx.new_nmod_poly()
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
        other = self.ctx.any_as_nmod_poly(other)
        if other is NotImplemented:
            raise TypeError("cannot convert input to nmod_poly")
        res = self.ctx.new_nmod_poly()
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
        g = self.ctx.any_as_nmod_poly(other)
        if g is NotImplemented:
            raise TypeError(f"cannot convert other = {other} to nmod_poly")

        h = self.ctx.any_as_nmod_poly(modulus)
        if h is NotImplemented:
            raise TypeError(f"cannot convert modulus = {modulus} to nmod_poly")

        if modulus.is_zero():
            raise ZeroDivisionError("cannot reduce modulo zero")

        res = self.ctx.new_nmod_poly()
        nmod_poly_compose_mod(res.val, self.val, (<nmod_poly>other).val, (<nmod_poly>modulus).val) 
        return res 

    def __call__(self, other):
        cdef nmod_poly r
        cdef mp_limb_t c
        if self.ctx.any_as_nmod(&c, other):
            v = nmod(0, self.modulus())
            (<nmod>v).val = nmod_poly_evaluate_nmod(self.val, c)
            return v
        t = self.ctx.any_as_nmod_poly(other)
        if t is not NotImplemented:
            r = self.ctx.new_nmod_poly()
            nmod_poly_compose(r.val, self.val, (<nmod_poly>t).val)
            return r
        raise TypeError("cannot call nmod_poly with input of type %s", type(other))

    def derivative(self):
        cdef nmod_poly res = self.ctx.new_nmod_poly()
        nmod_poly_derivative(res.val, self.val)
        return res

    def integral(self):
        cdef nmod_poly res = self.ctx.new_nmod_poly()
        nmod_poly_integral(res.val, self.val)
        return res

    def __pos__(self):
        return self

    def __neg__(self):
        cdef nmod_poly r = self.ctx.new_nmod_poly()
        nmod_poly_neg(r.val, self.val)
        return r

    def _add_(s, t):
        cdef nmod_poly r
        t = s.ctx.any_as_nmod_poly(t)
        if t is NotImplemented:
            return t
        if (<nmod_poly>s).val.mod.n != (<nmod_poly>t).val.mod.n:
            raise ValueError("cannot add nmod_polys with different moduli")
        r = s.ctx.new_nmod_poly()
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
        r = s.ctx.new_nmod_poly()
        nmod_poly_sub(r.val, (<nmod_poly>s).val, (<nmod_poly>t).val)
        return r

    def __sub__(s, t):
        t = s.ctx.any_as_nmod_poly(t)
        if t is NotImplemented:
            return t
        return s._sub_(t)

    def __rsub__(s, t):
        t = s.ctx.any_as_nmod_poly(t)
        if t is NotImplemented:
            return t
        return t._sub_(s)

    def _mul_(s, t):
        cdef nmod_poly r
        t = s.ctx.any_as_nmod_poly(t)
        if t is NotImplemented:
            return t
        if (<nmod_poly>s).val.mod.n != (<nmod_poly>t).val.mod.n:
            raise ValueError("cannot multiply nmod_polys with different moduli")
        r = s.ctx.new_nmod_poly()
        nmod_poly_mul(r.val, (<nmod_poly>s).val, (<nmod_poly>t).val)
        return r

    def __mul__(s, t):
        return s._mul_(t)

    def __rmul__(s, t):
        return s._mul_(t)

    def __truediv__(s, t):
        t = s.ctx.any_as_nmod_poly(t)
        if t is NotImplemented:
            return t
        res, r = s._divmod_(t)
        if not nmod_poly_is_zero((<nmod_poly>r).val):
            raise DomainError("nmod_poly inexact division")
        return res

    def __rtruediv__(s, t):
        t = s.ctx.any_as_nmod_poly(t)
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
        if not s.ctx._is_prime:
            raise DomainError("nmod_poly divmod: modulus {self.ctx.mod.n} is not prime")

        r = s.ctx.new_nmod_poly()
        nmod_poly_div(r.val, (<nmod_poly>s).val, (<nmod_poly>t).val)
        return r

    def __floordiv__(s, t):
        t = s.ctx.any_as_nmod_poly(t)
        if t is NotImplemented:
            return t
        return s._floordiv_(t)

    def __rfloordiv__(s, t):
        t = s.ctx.any_as_nmod_poly(t)
        if t is NotImplemented:
            return t
        return t._floordiv_(s)

    def _divmod_(s, t):
        cdef nmod_poly P, Q

        if (<nmod_poly>s).val.mod.n != (<nmod_poly>t).val.mod.n:
            raise ValueError("cannot divide nmod_polys with different moduli")
        if nmod_poly_is_zero((<nmod_poly>t).val):
            raise ZeroDivisionError("polynomial division by zero")
        if not s.ctx._is_prime:
            raise DomainError("nmod_poly divmod: modulus {self.ctx.mod.n} is not prime")

        P = s.ctx.new_nmod_poly()
        Q = s.ctx.new_nmod_poly()
        nmod_poly_divrem(P.val, Q.val, (<nmod_poly>s).val, (<nmod_poly>t).val)
        return P, Q

    def __divmod__(s, t):
        t = s.ctx.any_as_nmod_poly(t)
        if t is NotImplemented:
            return t
        return s._divmod_(t)

    def __rdivmod__(s, t):
        t = s.ctx.any_as_nmod_poly(t)
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
            raise ValueError("negative exponent")
        res = self.ctx.new_nmod_poly()
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

        modulus = self.ctx.any_as_nmod_poly(modulus)
        if modulus is NotImplemented:
            raise TypeError("cannot convert input to nmod_poly")

        # Output polynomial
        res = self.ctx.new_nmod_poly()
        
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
            mod_rev_inv = self.ctx.any_as_nmod_poly(mod_rev_inv)
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

        The modulus must be prime.
        """
        cdef nmod_poly res

        other = self.ctx.any_as_nmod_poly(other)
        if other is NotImplemented:
            raise TypeError("cannot convert input to nmod_poly")
        if self.val.mod.n != (<nmod_poly>other).val.mod.n:
            raise ValueError("moduli must be the same")
        if not self.ctx._is_prime:
            raise DomainError("nmod_poly gcd: modulus {self.ctx.mod.n} is not prime")

        res = self.ctx.new_nmod_poly()
        nmod_poly_gcd(res.val, self.val, (<nmod_poly>other).val)
        return res

    def xgcd(self, other):
        r"""
        Computes the extended gcd of self and other: (`G`, `S`, `T`)
        where `G` is the ``gcd(self, other)`` and `S`, `T` are such that:

        :math:`G = \textrm{self}*S +  \textrm{other}*T`

            >>> f = nmod_poly([143, 19, 37, 138, 102, 127, 95], 163)
            >>> g = nmod_poly([139, 9, 35, 154, 87, 120, 24], 163)
            >>> f.xgcd(g)
            (x^3 + 128*x^2 + 123*x + 91, 17*x^2 + 49*x + 104, 21*x^2 + 5*x + 25)

        The modulus must be prime.
        """
        cdef nmod_poly res1, res2, res3

        other = self.ctx.any_as_nmod_poly(other)
        if other is NotImplemented:
            raise TypeError("cannot convert input to fmpq_poly")
        if self.val.mod.n != (<nmod_poly>other).val.mod.n:
            raise ValueError("moduli must be the same")
        if not self.ctx._is_prime:
            raise DomainError("nmod_poly xgcd: modulus {self.ctx.mod.n} is not prime")

        res1 = self.ctx.new_nmod_poly()
        res2 = self.ctx.new_nmod_poly()
        res3 = self.ctx.new_nmod_poly()

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

        The modulus must be prime.
        """
        if algorithm is None:
            algorithm = 'irreducible'
        elif algorithm not in ('berlekamp', 'cantor-zassenhaus'):
            raise ValueError(f"unknown factorization algorithm: {algorithm}")
        if not self.ctx._is_prime:
            raise DomainError(f"nmod_poly factor: modulus {self.ctx.mod.n} is not prime")
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
        if not self.ctx._is_prime:
            raise DomainError(f"nmod_poly factor_squarefree: modulus {self.ctx.mod.n} is not prime")
        return self._factor('squarefree')

    def _factor(self, factor_type):
        cdef nmod_poly_factor_t fac
        cdef nmod_poly u
        cdef nmod c
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
            u = self.ctx.new_nmod_poly()
            nmod_poly_set(u.val, &fac.p[i])
            exp = fac.exp[i]
            res[i] = (u, exp)

        c = self.ctx.new_nmod()
        c.val = lead

        nmod_poly_factor_clear(fac)

        return c, res

    def sqrt(nmod_poly self):
        """Return exact square root or ``None``. """
        cdef nmod_poly

        if not self.ctx._is_prime:
            raise DomainError(f"nmod_poly sqrt: modulus {self.ctx.mod.n} is not prime")

        res = self.ctx.new_nmod_poly()

        if not nmod_poly_sqrt(res.val, self.val):
            raise DomainError(f"Cannot compute square root of {self}")

        return res

    def deflation(self):
        cdef nmod_poly v
        cdef ulong n
        if nmod_poly_is_zero(self.val):
            return self, 1
        n = nmod_poly_deflation(self.val)
        if n == 1:
            return self, int(n)
        else:
            v = self.ctx.new_nmod_poly()
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
