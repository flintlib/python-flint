from flint.pyflint cimport global_random_state
from flint.flintlib.functions.fmpz cimport(
    fmpz_t,
    fmpz_one,
    fmpz_zero,
    fmpz_set,
    fmpz_init,
    fmpz_clear,
    fmpz_equal,
    fmpz_is_probabprime,
    fmpz_mul,
    fmpz_invmod,
    fmpz_sqrtmod,
    fmpz_divexact,
    fmpz_gcd,
    fmpz_is_one,
    fmpz_is_zero,
    fmpz_randm
)
from flint.flintlib.functions.fmpz cimport fmpz_mod as fmpz_type_mod
from flint.flintlib.functions.fmpz_mod cimport *

from flint.utils.typecheck cimport typecheck
from flint.flint_base.flint_base cimport flint_scalar
from flint.types.fmpz cimport(
    fmpz,
    any_as_fmpz,
    fmpz_get_intlong
)
cimport libc.stdlib

from flint.utils.flint_exceptions import DomainError


cdef class fmpz_mod_ctx:
    r"""
    Context object for creating :class:`~.fmpz_mod` initialised
    with a modulus :math:`N`.

        >>> fmpz_mod_ctx(2**127 - 1)
        fmpz_mod_ctx(170141183460469231731687303715884105727)

    """
    def __cinit__(self):
        cdef fmpz one = fmpz.__new__(fmpz)
        fmpz_one(one.val)
        fmpz_mod_ctx_init(self.val, one.val)
        fmpz_mod_discrete_log_pohlig_hellman_clear(self.L)
        self._is_prime = 0

    def __dealloc__(self):
        fmpz_mod_ctx_clear(self.val)
        fmpz_mod_discrete_log_pohlig_hellman_clear(self.L)

    def __init__(self, mod):
        # Ensure modulus is fmpz type
        if not typecheck(mod, fmpz):
            mod = any_as_fmpz(mod)
            if mod is NotImplemented:
                raise TypeError(
                    "Context modulus must be able to be cast to an `fmpz` type"
                )

        # Ensure modulus is positive
        if mod < 1:
            raise ValueError("Modulus is expected to be positive")

        # Set the modulus
        fmpz_mod_ctx_set_modulus(self.val, (<fmpz>mod).val)

        # Check whether the modulus is prime
        # TODO: should we use a stronger test?
        self._is_prime = fmpz_is_probabprime(self.val.n)

    def modulus(self):
        """
        Return the modulus from the context as an fmpz
        type

            >>> mod_ctx = fmpz_mod_ctx(2**127 - 1)
            >>> mod_ctx.modulus()
            170141183460469231731687303715884105727

        """
        n = fmpz()
        fmpz_set(n.val, (<fmpz_t>self.val.n))
        return n

    def is_prime(self):
        """
        Return whether the modulus is prime

            >>> fmpz_mod_ctx(2**127).is_prime()
            False
            >>> fmpz_mod_ctx(2**127 - 1).is_prime()
            True
        """
        return self._is_prime == 1

    def zero(self):
        """
        Return the zero element

            >>> F = fmpz_mod_ctx(163)
            >>> F.zero()
            fmpz_mod(0, 163)
        """
        cdef fmpz_mod res
        res = fmpz_mod.__new__(fmpz_mod)
        fmpz_zero(res.val)
        res.ctx = self

        return res

    def one(self):
        """
        Return the one element

            >>> F = fmpz_mod_ctx(163)
            >>> F.one()
            fmpz_mod(1, 163)
        """
        cdef fmpz_mod res
        res = fmpz_mod.__new__(fmpz_mod)
        fmpz_one(res.val)
        res.ctx = self

        return res

    def random_element(self):
        r"""
        Return a random element in :math:`\mathbb{Z}/N\mathbb{Z}`
        """
        cdef fmpz_mod res
        res = fmpz_mod.__new__(fmpz_mod)
        res.ctx = self

        fmpz_randm(res.val, global_random_state, self.val.n)

        return res

    cdef discrete_log_pohlig_hellman_run(self, fmpz_t x, fmpz_t y):
        # First, Ensure that L has performed precomputations This generates a
        # base which is a primitive root, and used as the base in
        # fmpz_mod_discrete_log_pohlig_hellman_run
        if not self._init_L:
            fmpz_mod_discrete_log_pohlig_hellman_precompute_prime(self.L, self.val.n)
            self._init_L = True

        fmpz_mod_discrete_log_pohlig_hellman_run(x, self.L, y)

    cdef set_any_as_fmpz_mod(self, fmpz_t val, obj):
        # Try and convert obj to fmpz
        if typecheck(obj, fmpz_mod):
            if self != (<fmpz_mod>obj).ctx:
                raise ValueError("moduli must match")
            fmpz_set(val, (<fmpz_mod>obj).val)
            return 0

        # Try and convert obj to fmpz
        if not typecheck(obj, fmpz):
            obj = any_as_fmpz(obj)
            if obj is NotImplemented:
                return NotImplemented

        fmpz_mod_set_fmpz(val, (<fmpz>obj).val, self.val)

        return 0

    cdef any_as_fmpz_mod(self, obj):
        # If `obj` is an `fmpz_mod`, just check moduli
        # match
        if typecheck(obj, fmpz_mod):
            if self != (<fmpz_mod>obj).ctx:
                raise ValueError("moduli must match")
            return obj

        # We have been able to cast `obj` to an `fmpz` so
        # we create a new `fmpz_mod` and set the val
        cdef fmpz_mod res
        res = fmpz_mod.__new__(fmpz_mod)
        check = self.set_any_as_fmpz_mod(res.val, obj)
        if check is NotImplemented:
            return NotImplemented
        res.ctx = self

        return res

    def __eq__(self, other):
        # TODO:
        # If we could cache contexts, then we would ensure that only
        # the a is b check is needed for equality.

        # Most often, we expect both `fmpz_mod` to be pointing to the
        # same ctx, so this seems the fastest way to check
        if self is other:
            return True

        # If they're not the same object in memory, they may have the
        # same modulus, which is good enough
        if typecheck(other, fmpz_mod_ctx):
            return fmpz_equal(self.val.n, (<fmpz_mod_ctx>other).val.n) == 1
        return False

    def __hash__(self):
        return hash(self.modulus())

    def __str__(self):
        return f"Context for fmpz_mod with modulus: {self.modulus()}"

    def __repr__(self):
        return f"fmpz_mod_ctx({self.modulus()})"

    def __call__(self, val):
        return fmpz_mod(val, self)

cdef class fmpz_mod(flint_scalar):
    """
    The *fmpz_mod* type represents integer modulo an
    arbitrary-size modulus. For wordsize modulus, see
    :class:`~.nmod`.

    An *fmpz_mod* element is constructed from an :class:`~.fmpz_mod_ctx`
    either by passing it as an argument to the type, or
    by directly calling the context

        >>> fmpz_mod(-1, fmpz_mod_ctx(2**127 - 1))
        fmpz_mod(170141183460469231731687303715884105726, 170141183460469231731687303715884105727)
        >>> ZmodN = fmpz_mod_ctx(2**127 - 1)
        >>> ZmodN(-2)
        fmpz_mod(170141183460469231731687303715884105725, 170141183460469231731687303715884105727)

    """

    def __cinit__(self):
        fmpz_init(self.val)
        fmpz_init(self.x_g)

    def __dealloc__(self):
        fmpz_clear(self.val)
        fmpz_clear(self.x_g)

    def __init__(self, val, ctx):
        if not typecheck(ctx, fmpz_mod_ctx):
            raise TypeError
        self.ctx = ctx
        check = self.ctx.set_any_as_fmpz_mod(self.val, val)
        if check is NotImplemented:
            raise NotImplementedError(f"Cannot convert {val} to type 'fmpz_mod'")

    def is_zero(self):
        """
        Return whether an element is equal to zero

            >>> mod_ctx = fmpz_mod_ctx(163)
            >>> mod_ctx(0).is_zero()
            True
            >>> mod_ctx(1).is_zero()
            False
        """
        return self == 0

    def is_one(self):
        """
        Return whether an element is equal to one

            >>> mod_ctx = fmpz_mod_ctx(163)
            >>> mod_ctx(0).is_one()
            False
            >>> mod_ctx(1).is_one()
            True
        """

        cdef bint res
        res = fmpz_mod_is_one(self.val, self.ctx.val)
        return res == 1

    def is_unit(self):
        """
        Returns whether the element is invertible modulo `N`

            >>> from flint import *
            >>> F = fmpz_mod_ctx(10)
            >>> F(3).is_unit()
            True
            >>> F(2).is_unit()
            False
        """
        cdef fmpz_t g
        fmpz_init(g)
        fmpz_gcd(g, self.val, self.ctx.val.n)
        return 1 == fmpz_is_one(g)

    def inverse(self, check=True):
        r"""
        Computes :math:`a^{-1} \pmod N`

        When check=False, the solutions is assumed to exist and Flint will abort on
        failure.

            >>> mod_ctx = fmpz_mod_ctx(163)
            >>> mod_ctx(2).inverse()
            fmpz_mod(82, 163)
            >>> mod_ctx(2).inverse(check=False)
            fmpz_mod(82, 163)
        """
        cdef fmpz_mod res
        res = fmpz_mod.__new__(fmpz_mod)
        res.ctx = self.ctx

        if check is False:
            fmpz_mod_inv(res.val, self.val, self.ctx.val)
            return res

        cdef bint r
        cdef fmpz one = fmpz.__new__(fmpz)
        fmpz_one(one.val)

        r = fmpz_mod_divides(
            res.val, one.val, self.val, self.ctx.val
        )
        if r == 0:
            raise ZeroDivisionError(
                f"{self} is not invertible modulo {self.ctx.modulus()}"
            )

        return res

    def discrete_log(self, a):
        r"""
        Solve the discrete logarithm problem, using `self = g` as a base.
        Assumes a solution, :math:`a = g^x \pmod p` exists.

        NOTE: Requires that the context modulus is prime.

            >>> F = fmpz_mod_ctx(163)
            >>> g = F(2)
            >>> x = 123
            >>> a = g**123
            >>> g.discrete_log(a)
            123
        """
        # Ensure that the modulus is prime
        if not self.ctx.is_prime():
            raise NotImplementedError("algorithm assumes modulus is prime")

        # Then check the type of the input
        a = self.ctx.any_as_fmpz_mod(a)
        if a is NotImplemented:
            raise TypeError(f"Cannot solve the discrete log with {type(a)} as input")

        # Solve the discrete log for the chosen base and target
        # g = y^x_g and  a = y^x_a
        # We want to find x such that a = g^x =>
        # (y^x_a) = (y^x_g)^x => x = (x_a / x_g) mod (p-1)

        # For repeated calls to discrete_log, it's more efficient to
        # store x_g rather than keep computing it
        if fmpz_is_zero(self.x_g):
            self.ctx.discrete_log_pohlig_hellman_run(self.x_g, self.val)

        # Then we need to compute x_a which will be different for each call
        cdef fmpz_t x_a
        fmpz_init(x_a)
        self.ctx.discrete_log_pohlig_hellman_run(x_a, (<fmpz_mod>a).val)

        # If g is not a primitive root, then x_g and pm1 will share
        # a common factor. We can use this to compute the order of
        # g.
        cdef fmpz_t g, g_order, x_g
        fmpz_init(g)
        fmpz_init(g_order)
        fmpz_init(x_g)

        fmpz_gcd(g, self.x_g, self.ctx.L.pm1)
        if not fmpz_is_one(g):
            fmpz_divexact(x_g, self.x_g, g)
            fmpz_divexact(x_a, x_a, g)
            fmpz_divexact(g_order, self.ctx.L.pm1, g)
        else:
            fmpz_set(g_order, self.ctx.L.pm1)
            fmpz_set(x_g, self.x_g)

        # Finally, compute output exponent by computing
        # (x_a / x_g) mod g_order
        cdef fmpz x = fmpz.__new__(fmpz)
        fmpz_invmod(x.val, x_g, g_order)
        fmpz_mul(x.val, x.val, x_a)
        fmpz_type_mod(x.val, x.val, g_order)

        return x

    def __richcmp__(self, other, int op):
        cdef bint res
        if op != 2 and op != 3:
            raise TypeError("fmpz_mod cannot be ordered")

        other = self.ctx.any_as_fmpz_mod(other)
        if other is NotImplemented:
            return NotImplemented

        res = fmpz_equal(self.val, (<fmpz_mod>other).val) and \
            (self.ctx == (<fmpz_mod>other).ctx)
        return res if op == 2 else not res

    def __bool__(self):
        return not self.is_zero()

    def repr(self):
        return "fmpz_mod({}, {})".format(
            fmpz_get_intlong(self.val),
            self.ctx.modulus()
        )

    def __repr__(self):
        return self.repr()

    def __hash__(self):
        return hash((int(self)))

    def __int__(self):
        return fmpz_get_intlong(self.val)

    def str(self):
        return str(int(self))

    def __str__(self):
        return self.str()

    # ---------------- #
    #    Arithmetic    #
    # ---------------- #

    def __pos__(self):
        return self

    def __neg__(self):
        cdef fmpz_mod res
        res = fmpz_mod.__new__(fmpz_mod)
        res.ctx = self.ctx
        fmpz_mod_neg(res.val, self.val, self.ctx.val)
        return res

    def __add__(self, other):
        other = self.ctx.any_as_fmpz_mod(other)
        if other is NotImplemented:
            return NotImplemented

        cdef fmpz_mod res
        res = fmpz_mod.__new__(fmpz_mod)
        res.ctx = self.ctx
        fmpz_mod_add(
            res.val, self.val, (<fmpz_mod>other).val, self.ctx.val
        )
        return res

    def __radd__(self, other):
        return self.__add__(other)

    @staticmethod
    def _sub_(left, right):
        cdef fmpz_mod res
        res = fmpz_mod.__new__(fmpz_mod)

        # Case when left and right are already fmpz_mod
        if typecheck(left, fmpz_mod) and typecheck(right, fmpz_mod):
            if not (<fmpz_mod>left).ctx == (<fmpz_mod>right).ctx:
                raise ValueError("moduli must match")

        # Case when right is not fmpz_mod, try to convert to fmpz
        elif typecheck(left, fmpz_mod):
            right = (<fmpz_mod>left).ctx.any_as_fmpz_mod(right)
            if right is NotImplemented:
                return NotImplemented

        # Case when left is not fmpz_mod, try to convert to fmpz
        else:
            left = (<fmpz_mod>right).ctx.any_as_fmpz_mod(left)
            if left is NotImplemented:
                return NotImplemented

        res.ctx = (<fmpz_mod>left).ctx
        fmpz_mod_sub(
                res.val, (<fmpz_mod>left).val, (<fmpz_mod>right).val, res.ctx.val
        )
        return res

    def __sub__(s, t):
        return fmpz_mod._sub_(s, t)

    def __rsub__(s, t):
        return fmpz_mod._sub_(t, s)

    def __mul__(self, other):
        other = self.ctx.any_as_fmpz_mod(other)
        if other is NotImplemented:
            return NotImplemented

        cdef fmpz_mod res
        res = fmpz_mod.__new__(fmpz_mod)
        res.ctx = self.ctx

        fmpz_mod_mul(
            res.val, self.val, (<fmpz_mod>other).val, self.ctx.val
        )
        return res

    def __rmul__(self, other):
        return self.__mul__(other)

    @staticmethod
    def _div_(left, right):
        cdef bint check
        cdef fmpz_mod res
        res = fmpz_mod.__new__(fmpz_mod)

        # Case when left and right are already fmpz_mod
        if typecheck(left, fmpz_mod) and typecheck(right, fmpz_mod):
            if not (<fmpz_mod>left).ctx == (<fmpz_mod>right).ctx:
                raise ValueError("moduli must match")

        # Case when right is not fmpz_mod, try to convert to fmpz
        elif typecheck(left, fmpz_mod):
            right = (<fmpz_mod>left).ctx.any_as_fmpz_mod(right)
            if right is NotImplemented:
                return NotImplemented

        # Case when left is not fmpz_mod, try to convert to fmpz
        else:
            left = (<fmpz_mod>right).ctx.any_as_fmpz_mod(left)
            if left is NotImplemented:
                return NotImplemented

        res.ctx = (<fmpz_mod>left).ctx
        check = fmpz_mod_divides(
            res.val, (<fmpz_mod>left).val, (<fmpz_mod>right).val, res.ctx.val
        )
        if check == 0:
            raise ZeroDivisionError(
                f"{right} is not invertible modulo {res.ctx.modulus()}"
            )

        return res

    def __truediv__(s, t):
        return fmpz_mod._div_(s, t)

    def __rtruediv__(s, t):
        return fmpz_mod._div_(t, s)

    def __floordiv__(self, other):
        return NotImplemented

    def __invert__(self):
        return self.inverse()

    def __pow__(self, e):
        cdef bint check
        cdef fmpz_mod res
        res = fmpz_mod.__new__(fmpz_mod)
        res.ctx = self.ctx

        # Attempt to convert exponent to fmpz
        e = any_as_fmpz(e)
        if e is NotImplemented:
            raise NotImplementedError

        check = fmpz_mod_pow_fmpz(
            res.val, self.val, (<fmpz>e).val, self.ctx.val
        )

        if check == 0:
            raise ZeroDivisionError(
                f"{self} is not invertible modulo {self.ctx.modulus()}"
            )

        return res

    def sqrt(self):
        """
        Return the square root of this ``fmpz_mod`` or raise an exception.

            >>> ctx = fmpz_mod_ctx(13)
            >>> s = ctx(10).sqrt()
            >>> s
            fmpz_mod(6, 13)
            >>> s * s
            fmpz_mod(10, 13)
            >>> ctx(11).sqrt()
            Traceback (most recent call last):
                ...
            flint.utils.flint_exceptions.DomainError: no square root exists for 11 mod 13

        The modulus must be prime.

        """
        cdef fmpz_mod v

        v = fmpz_mod.__new__(fmpz_mod)
        v.ctx = self.ctx

        if fmpz_is_zero(self.val):
            return v

        if not fmpz_sqrtmod(v.val, self.val, self.ctx.val.n):
            raise DomainError("no square root exists for {} mod {}".format(self, self.ctx.modulus()))

        return v
