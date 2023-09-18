from flint.flintlib.fmpz cimport (
    fmpz_t,
    fmpz_one,
    fmpz_set,
    fmpz_init,
    fmpz_clear,
    fmpz_equal,
    fmpz_is_probabprime,
    fmpz_mul,
    fmpz_invmod,
    fmpz_divexact,
    fmpz_gcd,
    fmpz_is_one
)
from flint.flintlib.fmpz cimport fmpz_mod as fmpz_type_mod

from flint.flintlib.fmpz_mod cimport *

from flint.utils.typecheck cimport typecheck
from flint.flint_base.flint_base cimport flint_scalar
from flint.types.fmpz cimport (
    fmpz,
    any_as_fmpz,
    fmpz_get_intlong
)


cdef class fmpz_mod_ctx:
    r"""
    Context object for creating :class:`~.fmpz_mod` initalised 
    with a modulus :math:`N`.

        >>> fmpz_mod_ctx(2**127 - 1)
        fmpz_mod_ctx(170141183460469231731687303715884105727)

    """
    def __cinit__(self):
        cdef fmpz one = fmpz.__new__(fmpz)
        fmpz_one(one.val)
        fmpz_mod_ctx_init(self.val, one.val)

    def __dealloc__(self):
        fmpz_mod_ctx_clear(self.val)

    def __init__(self, mod):
        # Ensure modulus is fmpz type
        if not typecheck(mod, fmpz):
            mod = any_as_fmpz(mod)
            if mod is NotImplemented:
                raise TypeError("Context modulus must be able to be case to an `fmpz` type")

        # Ensure modulus is positive
        if mod < 1:
            raise ValueError("Modulus is expected to be positive")

        # Set the modulus
        fmpz_mod_ctx_set_modulus(self.val, (<fmpz>mod).val)
    
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

    cdef any_as_fmpz_mod(self, obj):
        # If `obj` is an `fmpz_mod`, just check moduli
        # match
        # TODO: we could allow conversion from one modulus to another?
        if typecheck(obj, fmpz_mod):
            if self != (<fmpz_mod>obj).ctx:
                raise ValueError("moduli must match")
            return obj
        
        # Try and convert obj to fmpz
        if not typecheck(obj, fmpz):
            obj = any_as_fmpz(obj)
            if obj is NotImplemented:
                return NotImplemented
        
        # We have been able to cast `obj` to an `fmpz` so 
        # we create a new `fmpz_mod` and set the val
        cdef fmpz_mod res
        res = fmpz_mod.__new__(fmpz_mod)
        res.ctx = self
        fmpz_mod_set_fmpz(res.val, (<fmpz>obj).val, self.val)
        
        return res

    def __eq__(self, other):
        if typecheck(other, fmpz_mod_ctx):
            return fmpz_equal(self.val.n, (<fmpz_mod_ctx>other).val.n)
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

    def __dealloc__(self):
        fmpz_clear(self.val)

    def __init__(self, val, ctx):
        if not typecheck(ctx, fmpz_mod_ctx):
            raise TypeError
        self.ctx = ctx

        # When the input is also an fmpz_mod we just need
        # moduli to match
        if typecheck(val, fmpz_mod):
            if self.ctx != (<fmpz_mod>val).ctx:
                raise ValueError("moduli must match")
            # fmpz_mod_set_fmpz(self.val, (<fmpz_mod>val).val, self.ctx.val)
            fmpz_set(self.val, (<fmpz_mod>val).val)
            return

        # For all other cases, the easiest is to first convert to
        # fmpz type and set this way
        if not typecheck(val, fmpz):
            val = any_as_fmpz(val)
            if val is NotImplemented:
                raise NotImplementedError
        fmpz_mod_set_fmpz(self.val, (<fmpz>val).val, self.ctx.val)

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
            >>> mod_ctx(1).is_zero()
            True
        """

        cdef bint res
        res = fmpz_mod_is_one(self.val, self.ctx.val)
        return res == 1

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
            raise ZeroDivisionError(f"{self} is not invertible modulo {self.ctx.modulus()}")

        return res

    def discrete_log(self, a, check=False):
        """
        Solve the discrete logarithm problem, using `self = g` as a base.
        Assumes a solution, :math:`a = g^x \pmod p` exists.
        
        NOTE: Requires that the context modulus is prime.

        TODO: This could instead be initalised as a class from a 
        given base and the precomputations could be stored to allow 
        faster computations for many discrete logs with the same base. 

            >>> F = fmpz_mod_ctx(163)
            >>> g = F(2)
            >>> x = 123
            >>> a = g**123
            >>> g.discrete_log(a)
            123
        """
        cdef fmpz_mod_discrete_log_pohlig_hellman_t L
        cdef bint is_prime

        print(f"[DEBUGGING!]")
        print("\t[DEBUG]: Starting discrete log")

        # Ensure that the modulus is prime
        if check:
            is_prime = fmpz_is_probabprime(self.ctx.val.n)
            if not is_prime:
                raise ValueError("modulus must be prime")

        print("\t[DEBUG]: Checked Primality")
    
        # Then check the type of the input
        if typecheck(a, fmpz_mod):
            if self.ctx != (<fmpz_mod>a).ctx:
                raise ValueError("moduli must match")
        else:
            a = self.ctx.any_as_fmpz_mod(a)
            if a is NotImplemented:
                raise TypeError

        print("\t[DEBUG]: Converted types")
        
        # Initalise the dlog data, all discrete logs are solved with an internally
        # chosen base `y`
        fmpz_mod_discrete_log_pohlig_hellman_init(L)

        print("\t[DEBUG]: Initialised L")

        fmpz_mod_discrete_log_pohlig_hellman_precompute_prime(L, self.ctx.val.n)

        print("\t[DEBUG]: Precomputed prime for L")

        # Solve the discrete log for the chosen base and target
        # g = y^x_g and  a = y^x_a
        # We want to find x such that a = g^x =>
        # (y^x_a) = (y^x_g)^x => x = (x_a / x_g) mod (p-1) 
        cdef fmpz_t x_a
        cdef fmpz_t x_g
        fmpz_init(x_a)
        fmpz_init(x_g)

        # TODO: should this value be stored for efficiency?
        fmpz_mod_discrete_log_pohlig_hellman_run(x_g, L, self.val)
        print("\t[DEBUG]: Solved dlog for base")
        fmpz_mod_discrete_log_pohlig_hellman_run(x_a, L, (<fmpz_mod>a).val)
        print("\t[DEBUG]: Solved dlog for input")

        # If g is not a primative root, then x_g and pm1 will share
        # a common factor. We can use this to compute the order of 
        # g.
        cdef fmpz_t g, g_order
        fmpz_init(g)
        fmpz_init(g_order)

        fmpz_gcd(g, x_g, L.pm1)
        if not fmpz_is_one(g):
            fmpz_divexact(x_g, x_g, g)
            fmpz_divexact(x_a, x_a, g)
            fmpz_divexact(g_order, L.pm1, g)
        else:
            fmpz_set(g_order, L.pm1)
        print("\t[DEBUG]: Fixed order")

        # Finally, compute output exponent
        cdef fmpz x = fmpz.__new__(fmpz)

        # Compute (x_a / x_g) mod g_order
        fmpz_invmod(x.val, x_g, g_order)
        print("\t[DEBUG]: Computed Inverse")

        fmpz_mul(x.val, x.val, x_a)
        print("\t[DEBUG]: Computed multiplication")

        fmpz_type_mod(x.val, x.val, g_order)
        print("\t[DEBUG]: Computed modular reduction")

        # Clear the dlog struct
        fmpz_mod_discrete_log_pohlig_hellman_clear(L)
        print("\t[DEBUG]: Cleared Struct")

        
        return x

    def __richcmp__(self, other, int op):
        cdef bint res
        if op != 2 and op != 3:
            raise TypeError("fmpz_mod cannot be ordered")

        if not typecheck(other, fmpz_mod):
            other = self.ctx.any_as_fmpz_mod(other)

        if typecheck(self, fmpz_mod) and typecheck(other, fmpz_mod):
            res = fmpz_equal(self.val, (<fmpz_mod>other).val) and \
                  (self.ctx == (<fmpz_mod>other).ctx)
            if op == 2:
                return res
            else:
                return not res
        else:
            return NotImplemented

    def __bool__(self):
        return not self.is_zero()

    def __repr__(self):
        return "fmpz_mod({}, {})".format(
            fmpz_get_intlong(self.val),
            self.ctx.modulus()
        )

    def __hash__(self):
        return  hash((int(self)))

    def __int__(self):
        return fmpz_get_intlong(self.val)

    def __str__(self):
        return str(int(self))

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
            raise ZeroDivisionError(f"{right} is not invertible modulo {res.ctx.modulus()}")

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
            raise ZeroDivisionError(f"{self} is not invertible modulo {self.ctx.modulus()}")

        return res
