from flint.flintlib.fmpz cimport (
    fmpz_t,
    fmpz_one,
    fmpz_set,
    fmpz_init,
    fmpz_clear,
    fmpz_equal
)
from flint.flintlib.fmpz_mod cimport *

from flint.utils.typecheck cimport typecheck
from flint.flint_base.flint_base cimport flint_scalar
from flint.types.fmpz cimport (
    fmpz,
    any_as_fmpz,
    fmpz_get_intlong
)


cdef class fmpz_mod_ctx:
    """
    Context object for *fmpz_mod* initalised with a modulus `n`

        >>> fmpz_mod_ctx(2**127 - 1)
        fmpz_mod_ctx(170141183460469231731687303715884105727)

    """
    def __cinit__(self):
        # TODO: is this the best method?
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

        # Init the context
        fmpz_mod_ctx_init(self.val, (<fmpz>mod).val)
    
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
    *nmod*.

        >>> fmpz_mod(-1, fmpz_mod_ctx(2**127 - 1))
        fmpz_mod(170141183460469231731687303715884105726, 170141183460469231731687303715884105727)

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

    def any_as_fmpz_mod(self, obj):
        try:
            return self.ctx(obj)
        except NotImplementedError:
            return NotImplemented

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

    def __richcmp__(self, other, int op):
        cdef bint res
        if op != 2 and op != 3:
            raise TypeError("fmpz_mod cannot be ordered")

        if not typecheck(other, fmpz_mod):
            other = self.any_as_fmpz_mod(other)

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
        other = self.any_as_fmpz_mod(other)
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

    def __sub__(self, other):
        other = self.any_as_fmpz_mod(other)
        if other is NotImplemented:
            return NotImplemented

        cdef fmpz_mod res
        res = fmpz_mod.__new__(fmpz_mod)
        res.ctx = self.ctx

        fmpz_mod_sub(
            res.val, self.val, (<fmpz_mod>other).val, self.ctx.val
        )
        return res

    def __rsub__(self, other):
        return self.__sub__(other).__neg__()

    def __mul__(self, other):
        other = self.any_as_fmpz_mod(other)
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
        
        # Division when left and right are fmpz_mod
        if typecheck(left, fmpz_mod) and typecheck(right, fmpz_mod):
            res.ctx = (<fmpz_mod>left).ctx
            if not (<fmpz_mod>left).ctx == (<fmpz_mod>right).ctx:
                raise ValueError("moduli must match")
            check = fmpz_mod_divides(
                res.val, (<fmpz_mod>left).val, (<fmpz_mod>right).val, res.ctx.val
            ) 
        
        # Case when only left is fmpz_mod
        elif typecheck(left, fmpz_mod):
            res.ctx = (<fmpz_mod>left).ctx
            right = any_as_fmpz(right)
            if right is NotImplemented:
                return NotImplemented   
            check = fmpz_mod_divides(
                res.val, (<fmpz_mod>left).val, (<fmpz>right).val, res.ctx.val
            ) 

        # Case when right is an fmpz_mod
        else:
            res.ctx = (<fmpz_mod>right).ctx
            left = any_as_fmpz(left)
            if left is NotImplemented:
                return NotImplemented   
            check = fmpz_mod_divides(
                res.val, (<fmpz>left).val, (<fmpz_mod>right).val, res.ctx.val
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

    def inverse(self, check=True):
        """
        Computes a^-1 mod N

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
