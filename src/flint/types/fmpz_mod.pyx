from flint.flintlib.flint cimport slong
from flint.flintlib.fmpz cimport (
    fmpz_t,
    fmpz_one,
    fmpz_set,
    COEFF_IS_MPZ,
    fmpz_get_str
)
from flint.flintlib.fmpz_mod cimport *

from flint.utils.typecheck cimport typecheck
from flint.utils.conversion cimport str_from_chars
from flint.flint_base.flint_base cimport flint_scalar
from flint.types.fmpz cimport (
    fmpz,
    any_as_fmpz,
)

cimport libc.stdlib

# TODO: import this from types.fmpz somehow, it's not good
# to have the function repeated here.
cdef fmpz_get_intlong(fmpz_t x):
    """
    Convert fmpz_t to a Python int or long.
    """
    cdef char * s
    if COEFF_IS_MPZ(x[0]):
        s = fmpz_get_str(NULL, 16, x)
        v = int(str_from_chars(s), 16)
        libc.stdlib.free(s)
        return v
    else:
        return <slong>x[0]

cdef class fmpz_mod_ctx:
    """
    """
    def __init__(self, mod):
        """
        """
        # Ensure modulus is fmpz type
        if not typecheck(mod, fmpz):
            mod = any_as_fmpz(mod)
            if mod is NotImplemented:
                raise NotImplementedError("TODO")

        # Ensure modulus is positive
        if mod < 1:
            raise ValueError("Modulus is expected to be positive")

        # Init the context
        fmpz_mod_ctx_init(self.val, (<fmpz>mod).val)
    
    # TODO: should this be cached if we make immutable?
    def modulus(self):
        """
        Return the modulus from the context as an fmpz
        type
        """
        n = fmpz()
        fmpz_set(n.val, (<fmpz_t>self.val.n))
        return n

    def __repr__(self):
        return "Context for fmpz_mod with modulus: {}".format(
            self.modulus()
        )

    def __str__(self):
        return "fmpz_mod_ctx({})".format(
            self.modulus()
        )

    def __call__(self, val):
        if not typecheck(val, fmpz):
            val = any_as_fmpz(val)
            if val is NotImplemented:
                raise NotImplementedError("TODO")

        return fmpz_mod(val, self)


cdef class fmpz_mod(flint_scalar):
    def __init__(self, val, ctx):
        self.ctx = ctx

        if not typecheck(val, fmpz):
            val = any_as_fmpz(val)
            if val is NotImplemented:
                raise NotImplementedError("TODO")

        fmpz_mod_set_fmpz(self.val, (<fmpz>val).val, self.ctx.val)

    def is_zero(self):
        return self == 0
    
    # TODO: kind of pointless, as we always ensure canonical on init?
    def is_canonical(self):
        cdef bint res
        res = fmpz_mod_is_canonical(self.val, self.ctx.val)
        return res == 1

    def is_one(self):
        cdef bint res
        res = fmpz_mod_is_one(self.val, self.ctx.val)
        return res == 1

    def __richcmp__(s, t, int op):
        cdef bint res
        if op != 2 and op != 3:
            raise TypeError("fmpz_mod cannot be ordered")
        # TODO: is this the best method for comparison?
        if typecheck(s, fmpz_mod) and typecheck(t, fmpz_mod):
            res = ((<fmpz_mod>s).val[0] == (<fmpz_mod>t).val[0]) and \
                  ((<fmpz_mod>s).ctx.val.n == (<fmpz_mod>t).ctx.val.n)
            if op == 2:
                return res
            else:
                return not res
        # TODO: is this the best method for comparison?
        # Seems like I'm doing too many type conversions?
        elif typecheck(s, fmpz_mod) and typecheck(t, int):
            res = int(s) == t % (<fmpz_mod>s).ctx.modulus()
            if op == 2:
                return res
            else:
                return not res
        return NotImplemented

    def __bool__(self):
        return not self.is_zero()

    def __repr__(self):
        return "fmpz_mod({}, {})".format(
            fmpz_get_intlong(self.val),
            self.ctx.modulus()
        )

    # TODO: seems ugly...
    def __hash__(self):
        return  hash(
            (int(self), int(self.ctx.modulus()))
        )

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

        fmpz_mod_neg(
            res.val, self.val, 
            (<fmpz_mod_ctx_t>self.ctx.val))
        return res

    def __add__(self, other):
        cdef fmpz_mod res
        res = fmpz_mod.__new__(fmpz_mod)
        res.ctx = self.ctx

        # Add two fmpz_mod if moduli match
        if typecheck(other, fmpz_mod):
            if not self.ctx.val.n == (<fmpz_mod>other).ctx.val.n:
                raise ValueError("moduli must match")

            fmpz_mod_add(
                res.val, self.val, (<fmpz_mod>other).val, self.ctx.val
            )
            return res

        # Attempt to add an fmpz to an fmpz_mod
        other = any_as_fmpz(other)
        if other is NotImplemented:
            raise NotImplementedError

        fmpz_mod_add_fmpz(
            res.val, self.val, (<fmpz>other).val, self.ctx.val
        )

        return res

    def __radd__(self, other):
        return self.__add__(other)

    def __iadd__(self, other):
        # Add two fmpz_mod if moduli match
        if typecheck(other, fmpz_mod):
            if not self.ctx.val.n == (<fmpz_mod>other).ctx.val.n:
                raise ValueError("moduli must match")

            fmpz_mod_add(
                self.val, self.val, (<fmpz_mod>other).val, self.ctx.val
            )
            return self

        # Add an fmpz to an fmpz_mod
        other = any_as_fmpz(other)
        if other is NotImplemented:
            raise NotImplementedError

        fmpz_mod_add_fmpz(
            self.val, self.val, (<fmpz>other).val, self.ctx.val
        )
        return self

    def __sub__(self, other):
        cdef fmpz_mod res
        res = fmpz_mod.__new__(fmpz_mod)
        res.ctx = self.ctx

        # Sub two fmpz_mod if moduli match
        if typecheck(other, fmpz_mod):
            if not self.ctx.val.n == (<fmpz_mod>other).ctx.val.n:
                raise ValueError("moduli must match")

            fmpz_mod_sub(
                res.val, self.val, (<fmpz_mod>other).val, self.ctx.val
            )
            return res

        # Attempt to sub an fmpz to an fmpz_mod
        other = any_as_fmpz(other)
        if other is NotImplemented:
            raise NotImplementedError

        fmpz_mod_fmpz_sub(
            res.val, self.val, (<fmpz>other).val, self.ctx.val
        )

        return res

    # TODO: is this bad? Should I just copy paste logic
    # above?
    def __rsub__(self, other):
        return self.__sub__(other).__neg__()

    def __isub__(self, other):
        if typecheck(other, fmpz_mod):
            if not self.ctx.val.n == (<fmpz_mod>other).ctx.val.n:
                raise ValueError("moduli must match")

            fmpz_mod_sub(
                self.val, self.val, (<fmpz_mod>other).val, self.ctx.val
            )
            return self

        other = any_as_fmpz(other)
        if other is NotImplemented:
            raise NotImplementedError

        fmpz_mod_sub_fmpz(
            self.val, self.val, (<fmpz>other).val, self.ctx.val
        )
        return self

    def __mul__(self, other):
        cdef fmpz_mod res
        res = fmpz_mod.__new__(fmpz_mod)
        res.ctx = self.ctx

        # Add two fmpz_mod if moduli match
        if typecheck(other, fmpz_mod):
            if not self.ctx.val.n == (<fmpz_mod>other).ctx.val.n:
                raise ValueError("moduli must match")

            fmpz_mod_mul(
                res.val, self.val, (<fmpz_mod>other).val, self.ctx.val
            )
            return res

        # Attempt to add an fmpz to an fmpz_mod
        other = any_as_fmpz(other)
        if other is NotImplemented:
            raise NotImplementedError

        # There's no `fmpz_mod_mul_fmpz`
        # So we need to make sure that other is canonical
        # TODO: should we check with `fmpz_mod_is_canoncial`
        # or just always reduce?
        fmpz_mod_set_fmpz(
            (<fmpz>other).val, (<fmpz>other).val, self.ctx.val
        )

        fmpz_mod_mul(
            res.val, self.val, (<fmpz>other).val, self.ctx.val
        )

        return res

    def __rmul__(self, other):
        return self.__mul__(other)

    def __imul__(self, other):
        # Add two fmpz_mod if moduli match
        if typecheck(other, fmpz_mod):
            if not self.ctx.val.n == (<fmpz_mod>other).ctx.val.n:
                raise ValueError("moduli must match")

            fmpz_mod_mul(
                self.val, self.val, (<fmpz_mod>other).val, self.ctx.val
            )
            return self

        # Attempt to add an fmpz to an fmpz_mod
        other = any_as_fmpz(other)
        if other is NotImplemented:
            raise NotImplementedError

        # There's no `fmpz_mod_mul_fmpz`
        # So we need to make sure that other is canonical
        # TODO: should we check with `fmpz_mod_is_canoncial`
        # or just always reduce?
        fmpz_mod_set_fmpz(
            (<fmpz>other).val, (<fmpz>other).val, self.ctx.val
        )

        fmpz_mod_mul(
            self.val, self.val, (<fmpz>other).val, self.ctx.val
        )

        return self

    @staticmethod
    def _div_(left, right):
        cdef bint check
        cdef fmpz_mod res
        res = fmpz_mod.__new__(fmpz_mod)
        
        # Division when left and right are fmpz_mod
        if typecheck(left, fmpz_mod) and typecheck(right, fmpz_mod):
            res.ctx = (<fmpz_mod>left).ctx
            if not (<fmpz_mod>left).ctx.val.n == (<fmpz_mod>right).ctx.val.n:
                raise ValueError("moduli must match")
            check = fmpz_mod_divides(
                res.val, (<fmpz_mod>left).val, (<fmpz_mod>right).val, res.ctx.val
            ) 
        
        # Case when only left is fmpz_mod
        elif typecheck(left, fmpz_mod):
            res.ctx = (<fmpz_mod>left).ctx
            right = any_as_fmpz(right)
            if right is NotImplemented:
                raise NotImplementedError   
            check = fmpz_mod_divides(
                res.val, (<fmpz_mod>left).val, (<fmpz>right).val, res.ctx.val
            ) 

        # Case when right is an fmpz_mod
        else:
            res.ctx = (<fmpz_mod>right).ctx
            left = any_as_fmpz(left)
            if left is NotImplemented:
                raise NotImplementedError        
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
        raise NotImplemented

    def inverse(self, check=True):
        cdef fmpz_mod res
        res = fmpz_mod.__new__(fmpz_mod)
        res.ctx = self.ctx

        # Warning! This will crash if there's no solution
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
