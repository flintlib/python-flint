from flint.flintlib.flint cimport slong
from flint.flintlib.fmpz cimport (
    fmpz_t,
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

    # TODO: should this be allowed, or should
    #       we make a ctx immutatble?
    def set_modulus(self, n):
        """
        """
        # Ensure modulus is fmpz type
        mod = any_as_fmpz(n)
        if mod is NotImplemented:
            raise NotImplementedError("TODO")

        # Ensure modulus is positive
        if mod < 1:
            raise ValueError("Modulus is expected to be positive")

        fmpz_mod_ctx_set_modulus(self.val, (<fmpz>mod).val)


cdef class fmpz_mod(flint_scalar):
    def __init__(self, val, ctx):
        self.ctx = ctx

        if not typecheck(val, fmpz):
            val = any_as_fmpz(val)
            if val is NotImplemented:
                raise NotImplementedError("TODO")

        fmpz_mod_set_fmpz(self.val, (<fmpz>val).val, (<fmpz_mod_ctx_t>self.ctx.val))

    def is_zero(self):
        return self == 0
    
    # TODO: kind of pointless, as we always ensure canonical on init?
    def is_canonical(self):
        cdef bint res
        res = fmpz_mod_is_canonical(self.val, (<fmpz_mod_ctx_t>self.ctx.val))
        return res == 1

    def is_one(self):
        cdef bint res
        res = fmpz_mod_is_one(self.val, (<fmpz_mod_ctx_t>self.ctx.val))
        return res == 1

    def __richcmp__(s, t, int op):
        cdef bint res
        if op != 2 and op != 3:
            raise TypeError("fmpz_mod cannot be ordered")
        # TODO: is this the best method for comparison?
        if typecheck(s, fmpz_mod) and typecheck(t, fmpz_mod):
            res = ((<fmpz_mod>s).val[0] == (<fmpz_mod>t).val[0]) and \
                  ((<fmpz_mod>s).ctx.modulus() == (<fmpz_mod>t).ctx.modulus())
            if op == 2:
                return res
            else:
                return not res
        # TODO: is this the best method for comparison?
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

    def __neg__(self):
        res = fmpz()
        fmpz_mod_neg(
            res.val, self.val, 
            (<fmpz_mod_ctx_t>self.ctx.val))
        return self.ctx(res)

    # TODO: proper type handing for the other...
    def __add__(self, other):
        res = fmpz()
        fmpz_mod_add(
            res.val, self.val, (<fmpz_mod>other).val, 
            (<fmpz_mod_ctx_t>self.ctx.val)
        )
        return self.ctx(res)

    def __radd__(self, other):
        return other + self

    def __iadd__(self, other):
        fmpz_mod_add(
            self.val, self.val, (<fmpz_mod>other).val, 
            (<fmpz_mod_ctx_t>self.ctx.val)
        )
        return self

    def __sub__(self, other):
        pass

    def __rsub__(self, other):
        pass

    def __isub__(self, other):
        pass

    def __mul__(self, other):
        pass

    def __rmul__(self, other):
        pass

    def __imul__(self, other):
        pass

    def __truediv__(self, other):
        pass

    def __floordiv__(self, other):
        raise NotImplemented
