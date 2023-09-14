from flint.flintlib.flint cimport slong
from flint.flintlib.fmpz cimport (
    fmpz_t,
    fmpz_set,
    COEFF_IS_MPZ,
    fmpz_get_str
)
from flint.flintlib.fmpz_mod cimport (
    fmpz_mod_ctx_t,
    fmpz_mod_ctx_init,
    fmpz_mod_ctx_set_modulus,
    fmpz_mod_set_fmpz
)
from flint.utils.typecheck cimport typecheck
from flint.utils.conversion cimport chars_from_str, str_from_chars, _str_trunc
from flint.flint_base.flint_base cimport flint_scalar
from flint.types.fmpz cimport (
    fmpz,
    any_as_fmpz,
)

cimport libc.stdlib

# TODO: import this from types.fmpz somehow
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
    def __init__(self, n):
        """
        """
        # Ensure modulus is fmpz type
        mod = any_as_fmpz(n)
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

# This is buggy and broken
cdef class fmpz_mod(flint_scalar):
    def __init__(self, val, ctx):
        self.ctx = ctx

        val_fmpz = any_as_fmpz(val)
        assert typecheck(val_fmpz, fmpz)
        fmpz_mod_set_fmpz(self.val, (<fmpz>val_fmpz).val, (<fmpz_mod_ctx_t>self.ctx.val))

    def __repr__(self):
        return "fmpz_mod({}, {})".format(
            fmpz_get_intlong(self.val),
            self.ctx.modulus()
        )

    def __str__(self):
        return str(fmpz_get_intlong(self.val))

    def __int__(self):
        return fmpz_get_intlong(self.val)
