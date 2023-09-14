from flint.flint_base.flint_base cimport flint_scalar
from flint.flintlib.fmpz_mod cimport (
    fmpz_mod_ctx_init,
    fmpz_mod_ctx_set_modulus,
)

from flint.types.fmpz cimport fmpz, any_as_fmpz

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

    def __repr__(self):
        # return "Context for fmpz_mod with modulus: %s" % self.ctx.n
        return "Stuff: (%s, %s, %s)" % (self.val.n, self.val.mod, self.val.n_limbs)

    def __str__(self):
        return "fmpz_mod_ctx(%s)" % self.val.mod

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
    pass
