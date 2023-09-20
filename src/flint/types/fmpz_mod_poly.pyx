from flint.flintlib.fmpz_mod_poly cimport *
from flint.flintlib.fmpz cimport(
    fmpz_equal,
    fmpz_set
)
from flint.types.fmpz cimport fmpz
from flint.types.fmpz_mod cimport fmpz_mod_ctx, fmpz_mod
from flint.types.fmpz_poly cimport any_as_fmpz_poly, fmpz_poly

from flint.flint_base.flint_base cimport flint_poly
from flint.utils.typecheck cimport typecheck

cdef class fmpz_mod_poly_ctx:
    """
    NOTE:

    Technically this could just be the same as `fmpz_mod_ctx`,
    however, the usage of fmpz_mod_ctx allows the creation of 
    `fmpz_mod` types by calling the context class. For symmetry
    we allow this to be the case here.
    """
    def __cinit__(self):
        pass

    def __dealloc__(self):
        pass

    def __init__(self, mod):
        # Allow context to be made from fmpz_mod_ctx
        if typecheck(mod, fmpz_mod_ctx):
            self.mod = mod
        else: # Otherwise attempt to create context from moduli
            self.mod = fmpz_mod_ctx(mod)

    def modulus(self):
        """
        Return the modulus from the context as an fmpz
        type

            >>> mod_ctx = fmpz_mod_poly_ctx(2**127 - 1)
            >>> mod_ctx.modulus()
            170141183460469231731687303715884105727

        """
        return self.mod.modulus()

    def __eq__(self, other):
        # Most often, we expect both `fmpz_mod_poly` to be pointing 
        # to the same ctx, so this seems the fastest way to check
        if self is other:
            return True
        
        # If they're not the same object in memory, they may have the
        # same modulus, which is good enough
        if typecheck(other, fmpz_mod_poly_ctx):
            return self.mod == (<fmpz_mod_poly_ctx>other).mod
        return False
    
    def __hash__(self):
        return hash(repr(self))

    def __str__(self):
        return f"Context for fmpz_mod_poly with modulus: {self.modulus()}"

    def __repr__(self):
        return f"fmpz_mod_poly_ctx({self.modulus()})"

    def __call__(self, val):
        return fmpz_mod_poly(val, self.mod)


cdef class fmpz_mod_poly(flint_poly):
    """
    """
    def __cinit__(self):
        fmpz_mod_poly_init(self.val, self.ctx.val)

    def __dealloc__(self):
        fmpz_mod_poly_clear(self.val, self.ctx.val)

    def __init__(self, val, ctx):
        if not typecheck(ctx, fmpz_mod_ctx):
            raise TypeError
        self.ctx = ctx

        val = any_as_fmpz_poly(val)
        if val is NotImplemented:
            raise TypeError

        fmpz_mod_poly_set_fmpz_poly(
            self.val, (<fmpz_poly>val).val, self.ctx.val
        )

    def __getitem__(self, long i):
        cdef fmpz_mod x
        x = fmpz_mod.__new__(fmpz_mod)
        x.ctx = self.ctx
        if i < 0:
            return x
        fmpz_mod_poly_get_coeff_fmpz(x.val, self.val, i, self.ctx.val)
        return x

    def __setitem__(self, long i, x):
        if i < 0:
            raise ValueError("cannot assign to index < 0 of polynomial")
        v = self.ctx.any_as_fmpz_mod(x)
        if v is NotImplemented:
            raise TypeError
        fmpz_mod_poly_set_coeff_fmpz(self.val, i, (<fmpz_mod>v).val, self.ctx.val)

    def repr(self):
        return "fmpz_mod_poly([%s])" % (", ".join(map(str, self.coeffs())))

    def __len__(self):
        return fmpz_mod_poly_length(self.val, self.ctx.val)

    cpdef long length(self):
        return fmpz_mod_poly_length(self.val, self.ctx.val)

    cpdef long degree(self):
        return fmpz_mod_poly_degree(self.val, self.ctx.val)
