from cpython.list cimport PyList_GET_SIZE
from flint.flintlib.fmpz_mod cimport fmpz_mod_set_fmpz
from flint.flintlib.fmpz_mod_poly cimport *
from flint.flintlib.fmpz_mod_poly_factor cimport *

from flint.flintlib.fmpz cimport(
    fmpz_equal,
    fmpz_set,
    fmpz_init,
    fmpz_clear,
    fmpz_set_si
)
from flint.types.fmpz cimport fmpz, any_as_fmpz
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

    def gen(self):
        """
        Return the generator of the polynomial `x`

        >>> mod_ctx = fmpz_mod_poly_ctx(2**127 - 1)
        >>> mod_ctx.gen()
        fmpz_mod_poly([0, 1], fmpz_mod_poly_ctx(170141183460469231731687303715884105727))
        >>> str(mod_ctx.gen())
        'x'
        """
        return self.any_as_fmpz_mod_poly([0, 1])

    cdef set_list_as_fmpz_mod_poly(self, fmpz_mod_poly_t poly, val):
        cdef long i, n
        cdef fmpz_t x

        n = PyList_GET_SIZE(val)
        fmpz_mod_poly_fit_length(poly, n, self.mod.val)
        
        # TODO: should we support conversion from nmod?
        fmpz_init(x)
        for i in range(n):
            if typecheck(val[i], fmpz_mod):
                if self.mod != (<fmpz_mod>(val[i])).ctx:
                    raise ValueError("moduli must match")
                fmpz_mod_poly_set_coeff_fmpz(
                    poly, i, (<fmpz_mod>(val[i])).val, self.mod.val
                )
            elif typecheck(val[i], fmpz):
                fmpz_mod_poly_set_coeff_fmpz(
                    poly, i, (<fmpz>(val[i])).val, self.mod.val
                )
            else:
                val_fmpz = any_as_fmpz(val[i])
                if val_fmpz is NotImplemented:
                    fmpz_clear(x)
                    raise TypeError("unsupported coefficient in list")
                fmpz_mod_poly_set_coeff_fmpz(
                    poly, i, (<fmpz>(val_fmpz)).val, self.mod.val
                )
        fmpz_clear(x)
        return 0

    cdef set_any_as_fmpz_mod_poly(self, fmpz_mod_poly_t poly, obj):
        if typecheck(obj, list):
            return self.set_list_as_fmpz_mod_poly(poly, obj)

        # Convert fmpz_mod to constant poly
        if typecheck(obj, fmpz_mod):
            if self.mod != (<fmpz_mod>obj).ctx:
                raise ValueError("moduli must match")
            fmpz_mod_poly_set_fmpz(
                poly, (<fmpz_mod>obj).val, self.mod.val
            )
            return 0

        # Reduced fmpz_poly modulo mod
        if typecheck(obj, fmpz_poly):
            fmpz_mod_poly_set_fmpz_poly(
                poly, (<fmpz_poly>obj).val, self.mod.val
            )
            return 0


        obj = any_as_fmpz(obj)
        if obj is NotImplemented:
            return NotImplemented
        fmpz_mod_poly_set_fmpz(
            poly, (<fmpz>obj).val, self.mod.val
        )
        return 0

    cdef any_as_fmpz_mod_poly(self, obj):
        # Convert fmpz_mod_poly
        if typecheck(obj, fmpz_mod_poly):
            if self.mod != (<fmpz_mod_poly>obj).ctx:
                raise ValueError("moduli must match")
            return obj
        
        cdef fmpz_mod_poly res
        res = fmpz_mod_poly.__new__(fmpz_mod_poly)
        check = self.set_any_as_fmpz_mod_poly(res.val, obj)
        if check is NotImplemented:
            return NotImplemented
        
        res.ctx = self
        return res

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
        return fmpz_mod_poly(val, self)


cdef class fmpz_mod_poly(flint_poly):
    """
    """
    def __cinit__(self):
        fmpz_mod_poly_init(self.val, self.ctx.mod.val)

    def __dealloc__(self):
        fmpz_mod_poly_clear(self.val, self.ctx.mod.val)

    def __init__(self, val, ctx):
        if not typecheck(ctx, fmpz_mod_poly_ctx):
            raise TypeError
        self.ctx = ctx
        self.ctx.set_any_as_fmpz_mod_poly(self.val, val)

    def __getitem__(self, long i):
        cdef fmpz_mod x
        x = fmpz_mod.__new__(fmpz_mod)
        x.ctx = self.ctx.mod
        if i < 0:
            return x
        fmpz_mod_poly_get_coeff_fmpz(
            x.val, self.val, i, self.ctx.mod.val
        )
        return x

    def __setitem__(self, long i, x):
        if i < 0:
            raise ValueError("cannot assign to index < 0 of polynomial")
        v = self.ctx.mod.any_as_fmpz_mod(x)
        if v is NotImplemented:
            raise TypeError
        fmpz_mod_poly_set_coeff_fmpz(
            self.val, i, (<fmpz_mod>v).val, self.ctx.mod.val
        )

    def __len__(self):
        return fmpz_mod_poly_length(self.val, self.ctx.mod.val)

    cpdef long length(self):
        return fmpz_mod_poly_length(self.val, self.ctx.mod.val)

    cpdef long degree(self):
        return fmpz_mod_poly_degree(self.val, self.ctx.mod.val)

    def is_zero(self):
        return 0 != fmpz_mod_poly_is_zero(self.val, self.ctx.mod.val)
    
    def is_one(self):
        return 0 != fmpz_mod_poly_is_one(self.val, self.ctx.mod.val)
    
    def is_gen(self):
        return 0 != fmpz_mod_poly_is_gen(self.val, self.ctx.mod.val)

    def __richcmp__(self, other, int op):
        cdef bint res
        if op != 2 and op != 3:
            raise TypeError("fmpz_mod_poly cannot be ordered")

        if not typecheck(other, fmpz_mod_poly):
            other = self.ctx.any_as_fmpz_mod_poly(other)

        if typecheck(other, fmpz_mod_poly):
            res = (self.ctx == (<fmpz_mod_poly>other).ctx) and \
                  fmpz_mod_poly_equal(self.val, (<fmpz_mod_poly>other).val, self.ctx.mod.val)
            if op == 2:
                return res
            else:
                return not res
        else:
            return NotImplemented


    def is_irreducible(self):
        pass

    def is_squarefree(self):
        pass

    def factor(self, algorithm=None):
        """
        Factors self into irreducible factors, returning a tuple
        (c, factors) where c is the content of the coefficients and
        factors is a list of (poly, exp) pairs.
        """
        cdef fmpz_mod_poly_factor_t fac
        cdef int i
        fmpz_mod_poly_factor_init(fac, self.ctx.mod.val)

        if algorithm == None:
            fmpz_mod_poly_factor(fac, self.val, self.ctx.mod.val)
        elif algorithm == "cantor_zassenhaus":
            fmpz_mod_poly_factor_cantor_zassenhaus(fac, self.val, self.ctx.mod.val)
        elif algorithm == "kaltofen_shoup":
            fmpz_mod_poly_factor_kaltofen_shoup(fac, self.val, self.ctx.mod.val)
        elif algorithm == "berlekamp":
            fmpz_mod_poly_factor_berlekamp(fac, self.val, self.ctx.mod.val)
        else:
            raise ValueError("unknown algorithm")
        res = [0] * fac.num

        cdef fmpz_mod_poly u
        for i in range(fac.num):
            u = fmpz_mod_poly.__new__(fmpz_mod_poly)
            u.ctx = self.ctx
            fmpz_mod_poly_set(u.val, &fac.poly[i], self.ctx.mod.val)
            exp = fac.exp[i]
            res[i] = (u, exp)
        return res
