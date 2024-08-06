from flint.types.fmpz cimport fmpz, any_as_fmpz
from flint.types.fmpz_mod_poly cimport fmpz_mod_poly, fmpz_mod_poly_ctx
from flint.types.nmod_poly cimport nmod_poly
from flint.utils.typecheck cimport typecheck

cdef class fq_default_ctx:
    r"""
    Context object for creating :class:`~.fq_default`.

    Finite fields can be initialized in one of two possible ways. The
    first is by providing characteristic and degree:

        >>> fq_default_ctx.from_order(fmpz(5), 2, 'y')
        fq_default_ctx.from_modulus(x^2 + 4*x + 2, b'y', 1)

    The second is by giving an irreducible polynomial of type
    :class:`~.nmod_poly` or :class:`~.fmpz_mod_poly`:

        >>> fq_default_ctx.from_modulus(fmpz_mod_poly([1,0,1], fmpz_mod_poly_ctx(11)), 'x')
        fq_default_ctx.from_modulus(x^2 + 1, b'x', 2)
    
    For more details, see the documentation of :method:`~.from_order`
    and :method:`~.from_modulus`.
    """
    def __dealloc__(self):
        if self.initialized is not None:
            fq_default_ctx_clear(self.val)

    def __init__(self):
        raise TypeError("This class cannot be instantiated directly. Use .from_order() or .from_modulus().")

    @staticmethod
    cdef fq_default_ctx c_from_order(fmpz p, int d, char *var,
                                     fq_default_type type=fq_default_type.DEFAULT):
        cdef fq_default_ctx ctx = fq_default_ctx.__new__(fq_default_ctx)
        ctx.var = var
        fq_default_ctx_init_type(ctx.val, p.val, d, ctx.var, type)
        ctx.initialized = True
        return ctx

    @staticmethod
    def from_order(p, d, var, type=fq_default_type.DEFAULT):
        """
        Construct a context for the finite field GF(p^d).

        `var` is a name for the ring generator of this field over GF(p).

        The optional parameter `type` select the implementation. The
        possible values are:

        - `fq_default_ctx.DEFAULT`: implementation automatically decided by Flint (default),
        - `fq_default_ctx.FQ_ZECH`: Use `fq_zech_t`,
        - `fq_default_ctx.FQ_NMOD`: Use `fq_nmod_t`,
        - `fq_default_ctx.FQ`: Use `fq_t`.

            >>> gf = fq_default_ctx.from_order(fmpz(5), 2, 'y')
            >>> gf
            fq_default_ctx.from_modulus(x^2 + 4*x + 2, b'y', 1)
            >>> gf.type
            <fq_default_type.FQ_ZECH: 1>
        
            >>> gf = fq_default_ctx.from_order(fmpz(5), 2, 'y', fq_default_type.FQ_NMOD)
            >>> gf
            fq_default_ctx.from_modulus(x^2 + 4*x + 2, b'y', 2)
            >>> gf.type
            <fq_default_type.FQ_NMOD: 2>
        

        """
        # c_from_order expects the characteristic to be fmpz type
        order = any_as_fmpz(p)
        if order is NotImplemented:
            raise TypeError(f"cannot coerce {p = } to type fmpz")
        
        # the degree must be strictly positive
        if d < 1:
            raise ValueError(f"the degree must be positive, got {d = }")

        # TODO: we should allow var to be None when d == 1
        if isinstance(var, str):
            var = var.encode()
        
        return fq_default_ctx.c_from_order(order, d, var, type)

    @staticmethod
    cdef fq_default_ctx c_from_modulus(modulus, char *var,
                                       fq_default_type type=fq_default_type.DEFAULT): 
        cdef fq_default_ctx ctx = fq_default_ctx.__new__(fq_default_ctx)

        ctx.var = var
        if typecheck(modulus, fmpz_mod_poly):
            fq_default_ctx_init_modulus_type(ctx.val, (<fmpz_mod_poly>modulus).val,
                                             (<fmpz_mod_poly>modulus).ctx.mod.val, ctx.var, type)
        elif typecheck(modulus, nmod_poly):
            fq_default_ctx_init_modulus_nmod_type(ctx.val, (<nmod_poly>modulus).val,
                                                  ctx.var, type)
        else:
            raise TypeError(f"modulus must be fmpz_mod_poly or nmod_poly, got {modulus!r}")
        ctx.initialized = True

        return ctx

    @staticmethod
    def from_modulus(modulus, var, type=fq_default_type.DEFAULT):
        """
        Construct a context for a finite field from an irreducible polynomial.

        `modulus` may be of type :class:`~.fmpz_mod_poly` or :class:`~.nmod_poly`.
        
        `var` is a name for the ring generator of this field over the prime field.

        The optional parameter `type` select the implementation. The
        possible values are:

        - `fq_default_ctx.DEFAULT`: implementation automatically decided by Flint (default),
        - `fq_default_ctx.FQ_ZECH`: Use `fq_zech_t`,
        - `fq_default_ctx.FQ_NMOD`: Use `fq_nmod_t`,
        - `fq_default_ctx.FQ`: Use `fq_t`.

            >>> gf = fq_default_ctx.from_modulus(fmpz_mod_poly([1,0,1], fmpz_mod_poly_ctx(11)), 'x')
            >>> gf
            fq_default_ctx.from_modulus(x^2 + 1, b'x', 2)
            >>> gf.type
            <fq_default_type.FQ_NMOD: 2>
        
            >>> gf = fq_default_ctx.from_modulus(fmpz_mod_poly([1,0,1], fmpz_mod_poly_ctx(11)), 'x', fq_default_type.FQ)
            >>> gf
            fq_default_ctx.from_modulus(x^2 + 1, b'x', 2)
            >>> gf.type
            <fq_default_type.FQ: 3>

        """
        if isinstance(var, str):
            var = var.encode()
        return fq_default_ctx.c_from_modulus(modulus, var, type)
    
    @property
    def type(self):
        """
        Return the implementation of this context. It is one of:
        
        - `fq_default_ctx.FQ_ZECH`: Using `fq_zech_t`,
        - `fq_default_ctx.FQ_NMOD`: Using `fq_nmod_t`,
        - `fq_default_ctx.FQ`: Using `fq_t`.
        """
        return fq_default_type(fq_default_ctx_type(self.val))

    def degree(self):
        """
            >>> gf = fq_default_ctx.from_order(fmpz(5), 2, 'y')
            >>> gf.degree()
            2
        """
        return fq_default_ctx_degree(self.val)

    def prime(self):
        """
            >>> gf = fq_default_ctx.from_order(fmpz(5), 2, 'y')
            >>> gf.prime()
            5
        """
        cdef fmpz p
        p = fmpz.__new__(fmpz)
        fq_default_ctx_prime(p.val, self.val)
        return p
    
    def order(self):
        """
            >>> gf = fq_default_ctx.from_order(fmpz(5), 2, 'y')
            >>> gf.prime()
            25
        """
        cdef fmpz q
        q = fmpz.__new__(fmpz)
        fq_default_ctx_order(q.val, self.val)
        return q

    def modulus(self):
        """
        Return the modulus from the context as an fmpz_mod_poly type

            >>> gf = fq_default_ctx.from_order(fmpz(5), 2, 'y')
            >>> gf.modulus()
            x^2 + 4*x + 2

        """
        cdef fmpz_mod_poly_ctx ctx
        cdef fmpz_mod_poly pol
        ctx = fmpz_mod_poly_ctx(self.prime())
        pol = ctx.new_ctype_poly()
        fq_default_ctx_modulus(pol.val, self.val)
        return pol
    
    def zero(self):
        """
        Return the zero element

            >>> F = fmpz_mod_ctx(163)
            >>> F.zero()
            fmpz_mod(0, 163)
        """
        pass
    
    def one(self):
        """
        Return the one element

            >>> F = fmpz_mod_ctx(163)
            >>> F.one()
            fmpz_mod(1, 163)
        """
        pass

    def random_element(self):
        r"""
        Return a random element in :math:`\mathbb{Z}/N\mathbb{Z}`
        """
        pass

    cdef set_any_as_fq_default(self, fq_default_t val, obj):
        pass

    cdef any_as_fmpz_mod(self, obj):
        pass

    def __eq__(self, other):
        """
        Two finite field context compare equal if they have same
        characteristic, modulus, type and variable
        
            >>> gf = fq_default_ctx.from_order(fmpz(5), 2, 'y')
            >>> gf2 = fq_default_ctx.from_modulus(fmpz_mod_poly([2,4,1], fmpz_mod_poly_ctx(5)), 'y', 1)
            >>> gf2 == gf
            True
            >>> gf3 = fq_default_ctx.from_modulus(fmpz_mod_poly([2,4,1], fmpz_mod_poly_ctx(5)), 'x', 1)
            >>> gf3 == gf
            False

        """
        if self is other:
            return True
        
        if typecheck(other, fq_default_ctx):
            return (self.type == other.type
                    and self.var == other.var
                    and self.prime() == other.prime()
                    and self.modulus() == other.modulus())
        else:
            raise TypeError(f"Cannot compare fq_default_ctx with {type(other)}")

    def __hash__(self):
        return hash((self.type, self.var, self.prime(), self.modulus()))

    def __str__(self):
        return f"Context for fq_default in GF({self.prime()}^{self.degree()})[{self.var.decode()}]/({self.modulus().str(var=self.var.decode())})"

    def __repr__(self):
        return f"fq_default_ctx.from_modulus({self.modulus()!r}, {self.var.encode()}, {self.type})"

    def __call__(self, val):
        return fq_default(val, self)


cdef class fq_default(flint_scalar):
    pass
