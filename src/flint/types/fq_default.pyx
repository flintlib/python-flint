from flint.pyflint cimport global_random_state
from flint.types.fmpz cimport fmpz, any_as_fmpz
from flint.types.fmpz_poly cimport fmpz_poly
from flint.types.fmpz_mod_poly cimport fmpz_mod_poly, fmpz_mod_poly_ctx
from flint.types.nmod_poly cimport nmod_poly
from flint.utils.typecheck cimport typecheck

cdef class fq_default_ctx:
    r"""
    Context object for creating :class:`~.fq_default`.

    Finite fields can be initialized in one of two possible ways. The
    first is by providing characteristic and degree:

        >>> fq_default_ctx.from_order(5, 2, 'y')
        fq_default_ctx.from_modulus(x^2 + 4*x + 2, 'y', 1)

    The second is by giving an irreducible polynomial of type
    :class:`~.nmod_poly` or :class:`~.fmpz_mod_poly`:

        >>> from flint import fmpz_mod_poly_ctx
        >>> modulus = fmpz_mod_poly_ctx(11)([1,0,1])
        >>> fq_default_ctx.from_modulus(modulus, 'x')
        fq_default_ctx.from_modulus(x^2 + 1, 'x', 2)
    
    For more details, see the documentation of :method:`~.from_order`
    and :method:`~.from_modulus`.
    """
    def __cinit__(self):
        pass
    
    def __dealloc__(self):
        if self._initialized:
            fq_default_ctx_clear(self.val)

    def __init__(self, *args, **kwargs):
        raise TypeError("This class cannot be instantiated directly. Use .from_order() or .from_modulus().")

    @staticmethod
    cdef fq_default_ctx c_from_order(fmpz p, int d, char *var,
                                     fq_default_type type=fq_default_type.DEFAULT):
        cdef fq_default_ctx ctx = fq_default_ctx.__new__(fq_default_ctx)
        ctx.var = var
        fq_default_ctx_init_type(ctx.val, p.val, d, ctx.var, type)
        ctx._initialized = True
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

            >>> gf = fq_default_ctx.from_order(5, 2, 'y')
            >>> gf
            fq_default_ctx.from_modulus(x^2 + 4*x + 2, 'y', 1)
            >>> gf.type
            <fq_default_type.FQ_ZECH: 1>
        
            >>> gf = fq_default_ctx.from_order(5, 2, 'y', fq_default_type.FQ_NMOD)
            >>> gf
            fq_default_ctx.from_modulus(x^2 + 4*x + 2, 'y', 2)
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
        ctx._initialized = True

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

            >>> from flint import fmpz_mod_poly_ctx
            >>> modulus = fmpz_mod_poly_ctx(11)([1,0,1])
            >>> gf = fq_default_ctx.from_modulus(modulus, 'x')
            >>> gf
            fq_default_ctx.from_modulus(x^2 + 1, 'x', 2)
            >>> gf.type
            <fq_default_type.FQ_NMOD: 2>
        
            >>> gf = fq_default_ctx.from_modulus(modulus, 'x', fq_default_type.FQ)
            >>> gf
            fq_default_ctx.from_modulus(x^2 + 1, 'x', 3)
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
            >>> gf = fq_default_ctx.from_order(5, 2, 'y')
            >>> gf.degree()
            2
        """
        return fq_default_ctx_degree(self.val)

    def prime(self):
        """
            >>> gf = fq_default_ctx.from_order(5, 2, 'y')
            >>> gf.prime()
            5
        """
        cdef fmpz p
        p = fmpz.__new__(fmpz)
        fq_default_ctx_prime(p.val, self.val)
        return p

    characteristic = prime
    
    def order(self):
        """
            >>> gf = fq_default_ctx.from_order(5, 2, 'y')
            >>> gf.order()
            25
        """
        cdef fmpz q
        q = fmpz.__new__(fmpz)
        fq_default_ctx_order(q.val, self.val)
        return q

    def modulus(self):
        """
        Return the modulus from the context as an fmpz_mod_poly type

            >>> gf = fq_default_ctx.from_order(5, 2, 'y')
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

            >>> gf = fq_default_ctx.from_order(5, 1, "x")
            >>> gf.zero()
            0
        """
        cdef fq_default res
        res = self.new_ctype_fq_default()
        res.ctx = self
        fq_default_zero(res.val, self.val)
        return res
    
    def one(self):
        """
        Return the one element

            >>> gf = fq_default_ctx.from_order(5, 1, "x")
            >>> gf.one()
            1
        """
        cdef fq_default res
        res = self.new_ctype_fq_default()
        res.ctx = self
        fq_default_one(res.val, self.val)
        return res

    def random_element(self, not_zero=False):
        r"""
        Return a random element of the finite field
        """
        cdef fq_default res
        res = self.new_ctype_fq_default()
        res.ctx = self
        if not_zero:
            fq_default_rand_not_zero(res.val, global_random_state, self.val)
        else:
            fq_default_rand(res.val, global_random_state, self.val)
        return res

    cdef new_ctype_fq_default(self):
        return fq_default.__new__(fq_default, None, self)

    cdef set_any_as_fq_default(self, fq_default_t val, obj):
        if typecheck(obj, fmpz):
            fq_default_set_fmpz(val, (<fmpz>obj).val, self.val)
            return 0
        return NotImplemented

    def __eq__(self, other):
        """
        Two finite field context compare equal if they have same
        characteristic, modulus, type and variable
        
            >>> from flint import fmpz_mod_poly_ctx
            >>> R = fmpz_mod_poly_ctx(5)
            >>> modulus = R([2,4,1])
            >>> gf = fq_default_ctx.from_order(5, 2, 'y')
            >>> gf2 = fq_default_ctx.from_modulus(modulus, 'y', 1)
            >>> gf2 == gf
            True
            >>> gf3 = fq_default_ctx.from_modulus(modulus, 'x', 1)
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
        return f"fq_default_ctx.from_modulus({self.modulus()!r}, '{self.var.decode()}', {self.type})"

    def __call__(self, val):
        return fq_default(val, self)


cdef class fq_default(flint_scalar):
    def __cinit__(self, val, ctx):
        if not typecheck(ctx, fq_default_ctx):
            raise TypeError
        self.ctx = ctx
        fq_default_init(self.val, self.ctx.val)

    def __dealloc__(self):
        if self.ctx is not None:
            fq_default_clear(self.val, self.ctx.val)

    def __init__(self, val, ctx):
        if not typecheck(ctx, fq_default_ctx):
            raise TypeError
        self.ctx = ctx

        check = self.ctx.set_any_as_fq_default(self.val, val)
        if check is NotImplemented:
            raise TypeError

    def __int__(self):
        """
        Attempts to lift self to an integer of type fmpz in [0, p-1]
        """
        cdef fmpz x = fmpz.__new__(fmpz)
        res = fq_default_get_fmpz(x.val, self.val, self.ctx.val)
        if res == 1:
            return int(x)
        raise ValueError("fq element has no lift to the integers")

    def polynomial(self):
        """
        Returns a representative of ``self`` as a polynomial in `(Z/pZ)[x] / h(x)`
        where `h(x)` is the defining polynomial of the finite field.
        """
        cdef fmpz_mod_poly_ctx ctx
        cdef fmpz_mod_poly pol

        ring_ctx = fmpz_mod_poly_ctx(self.ctx.prime())
        pol = ring_ctx.new_ctype_poly()
        fq_default_get_fmpz_mod_poly((<fmpz_mod_poly>pol).val, self.val, self.ctx.val)

        return pol

    def str(self):
        return self.polynomial().str(var=self.ctx.var)

    def __repr__(self):
        # TODO: what do we want here?
        return str(self)

    # =================================================
    # Comparisons
    # =================================================  
    def is_zero(self):
        return 1 == fq_default_is_zero(self.val, self.ctx.val)

    def is_one(self):
        return 1 == fq_default_is_zero(self.val, self.ctx.val)

    # =================================================
    # Generic arithmetic required by flint_scalar
    # =================================================   

    def _neg_(self):
        cdef fq_default res
        res = self.ctx.new_ctype_fq_default()
        fq_default_neg(res.val, self.val, self.ctx.val)
        return res

    def _add_(self, other):
        return NotImplemented

    @staticmethod
    def _sub_(left, right):
        return NotImplemented

    def _mul_(self, other):
        return NotImplemented

    @staticmethod
    def _div_(left, right):
        return NotImplemented

    def _invert_(self):
        cdef fq_default res
        res = self.ctx.new_ctype_fq_default()
        fq_default_inv(res.val, self.val, self.ctx.val)
        return res

    # =================================================
    # Additional arithmetic
    # ================================================= 

    def square(self):
        cdef fq_default res
        res = self.ctx.new_ctype_fq_default()
        fq_default_sqr(res.val, self.val, self.ctx.val)
        return res

    def __pow__(self, e):
        return NotImplemented

    def sqrt(self):
        cdef fq_default res
        res = self.ctx.new_ctype_fq_default()
        check = fq_default_sqrt(res.val, self.val, self.ctx.val)
        if check:
            return res
        raise ValueError("element is not a square")

    def is_square(self):
        return 1 == fq_default_is_square(self.val, self.ctx.val)

    def pth_root(self):
        cdef fq_default res
        res = self.ctx.new_ctype_fq_default()
        fq_default_pth_root(res.val, self.val, self.ctx.val)
        return res

    # =================================================
    # Special functions
    # ================================================= 

    def trace(self):
        return NotImplemented

    def norm(self):
        return NotImplemented

    def frobenius(self):
        return NotImplemented
