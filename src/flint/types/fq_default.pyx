from flint.pyflint cimport global_random_state
from flint.types.fmpz cimport fmpz, any_as_fmpz
from flint.types.fmpz_poly cimport fmpz_poly, any_as_fmpz_poly, fmpz_poly_set_list
from flint.types.fmpz_mod_poly cimport fmpz_mod_poly, fmpz_mod_poly_ctx
from flint.types.nmod_poly cimport nmod_poly
from flint.utils.typecheck cimport typecheck

cdef class fq_default_ctx:
    r"""
    Context object for creating :class:`~.fq_default`.

    Finite fields can be initialized in one of two possible ways. The
    first is by providing characteristic and degree:

        >>> fq_default_ctx(5, 2, 'y')
        fq_default_ctx(5, 2, 'y', x^2 + 4*x + 2, 1)

    The second is by giving an irreducible polynomial of type
    :class:`~.nmod_poly` or :class:`~.fmpz_mod_poly`:

        >>> from flint import fmpz_mod_poly_ctx
        >>> mod = fmpz_mod_poly_ctx(11)([1,0,1])
        >>> fq_default_ctx(modulus=mod)
        fq_default_ctx(11, 2, 'x', x^2 + 1, 2)
    
    For more details, see the documentation of :method:`~.from_order`
    and :method:`~.from_modulus`.
    """
    def __cinit__(self):
        pass
    
    def __dealloc__(self):
        if self._initialized:
            fq_default_ctx_clear(self.val)

    def __init__(self, p=None, degree=None, var=None, modulus=None, type=fq_default_type.DEFAULT):
        # If no variable is given, use x
        if var is None:
            var = "x"
        
        # If a modulus is given, attempt to construct from this
        if modulus is not None:
            self.set_from_modulus(modulus, var, type)
            return

        # If there's no modulus and no prime, we can't continue
        if p is None:
            raise ValueError("either a prime or modulus must be passed for construction")
        
        # If we're not given a degree, construct GF(p)
        if degree is None:
            degree = 1

        # Construct the field from the prime and degree GF(p^d)
        self.set_from_order(p, degree, var, type)

    cdef _c_set_from_order(self, fmpz p, int d, char *var, fq_default_type type=fq_default_type.DEFAULT):
        self.var = var
        fq_default_ctx_init_type(self.val, p.val, d, self.var, type)
        self._initialized = True

    def set_from_order(self, p, d, var, type=fq_default_type.DEFAULT, check_prime=True):
        """
        Construct a context for the finite field GF(p^d).

        `var` is a name for the ring generator of this field over GF(p).

        The optional parameter `type` select the implementation. The
        possible values are:

        - `fq_default_ctx.DEFAULT`: implementation automatically decided by Flint (default),
        - `fq_default_ctx.FQ_ZECH`: Use `fq_zech_t`,
        - `fq_default_ctx.FQ_NMOD`: Use `fq_nmod_t`,
        - `fq_default_ctx.FQ`: Use `fq_t`.
        """
        # c_from_order expects the characteristic to be fmpz type
        prime = any_as_fmpz(p)
        if prime is NotImplemented:
            raise TypeError(f"cannot coerce {p = } to type fmpz")
        
        if check_prime and not prime.is_prime():
            raise ValueError("characteristic is not prime")

        # the degree must be strictly positive
        if d < 1:
            raise ValueError(f"the degree must be positive, got {d = }")

        if isinstance(var, str):
            var = var.encode()

        # TODO: Flint only wants one-character inputs
        if len(var) > 1:
            raise ValueError
        
        # Cython type conversion and context initalisation
        self._c_set_from_order(prime, d, var, type)

    cdef _c_set_from_modulus(self, modulus, char *var, fq_default_type type=fq_default_type.DEFAULT): 
        self.var = var
        if typecheck(modulus, fmpz_mod_poly):
            fq_default_ctx_init_modulus_type(self.val, (<fmpz_mod_poly>modulus).val,
                                             (<fmpz_mod_poly>modulus).ctx.mod.val, self.var, type)
        elif typecheck(modulus, nmod_poly):
            fq_default_ctx_init_modulus_nmod_type(self.val, (<nmod_poly>modulus).val,
                                                  self.var, type)
        else:
            raise TypeError(f"modulus must be fmpz_mod_poly or nmod_poly, got {modulus!r}")
        self._initialized = True

    def set_from_modulus(self, modulus, var, type=fq_default_type.DEFAULT, check_modulus=True):
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
        """
        if isinstance(var, str):
            var = var.encode()

        if check_modulus and not modulus.is_irreducible():
            raise ValueError("modulus must be irreducible")

        # Cython type conversion and context initalisation
        self._c_set_from_modulus(modulus, var, type)
    
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
        The extension degree of the finite field

            >>> gf = fq_default_ctx(5, 2)
            >>> gf.degree()
            2
        """
        return fq_default_ctx_degree(self.val)

    def characteristic(self):
        """
        Return the characteristic of the finite field

            >>> gf = fq_default_ctx(5, 2)
            >>> gf.characteristic()
            5
        """
        cdef fmpz p
        p = fmpz.__new__(fmpz)
        fq_default_ctx_prime(p.val, self.val)
        return p

    prime = characteristic
    
    def order(self):
        """
        Return the order of the finite field

            >>> gf = fq_default_ctx(5, 2)
            >>> gf.order()
            25
        """
        cdef fmpz q
        q = fmpz.__new__(fmpz)
        fq_default_ctx_order(q.val, self.val)
        return q

    def multiplicative_order(self):
        """"
        Return the multiplicative order of the finite field

            >>> gf = fq_default_ctx(5, 2)
            >>> gf.multiplicative_order()
            24
        """
        return self.order() - 1

    def modulus(self):
        """
        Return the modulus from the context as an fmpz_mod_poly type

            >>> gf = fq_default_ctx(5, 2)
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

            >>> gf = fq_default_ctx(5)
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

            >>> gf = fq_default_ctx(5)
            >>> gf.one()
            1
        """
        cdef fq_default res
        res = self.new_ctype_fq_default()
        res.ctx = self
        fq_default_one(res.val, self.val)
        return res

    def gen(self):
        """
        Return the one element

            >>> gf = fq_default_ctx(5, 2, var="w")
            >>> gf.gen()
            w
        """
        cdef fq_default res
        res = self.new_ctype_fq_default()
        res.ctx = self
        fq_default_gen(res.val, self.val)
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

    cdef set_list_as_fq_default(self, fq_default_t fq_ele, obj):
        cdef fmpz_poly poly
        poly = fmpz_poly.__new__(fmpz_poly)
        fmpz_poly_set_list(poly.val, obj)
        
        # Now set the value from the fmpz_poly
        fq_default_set_fmpz_poly(fq_ele, poly.val, self.val)

        return 0

    cdef set_any_as_fq_default(self, fq_default_t fq_ele, obj):
        # Converts the list to an fmpz_poly and then sets from this
        if typecheck(obj, list):
            return self.set_list_as_fq_default(fq_ele, obj)

        # For small integers we can convert directly
        if typecheck(obj, int) and obj.bit_length() < 32:
            fq_default_set_si(fq_ele, <slong>obj, self.val)
            return 0        

        # Assumes that the modulus of the polynomial matches
        # the context for the fq_default
        if typecheck(obj, fmpz_mod_poly):
            fq_default_set_fmpz_mod_poly(fq_ele, (<fmpz_mod_poly>obj).val, self.val)
            return 0

        # Assumes that the modulus of the polynomial matches
        # the context for the fq_default
        if typecheck(obj, nmod_poly):
            fq_default_set_nmod_poly(fq_ele, (<nmod_poly>obj).val, self.val)
            return 0

        # If the input is not fmpz_mod_poly or nmod_poly or a list, we cast the
        # input to an fmpz_poly and then set from this
        poly = any_as_fmpz_poly(obj)
        if poly is NotImplemented:
            return NotImplemented

        fq_default_set_fmpz_poly(fq_ele, (<fmpz_poly>poly).val, self.val)
        return 0

    cdef any_as_fq_default(self, obj):
        # convert from fq_default
        if typecheck(obj, fq_default):
            if self != (<fq_default>obj).ctx:
                raise ValueError("contexts dont match")
            return obj

        cdef fq_default res
        res = self.new_ctype_fq_default()
        check = self.set_any_as_fq_default(res.val, obj)
        if check is NotImplemented:
            return NotImplemented        
        return res

    def __eq__(self, other):
        """
        Two finite field context compare equal if they have same
        characteristic, modulus, type and variable
        
            >>> from flint import fmpz_mod_poly_ctx
            >>> modulus = fmpz_mod_poly_ctx(5)([2,4,1])
            >>> gf = fq_default_ctx(5, 2)
            >>> gf2 = fq_default_ctx(modulus=modulus)
            >>> gf2 == gf
            True
            >>> gf3 = fq_default_ctx(modulus=modulus, var="y")
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
        if self.degree() == 1:
            return f"Context for fq_default in GF({self.prime()})"
        return f"Context for fq_default in GF({self.prime()}^{self.degree()})[{self.var.decode()}]/({self.modulus().str(var=self.var.decode())})"

    def __repr__(self):
        if self.degree() == 1:
            return f"fq_default_ctx({self.prime()}, var='{self.var.decode()}' type={self.type})"
        return f"fq_default_ctx({self.prime()}, {self.degree()}, '{self.var.decode()}', {self.modulus()!r}, {self.type})"

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
        return self.polynomial().str(var=self.ctx.var.decode())

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

    def __richcmp__(self, other, int op):
        cdef bint res
        if op != 2 and op != 3:
            raise TypeError("fq_default cannot be ordered")

        if not typecheck(other, fq_default):
            other = self.ctx.any_as_fq_default(other)

        if typecheck(other, fq_default):
            res = (self.ctx == (<fq_default>other).ctx) and \
                  fq_default_equal(self.val, (<fq_default>other).val, self.ctx.val)
            if op == 2:
                return res
            else:
                return not res
        else:
            return NotImplemented

    # =================================================
    # Generic arithmetic required by flint_scalar
    # =================================================   

    def _neg_(self):
        cdef fq_default res
        res = self.ctx.new_ctype_fq_default()
        fq_default_neg(res.val, self.val, self.ctx.val)
        return res

    def _add_(self, other):
        other = self.ctx.any_as_fq_default(other)
        if other is NotImplemented:
            return NotImplemented

        cdef fq_default res
        res = self.ctx.new_ctype_fq_default()
        fq_default_add(res.val, self.val, (<fq_default>other).val, self.ctx.val)
        return res

    @staticmethod
    def _sub_(left, right):
        # Case when left and right are already fq_default
        if typecheck(left, fq_default) and typecheck(right, fq_default):
            if not (<fq_default>left).ctx == (<fq_default>right).ctx:
                raise ValueError("moduli must match")

        # Otherwise, exactly one of these values is fq_default, if it's left
        # convert the right
        elif typecheck(left, fq_default):
            right = (<fq_default>left).ctx.any_as_fq_default(right)
            if right is NotImplemented:
                return NotImplemented
        # Otherwise convert the left
        else:
            left = (<fq_default>right).ctx.any_as_fq_default(left)
            if left is NotImplemented:
                return NotImplemented

        # Now perform subtraction
        cdef fq_default res
        res = (<fq_default>left).ctx.new_ctype_fq_default()
        fq_default_sub(res.val, (<fq_default>left).val, (<fq_default>right).val, res.ctx.val)
        return res

    def _mul_(self, other):
        other = self.ctx.any_as_fq_default(other)
        if other is NotImplemented:
            return NotImplemented

        cdef fq_default res
        res = self.ctx.new_ctype_fq_default()
        fq_default_mul(res.val, self.val, (<fq_default>other).val, self.ctx.val)
        return res

    @staticmethod
    def _div_(left, right):
        # Case when left and right are already fq_default
        if typecheck(left, fq_default) and typecheck(right, fq_default):
            if not (<fq_default>left).ctx == (<fq_default>right).ctx:
                raise ValueError("moduli must match")

        # Otherwise, exactly one of these values is fq_default, if it's left
        # convert the right
        elif typecheck(left, fq_default):
            right = (<fq_default>left).ctx.any_as_fq_default(right)
            if right is NotImplemented:
                return NotImplemented
        # Otherwise convert the left
        else:
            left = (<fq_default>right).ctx.any_as_fq_default(left)
            if left is NotImplemented:
                return NotImplemented

        # Now perform division
        cdef fq_default res
        res = (<fq_default>left).ctx.new_ctype_fq_default()
        fq_default_div(res.val, (<fq_default>left).val, (<fq_default>right).val, res.ctx.val)
        return res

    def _invert_(self):
        cdef fq_default res
        res = self.ctx.new_ctype_fq_default()
        fq_default_inv(res.val, self.val, self.ctx.val)
        return res

    # =================================================
    # Additional arithmetic
    # ================================================= 

    def square(self):
        """
        Computes the square of ``self``
        """
        cdef fq_default res
        res = self.ctx.new_ctype_fq_default()
        fq_default_sqr(res.val, self.val, self.ctx.val)
        return res

    def __pow__(self, e):
        """
        """
        cdef fq_default res
        res = self.ctx.new_ctype_fq_default()

        # If e is negative, we invert first
        if e < 0:
            if self.is_zero():
                raise ZeroDivisionError
            e = -e
            fq_default_inv(res.val, self.val, self.ctx.val)
        else:
            fq_default_set(res.val, self.val, self.ctx.val)

        # For small e we dont need to make an fmpz
        if e.bit_length() < 32:
            fq_default_pow_ui(res.val, res.val, <ulong>e, self.ctx.val)
            return res

        # Attempt to cast the exponent to an fmpz type then exponentiate
        e_fmpz = any_as_fmpz(e)
        if e_fmpz is NotImplemented:
            raise TypeError(f"exponent {e = } cannot be cast to fmpz")

        fq_default_pow(res.val, res.val, (<fmpz>e_fmpz).val, self.ctx.val)

        return res

    def sqrt(self):
        """
        Returns the square root of the element, if not square root exists,
        throws a ValueError.
        """
        cdef fq_default res
        res = self.ctx.new_ctype_fq_default()
        check = fq_default_sqrt(res.val, self.val, self.ctx.val)
        if check:
            return res
        raise ValueError("element is not a square")

    def is_square(self):
        """
        Returns if the element is a square in the field
        """
        return 1 == fq_default_is_square(self.val, self.ctx.val)

    def pth_root(self):
        """
        Returns the pth root of the element.
        
        This is computed by  raising ``self`` to the `p^(d-1)` power, 
        `p` is the characteristic of the field and `d` is the degree 
        of the extension.
        """
        cdef fq_default res
        res = self.ctx.new_ctype_fq_default()
        fq_default_pth_root(res.val, self.val, self.ctx.val)
        return res

    # =================================================
    # Special functions
    # ================================================= 

    def trace(self):
        """
        Returns the trace of self
        """
        cdef fmpz tr = fmpz.__new__(fmpz)
        fq_default_trace(tr.val, self.val, self.ctx.val)
        return tr

    def norm(self):
        """
        Returns the norm of self
        """
        cdef fmpz nrm = fmpz.__new__(fmpz)
        fq_default_norm(nrm.val, self.val, self.ctx.val)
        return nrm

    def frobenius(self, e):
        """
        Evaluates the homomorphism `\Sigma^e` on ``self``.
        """
        cdef fq_default res
        res = self.ctx.new_ctype_fq_default()
        fq_default_frobenius(res.val, self.val, <slong>e, self.ctx.val)
        return res
