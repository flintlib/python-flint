from cpython.list cimport PyList_GET_SIZE
from flint.flintlib.fmpz_mod_poly cimport *
from flint.flintlib.fmpz_mod_poly_factor cimport *

from flint.flintlib.fmpz cimport(
    fmpz_init,
    fmpz_clear,
    fmpz_is_one
)
from flint.types.fmpz cimport fmpz, any_as_fmpz
from flint.types.fmpz_mod cimport fmpz_mod_ctx, fmpz_mod
from flint.types.fmpz_poly cimport fmpz_poly

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

    def is_prime(self):
        """
        Return whether the modulus is prime

            >>> fmpz_mod_poly_ctx(2**127).is_prime()
            False
            >>> fmpz_mod_poly_ctx(2**127 - 1).is_prime()
            True
        """
        return self.mod.is_prime()

    def zero(self):
        """
        Return the zero element of this polynomial ring

            >>> R = fmpz_mod_poly_ctx(163)
            >>> R.zero()
            0
        """
        cdef fmpz_mod_poly res
        res = fmpz_mod_poly.__new__(fmpz_mod_poly)
        res.ctx = self
        fmpz_mod_poly_zero(res.val, self.mod.val)
        return res

    def one(self):
        """
        Return the one element of this polynomial ring

            >>> R = fmpz_mod_poly_ctx(163)
            >>> R.one()
            1
        """
        cdef fmpz_mod_poly res
        res = fmpz_mod_poly.__new__(fmpz_mod_poly)
        res.ctx = self
        fmpz_mod_poly_one(res.val, self.mod.val)
        return res

    def gen(self):
        """
        Return the generator of the polynomial `x`

            >>> R = fmpz_mod_poly_ctx(2**127 - 1)
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
            if self != (<fmpz_mod_poly>obj).ctx:
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
        return hash(self.modulus())

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

    # ---------------- #
    #    Arithmetic    #
    # ---------------- #

    def __pos__(self):
        return self

    def __neg__(self):
        cdef fmpz_mod_poly res
        res = fmpz_mod_poly.__new__(fmpz_mod_poly)
        res.ctx = self.ctx
        fmpz_mod_poly_neg(res.val, self.val, self.ctx.mod.val)
        return res

    def __add__(self, other):
        cdef fmpz_mod_poly res
        other = self.ctx.any_as_fmpz_mod_poly(other)
        if other is NotImplemented:
            return other

        res = fmpz_mod_poly.__new__(fmpz_mod_poly)
        fmpz_mod_poly_add(
            res.val, self.val, (<fmpz_mod_poly>other).val, self.ctx.mod.val
        )
        res.ctx = self.ctx
        return res

    def __radd__(self, other):
        return self.__add__(other)

    @staticmethod
    def _sub_(left, right):
        cdef fmpz_mod_poly res
        res = fmpz_mod_poly.__new__(fmpz_mod_poly)

        # Case when left and right are already fmpz_mod_poly
        if typecheck(left, fmpz_mod_poly) and typecheck(right, fmpz_mod_poly):
            if not (<fmpz_mod_poly>left).ctx == (<fmpz_mod_poly>right).ctx:
                raise ValueError("moduli must match")

        # Case when right is not fmpz_mod_poly, try to convert to fmpz
        elif typecheck(left, fmpz_mod_poly):
            right = (<fmpz_mod_poly>left).ctx.any_as_fmpz_mod_poly(right)
            if right is NotImplemented:
                return NotImplemented

        # Case when left is not fmpz_mod_poly, try to convert to fmpz
        else:
            left = (<fmpz_mod_poly>right).ctx.any_as_fmpz_mod_poly(left)
            if left is NotImplemented:
                return NotImplemented

        res.ctx = (<fmpz_mod_poly>left).ctx
        fmpz_mod_poly_sub(
                res.val, (<fmpz_mod_poly>left).val, (<fmpz_mod_poly>right).val, res.ctx.mod.val
        )
        return res

    def __sub__(s, t):
        return fmpz_mod_poly._sub_(s, t)

    def __rsub__(s, t):
        return fmpz_mod_poly._sub_(t, s)

    def __mul__(self, other):
        # TODO:
        # Allow scalar multiplication for efficiency, rather 
        # than casting `other` to a polynomial?
        cdef fmpz_mod_poly res
        other = self.ctx.any_as_fmpz_mod_poly(other)
        if other is NotImplemented:
            return other

        res = fmpz_mod_poly.__new__(fmpz_mod_poly)
        fmpz_mod_poly_mul(
            res.val, self.val, (<fmpz_mod_poly>other).val, self.ctx.mod.val
        )
        res.ctx = self.ctx
        return res

    def __rmul__(self, other):
        return self.__mul__(other)

    @staticmethod
    def _div_(left, right):
        # TODO:
        # Allow scalar division for efficiency, rather 
        # than casting `other` to a polynomial?
        cdef bint check
        cdef fmpz_mod_poly res

        # Case when left and right are already fmpz_mod_poly
        if typecheck(left, fmpz_mod_poly) and typecheck(right, fmpz_mod_poly):
            if not (<fmpz_mod_poly>left).ctx == (<fmpz_mod_poly>right).ctx:
                raise ValueError("moduli must match")

        # Case when right is not fmpz_mod_poly, try to convert to fmpz
        elif typecheck(left, fmpz_mod_poly):
            right = (<fmpz_mod_poly>left).ctx.any_as_fmpz_mod_poly(right)
            if right is NotImplemented:
                return NotImplemented

        # Case when left is not fmpz_mod_poly, try to convert to fmpz
        else:
            left = (<fmpz_mod_poly>right).ctx.any_as_fmpz_mod_poly(left)
            if left is NotImplemented:
                return NotImplemented

        res = fmpz_mod_poly.__new__(fmpz_mod_poly)
        res.ctx = (<fmpz_mod_poly>left).ctx
        check = fmpz_mod_poly_divides(
            res.val, (<fmpz_mod_poly>left).val, (<fmpz_mod_poly>right).val, res.ctx.mod.val
        )
        if check == 0:
            raise ValueError(
                f"{right} does not divide {left}"
            )

        return res

    def __truediv__(s, t):
        return fmpz_mod_poly._div_(s, t)

    def __rtruediv__(s, t):
        return fmpz_mod_poly._div_(t, s)

    def __floordiv__(self, other):
        return NotImplemented

    def __pow__(self, e):
        cdef fmpz_mod_poly res
        if e < 0:
            raise ValueError("Exponent must be non-negative")

        res = fmpz_mod_poly.__new__(fmpz_mod_poly)
        res.ctx = self.ctx
        fmpz_mod_poly_pow(
            res.val, self.val, (<ulong>e), self.ctx.mod.val
        )
        return res

    @staticmethod
    def _mod_(s, t):
        pass

    def __mod__(s, t):
        return fmpz_mod_poly._mod_(s, t)

    def __rmod__(s, t):
        return fmpz_mod_poly._mod_(t, s)    

    # =
    # Other Magic Methods
    # =

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

    def __hash__(self):
        return hash(map(int, self.coeffs()))

    cpdef long length(self):
        return fmpz_mod_poly_length(self.val, self.ctx.mod.val)

    cpdef long degree(self):
        return fmpz_mod_poly_degree(self.val, self.ctx.mod.val)

    def is_zero(self):
        """
        Return `True` if the polynomial is the zero polynomial
        and `False` otherwise

            >>> R = fmpz_mod_poly_ctx(163)
            >>> f = R(0)
            >>> f.is_zero()
            True
        """
        return 0 != fmpz_mod_poly_is_zero(self.val, self.ctx.mod.val)
    
    def is_one(self):
        """
        Return `True` if the polynomial is the zero polynomial
        and `False` otherwise

            >>> R = fmpz_mod_poly_ctx(163)
            >>> f = R(1)
            >>> f.is_one()
            True
        """
        return 0 != fmpz_mod_poly_is_one(self.val, self.ctx.mod.val)
    
    def is_gen(self):
        """
        Return `True` if the polynomial is the zero polynomial
        and `False` otherwise

            >>> R = fmpz_mod_poly_ctx(163)
            >>> f = R([0,1])
            >>> f.is_gen()
            True
        """
        return 0 != fmpz_mod_poly_is_gen(self.val, self.ctx.mod.val)

    def is_constant(self):
        """
        Return True if this is a constant polynomial.

            >>> R = fmpz_mod_poly_ctx(163)
            >>> x = R.gen()
            >>> x.is_constant()
            False
            >>> R(123).is_constant()
            True
        """
        return self.degree() <= 0

    def constant_coefficient(self):
        """
        Return the leading coefficient of this polynomial.

            >>> R = fmpz_mod_poly_ctx(163)
            >>> f = R([1,2,3])
            >>> f.leading_coefficient()
            fmpz_mod(1, 163)
        """
        return self[0]

    def leading_coefficient(self):
        """
        Return the leading coefficient of this polynomial.

            >>> R = fmpz_mod_poly_ctx(163)
            >>> f = R([1,2,3])
            >>> f.leading_coefficient()
            fmpz_mod(3, 163)
        """
        return self[self.degree()]

    def reverse(self, degree=None):
        """
        Return a polynomial with the coefficients of this polynomial
        reversed.

        If `degree` is not None, the output polynomial will be zero-padded
        or truncated before being reversed. Note: degree must be non-negative.

            >>> R = fmpz_mod_poly_ctx(163)
            >>> f = R([1,2,3,4,5])
            >>> f.reverse()
            x^4 + 2*x^3 + 3*x^2 + 4*x + 5
            >>> f.reverse(degree=1)
            x + 2
            >>> f.reverse(degree=100)
            x^100 + 2*x^99 + 3*x^98 + 4*x^97 + 5*x^96
        """
        cdef fmpz_mod_poly res
        cdef slong d

        res =  fmpz_mod_poly.__new__(fmpz_mod_poly)
        res.ctx = self.ctx

        if degree is not None:
            d = degree
            if d != degree or d < 0:
                raise ValueError(f"degree argument must be a non-negative integer, got {degree}")
        else:
            d = fmpz_mod_poly_degree(self.val, self.ctx.mod.val)

        length = d + 1
        fmpz_mod_poly_reverse(res.val, self.val, length, self.ctx.mod.val)
        return res

    def monic(self, check=True):
        """
        Return this polynomial divided by its leading coefficient.

        If `check` is True, raises ValueError if the leading coefficient 
        is not invertible modulo N.

            >>> R = fmpz_mod_poly_ctx(163)
            >>> f = R([1,2,3])
            >>> f.monic()
            x^2 + 55*x + 109
        """
        cdef fmpz_mod_poly res
        cdef fmpz_t f

        res =  fmpz_mod_poly.__new__(fmpz_mod_poly)
        if not check:
            fmpz_mod_poly_make_monic(
                res.val, self.val, self.ctx.mod.val
            )
        else:
            fmpz_init(f) 
            fmpz_mod_poly_make_monic_f(
                f, res.val, self.val, self.ctx.mod.val
            )
            if not fmpz_is_one(f):
                raise ValueError(f"Leading coefficient is not invertible")
        res.ctx = self.ctx
        return res

    def is_irreducible(self):
        """
        Return whether this polynomial is irreducible.

            >>> R = fmpz_mod_poly_ctx(163)
            >>> f = R([1,2,3])
            >>> f = (x**2 + 5*x + 3)
            >>> f.is_irreducible()
            True
            >>> f = (x**2 + x + 3)
            >>> f.is_irreducible()
            False
        """
        return 1 == fmpz_mod_poly_is_irreducible(self.val, self.ctx.mod.val)

    def is_squarefree(self):
        """
        Return whether this polynomial is squarefree.

            >>> R = fmpz_mod_poly_ctx(163)
            >>> x = R.gen()
            >>> f = (x + 1)**2 * (x + 3)
            >>> f.is_squarefree()
            False
            >>> f = (x + 1) * (x + 3)
            >>> f.is_squarefree()
            True

        """
        return 1 == fmpz_mod_poly_is_squarefree(self.val, self.ctx.mod.val)

    def gcd(self, other):
        """
        Return the greatest common divisor of self and other.

            >>> R = fmpz_mod_poly_ctx(163)
            >>> x = R.gen()
            >>> f = x*(x + 1)
            >>> f.gcd(x+1)
            x + 1
            >>> f.gcd(x*x)
            x

        """
        cdef fmpz_mod_poly res

        if not self.ctx.is_prime():
            raise NotImplementedError("gcd algorithm assumes that the modulus is prime")
                 
        other = self.ctx.any_as_fmpz_mod_poly(other)
        if other is NotImplemented:
            return other

        res = fmpz_mod_poly.__new__(fmpz_mod_poly)
        res.ctx = self.ctx
        fmpz_mod_poly_gcd(
            res.val, self.val, (<fmpz_mod_poly>other).val, self.ctx.mod.val
        )
        return res

    def derivative(self):
        """
        The formal derivative of this polynomial

            >>> R = fmpz_mod_poly_ctx(163)
            >>> x = R.gen()
            >>> f = 111*x**4 + 58*x**3 + 98*x**2 + 117*x + 7
            >>> f.derivative()
            118*x^3 + 11*x^2 + 33*x + 117

        """
        cdef fmpz_mod_poly res
        res = fmpz_mod_poly.__new__(fmpz_mod_poly)
        res.ctx = self.ctx
        fmpz_mod_poly_derivative(
            res.val, self.val, self.ctx.mod.val
        )
        return res

    def discriminant(self):
        """
        Return the discriminant of self.

            >>> R = fmpz_mod_poly_ctx(163)
            >>> x = R.gen()
            >>> f = 6*x**4 + 7*x**3 + 7*x**2 + 8*x + 6
            >>> f.discriminant()
            fmpz_mod(50, 163)

        """
        cdef fmpz_mod res

        if not self.ctx.is_prime():
            raise NotImplementedError("discriminant algorithm assumes that the base is a field")

        res = fmpz_mod.__new__(fmpz_mod)
        res.ctx = self.ctx.mod
        fmpz_mod_poly_discriminant(
            res.val, self.val, self.ctx.mod.val
        )
        return res

    def radical(self):
        """
        Return the radical of self, the product of the irreducible
        factors of the polynomial. This is also referred to as the
        square-free part of the polynomial.

            >>> R = fmpz_mod_poly_ctx(163)
            >>> x = R.gen()
            >>> f = (x + 1)**3 * (x + 2)
            >>> f.radical()
            x^2 + 3*x + 2

        """
        if not self.ctx.is_prime():
            raise NotImplementedError("radical algorithm assumes that the base is a field")

        return self / self.gcd(self.derivative())    

    # TODO: we could make a factorisation class which we could then
    # implement the factor methods such as pow and concat. I think
    # sage does something like this with `Factorisation` classes.

    def factor_squarefree(self):
        """
        Factors self into irreducible factors, returning a tuple
        (c, factors) where c is the content of the coefficients and
        factors is a list of (poly, exp) pairs.

            >>> R = fmpz_mod_poly_ctx(163)
            >>> x = R.gen()
            >>> f = (x + 1) * (x + 2)
            >>> f.factor_squarefree()
            (fmpz_mod(1, 163), [(x^2 + 3*x + 2, 1)])
            >>> f = (x + 1) * (x + 2)**5
            >>> f.factor_squarefree()
            (fmpz_mod(1, 163), [(x + 1, 1), (x + 2, 5)])
        """
        cdef fmpz_mod_poly_factor_t fac
        cdef int i

        if not self.ctx.is_prime():
            raise NotImplementedError("factor_squarefree algorithm assumes that the modulus is prime")

        fmpz_mod_poly_factor_init(fac, self.ctx.mod.val)
        fmpz_mod_poly_factor_squarefree(fac, self.val, self.ctx.mod.val)

        res = [0] * fac.num

        cdef fmpz_mod_poly u
        for i in range(fac.num):
            u = fmpz_mod_poly.__new__(fmpz_mod_poly)
            u.ctx = self.ctx
            fmpz_mod_poly_set(u.val, &fac.poly[i], self.ctx.mod.val)
            exp = fac.exp[i]
            res[i] = (u, exp)
        return self.leading_coefficient(), res

    def factor(self, algorithm=None):
        """
        Factors self into irreducible factors, returning a tuple
        (c, factors) where c is the content of the coefficients and
        factors is a list of (poly, exp) pairs.

            >>> R = fmpz_mod_poly_ctx(163)
            >>> x = R.gen()
            >>> f = 6*x**4 + 7*x**3 + 7*x**2 + 8*x + 6
            >>> f.factor()
            (fmpz_mod(6, 163), [(x^4 + 137*x^3 + 137*x^2 + 110*x + 1, 1)])
            >>> f = (x + 1)**3 * (x + 2)
            >>> f.factor()
            (fmpz_mod(1, 163), [(x + 1, 3), (x + 2, 1)])
        """
        cdef fmpz_mod_poly_factor_t fac
        cdef int i

        if not self.ctx.is_prime():
            raise NotImplementedError("factor algorithm assumes that the modulus is prime")
                 
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
        return self.leading_coefficient(), res

    def roots(self):
        return NotImplemented

    def complex_roots(self):
        return NotImplemented
