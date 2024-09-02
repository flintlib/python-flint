from flint.pyflint cimport global_random_state
from flint.types.fmpz cimport fmpz, any_as_fmpz
from flint.types.nmod cimport nmod
from flint.types.fmpz_mod cimport fmpz_mod
from flint.types.fmpz_poly cimport fmpz_poly, any_as_fmpz_poly, fmpz_poly_set_list
from flint.types.fmpz_mod_poly cimport fmpz_mod_poly, fmpz_mod_poly_ctx
from flint.types.nmod_poly cimport nmod_poly
from flint.utils.typecheck cimport typecheck

from flint.utils.flint_exceptions import DomainError

# Allow the type to be denoted by strings or integers
FQ_TYPES = {
    "FQ_ZECH" : 1,
    "FQ_NMOD" : 2,
    "FQ" : 3,
    "NMOD" : 4,
    "FMPZ_MOD" : 5
}

cdef class fq_default_ctx:
    r"""
    Context object for creating :class:`~.fq_default`.

    Finite fields can be initialized in one of two possible ways. The
    first is by providing characteristic and degree:

        >>> fq_default_ctx(5, 2, 'y', fq_type='FQ_ZECH')
        fq_default_ctx(5, 2, 'y', x^2 + 4*x + 2, 'FQ_ZECH')

    The second is by giving an irreducible polynomial of type
    :class:`~.nmod_poly` or :class:`~.fmpz_mod_poly`:

        >>> from flint import fmpz_mod_poly_ctx
        >>> mod = fmpz_mod_poly_ctx(11)([1,0,1])
        >>> fq_default_ctx(modulus=mod, fq_type='FQ_NMOD')
        fq_default_ctx(11, 2, 'z', x^2 + 1, 'FQ_NMOD')

    For more details, see the documentation of :method:`~._set_from_order`
    and :method:`~._set_from_modulus`.
    """
    def __cinit__(self):
        pass

    def __dealloc__(self):
        if self._initialized:
            fq_default_ctx_clear(self.val)
            self._initialized = False

    @staticmethod
    def _parse_input_fq_type(fq_type):
        if typecheck(fq_type, str):
            fq_type = FQ_TYPES.get(fq_type, None)
            if fq_type is None:
                raise ValueError("invalid fq_type string")

        # Now fq_type should be an int between 0, 5
        if not typecheck(fq_type, int):
            raise TypeError(f"fq_type = {fq_type} is invalid")
        if fq_type < 0 or fq_type > 5:
            raise ValueError(f"fq_type = {fq_type} should be between 0 and 5")

        return fq_type

    @staticmethod
    def _parse_input_var(var):
        # If no variable is given, use x
        if var is None:
            var = b"z"

        # Encode to bytes for cython to parse
        if isinstance(var, str):
            var = var.encode()

        # TODO: Flint only wants one-character inputs
        if len(var) > 1:
            raise ValueError("variable for GF(p^k) generator can only be one character")

        return var

    def __init__(self, p=None, degree=None, var=None, modulus=None, fq_type=fq_default_type.DEFAULT,
                 check_prime=True, check_modulus=True):
        # Ensure the var used for the generator of GF(p^d) is a single byte
        var = self._parse_input_var(var)
        self.var = var

        # Ensure the fq_type is an integer between 0, 5 -- we allow users to
        # input a string which is converted as an enum
        fq_type = self._parse_input_fq_type(fq_type)

        # If a modulus is given, attempt to construct from this
        if modulus is not None:
            # If the polynomial has no known characteristic, we can try and create one
            # using the supplied prime
            if not typecheck(modulus, fmpz_mod_poly):
                if p is None:
                    raise ValueError("cannot create from modulus if no characteristic is known")
                ring_ctx = fmpz_mod_poly_ctx(p)
                if not ring_ctx.is_prime():
                    raise ValueError("characteristic is not prime")
                modulus = ring_ctx.any_as_fmpz_mod_poly(modulus)
                if modulus is NotImplemented:
                    raise TypeError("modulus cannot be cast to fmpz_mod_poly")

            self._set_from_modulus(modulus, var, fq_type, check_prime=check_prime, check_modulus=check_modulus)
            return

        # If there's no modulus and no prime, we can't continue
        if p is None:
            raise ValueError("either a prime or modulus must be passed for construction")

        # If we're not given a degree, construct GF(p)
        if degree is None:
            degree = 1

        # Construct the field from the prime and degree GF(p^d)
        self._set_from_order(p, degree, var, fq_type, check_prime=check_prime)

    cdef _set_from_order(self, p, d, var, fq_type=fq_default_type.DEFAULT, check_prime=True):
        """
        Construct a context for the finite field GF(p^d).

        `var` is a name for the ring generator of this field over GF(p).

        The optional parameter `type` select the implementation. For more
        information about the types available, see :class:`~.fq_default_type`
        for possible types.
        """
        # c_from_order expects the characteristic to be fmpz type
        prime = any_as_fmpz(p)
        if prime is NotImplemented:
            raise TypeError(f"cannot coerce p = {p} to type fmpz")

        if check_prime and not prime.is_prime():
            raise ValueError("characteristic is not prime")

        # the degree must be strictly positive
        if d < 1:
            raise ValueError(f"the degree must be positive, got d = {d}")

        fq_default_ctx_init_type(self.val, (<fmpz>prime).val, d, self.var, <fq_default_type>fq_type)
        self._initialized = True

    cdef _set_from_modulus(self, modulus, var, fq_type=fq_default_type.DEFAULT,
                           check_prime=True, check_modulus=True):
        """
        Construct a context for a finite field from an irreducible polynomial.

        `modulus` may be of type :class:`~.fmpz_mod_poly` or :class:`~.nmod_poly`.

        `var` is a name for the ring generator of this field over the prime field.

        The optional parameter `type` select the implementation. For more
        information about the types available, see :class:`~.fq_default_type`
        for possible types.
        """
        if check_prime and not (<fmpz_mod_poly>modulus).ctx.is_prime():
            raise ValueError("characteristic is not prime")

        if check_modulus and not modulus.is_irreducible():
            raise ValueError("modulus must be irreducible")

        fq_default_ctx_init_modulus_type(self.val, (<fmpz_mod_poly>modulus).val,
                                         (<fmpz_mod_poly>modulus).ctx.mod.val, self.var, <fq_default_type>fq_type)
        self._initialized = True

    @property
    def fq_type(self):
        """
        Return the implementation of this context
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
        Return the unit element

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

            >>> gf = fq_default_ctx(163, 3)
            >>> a = gf.random_element()
            >>> type(a) is fq_default
            True
            >>> a = gf.random_element(not_zero=True)
            >>> not a.is_zero()
            True
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

    cdef set_any_scalar_as_fq_default(self, fq_default_t fq_ele, obj):
        cdef slong i
        if typecheck(obj, int):
            # For small integers we can convert directly
            try:
                i = obj
                fq_default_set_si(fq_ele, i, self.val)
                return 0
            # For larger integers fall through to conversion to fmpz
            except OverflowError:
                pass

            obj_fmpz = any_as_fmpz(obj)
            fq_default_set_fmpz(fq_ele, (<fmpz>obj_fmpz).val, self.val)
            return 0

        # For fmpz we can also convert directly
        if typecheck(obj, fmpz):
            fq_default_set_fmpz(fq_ele, (<fmpz>obj).val, self.val)
            return 0

        # For nmod we can convert by taking it as an int
        if typecheck(obj, nmod) and self.prime() == obj.modulus():
            fq_default_set_ui(fq_ele, <ulong>(<nmod>obj).val, self.val)
            return 0

        # For fmpz_mod we can also convert directly
        if typecheck(obj, fmpz_mod) and self.prime() == (<fmpz_mod>obj).ctx.modulus():
            fq_default_set_fmpz(fq_ele, (<fmpz_mod>obj).val, self.val)
            return 0

        # Otherwise the object wasn't a scalar we convert from
        return NotImplemented

    cdef set_any_as_fq_default(self, fq_default_t fq_ele, obj):
        # First try and convert from scalars
        check = self.set_any_scalar_as_fq_default(fq_ele, obj)
        if check is not NotImplemented:
            return 0

        if typecheck(obj, fmpz_mod_poly) and self.prime() == (<fmpz_mod_poly>obj).ctx.mod.modulus():
            fq_default_set_fmpz_mod_poly(fq_ele, (<fmpz_mod_poly>obj).val, self.val)
            return 0

        if typecheck(obj, nmod_poly) and self.prime() == obj.modulus():
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
                raise ValueError("fields must match")
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
            return (self.fq_type == other.fq_type
                    and self.var == other.var
                    and self.prime() == other.prime()
                    and self.modulus() == other.modulus())
        return False

    def __hash__(self):
        return hash((self.fq_type, self.var, self.prime(), self.modulus()))

    def __str__(self):
        if self.degree() == 1:
            return f"Context for fq_default in GF({self.prime()})"
        return f"Context for fq_default in GF({self.prime()}^{self.degree()})[{self.var.decode()}]/({self.modulus().str(var=self.var.decode())})"

    def __repr__(self):
        if self.degree() == 1:
            return f"fq_default_ctx({self.prime()}, var='{self.var.decode()}' type='{self.fq_type._name_}')"
        return f"fq_default_ctx({self.prime()}, {self.degree()}, '{self.var.decode()}', {self.modulus()!r}, '{self.fq_type._name_}')"

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

        # Converts the list to an fmpz_poly and then sets from this
        if typecheck(val, list):
            self.ctx.set_list_as_fq_default(self.val, val)
            return

        # Otherwise cascades through types to convert
        check = self.ctx.set_any_as_fq_default(self.val, val)
        if check is NotImplemented:
            raise TypeError

    def __int__(self):
        """
        Attempts to lift self to an integer of type fmpz in [0, p-1]

            >>> gf = fq_default_ctx(163, 3)
            >>> a = gf([1,2,3])
            >>> int(gf(123))
            123
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

            >>> gf = fq_default_ctx(163, 3)
            >>> a = gf([1,2,3])
            >>> a.polynomial()
            3*x^2 + 2*x + 1
            >>> gf(123).polynomial()
            123
        """
        cdef fmpz_mod_poly pol

        ring_ctx = fmpz_mod_poly_ctx(self.ctx.prime())
        pol = ring_ctx.new_ctype_poly()
        fq_default_get_fmpz_mod_poly((<fmpz_mod_poly>pol).val, self.val, self.ctx.val)

        return pol

    def to_list(self):
        """
        Returns self as a list of fmpz types corresponding to a
        list of coefficients

            >>> gf = fq_default_ctx(163, 3)
            >>> gf([-1,2,1]).to_list()
            [162, 2, 1]
            >>> gf.one().to_list()
            [1, 0, 0]
        """
        coeffs = [None for _ in range(self.ctx.degree())]
        for i in range(self.ctx.degree()):
            c = fmpz.__new__(fmpz)
            fq_default_get_coeff_fmpz((<fmpz>c).val, self.val, <slong>i, self.ctx.val)
            coeffs[i] = c
        return coeffs

    def str(self):
        return self.polynomial().str(var=self.ctx.var.decode())

    def __hash__(self):
        return hash((self.polynomial(), hash(self.ctx)))

    # =================================================
    # Comparisons
    # =================================================
    def is_zero(self):
        """
        Returns true is self is zero and false otherwise

            >>> gf = fq_default_ctx(163, 3)
            >>> gf(0).is_zero()
            True
            >>> gf(-1).is_zero()
            False
        """
        return 1 == fq_default_is_zero(self.val, self.ctx.val)

    def is_one(self):
        """
        Returns true is self is one and false otherwise

            >>> gf = fq_default_ctx(163, 3)
            >>> gf(-1).is_one()
            False
            >>> gf(1).is_one()
            True
        """
        return 1 == fq_default_is_one(self.val, self.ctx.val)

    def __richcmp__(self, other, int op):
        cdef bint res
        if op != 2 and op != 3:
            raise TypeError("fq_default cannot be ordered")

        # If other is not an fq_default element, we attempt to convert to fq_default
        if not typecheck(other, fq_default):
            # For nmod and fmpz_mod if the modulus does not match the characteristic
            # then we return false.
            if typecheck(other, nmod) and self.ctx.characteristic() != (<nmod>other).modulus():
                res = False

            elif typecheck(other, fmpz_mod) and self.ctx.characteristic() != (<fmpz_mod>other).ctx.modulus():
                res = False

            else:
                # Convert from int, fmpz, fmpz_mod and nmod to fq_default
                cmp = self.ctx.new_ctype_fq_default()
                check = self.ctx.set_any_scalar_as_fq_default((<fq_default>cmp).val, other)

                # Conversion failed, element was not a compatible scalar
                if check is NotImplemented:
                    return NotImplemented

                # We now have an fq_default element with the same context, so compare value only
                res = fq_default_equal(self.val, (<fq_default>cmp).val, self.ctx.val)

            # Flip the result of res if we're doing not equals
            if op == 2:
                return res
            return not res

        # Otherwise we're in the case where other is also an fq_default but may have a
        # different context.

        # If the contexts match exactly, only check values
        if self.ctx == (<fq_default>other).ctx:
            res = fq_default_equal(self.val, (<fq_default>other).val, self.ctx.val)

        # Otherwise if both contexts have the same characteristic check
        # if both fq_default lift to the same integer
        elif self.ctx.characteristic() == (<fq_default>other).ctx.characteristic():
            try:
                res = int(self) == int(other)
            # One or both lifts failed, so values are not equal
            except ValueError:
                res = False

        # Otherwise, values are considered not equal
        else:
            res = False

        # Flip the result of res if we're doing not equals
        if op == 2:
            return res
        return not res

    # =================================================
    # Generic arithmetic required by flint_scalar
    # =================================================

    def _any_as_self(self, other):
        return self.ctx.any_as_fq_default(other)

    cpdef fq_default _neg_(fq_default self):
        cdef fq_default res = self.ctx.new_ctype_fq_default()
        fq_default_neg(res.val, self.val, self.ctx.val)
        return res

    cpdef fq_default _add_(fq_default self, fq_default other):
        """
        Assumes that __add__() has ensured other is of type self
        """
        cdef fq_default res = self.ctx.new_ctype_fq_default()
        fq_default_add(res.val, self.val, other.val, self.ctx.val)
        return res

    cpdef fq_default _sub_(fq_default self, fq_default other):
        """
        Assumes that __sub__() has ensured other is of type self
        """
        cdef fq_default res = self.ctx.new_ctype_fq_default()
        fq_default_sub(res.val, self.val, other.val, res.ctx.val)
        return res

    cpdef fq_default _rsub_(fq_default self, fq_default other):
        """
        Assumes that __rsub__() has ensured other is of type self
        """
        cdef fq_default res = self.ctx.new_ctype_fq_default()
        fq_default_sub(res.val, other.val, self.val, res.ctx.val)
        return res

    cpdef fq_default _mul_(fq_default self, fq_default other):
        """
        Assumes that __mul__() has ensured other is of type self

        TODO: this could be optimised by using mul_si and others
        """
        cdef fq_default res = self.ctx.new_ctype_fq_default()
        fq_default_mul(res.val, self.val, other.val, self.ctx.val)
        return res

    cpdef fq_default _div_(fq_default self, fq_default other):
        """
        Assumes that __div__() has ensured other is of type self
        """
        cdef fq_default res = self.ctx.new_ctype_fq_default()
        fq_default_div(res.val, self.val, other.val, res.ctx.val)
        return res

    cpdef fq_default _rdiv_(fq_default self, fq_default other):
        """
        Assumes that __div__() has ensured other is of type self
        """
        cdef fq_default res = self.ctx.new_ctype_fq_default()
        fq_default_div(res.val, other.val, self.val, res.ctx.val)
        return res

    cpdef fq_default _invert_(fq_default self):
        cdef fq_default res = self.ctx.new_ctype_fq_default()
        fq_default_inv(res.val, self.val, self.ctx.val)
        return res

    # =================================================
    # Additional arithmetic
    # =================================================

    def inverse(self):
        """
        Computes the inverse of self

            >>> gf = fq_default_ctx(163, 3)
            >>> a = gf([1,2,3])
            >>> b = a.inverse()
            >>> b
            68*z^2 + 17*z + 116
            >>> a*b == gf.one()
            True
        """
        return self._invert_()

    def square(self):
        """
        Computes the square of ``self``

            >>> gf = fq_default_ctx(163, 3)
            >>> a = gf([1,2,3])
            >>> a.square() == a*a
            True
            >>> a.square()
            110*z^2 + 101*z + 25
        """
        cdef fq_default res
        res = self.ctx.new_ctype_fq_default()
        fq_default_sqr(res.val, self.val, self.ctx.val)
        return res

    def __pow__(self, e):
        """
        Compute `a^e` for `a` equal to ``self``.

            >>> gf = fq_default_ctx(163, 3)
            >>> a = gf([1,2,3])
            >>> pow(a, -1) == 1/a
            True
            >>> pow(a, 2) == a * a
            True
            >>> pow(a, 2**128) == pow(a, 2**128 % (163**3 - 1))
            True
            >>> pow(a, 123)
            46*z^2 + 110*z + 155

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

        # Cast the exponent to fmpz type
        e_fmpz = any_as_fmpz(e)
        fq_default_pow(res.val, res.val, (<fmpz>e_fmpz).val, self.ctx.val)

        return res

    def sqrt(self):
        """
        Returns the square root of the element, if not square root exists,
        throws a ValueError.

            >>> gf = fq_default_ctx(163, 3)
            >>> a = gf([3,2,1])
            >>> a.is_square()
            True
            >>> b = a.sqrt()
            >>> b
            95*z^2 + 36*z + 34
            >>> b**2 in [a, -a]
            True
        """
        cdef fq_default res
        res = self.ctx.new_ctype_fq_default()
        check = fq_default_sqrt(res.val, self.val, self.ctx.val)
        if check:
            return res
        raise DomainError("element is not a square")

    def is_square(self):
        """
        Returns if the element is a square in the field

            >>> gf = fq_default_ctx(163, 3)
            >>> a = gf([1,2,3])
            >>> a.is_square()
            False
            >>> (a*a).is_square()
            True
        """
        return 1 == fq_default_is_square(self.val, self.ctx.val)

    def pth_root(self):
        """
        Returns the pth root of the element.

        This is computed by  raising ``self`` to the `p^(d-1)` power,
        `p` is the characteristic of the field and `d` is the degree
        of the extension.

            >>> gf = fq_default_ctx(163, 3)
            >>> a = gf([1,2,3])
            >>> a.pth_root()
            5*z^2 + 152*z + 119
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

            >>> gf = fq_default_ctx(163, 3)
            >>> a = gf([1,2,3])
            >>> a.trace()
            124
        """
        cdef fmpz tr = fmpz.__new__(fmpz)
        fq_default_trace(tr.val, self.val, self.ctx.val)
        return tr

    def norm(self):
        """
        Returns the norm of self

            >>> gf = fq_default_ctx(163, 3)
            >>> a = gf([1,2,3])
            >>> a.norm()
            116
        """
        cdef fmpz nrm = fmpz.__new__(fmpz)
        fq_default_norm(nrm.val, self.val, self.ctx.val)
        return nrm

    def frobenius(self, e=1):
        r"""
        Evaluates the homomorphism `\Sigma^e` on ``self``.

            >>> gf = fq_default_ctx(163, 3)
            >>> a = gf([1,2,3])
            >>> a.frobenius()
            155*z^2 + 9*z + 4
            >>> a.frobenius(2)
            5*z^2 + 152*z + 119
            >>> a == a.frobenius(3)
            True
            >>> a.frobenius(2) == a.frobenius(-1)
            True
        """
        cdef fq_default res
        res = self.ctx.new_ctype_fq_default()
        fq_default_frobenius(res.val, self.val, <slong>e, self.ctx.val)
        return res
