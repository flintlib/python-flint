from flint.flint_base.flint_base cimport (
    flint_mpoly,
    flint_mpoly_context,
    ordering_py_to_c,
    ordering_c_to_py,
)
from flint.flint_base.flint_base import Ordering

from flint.utils.typecheck cimport typecheck
from flint.utils.flint_exceptions import DomainError, IncompatibleContextError

from flint.types.fmpz cimport any_as_fmpz, fmpz
from flint.types.fmpz_vec cimport fmpz_vec
from flint.types.fmpz_mod cimport fmpz_mod

from flint.flintlib.fmpz cimport fmpz_set
from flint.flintlib.fmpz_mod_mpoly cimport (
    fmpz_mod_mpoly_add,
    fmpz_mod_mpoly_add_fmpz,
    fmpz_mod_mpoly_clear,
    fmpz_mod_mpoly_compose_fmpz_mod_mpoly,
    fmpz_mod_mpoly_ctx_get_modulus,
    fmpz_mod_mpoly_ctx_init,
    fmpz_mod_mpoly_degrees_fmpz,
    fmpz_mod_mpoly_derivative,
    fmpz_mod_mpoly_div,
    fmpz_mod_mpoly_divides,
    fmpz_mod_mpoly_divrem,
    fmpz_mod_mpoly_equal,
    fmpz_mod_mpoly_equal_fmpz,
    fmpz_mod_mpoly_evaluate_all_fmpz,
    fmpz_mod_mpoly_evaluate_one_fmpz,
    fmpz_mod_mpoly_gcd,
    fmpz_mod_mpoly_gen,
    fmpz_mod_mpoly_get_coeff_fmpz_fmpz,
    fmpz_mod_mpoly_get_str_pretty,
    fmpz_mod_mpoly_get_term,
    fmpz_mod_mpoly_get_term_coeff_fmpz,
    fmpz_mod_mpoly_get_term_exp_fmpz,
    fmpz_mod_mpoly_is_one,
    fmpz_mod_mpoly_is_zero,
    fmpz_mod_mpoly_length,
    fmpz_mod_mpoly_mul,
    fmpz_mod_mpoly_neg,
    fmpz_mod_mpoly_pow_fmpz,
    fmpz_mod_mpoly_push_term_fmpz_ffmpz,
    fmpz_mod_mpoly_scalar_mul_fmpz,
    fmpz_mod_mpoly_set,
    fmpz_mod_mpoly_set_coeff_fmpz_fmpz,
    fmpz_mod_mpoly_set_fmpz,
    fmpz_mod_mpoly_set_str_pretty,
    fmpz_mod_mpoly_sort_terms,
    fmpz_mod_mpoly_sub,
    fmpz_mod_mpoly_sub_fmpz,
    fmpz_mod_mpoly_total_degree_fmpz,
    fmpz_mod_mpoly_sqrt,
)
from flint.flintlib.fmpz_mod_mpoly_factor cimport (
    fmpz_mod_mpoly_factor,
    fmpz_mod_mpoly_factor_clear,
    fmpz_mod_mpoly_factor_init,
    fmpz_mod_mpoly_factor_squarefree,
    fmpz_mod_mpoly_factor_t,
)

from cpython.object cimport Py_EQ, Py_NE
cimport libc.stdlib

cdef dict _fmpz_mod_mpoly_ctx_cache = {}


cdef class fmpz_mod_mpoly_ctx(flint_mpoly_context):
    """
    A class for storing the polynomial context

    :param nvars: The number of variables in the ring
    :param ordering:  The term order for the ring
    :param names:  A tuple containing the names of the variables of the ring.

    Do not construct one of these directly, use `fmpz_mod_mpoly_ctx.get_context`.
    """

    _ctx_cache = _fmpz_mod_mpoly_ctx_cache

    def __init__(self, slong nvars, ordering, names, modulus):
        cdef fmpz m
        if not typecheck(modulus, fmpz):
            m = any_as_fmpz(modulus)
            if m is NotImplemented:
                raise TypeError(f"modulus ({modulus}) is not coercible to fmpz")
        else:
            m = modulus
        fmpz_mod_mpoly_ctx_init(self.val, nvars, ordering_py_to_c(ordering), m.val)
        self.__prime_modulus = None
        super().__init__(nvars, names)

    @classmethod
    def create_context_key(
            cls,
            slong nvars=1,
            ordering=Ordering.lex,
            modulus = None,
            names: Optional[str] = "x",
            nametup: Optional[tuple] = None,
    ):
        """
        Create a key for the context cache via the number of variables, the ordering, the modulus, and either a
        variable name string, or a tuple of variable names.
        """
        # A type hint of `ordering: Ordering` results in the error "TypeError: an integer is required" if a Ordering
        # object is not provided. This is pretty obtuse so we check its type ourselves
        if not isinstance(ordering, Ordering):
            raise TypeError(f"`ordering` ('{ordering}') is not an instance of flint.Ordering")
        elif not typecheck(modulus, fmpz):
            m = any_as_fmpz(modulus)
            if m is NotImplemented:
                raise TypeError(f"`modulus` ('{modulus}') is not coercible to fmpz")
            else:
                modulus = m

        if nametup is not None:
            key = nvars, ordering, nametup, modulus
        elif nametup is None and names is not None:
            key = nvars, ordering, cls.create_variable_names(nvars, names), modulus
        else:
            raise ValueError("must provide either `names` or `nametup`")
        return key

    def nvars(self):
        """
        Return the number of variables in the context

            >>> from flint import Ordering
            >>> ctx = fmpz_mod_mpoly_ctx.get_context(4, Ordering.lex, 11, 'x')
            >>> ctx.nvars()
            4
        """
        return self.val.minfo.nvars

    def ordering(self):
        """
        Return the term order of the context object.

            >>> from flint import Ordering
            >>> ctx = fmpz_mod_mpoly_ctx.get_context(4, Ordering.deglex, 11, 'w')
            >>> ctx.ordering()
            <Ordering.deglex: 1>
        """
        return ordering_c_to_py(self.val.minfo.ord)

    def modulus(self):
        """
        Return the modulus of the context object.

            >>> from flint import Ordering
            >>> ctx = fmpz_mod_mpoly_ctx.get_context(4, Ordering.deglex, 2, 'w')
            >>> ctx.modulus()
            2

        """
        cdef fmpz m = fmpz.__new__(fmpz)
        fmpz_mod_mpoly_ctx_get_modulus(m.val, self.val)
        return m

    def is_prime(self):
        """
        Return whether the modulus is prime

            >>> from flint import Ordering
            >>> ctx = fmpz_mod_mpoly_ctx.get_context(4, Ordering.degrevlex, 2**127, 'z')
            >>> ctx.is_prime()
            False
            >>> ctx = fmpz_mod_mpoly_ctx.get_context(4, Ordering.degrevlex, 2**127 - 1, 'z')
            >>> ctx.is_prime()
            True
        """
        if self.__prime_modulus is None:
            self.__prime_modulus = self.modulus().is_prime()
        return self.__prime_modulus

    def gen(self, slong i):
        """
        Return the `i`th generator of the polynomial ring

            >>> from flint import Ordering
            >>> ctx = fmpz_mod_mpoly_ctx.get_context(3, Ordering.degrevlex, 11, 'z')
            >>> ctx.gen(1)
            z1
        """
        cdef fmpz_mod_mpoly res
        if not 0 <= i < self.val.minfo.nvars:
            raise IndexError("generator index out of range")
        res = create_fmpz_mod_mpoly(self)
        fmpz_mod_mpoly_gen(res.val, i, res.ctx.val)
        return res

    def constant(self, z):
        """
        Create a constant polynomial in this context
        """
        cdef fmpz_mod_mpoly res
        z = any_as_fmpz(z)
        if z is NotImplemented:
            raise TypeError("cannot coerce argument to fmpz")
        res = create_fmpz_mod_mpoly(self)
        fmpz_mod_mpoly_set_fmpz(res.val, (<fmpz>z).val, res.ctx.val)
        return res

    def from_dict(self, d):
        """
        Create a fmpz_mod_mpoly from a dictionary in this context.

        The dictionary's keys are tuples of ints (or anything that implicitly converts
        to fmpz) representing exponents, and corresponding values of fmpz.

            >>> from flint import Ordering
            >>> ctx = fmpz_mod_mpoly_ctx.get_context(2, Ordering.lex, 11, 'x,y')
            >>> ctx.from_dict({(1,0):2, (1,1):3, (0,1):1})
            3*x*y + 2*x + y
        """
        cdef:
            fmpz_vec exp_vec
            slong i, nvars = self.nvars()
            fmpz_mod_mpoly res

        if not isinstance(d, dict):
            raise ValueError("expected a dictionary")

        res = create_fmpz_mod_mpoly(self)

        for i, (exps, coeff) in enumerate(d.items()):
            if len(exps) != nvars:
                raise ValueError(f"expected {nvars} exponents, got {len(exps)}")
            elif not coeff:
                continue

            exp_vec = fmpz_vec(exps)

            if isinstance(coeff, int) or typecheck(coeff, fmpz):
                coeff_fmpz = any_as_fmpz(coeff)
                if coeff_fmpz is NotImplemented:
                    raise TypeError(f"cannot coerce '{repr(coeff)}' to fmpz")

                fmpz_mod_mpoly_push_term_fmpz_ffmpz(
                    res.val,
                    (<fmpz>coeff_fmpz).val,
                    exp_vec.val,
                    self.val
                )
            elif typecheck(coeff, fmpz_mod):
                fmpz_mod_mpoly_push_term_fmpz_ffmpz(
                    res.val,
                    (<fmpz_mod>coeff).val,
                    exp_vec.val,
                    self.val
                )
            else:
                raise TypeError(f"cannot coerce {repr(coeff)} to nmod_mpoly coefficient")

        fmpz_mod_mpoly_sort_terms(res.val, self.val)
        return res


cdef class fmpz_mod_mpoly(flint_mpoly):
    """
    The *fmpz_mod_mpoly* type represents sparse multivariate polynomials over
    the integers modulo `n`, for large `n`.
    """

    def __cinit__(self):
        self._init = False

    def __dealloc__(self):
        if self._init:
            fmpz_mod_mpoly_clear(self.val, self.ctx.val)
            self._init = False

    def __init__(self, val=0, ctx=None):
        if typecheck(val, fmpz_mod_mpoly):
            if ctx is None or ctx == (<fmpz_mod_mpoly>val).ctx:
                init_fmpz_mod_mpoly(self, (<fmpz_mod_mpoly>val).ctx)
                fmpz_mod_mpoly_set(self.val, (<fmpz_mod_mpoly>val).val, self.ctx.val)
            else:
                raise IncompatibleContextError(f"{ctx} is not {(<fmpz_mod_mpoly>val).ctx}")
        elif isinstance(val, dict):
            if ctx is None:
                raise ValueError("a context is required to create a fmpz_mod_mpoly from a dict")
            x = ctx.from_dict(val)
            # XXX: this copy is silly, have a ctx function that assigns an fmpz_mod_mpoly_t
            init_fmpz_mod_mpoly(self, ctx)
            fmpz_mod_mpoly_set(self.val, (<fmpz_mod_mpoly>x).val, self.ctx.val)
        elif isinstance(val, str):
            if ctx is None:
                raise ValueError("cannot parse a polynomial without context")
            val = bytes(val, 'utf-8')
            init_fmpz_mod_mpoly(self, ctx)
            if fmpz_mod_mpoly_set_str_pretty(self.val, val, self.ctx.c_names, self.ctx.val) == -1:
                raise ValueError("unable to parse fmpz_mod_mpoly from string")
            fmpz_mod_mpoly_sort_terms(self.val, self.ctx.val)
        elif typecheck(val, fmpz_mod):
            if ctx is None:
                raise ValueError("need context to convert fmpz_mod to fmpz_mod_mpoly")
            init_fmpz_mod_mpoly(self, ctx)
            fmpz_mod_mpoly_set_fmpz(self.val, (<fmpz_mod>val).val, self.ctx.val)
        else:
            v = any_as_fmpz(val)
            if v is NotImplemented:
                raise TypeError("cannot create fmpz_mod_mpoly from type %s" % type(val))
            elif ctx is None:
                raise ValueError("need context to convert fmpz to fmpz_mod_mpoly")
            init_fmpz_mod_mpoly(self, ctx)
            fmpz_mod_mpoly_set_fmpz(self.val, (<fmpz>v).val, self.ctx.val)

    def __bool__(self):
        return not fmpz_mod_mpoly_is_zero(self.val, self.ctx.val)

    def is_zero(self):
        return <bint>fmpz_mod_mpoly_is_zero(self.val, self.ctx.val)

    def is_one(self):
        return <bint>fmpz_mod_mpoly_is_one(self.val, self.ctx.val)

    def __richcmp__(self, other, int op):
        if not (op == Py_EQ or op == Py_NE):
            return NotImplemented
        elif other is None:
            return op == Py_NE
        elif typecheck(self, fmpz_mod_mpoly) and typecheck(other, fmpz_mod_mpoly):
            if (<fmpz_mod_mpoly>self).ctx is (<fmpz_mod_mpoly>other).ctx:
                return (op == Py_NE) ^ bool(
                    fmpz_mod_mpoly_equal((<fmpz_mod_mpoly>self).val, (<fmpz_mod_mpoly>other).val, (<fmpz_mod_mpoly>self).ctx.val)
                )
            else:
                return op == Py_NE
        elif typecheck(other, fmpz):
            return (op == Py_NE) ^ <bint>fmpz_mod_mpoly_equal_fmpz(self.val, (<fmpz>other).val, self.ctx.val)
        elif typecheck(other, fmpz_mod):
            return (op == Py_NE) ^ <bint>fmpz_mod_mpoly_equal_fmpz(self.val, (<fmpz_mod>other).val, self.ctx.val)
        elif isinstance(other, int):
            other = any_as_fmpz(other)
            if other is NotImplemented:
                return NotImplemented
            return (op == Py_NE) ^ <bint>fmpz_mod_mpoly_equal_fmpz(self.val, (<fmpz>other).val, self.ctx.val)
        else:
            return NotImplemented

    def __len__(self):
        return fmpz_mod_mpoly_length(self.val, self.ctx.val)

    def __getitem__(self, x):
        """
        Return the coefficient of the term with the exponent vector `x`.
        Always returns a value, missing keys will return `0`.
        Negative exponents are made positive.

            >>> from flint import Ordering
            >>> ctx = fmpz_mod_mpoly_ctx.get_context(2, Ordering.lex, 11, 'x')
            >>> p = ctx.from_dict({(0, 1): 2, (1, 1): 3})
            >>> p[1, 1]
            3

        """
        cdef:
            slong nvars = self.ctx.nvars()

        if not isinstance(x, tuple):
            raise TypeError("exponent vector index is not a tuple")
        elif len(x) != nvars:
            raise ValueError("exponent vector provided does not match number of variables")

        res = fmpz()
        exp_vec = fmpz_vec(x, double_indirect=True)
        fmpz_mod_mpoly_get_coeff_fmpz_fmpz((<fmpz>res).val, self.val, exp_vec.double_indirect, self.ctx.val)
        return res

    def __setitem__(self, x, y):
        """
        Set the coefficient of the term with the exponent vector `x` to `y`.
        Will always set a value, missing keys will create a new term.
        Negative exponents are made positive.

            >>> from flint import Ordering
            >>> ctx = fmpz_mod_mpoly_ctx.get_context(2, Ordering.lex, 11, 'x')
            >>> p = ctx.from_dict({(0, 1): 2, (1, 1): 3})
            >>> p[1, 1] = 20
            >>> p
            20*x0*x1 + 2*x1

        """
        cdef:
            slong nvars = self.ctx.nvars()

        if not isinstance(x, tuple):
            raise TypeError("exponent vector index is not a tuple")
        elif len(x) != nvars:
            raise ValueError("exponent vector provided does not match number of variables")
        exp_vec = fmpz_vec(x, double_indirect=True)

        if typecheck(y, fmpz_mod):
            if self.ctx.modulus() != (<fmpz_mod>y).ctx.modulus():
                raise ValueError(
                    f"cannot set a coefficient with different modulus {(<fmpz_mod>y).ctx.modulus()} vs {self.ctx.modulus()}"
                )
            fmpz_mod_mpoly_set_coeff_fmpz_fmpz(self.val, (<fmpz_mod>y).val, exp_vec.double_indirect, self.ctx.val)
        else:
            coeff = any_as_fmpz(y)
            if coeff is NotImplemented:
                raise TypeError("provided coefficient not coercible to fmpz")
            fmpz_mod_mpoly_set_coeff_fmpz_fmpz(self.val, (<fmpz>coeff).val, exp_vec.double_indirect, self.ctx.val)

    def __neg__(self):
        cdef fmpz_mod_mpoly res
        res = create_fmpz_mod_mpoly(self.ctx)
        fmpz_mod_mpoly_neg(res.val, (<fmpz_mod_mpoly>self).val, res.ctx.val)
        return res

    def __add__(self, other):
        cdef fmpz_mod_mpoly res
        if typecheck(other, fmpz_mod_mpoly):
            if (<fmpz_mod_mpoly>self).ctx is not (<fmpz_mod_mpoly>other).ctx:
                raise IncompatibleContextError(f"{(<fmpz_mod_mpoly>self).ctx} is not {(<fmpz_mod_mpoly>other).ctx}")
            res = create_fmpz_mod_mpoly(self.ctx)
            fmpz_mod_mpoly_add(res.val, (<fmpz_mod_mpoly>self).val, (<fmpz_mod_mpoly>other).val, res.ctx.val)
            return res
        elif typecheck(other, fmpz_mod):
            if self.ctx.modulus() != (<fmpz_mod>other).ctx.modulus():
                raise ValueError(
                    f"cannot add with different modulus {(<fmpz_mod>other).ctx.modulus()} vs {self.ctx.modulus()}"
                )
            res = create_fmpz_mod_mpoly(self.ctx)
            fmpz_mod_mpoly_add_fmpz(res.val, (<fmpz_mod_mpoly>self).val, (<fmpz_mod>other).val, self.ctx.val)
            return res
        else:
            other = any_as_fmpz(other)
            if other is not NotImplemented:
                res = create_fmpz_mod_mpoly(self.ctx)
                fmpz_mod_mpoly_add_fmpz(res.val, (<fmpz_mod_mpoly>self).val, (<fmpz>other).val, self.ctx.val)
                return res
        return NotImplemented

    def __radd__(self, other):
        cdef fmpz_mod_mpoly res
        if typecheck(other, fmpz_mod):
            if self.ctx.modulus() != (<fmpz_mod>other).ctx.modulus():
                raise ValueError(
                    f"cannot add with different modulus {(<fmpz_mod>other).ctx.modulus()} vs {self.ctx.modulus()}"
                )
            res = create_fmpz_mod_mpoly(self.ctx)
            fmpz_mod_mpoly_add_fmpz(res.val, (<fmpz_mod_mpoly>self).val, (<fmpz_mod>other).val, self.ctx.val)
            return res
        else:
            other = any_as_fmpz(other)
            if other is not NotImplemented:
                res = create_fmpz_mod_mpoly(self.ctx)
                fmpz_mod_mpoly_add_fmpz(res.val, (<fmpz_mod_mpoly>self).val, (<fmpz>other).val, self.ctx.val)
                return res
        return NotImplemented

    def __sub__(self, other):
        cdef fmpz_mod_mpoly res
        if typecheck(other, fmpz_mod_mpoly):
            if (<fmpz_mod_mpoly>self).ctx is not (<fmpz_mod_mpoly>other).ctx:
                raise IncompatibleContextError(f"{(<fmpz_mod_mpoly>self).ctx} is not {(<fmpz_mod_mpoly>other).ctx}")
            res = create_fmpz_mod_mpoly(self.ctx)
            fmpz_mod_mpoly_sub(res.val, (<fmpz_mod_mpoly>self).val, (<fmpz_mod_mpoly>other).val, res.ctx.val)
            return res
        elif typecheck(other, fmpz_mod):
            if self.ctx.modulus() != (<fmpz_mod>other).ctx.modulus():
                raise ValueError(
                    f"cannot subtract with different modulus {(<fmpz_mod>other).ctx.modulus()} vs {self.ctx.modulus()}"
                )
            res = create_fmpz_mod_mpoly(self.ctx)
            fmpz_mod_mpoly_sub_fmpz(res.val, (<fmpz_mod_mpoly>self).val, (<fmpz_mod>other).val, self.ctx.val)
            return res
        else:
            other = any_as_fmpz(other)
            if other is not NotImplemented:
                res = create_fmpz_mod_mpoly(self.ctx)
                fmpz_mod_mpoly_sub_fmpz(res.val, (<fmpz_mod_mpoly>self).val, (<fmpz>other).val, self.ctx.val)
                return res
        return NotImplemented

    def __rsub__(self, other):
        cdef fmpz_mod_mpoly res
        if typecheck(other, fmpz_mod):
            if self.ctx.modulus() != (<fmpz_mod>other).ctx.modulus():
                raise ValueError(
                    f"cannot subtract with different modulus {(<fmpz_mod>other).ctx.modulus()} vs {self.ctx.modulus()}"
                )
            res = create_fmpz_mod_mpoly(self.ctx)
            fmpz_mod_mpoly_sub_fmpz(res.val, (<fmpz_mod_mpoly>self).val, (<fmpz_mod>other).val, self.ctx.val)
            return -res
        else:
            other = any_as_fmpz(other)
            if other is not NotImplemented:
                res = create_fmpz_mod_mpoly(self.ctx)
                fmpz_mod_mpoly_sub_fmpz(res.val, (<fmpz_mod_mpoly>self).val, (<fmpz>other).val, res.ctx.val)
                return -res
        return NotImplemented

    def __mul__(self, other):
        cdef fmpz_mod_mpoly res
        if typecheck(other, fmpz_mod_mpoly):
            if (<fmpz_mod_mpoly>self).ctx is not (<fmpz_mod_mpoly>other).ctx:
                raise IncompatibleContextError(f"{(<fmpz_mod_mpoly>self).ctx} is not {(<fmpz_mod_mpoly>other).ctx}")
            res = create_fmpz_mod_mpoly(self.ctx)
            fmpz_mod_mpoly_mul(res.val, (<fmpz_mod_mpoly>self).val, (<fmpz_mod_mpoly>other).val, res.ctx.val)
            return res
        elif typecheck(other, fmpz_mod):
            if self.ctx.modulus() != (<fmpz_mod>other).ctx.modulus():
                raise ValueError(
                    f"cannot multiply with different modulus {(<fmpz_mod>other).ctx.modulus()} vs {self.ctx.modulus()}"
                )
            res = create_fmpz_mod_mpoly(self.ctx)
            fmpz_mod_mpoly_scalar_mul_fmpz(res.val, (<fmpz_mod_mpoly>self).val, (<fmpz_mod>other).val, self.ctx.val)
            return res
        else:
            other = any_as_fmpz(other)
            if other is not NotImplemented:
                res = create_fmpz_mod_mpoly(self.ctx)
                fmpz_mod_mpoly_scalar_mul_fmpz(res.val, (<fmpz_mod_mpoly>self).val, (<fmpz>other).val, res.ctx.val)
                return res
        return NotImplemented

    def __rmul__(self, other):
        cdef fmpz_mod_mpoly res
        if typecheck(other, fmpz_mod):
            if self.ctx.modulus() != (<fmpz_mod>other).ctx.modulus():
                raise ValueError(
                    f"cannot multiply with different modulus {(<fmpz_mod>other).ctx.modulus()} vs {self.ctx.modulus()}"
                )
            res = create_fmpz_mod_mpoly(self.ctx)
            fmpz_mod_mpoly_scalar_mul_fmpz(res.val, (<fmpz_mod_mpoly>self).val, (<fmpz_mod>other).val, self.ctx.val)
            return res
        else:
            other = any_as_fmpz(other)
            if other is not NotImplemented:
                res = create_fmpz_mod_mpoly(self.ctx)
                fmpz_mod_mpoly_scalar_mul_fmpz(res.val, (<fmpz_mod_mpoly>self).val, (<fmpz>other).val, res.ctx.val)
                return res
        return NotImplemented

    def __pow__(self, other, modulus):
        cdef fmpz_mod_mpoly res
        if modulus is not None:
            raise NotImplementedError
        other = any_as_fmpz(other)
        if other is NotImplemented:
            return other
        if other < 0:
            raise ValueError("cannot raise fmpz_mod_mpoly to negative power")
        res = create_fmpz_mod_mpoly(self.ctx)
        if fmpz_mod_mpoly_pow_fmpz(res.val, (<fmpz_mod_mpoly>self).val, (<fmpz>other).val, res.ctx.val) == 0:
            raise ValueError("unreasonably large polynomial")  # pragma: no cover
        return res

    def __divmod__(self, other):
        cdef fmpz_mod_mpoly res, res2

        if typecheck(other, fmpz_mod_mpoly):
            if not other:
                raise ZeroDivisionError("fmpz_mod_mpoly division by zero")
            elif not self.ctx.is_prime():
                raise DomainError("division with non-prime modulus is not supported")
            elif (<fmpz_mod_mpoly>self).ctx is not (<fmpz_mod_mpoly>other).ctx:
                raise IncompatibleContextError(f"{(<fmpz_mod_mpoly>self).ctx} is not {(<fmpz_mod_mpoly>other).ctx}")
            res = create_fmpz_mod_mpoly(self.ctx)
            res2 = create_fmpz_mod_mpoly(self.ctx)
            fmpz_mod_mpoly_divrem(res.val, res2.val, (<fmpz_mod_mpoly>self).val, (<fmpz_mod_mpoly>other).val, res.ctx.val)
            return (res, res2)
        elif typecheck(other, fmpz) or typecheck(other, fmpz_mod) or isinstance(other, int):
            if not other:
                raise ZeroDivisionError("fmpz_mod_mpoly division by zero")
            elif not self.ctx.is_prime():
                raise DomainError("division with non-prime modulus is not supported")
            elif typecheck(other, fmpz_mod) and self.ctx.modulus() != (<fmpz_mod>other).ctx.modulus():
                raise ValueError(
                    f"cannot divide with different modulus {(<fmpz_mod>other).ctx.modulus()} vs {self.ctx.modulus()}"
                )

            other = fmpz_mod_mpoly(other, self.ctx)
            res = create_fmpz_mod_mpoly(self.ctx)
            res2 = create_fmpz_mod_mpoly(self.ctx)
            fmpz_mod_mpoly_divrem(res.val, res2.val, (<fmpz_mod_mpoly>self).val, (<fmpz_mod_mpoly>other).val, res.ctx.val)
            return (res, res2)
        return NotImplemented

    def __rdivmod__(self, other):
        cdef fmpz_mod_mpoly res, res2
        if typecheck(other, fmpz) or typecheck(other, fmpz_mod) or isinstance(other, int):
            if not self:
                raise ZeroDivisionError("fmpz_mod_mpoly division by zero")
            elif not self.ctx.is_prime():
                raise DomainError("division with non-prime modulus is not supported")
            elif typecheck(other, fmpz_mod) and self.ctx.modulus() != (<fmpz_mod>other).ctx.modulus():
                raise ValueError(
                    f"cannot divide with different modulus {(<fmpz_mod>other).ctx.modulus()} vs {self.ctx.modulus()}"
                )

            other = fmpz_mod_mpoly(other, self.ctx)
            res = create_fmpz_mod_mpoly(self.ctx)
            res2 = create_fmpz_mod_mpoly(self.ctx)
            fmpz_mod_mpoly_divrem(res.val, res2.val, (<fmpz_mod_mpoly>other).val, (<fmpz_mod_mpoly>self).val, res.ctx.val)
            return (res, res2)
        return NotImplemented

    def __floordiv__(self, other):
        cdef fmpz_mod_mpoly res
        if typecheck(other, fmpz_mod_mpoly):
            if not other:
                raise ZeroDivisionError("fmpz_mod_mpoly division by zero")
            elif not self.ctx.is_prime():
                raise DomainError("division with non-prime modulus is not supported")
            elif (<fmpz_mod_mpoly>self).ctx is not (<fmpz_mod_mpoly>other).ctx:
                raise IncompatibleContextError(f"{(<fmpz_mod_mpoly>self).ctx} is not {(<fmpz_mod_mpoly>other).ctx}")
            res = create_fmpz_mod_mpoly(self.ctx)
            fmpz_mod_mpoly_div(res.val, (<fmpz_mod_mpoly>self).val, (<fmpz_mod_mpoly>other).val, res.ctx.val)
            return res
        elif typecheck(other, fmpz) or typecheck(other, fmpz_mod) or isinstance(other, int):
            if not other:
                raise ZeroDivisionError("fmpz_mod_mpoly division by zero")
            elif not self.ctx.is_prime():
                raise DomainError("division with non-prime modulus is not supported")
            elif typecheck(other, fmpz_mod) and self.ctx.modulus() != (<fmpz_mod>other).ctx.modulus():
                raise ValueError(
                    f"cannot divide with different modulus {(<fmpz_mod>other).ctx.modulus()} vs {self.ctx.modulus()}"
                )

            other = fmpz_mod_mpoly(other, self.ctx)
            res = create_fmpz_mod_mpoly(self.ctx)
            fmpz_mod_mpoly_div(res.val, (<fmpz_mod_mpoly>self).val, (<fmpz_mod_mpoly>other).val, res.ctx.val)
            return res
        return NotImplemented

    def __rfloordiv__(self, other):
        cdef fmpz_mod_mpoly res
        if typecheck(other, fmpz) or typecheck(other, fmpz_mod) or isinstance(other, int):
            if not self:
                raise ZeroDivisionError("fmpz_mod_mpoly division by zero")
            elif not self.ctx.is_prime():
                raise DomainError("division with non-prime modulus is not supported")
            elif typecheck(other, fmpz_mod) and self.ctx.modulus() != (<fmpz_mod>other).ctx.modulus():
                raise ValueError(
                    f"cannot divide with different modulus {(<fmpz_mod>other).ctx.modulus()} vs {self.ctx.modulus()}"
                )

            other = fmpz_mod_mpoly(other, self.ctx)
            res = create_fmpz_mod_mpoly(self.ctx)
            fmpz_mod_mpoly_div(res.val, (<fmpz_mod_mpoly>other).val,  self.val, res.ctx.val)
            return res
        return NotImplemented

    def __truediv__(self, other):
        cdef:
            fmpz_mod_mpoly res, div

        if typecheck(other, fmpz_mod_mpoly):
            if not other:
                raise ZeroDivisionError("fmpz_mod_mpoly division by zero")
            elif not self.ctx.is_prime():
                raise DomainError("division with non-prime modulus is not supported")
            elif self.ctx is not (<fmpz_mod_mpoly>other).ctx:
                raise IncompatibleContextError(f"{self.ctx} is not {(<fmpz_mod_mpoly>other).ctx}")

            res = create_fmpz_mod_mpoly(self.ctx)
            if fmpz_mod_mpoly_divides(res.val, self.val, (<fmpz_mod_mpoly>other).val, self.ctx.val):
                return res
            else:
                raise DomainError("fmpz_mod_mpoly division is not exact")
        elif typecheck(other, fmpz) or typecheck(other, fmpz_mod) or isinstance(other, int):
            if not other:
                raise ZeroDivisionError("fmpz_mod_mpoly division by zero")
            elif not self.ctx.is_prime():
                raise DomainError("division with non-prime modulus is not supported")
            elif typecheck(other, fmpz_mod) and self.ctx.modulus() != (<fmpz_mod>other).ctx.modulus():
                raise ValueError(
                    f"cannot divide with different modulus {(<fmpz_mod>other).ctx.modulus()} vs {self.ctx.modulus()}"
                )

            div = fmpz_mod_mpoly(other, self.ctx)
            res = create_fmpz_mod_mpoly(self.ctx)
            if fmpz_mod_mpoly_divides(res.val, self.val, div.val, self.ctx.val):
                return res
            else:
                raise DomainError("fmpz_mod_mpoly division is not exact")
        return NotImplemented

    def __rtruediv__(self, other):
        cdef fmpz_mod_mpoly res
        if typecheck(other, fmpz) or typecheck(other, fmpz_mod) or isinstance(other, int):
            if not self:
                raise ZeroDivisionError("fmpz_mod_mpoly division by zero")
            elif not self.ctx.is_prime():
                raise DomainError("division with non-prime modulus is not supported")
            elif typecheck(other, fmpz_mod) and self.ctx.modulus() != (<fmpz_mod>other).ctx.modulus():
                raise ValueError(
                    f"cannot divide with different modulus {(<fmpz_mod>other).ctx.modulus()} vs {self.ctx.modulus()}"
                )
            other = fmpz_mod_mpoly(other, self.ctx)
            res = create_fmpz_mod_mpoly(self.ctx)
            if fmpz_mod_mpoly_divides(res.val, (<fmpz_mod_mpoly>other).val, self.val, self.ctx.val):
                return res
            else:
                raise DomainError("fmpz_mod_mpoly division is not exact")
        return NotImplemented

    def __mod__(self, other):
        return divmod(self, other)[1]

    def __rmod__(self, other):
        return divmod(other, self)[1]

    def __call__(self, *args) -> fmpz:
        cdef:
            fmpz_vec V
            fmpz vres
            slong nvars = self.ctx.nvars(), nargs = len(args)

        if nargs != nvars:
            raise ValueError("number of generators does not match number of arguments")  # TODO: fix for all of them

        V = fmpz_vec(args, double_indirect=True)
        vres = fmpz.__new__(fmpz)
        fmpz_mod_mpoly_evaluate_all_fmpz(vres.val, self.val, V.double_indirect, self.ctx.val)
        return vres

    def iadd(self, other):
        """
        In-place addition, mutates self.

            >>> from flint import Ordering
            >>> ctx = fmpz_mod_mpoly_ctx.get_context(2, Ordering.lex, 11, 'x')
            >>> f = ctx.from_dict({(1, 0): 2, (0, 1): 3, (1, 1): 4})
            >>> f
            4*x0*x1 + 2*x0 + 3*x1
            >>> f.iadd(5)
            >>> f
            4*x0*x1 + 2*x0 + 3*x1 + 5

        """
        if typecheck(other, fmpz_mod_mpoly):
            if (<fmpz_mod_mpoly>self).ctx is not (<fmpz_mod_mpoly>other).ctx:
                raise IncompatibleContextError(f"{(<fmpz_mod_mpoly>self).ctx} is not {(<fmpz_mod_mpoly>other).ctx}")
            fmpz_mod_mpoly_add((<fmpz_mod_mpoly>self).val, (<fmpz_mod_mpoly>self).val, (<fmpz_mod_mpoly>other).val, self.ctx.val)
            return
        elif typecheck(other, fmpz_mod):
            if self.ctx.modulus() != (<fmpz_mod>other).ctx.modulus():
                raise ValueError(
                    f"cannot add with different modulus {(<fmpz_mod>other).ctx.modulus()} vs {self.ctx.modulus()}"
                )
            fmpz_mod_mpoly_add_fmpz(self.val, self.val, (<fmpz_mod>other).val, self.ctx.val)
            return
        else:
            zval = any_as_fmpz(other)
            if zval is not NotImplemented:
                fmpz_mod_mpoly_add_fmpz((<fmpz_mod_mpoly>self).val, (<fmpz_mod_mpoly>self).val, (<fmpz>zval).val, self.ctx.val)
                return
        raise NotImplementedError()

    def isub(self, other):
        """
        In-place subtraction, mutates self.

            >>> from flint import Ordering
            >>> ctx = fmpz_mod_mpoly_ctx.get_context(2, Ordering.lex, 11, 'x')
            >>> f = ctx.from_dict({(1, 0): 2, (0, 1): 3, (1, 1): 4})
            >>> f
            4*x0*x1 + 2*x0 + 3*x1
            >>> f.isub(5)
            >>> f
            4*x0*x1 + 2*x0 + 3*x1 - 5

        """
        if typecheck(other, fmpz_mod_mpoly):
            if (<fmpz_mod_mpoly>self).ctx is not (<fmpz_mod_mpoly>other).ctx:
                raise IncompatibleContextError(f"{(<fmpz_mod_mpoly>self).ctx} is not {(<fmpz_mod_mpoly>other).ctx}")
            fmpz_mod_mpoly_sub((<fmpz_mod_mpoly>self).val, (<fmpz_mod_mpoly>self).val, (<fmpz_mod_mpoly>other).val, self.ctx.val)
            return
        elif typecheck(other, fmpz_mod):
            if self.ctx.modulus() != (<fmpz_mod>other).ctx.modulus():
                raise ValueError(
                    f"cannot subtract with different modulus {(<fmpz_mod>other).ctx.modulus()} vs {self.ctx.modulus()}"
                )
            fmpz_mod_mpoly_sub_fmpz(self.val, self.val, (<fmpz_mod>other).val, self.ctx.val)
            return
        else:
            other = any_as_fmpz(other)
            if other is not NotImplemented:
                fmpz_mod_mpoly_sub_fmpz((<fmpz_mod_mpoly>self).val, (<fmpz_mod_mpoly>self).val, (<fmpz>other).val, self.ctx.val)
                return
        raise NotImplementedError()

    def imul(self, other):
        """
        In-place multiplication, mutates self.

            >>> from flint import Ordering
            >>> ctx = fmpz_mod_mpoly_ctx.get_context(2, Ordering.lex, 11, 'x')
            >>> f = ctx.from_dict({(1, 0): 2, (0, 1): 3, (1, 1): 4})
            >>> f
            4*x0*x1 + 2*x0 + 3*x1
            >>> f.imul(2)
            >>> f
            8*x0*x1 + 4*x0 + 6*x1

        """
        if typecheck(other, fmpz_mod_mpoly):
            if (<fmpz_mod_mpoly>self).ctx is not (<fmpz_mod_mpoly>other).ctx:
                raise IncompatibleContextError(f"{(<fmpz_mod_mpoly>self).ctx} is not {(<fmpz_mod_mpoly>other).ctx}")
            fmpz_mod_mpoly_mul((<fmpz_mod_mpoly>self).val, (<fmpz_mod_mpoly>self).val, (<fmpz_mod_mpoly>other).val, self.ctx.val)
            return
        elif typecheck(other, fmpz_mod):
            if self.ctx.modulus() != (<fmpz_mod>other).ctx.modulus():
                raise ValueError(
                    f"cannot multiply with different modulus {(<fmpz_mod>other).ctx.modulus()} vs {self.ctx.modulus()}"
                )
            fmpz_mod_mpoly_scalar_mul_fmpz(self.val, self.val, (<fmpz_mod>other).val, self.ctx.val)
            return
        else:
            other = any_as_fmpz(other)
            if other is not NotImplemented:
                fmpz_mod_mpoly_scalar_mul_fmpz(self.val, (<fmpz_mod_mpoly>self).val, (<fmpz>other).val, self.ctx.val)
                return
        raise NotImplementedError()

    def monoms(self):
        """
        Return the exponent vectors of each term as a tuple of fmpz.

            >>> from flint import Ordering
            >>> ctx = fmpz_mod_mpoly_ctx.get_context(2, Ordering.lex, 11, 'x')
            >>> f = ctx.from_dict({(0, 0): 1, (1, 0): 2, (0, 1): 3, (1, 1): 4})
            >>> f.monoms()
            [(1, 1), (1, 0), (0, 1), (0, 0)]

        """
        cdef:
            slong i, nvars = self.ctx.nvars()
            fmpz_vec vec = fmpz_vec(nvars, double_indirect=True)

        res = []
        for i in range(len(self)):
            fmpz_mod_mpoly_get_term_exp_fmpz(vec.double_indirect, self.val, i, self.ctx.val)
            res.append(vec.to_tuple())

        return res

    def coeffs(self):
        """
        Return the coefficients of each term as a fmpz

            >>> from flint import Ordering
            >>> ctx = fmpz_mod_mpoly_ctx.get_context(2, Ordering.lex, 11, 'x')
            >>> f = ctx.from_dict({(0, 0): 1, (1, 0): 2, (0, 1): 3, (1, 1): 4})
            >>> f.coeffs()
            [4, 2, 3, 1]

        """
        cdef:
            fmpz coeff
            slong i

        res = []
        for i in range(len(self)):
            coeff = fmpz.__new__(fmpz)
            fmpz_mod_mpoly_get_term_coeff_fmpz(coeff.val, self.val, i, self.ctx.val)
            res.append(coeff)

        return res

    # def terms(self):
    #     """
    #     Return the terms of this polynomial as a list of fmpz_mod_mpolys.

    #         >>> from flint import Ordering
    #         >>> ctx = fmpz_mod_mpoly_ctx.get_context(2, Ordering.lex, 11, 'x')
    #         >>> f = ctx.from_dict({(0, 0): 1, (1, 0): 2, (0, 1): 3, (1, 1): 4})
    #         >>> f.terms()
    #         [4*x0*x1, 2*x0, 3*x1, 1]

    #     """
    #     cdef:
    #         fmpz_mod_mpoly term
    #         slong i

    #     res = []
    #     for i in range(len(self)):
    #         term = create_fmpz_mod_mpoly(self.ctx)
    #         fmpz_mod_mpoly_get_term(term.val, self.val, i, self.ctx.val)
    #         res.append(term)

    #     return res

    def subs(self, dict_args) -> fmpz_mod_mpoly:
        """
        Partial evaluate this polynomial with select constants. Keys must be generator names or generator indices,
        all values must be fmpz.

            >>> from flint import Ordering
            >>> ctx = fmpz_mod_mpoly_ctx.get_context(2, Ordering.lex, 11, 'x')
            >>> f = ctx.from_dict({(0, 0): 1, (1, 0): 2, (0, 1): 3, (1, 1): 4})
            >>> f.subs({"x1": 0})
            2*x0 + 1

        """
        cdef:
            fmpz_mod_mpoly res
            slong i

        args = tuple((self.ctx.variable_to_index(k), v) for k, v in dict_args.items())

        # Partial application with args in Z. We evaluate the polynomial one variable at a time
        res = create_fmpz_mod_mpoly(self.ctx)

        fmpz_mod_mpoly_set(res.val, self.val, self.ctx.val)
        for i, arg in args:
            if typecheck(arg, fmpz_mod):
                fmpz_mod_mpoly_evaluate_one_fmpz(res.val, res.val, i, (<fmpz_mod>arg).val, self.ctx.val)
            else:
                arg_fmpz = any_as_fmpz(arg)
                if arg_fmpz is NotImplemented:
                    raise TypeError(f"cannot coerce {type(arg)} to fmpz")
                fmpz_mod_mpoly_evaluate_one_fmpz(res.val, res.val, i, (<fmpz>arg).val, self.ctx.val)

        return res

    def compose(self, *args, ctx=None) -> fmpz_mod_mpoly:
        """
        Compose this polynomial with other fmpz_mod_mpolys. All arguments must share the same context, it may different
        from this polynomials context.

            >>> from flint import Ordering
            >>> ctx = fmpz_mod_mpoly_ctx.get_context(1, Ordering.lex, 11, 'x')
            >>> ctx1 = fmpz_mod_mpoly_ctx.get_context(2, Ordering.lex, 11, 'y')
            >>> f = ctx.from_dict({(2,): 1})
            >>> g = ctx1.from_dict({(1, 0): 1, (0, 1): 1})
            >>> f
            x^2
            >>> g
            y0 + y1
            >>> f.compose(g)
            y0^2 + 2*y0*y1 + y1^2

        """
        cdef:
            fmpz_mod_mpoly res
            fmpz_mod_mpoly_ctx res_ctx
            fmpz_mod_mpoly_vec C
            slong i, nvars = self.ctx.nvars(), nargs = len(args)

        if nargs < nvars:
            raise ValueError("not enough arguments provided")
        elif nargs > nvars:
            raise ValueError("more arguments provided than variables")
        elif self.ctx.nvars() == 0 and ctx is None:
            raise ValueError("a context must be provided when composing a polynomial with no generators")
        elif not all(typecheck(arg, fmpz_mod_mpoly) for arg in args):
            raise TypeError("all arguments must be fmpz_mod_mpolys")

        if ctx is None:
            res_ctx = (<fmpz_mod_mpoly> args[0]).ctx
        elif typecheck(ctx, fmpz_mod_mpoly_ctx):
            res_ctx = <fmpz_mod_mpoly_ctx>ctx
        else:
            raise TypeError(f"provided context ({ctx}) is not a fmpz_mod_mpoly_ctx")

        if not all((<fmpz_mod_mpoly> arg).ctx is res_ctx for arg in args):
            raise IncompatibleContextError(
                "all arguments must share the " + ("same" if ctx is not None else "provided") + " context"
            )

        C = fmpz_mod_mpoly_vec(args, res_ctx, double_indirect=True)
        res = create_fmpz_mod_mpoly(res_ctx)
        if fmpz_mod_mpoly_compose_fmpz_mod_mpoly(res.val, self.val, C.double_indirect, self.ctx.val, res_ctx.val) == 0:
            raise ValueError("unreasonably large polynomial")  # pragma: no cover
        return res

    def context(self):
        """
        Return the context object for this polynomials.

            >>> from flint import Ordering
            >>> ctx = fmpz_mod_mpoly_ctx.get_context(2, Ordering.lex, 11, 'x')
            >>> p = ctx.from_dict({(0, 1): 2})
            >>> ctx is p.context()
            True
        """
        return self.ctx

    def coefficient(self, slong i):
        """
        Return the coefficient at index `i`.

            >>> from flint import Ordering
            >>> ctx = fmpz_mod_mpoly_ctx.get_context(2, Ordering.lex, 11, 'x')
            >>> p = ctx.from_dict({(0, 1): 2, (1, 1): 3})
            >>> p.coefficient(1)
            2
        """
        cdef fmpz v
        if not 0 <= i < fmpz_mod_mpoly_length(self.val, self.ctx.val):
            raise IndexError("term index out of range")
        else:
            v = fmpz.__new__(fmpz)
            fmpz_mod_mpoly_get_term_coeff_fmpz(v.val, self.val, i, self.ctx.val)
            return v

    def monomial(self, slong i):
        """
        Return the exponent vector at index `i` as a tuple.

            >>> from flint import Ordering
            >>> ctx = fmpz_mod_mpoly_ctx.get_context(2, Ordering.lex, 11, 'x')
            >>> p = ctx.from_dict({(0, 1): 2, (1, 1): 3})
            >>> p.monomial(1)
            (0, 1)
        """
        cdef:
            slong nvars = self.ctx.nvars()

        if not 0 <= i < fmpz_mod_mpoly_length(self.val, self.ctx.val):
            raise IndexError("term index out of range")
        res = fmpz_vec(nvars, double_indirect=True)
        fmpz_mod_mpoly_get_term_exp_fmpz(res.double_indirect, self.val, i, self.ctx.val)
        return res.to_tuple()

    def degrees(self):
        """
        Return a dictionary of variable name to degree.

            >>> from flint import Ordering
            >>> ctx = fmpz_mod_mpoly_ctx.get_context(4, Ordering.lex, 11, 'x')
            >>> p = ctx.from_dict({(1, 0, 0, 0): 1, (0, 2, 0, 0): 2, (0, 0, 3, 0): 3})
            >>> p.degrees()
            (1, 2, 3, 0)
        """
        cdef:
            slong nvars = self.ctx.nvars()

        res = fmpz_vec(nvars, double_indirect=True)
        fmpz_mod_mpoly_degrees_fmpz(res.double_indirect, self.val, self.ctx.val)
        return res.to_tuple()

    def total_degree(self):
        """
        Return the total degree.

            >>> from flint import Ordering
            >>> ctx = fmpz_mod_mpoly_ctx.get_context(4, Ordering.lex, 11, 'x')
            >>> p = ctx.from_dict({(1, 0, 0, 0): 1, (0, 2, 0, 0): 2, (0, 0, 3, 0): 3})
            >>> p.total_degree()
            3
        """
        cdef fmpz res = fmpz()
        fmpz_mod_mpoly_total_degree_fmpz((<fmpz> res).val, self.val, self.ctx.val)
        return res

    def repr(self):
        return f"{self.ctx}.from_dict({self.to_dict()})"

    def str(self):
        cdef bytes s = <bytes> fmpz_mod_mpoly_get_str_pretty(self.val, self.ctx.c_names, self.ctx.val)
        res = s.decode().replace("+", " + ").replace("-", " - ")
        if res.startswith(" - "):
            res = "-" + res[3:]
        return res

    def gcd(self, other):
        """
        Return the gcd of self and other.

            >>> from flint import Ordering
            >>> ctx = fmpz_mod_mpoly_ctx.get_context(2, Ordering.lex, 11, 'x')
            >>> f = ctx.from_dict({(1, 1): 4, (0, 0): 1})
            >>> g = ctx.from_dict({(0, 1): 2, (1, 0): 2})
            >>> (f * g).gcd(f)
            4*x0*x1 + 1
        """
        cdef fmpz_mod_mpoly res
        if not self.ctx.is_prime():
            raise DomainError("gcd with non-prime modulus is not supported")
        elif not typecheck(other, fmpz_mod_mpoly):
            raise TypeError("argument must be a fmpz_mod_mpoly")
        elif (<fmpz_mod_mpoly>self).ctx is not (<fmpz_mod_mpoly>other).ctx:
            raise IncompatibleContextError(f"{(<fmpz_mod_mpoly>self).ctx} is not {(<fmpz_mod_mpoly>other).ctx}")
        res = create_fmpz_mod_mpoly(self.ctx)
        fmpz_mod_mpoly_gcd(res.val, (<fmpz_mod_mpoly>self).val, (<fmpz_mod_mpoly>other).val, res.ctx.val)
        return res

    def sqrt(self):
        """
        Return the square root of self.

            >>> from flint import Ordering
            >>> ctx = fmpz_mod_mpoly_ctx.get_context(2, Ordering.lex, 11, 'x')
            >>> f = ctx.from_dict({(1, 1): 4, (0, 0): 1})
            >>> (f * f).sqrt()
            4*x0*x1 + 1
        """
        cdef fmpz_mod_mpoly res
        if not self.ctx.is_prime():
            raise DomainError("square root with non-prime modulus is not supported")

        res = create_fmpz_mod_mpoly(self.ctx)

        if fmpz_mod_mpoly_sqrt(res.val, self.val, self.ctx.val):
            return res
        else:
            raise ValueError("polynomial is not a perfect square")

    def factor(self):
        """
        Factors self into irreducible factors, returning a tuple
        (c, factors) where c is the content of the coefficients and
        factors is a list of (poly, exp) pairs.

            >>> from flint import Ordering
            >>> Zm = fmpz_mod_mpoly
            >>> ctx = fmpz_mod_mpoly_ctx.get_context(3, Ordering.lex, 11, 'x,y,z')
            >>> p1 = Zm("2*x + 4", ctx)
            >>> p2 = Zm("3*x*z + 3*x + 3*z + 3", ctx)
            >>> (p1 * p2).factor()
            (6, [(z + 1, 1), (x + 2, 1), (x + 1, 1)])
            >>> (p2 * p1 * p2).factor()
            (18, [(z + 1, 2), (x + 2, 1), (x + 1, 2)])
        """
        cdef:
            fmpz_mod_mpoly_factor_t fac
            fmpz c
            fmpz_mod_mpoly u

        fmpz_mod_mpoly_factor_init(fac, self.ctx.val)
        if not fmpz_mod_mpoly_factor(fac, self.val, self.ctx.val):
            raise RuntimeError("factorisation failed")
        res = [0] * fac.num

        for i in range(fac.num):
            u = create_fmpz_mod_mpoly(self.ctx)
            fmpz_mod_mpoly_set((<fmpz_mod_mpoly>u).val, &fac.poly[i], self.ctx.val)

            c = fmpz.__new__(fmpz)
            fmpz_set((<fmpz>c).val, &fac.exp[i])

            res[i] = (u, c)

        c = fmpz.__new__(fmpz)
        fmpz_set((<fmpz>c).val, fac.constant)
        fmpz_mod_mpoly_factor_clear(fac, self.ctx.val)
        return c, res

    def factor_squarefree(self):
        """
        Factors self into irreducible factors, returning a tuple
        (c, factors) where c is the content of the coefficients and
        factors is a list of (poly, exp) pairs.

            >>> from flint import Ordering
            >>> Zm = fmpz_mod_mpoly
            >>> ctx = fmpz_mod_mpoly_ctx.get_context(3, Ordering.lex, 11, 'x,y,z')
            >>> p1 = Zm("2*x + 4", ctx)
            >>> p2 = Zm("3*x*y + 3*x + 3*y + 3", ctx)
            >>> (p1 * p2).factor_squarefree()
            (6, [(y + 1, 1), (x^2 + 3*x + 2, 1)])
            >>> (p1 * p2 * p1).factor_squarefree()
            (12, [(y + 1, 1), (x + 1, 1), (x + 2, 2)])
        """
        cdef:
            fmpz_mod_mpoly_factor_t fac
            fmpz c
            fmpz_mod_mpoly u

        fmpz_mod_mpoly_factor_init(fac, self.ctx.val)
        if not fmpz_mod_mpoly_factor_squarefree(fac, self.val, self.ctx.val):
            raise RuntimeError("factorisation failed")
        res = [0] * fac.num

        for i in range(fac.num):
            u = create_fmpz_mod_mpoly(self.ctx)
            fmpz_mod_mpoly_init(u.val, u.ctx.val)
            fmpz_mod_mpoly_set((<fmpz_mod_mpoly>u).val, &fac.poly[i], self.ctx.val)

            c = fmpz.__new__(fmpz)
            fmpz_set((<fmpz>c).val, &fac.exp[i])

            res[i] = (u, c)

        c = fmpz.__new__(fmpz)
        fmpz_set((<fmpz>c).val, fac.constant)
        fmpz_mod_mpoly_factor_clear(fac, self.ctx.val)
        return c, res

    # TODO: Rethink context conversions, particularly the proposed methods in #132
    # def coerce_to_context(self, ctx):
    #     cdef:
    #         fmpz_mod_mpoly res
    #         slong *C
    #         slong i

    #     if not typecheck(ctx, fmpz_mod_mpoly_ctx):
    #         raise ValueError("provided context is not a fmpz_mod_mpoly_ctx")

    #     if self.ctx is ctx:
    #         return self

    #     C = <slong *> libc.stdlib.malloc(self.ctx.val.minfo.nvars * sizeof(slong *))
    #     if C is NULL:
    #         raise MemoryError("malloc returned a null pointer")
    #     res = create_fmpz_mod_mpoly(self.ctx)

    #     vars = {x: i for i, x in enumerate(ctx.py_names)}
    #     for i, var in enumerate(self.ctx.py_names):
    #         C[i] = <slong>vars[var]

    #     fmpz_mod_mpoly_compose_fmpz_mod_mpoly_gen(res.val, self.val, C, self.ctx.val, (<fmpz_mod_mpoly_ctx>ctx).val)

    #     libc.stdlib.free(C)
    #     return res

    def derivative(self, var):
        """
        Return the derivative of this polynomial with respect to the provided variable.
        The argument can either be the variable as a string, or the index of the
        variable in the context.

            >>> from flint import Ordering
            >>> ctx = fmpz_mod_mpoly_ctx.get_context(2, Ordering.lex, 11, 'x')
            >>> p = ctx.from_dict({(0, 3): 2, (2, 1): 3})
            >>> p
            3*x0^2*x1 + 2*x1^3
            >>> p.derivative("x0")
            6*x0*x1
            >>> p.derivative(1)
            3*x0^2 + 6*x1^2

        """
        cdef:
            fmpz_mod_mpoly res
            slong i = self.ctx.variable_to_index(var)

        res = create_fmpz_mod_mpoly(self.ctx)

        fmpz_mod_mpoly_derivative(res.val, self.val, i, self.ctx.val)
        return res


cdef class fmpz_mod_mpoly_vec:
    """
    A class representing a vector of fmpz_mod_mpolys.
    """

    def __cinit__(self, iterable_or_len, fmpz_mod_mpoly_ctx ctx, bint double_indirect = False):
        if isinstance(iterable_or_len, int):
            self.length = iterable_or_len
        else:
            self.length = len(iterable_or_len)

        self.ctx = ctx

        self.val = <fmpz_mod_mpoly_struct *> libc.stdlib.malloc(self.length * sizeof(fmpz_mod_mpoly_struct))
        for i in range(self.length):
            fmpz_mod_mpoly_init(&self.val[i], self.ctx.val)

        if double_indirect:
            self.double_indirect = <fmpz_mod_mpoly_struct **> libc.stdlib.malloc(self.length * sizeof(fmpz_mod_mpoly_struct *))
            if self.double_indirect is NULL:
                raise MemoryError("malloc returned a null pointer")  # pragma: no cover

            for i in range(self.length):
                self.double_indirect[i] = &self.val[i]
        else:
            self.double_indirect = NULL

    def __init__(self, iterable_or_len, _, double_indirect: bool = False):
        if not isinstance(iterable_or_len, int):
            for i, x in enumerate(iterable_or_len):
                self[i] = x

    def __dealloc__(self):
        libc.stdlib.free(self.double_indirect)
        for i in range(self.length):
            fmpz_mod_mpoly_clear(&self.val[i], self.ctx.val)
        libc.stdlib.free(self.val)

    def __getitem__(self, x):
        if not isinstance(x, int):
            raise TypeError("index is not integer")
        elif not 0 <= x < self.length:
            raise IndexError("index out of range")

        cdef fmpz_mod_mpoly z = create_fmpz_mod_mpoly(self.ctx)
        fmpz_mod_mpoly_set(z.val, &self.val[x], self.ctx.val)
        return z

    def __setitem__(self, x, y):
        if not isinstance(x, int):
            raise TypeError("index is not integer")
        elif not 0 <= x < self.length:
            raise IndexError("index out of range")
        elif not typecheck(y, fmpz_mod_mpoly):
            raise TypeError("argument is not fmpz_mod_mpoly")
        elif (<fmpz_mod_mpoly>y).ctx is not self.ctx:
            raise IncompatibleContextError(f"{(<fmpz_mod_mpoly>y).ctx} is not {self.ctx}")

        fmpz_mod_mpoly_set(&self.val[x], (<fmpz_mod_mpoly>y).val, self.ctx.val)

    def __len__(self):
        return self.val.length

    def __str__(self):
        s = [None] * self.length
        for i in range(self.length):
            x = create_fmpz_mod_mpoly(self.ctx)
            fmpz_mod_mpoly_set(x.val, &self.val[x], self.ctx.val)
            s[i] = str(x)
        return f"[{', '.join(s)}]"

    def __repr__(self):
        return f"fmpz_mod_mpoly_vec({self}, ctx={self.ctx})"

    def to_tuple(self):
        return tuple(self[i] for i in range(self.val.length))
