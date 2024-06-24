from flint.flint_base.flint_base cimport (
    flint_mpoly,
    flint_mpoly_context,
    Ordering,
    ordering_py_to_c,
    ordering_c_to_py,
)

from flint.utils.typecheck cimport typecheck
from flint.utils.flint_exceptions import DomainError, IncompatibleContextError

from flint.types.fmpq cimport fmpq, any_as_fmpq
from flint.types.fmpz_vec cimport fmpz_vec
from flint.types.fmpq_vec cimport fmpq_vec

from flint.types.fmpz cimport fmpz, any_as_fmpz
from flint.types.fmpz_mpoly cimport fmpz_mpoly

from flint.flintlib.fmpq cimport fmpq_set, fmpq_one
from flint.flintlib.mpoly cimport ordering_t
from flint.flintlib.fmpq_mpoly cimport (
    fmpq_mpoly_add,
    fmpq_mpoly_add_fmpq,
    fmpq_mpoly_clear,
    fmpq_mpoly_compose_fmpq_mpoly,
    fmpq_mpoly_ctx_init,
    fmpq_mpoly_degrees_fmpz,
    fmpq_mpoly_derivative,
    fmpq_mpoly_div,
    fmpq_mpoly_divides,
    fmpq_mpoly_divrem,
    fmpq_mpoly_equal,
    fmpq_mpoly_evaluate_all_fmpq,
    fmpq_mpoly_evaluate_one_fmpq,
    fmpq_mpoly_gcd,
    fmpq_mpoly_gen,
    fmpq_mpoly_get_coeff_fmpq_fmpz,
    fmpq_mpoly_get_str_pretty,
    fmpq_mpoly_get_term,
    fmpq_mpoly_get_term_coeff_fmpq,
    fmpq_mpoly_get_term_exp_fmpz,
    fmpq_mpoly_integral,
    fmpq_mpoly_is_one,
    fmpq_mpoly_is_zero,
    fmpq_mpoly_length,
    fmpq_mpoly_mul,
    fmpq_mpoly_neg,
    fmpq_mpoly_pow_fmpz,
    fmpq_mpoly_push_term_fmpq_ffmpz,
    fmpq_mpoly_reduce,
    fmpq_mpoly_scalar_div_fmpq,
    fmpq_mpoly_scalar_mul_fmpq,
    fmpq_mpoly_set,
    fmpq_mpoly_set_coeff_fmpq_fmpz,
    fmpq_mpoly_set_fmpq,
    fmpq_mpoly_set_str_pretty,
    fmpq_mpoly_set_term_coeff_fmpq,
    fmpq_mpoly_sort_terms,
    fmpq_mpoly_sqrt,
    fmpq_mpoly_sub,
    fmpq_mpoly_sub_fmpq,
    fmpq_mpoly_total_degree_fmpz,
)
from flint.flintlib.fmpq_mpoly_factor cimport (
    fmpq_mpoly_factor,
    fmpq_mpoly_factor_clear,
    fmpq_mpoly_factor_init,
    fmpq_mpoly_factor_squarefree,
    fmpq_mpoly_factor_t,
)

from flint.flintlib.fmpz cimport fmpz_init_set
from flint.flintlib.fmpz_mpoly cimport fmpz_mpoly_set

from cpython.object cimport Py_EQ, Py_NE
cimport libc.stdlib

cdef dict _fmpq_mpoly_ctx_cache = {}


cdef class fmpq_mpoly_ctx(flint_mpoly_context):
    """
    A class for storing the polynomial context

    :param nvars: The number of variables in the ring
    :param ordering:  The term order for the ring
    :param names:  A tuple containing the names of the variables of the ring.

    Do not construct one of these directly, use `fmpz_mpoly_ctx.get_context`.
    """

    _ctx_cache = _fmpq_mpoly_ctx_cache

    def __init__(self, slong nvars, ordering, names):
        fmpq_mpoly_ctx_init(self.val, nvars, ordering_py_to_c(ordering))
        super().__init__(nvars, names)

    def nvars(self):
        """
        Return the number of variables in the context

            >>> from flint import Ordering
            >>> ctx = fmpq_mpoly_ctx.get_context(4, Ordering.lex, 'x')
            >>> ctx.nvars()
            4
        """
        return self.val.zctx.minfo.nvars

    def ordering(self):
        """
        Return the term order of the context object.

            >>> from flint import Ordering
            >>> ctx = fmpq_mpoly_ctx.get_context(4, Ordering.deglex, 'w')
            >>> ctx.ordering()
            <Ordering.deglex: 1>
        """
        return ordering_c_to_py(self.val.zctx.minfo.ord)

    def gen(self, slong i):
        """
        Return the `i`th generator of the polynomial ring

            >>> from flint import Ordering
            >>> ctx = fmpq_mpoly_ctx.get_context(3, Ordering.degrevlex, 'z')
            >>> ctx.gen(1)
            z1
        """
        cdef fmpq_mpoly res
        if not 0 <= i < self.val.zctx.minfo.nvars:
            raise IndexError("generator index out of range")
        res = create_fmpq_mpoly(self)
        fmpq_mpoly_gen(res.val, i, res.ctx.val)
        return res

    def constant(self, z):
        """
        Create a constant polynomial in this context
        """
        cdef fmpq_mpoly res
        z = any_as_fmpq(z)
        if z is NotImplemented:
            raise TypeError("cannot coerce argument to fmpq")
        res = create_fmpq_mpoly(self)
        fmpq_mpoly_set_fmpq(res.val, (<fmpq>z).val, res.ctx.val)
        return res

    def from_dict(self, d):
        """
        Create a fmpq_mpoly from a dictionary.

        The dictionary's keys are tuples of ints (or anything that implicitly converts
        to fmpz) representing exponents, and corresponding coefficient values of fmpq.

            >>> from flint import Ordering
            >>> ctx = fmpq_mpoly_ctx.get_context(2, Ordering.lex, 'x,y')
            >>> ctx.from_dict({(1,0):2, (1,1):3, (0,1):1})
            3*x*y + 2*x + y
        """
        cdef:
            fmpz_vec exp_vec
            slong nvars = self.nvars()
            fmpq_mpoly res

        if not isinstance(d, dict):
            raise ValueError("Expected a dictionary")

        res = create_fmpq_mpoly(self)

        for k, v in d.items():
            o = any_as_fmpq(v)
            if o is NotImplemented:
                raise TypeError(f"Cannot coerce coefficient '{v}' to fmpq")
            elif len(k) != nvars:
                raise ValueError(f"Expected {nvars} exponents, got {len(k)}")

            exp_vec = fmpz_vec(k)

            if o:
                fmpq_mpoly_push_term_fmpq_ffmpz(res.val, (<fmpq>o).val, exp_vec.val, self.val)

        fmpq_mpoly_sort_terms(res.val, self.val)
        fmpq_mpoly_reduce(res.val, self.val)
        return res


cdef class fmpq_mpoly(flint_mpoly):
    """
    The *fmpq_mpoly* type represents sparse multivariate polynomials over
    the integers.
    """

    def __cinit__(self):
        self._init = False

    def __dealloc__(self):
        if self._init:
            fmpq_mpoly_clear(self.val, self.ctx.val)
            self._init = False

    def __init__(self, val=0, ctx=None):
        if typecheck(val, fmpq_mpoly):
            if ctx is None or ctx == (<fmpq_mpoly>val).ctx:
                init_fmpq_mpoly(self, (<fmpq_mpoly>val).ctx)
                fmpq_mpoly_set(self.val, (<fmpq_mpoly>val).val, self.ctx.val)
            else:
                raise IncompatibleContextError(f"{ctx} is not {(<fmpq_mpoly>val).ctx}")
        elif typecheck(val, fmpz_mpoly):
            if ctx is None:
                ctx = fmpq_mpoly_ctx.from_context((<fmpz_mpoly>val).ctx)
            elif not typecheck(ctx, fmpq_mpoly_ctx) and not typecheck(ctx, fmpz_mpoly_ctx):
                raise TypeError(f"{ctx} is not a fmpq_mpoly_ctx or fmpz_mpoly_ctx")
            elif ctx.nvars() != val.context().nvars():
                raise ValueError(
                    f"Provided context ('{ctx}') and provided fmpz_mpoly ('{val}') don't share the same number of variables"
                )
            init_fmpq_mpoly(self, ctx)
            fmpz_mpoly_set(self.val.zpoly, (<fmpz_mpoly>val).val, (<fmpz_mpoly> val).ctx.val)
            fmpq_one(self.val.content)
            fmpq_mpoly_reduce(self.val, self.ctx.val)
        elif isinstance(val, dict):
            if ctx is None:
                raise ValueError("A context is required to create a fmpq_mpoly from a dict")
            x = ctx.from_dict(val)
            # XXX: this copy is silly, have a ctx function that assigns an fmpz_mpoly_t
            init_fmpq_mpoly(self, ctx)
            fmpq_mpoly_set(self.val, (<fmpq_mpoly>x).val, self.ctx.val)
        elif isinstance(val, str):
            if ctx is None:
                raise ValueError("Cannot parse a polynomial without context")
            val = val.encode("ascii")
            init_fmpq_mpoly(self, ctx)
            fmpq_mpoly_set_str_pretty(self.val, val, self.ctx.c_names, self.ctx.val)
            fmpq_mpoly_sort_terms(self.val, self.ctx.val)
        else:
            v = any_as_fmpq(val)
            if v is NotImplemented:
                raise TypeError("cannot create fmpz_mpoly from type %s" % type(val))
            if ctx is None:
                raise ValueError("Need context to convert  fmpz to fmpq_mpoly")
            init_fmpq_mpoly(self, ctx)
            fmpq_mpoly_set_fmpq(self.val, (<fmpq>v).val, self.ctx.val)

    def __bool__(self):
        return not fmpq_mpoly_is_zero(self.val, self.ctx.val)

    def __richcmp__(self, other, int op):
        if not (op == Py_EQ or op == Py_NE):
            return NotImplemented
        elif other is None:
            return op == Py_NE
        elif typecheck(self, fmpq_mpoly) and typecheck(other, fmpq_mpoly):
            if (<fmpq_mpoly>self).ctx is (<fmpq_mpoly>other).ctx:
                return (op == Py_NE) ^ bool(
                    fmpq_mpoly_equal((<fmpq_mpoly>self).val, (<fmpq_mpoly>other).val, (<fmpq_mpoly>self).ctx.val)
                )
            else:
                return op == Py_NE
        else:
            return NotImplemented

    def __len__(self):
        return fmpq_mpoly_length(self.val, self.ctx.val)

    def __getitem__(self, x):
        """
        Return the term at index `x` if `x` is an `int`, or the term with the exponent
        vector `x` if `x` is a tuple of `int`s or `fmpq`s.

            >>> from flint import Ordering
            >>> ctx = fmpq_mpoly_ctx.get_context(2, Ordering.lex, 'x')
            >>> p = ctx.from_dict({(0, 1): 2, (1, 1): 3})
            >>> p[1]
            2*x1
            >>> p[1, 1]
            3

        """
        cdef:
            slong nvars = self.ctx.nvars()

        if isinstance(x, int):
            if not 0 <= x < fmpq_mpoly_length(self.val, self.ctx.val):
                raise IndexError("term index out of range")
            res = create_fmpq_mpoly(self.ctx)
            fmpq_mpoly_get_term((<fmpq_mpoly> res).val, self.val, x, self.ctx.val)
            return res
        elif isinstance(x, tuple):
            if len(x) != nvars:
                raise ValueError("exponent vector provided does not match number of variables")
            res = fmpq()
            exp_vec = fmpz_vec(x, double_indirect=True)
            fmpq_mpoly_get_coeff_fmpq_fmpz((<fmpq>res).val, self.val, exp_vec.double_indirect, self.ctx.val)
            return res
        else:
            raise TypeError("index is not integer or tuple")

    def __setitem__(self, x, y):
        """
        Set the coefficient at index `x` to `y` if `x` is an `int`, or the term with
        the exponent vector `x` if `x` is a tuple of `int`s or `fmpq`s.

            >>> from flint import Ordering
            >>> ctx = fmpq_mpoly_ctx.get_context(2, Ordering.lex, 'x')
            >>> p = ctx.from_dict({(0, 1): 2, (1, 1): 3})
            >>> p[1] = 10
            >>> p[1, 1] = 20
            >>> p
            20*x0*x1 + 10*x1

        """
        cdef:
            slong nvars = self.ctx.nvars()

        coeff = any_as_fmpq(y)
        if coeff is NotImplemented:
            raise TypeError("provided coefficient not coercible to fmpq")

        if isinstance(x, int):
            if not 0 <= x < fmpq_mpoly_length(self.val, self.ctx.val):
                raise IndexError("term index out of range")
            fmpq_mpoly_set_term_coeff_fmpq(self.val, x, (<fmpq>coeff).val, self.ctx.val)
        elif isinstance(x, tuple):
            if len(x) != nvars:
                raise ValueError("exponent vector provided does not match number of variables")
            exp_vec = fmpz_vec(x, double_indirect=True)
            fmpq_mpoly_set_coeff_fmpq_fmpz(self.val, (<fmpq>coeff).val, exp_vec.double_indirect, self.ctx.val)
        else:
            raise TypeError("index is not integer or tuple")

    def __pos__(self):
        return self

    def __neg__(self):
        cdef fmpq_mpoly res
        res = create_fmpq_mpoly(self.ctx)
        fmpq_mpoly_neg(res.val, (<fmpq_mpoly>self).val, res.ctx.val)
        return res

    def __add__(self, other):
        cdef fmpq_mpoly res
        if typecheck(other, fmpq_mpoly):
            if (<fmpq_mpoly>self).ctx is not (<fmpq_mpoly>other).ctx:
                raise IncompatibleContextError(f"{(<fmpq_mpoly>self).ctx} is not {(<fmpq_mpoly>other).ctx}")
            res = create_fmpq_mpoly(self.ctx)
            fmpq_mpoly_add(res.val, (<fmpq_mpoly>self).val, (<fmpq_mpoly>other).val, res.ctx.val)
            return res
        else:
            other = any_as_fmpq(other)
            if other is not NotImplemented:
                res = create_fmpq_mpoly(self.ctx)
                fmpq_mpoly_add_fmpq(res.val, (<fmpq_mpoly>self).val, (<fmpq>other).val, res.ctx.val)
                return res
        return NotImplemented

    def __radd__(self, other):
        cdef fmpq_mpoly res
        other = any_as_fmpq(other)
        if other is not NotImplemented:
            res = create_fmpq_mpoly(self.ctx)
            fmpq_mpoly_add_fmpq(res.val, (<fmpq_mpoly>self).val, (<fmpq>other).val, res.ctx.val)
            return res
        return NotImplemented

    def __iadd__(self, other):
        if typecheck(other, fmpq_mpoly):
            if (<fmpq_mpoly>self).ctx is not (<fmpq_mpoly>other).ctx:
                raise IncompatibleContextError(f"{(<fmpq_mpoly>self).ctx} is not {(<fmpq_mpoly>other).ctx}")
            fmpq_mpoly_add((<fmpq_mpoly>self).val, (<fmpq_mpoly>self).val, (<fmpq_mpoly>other).val, self.ctx.val)
            return self
        else:
            other = any_as_fmpq(other)
            if other is not NotImplemented:
                fmpq_mpoly_add_fmpq((<fmpq_mpoly>self).val, (<fmpq_mpoly>self).val, (<fmpq>other).val, self.ctx.val)
                return self
        return NotImplemented

    def __sub__(self, other):
        cdef fmpq_mpoly res
        if typecheck(other, fmpq_mpoly):
            if (<fmpq_mpoly>self).ctx is not (<fmpq_mpoly>other).ctx:
                raise IncompatibleContextError(f"{(<fmpq_mpoly>self).ctx} is not {(<fmpq_mpoly>other).ctx}")
            res = create_fmpq_mpoly(self.ctx)
            fmpq_mpoly_sub(res.val, (<fmpq_mpoly>self).val, (<fmpq_mpoly>other).val, res.ctx.val)
            return res
        else:
            other = any_as_fmpq(other)
            if other is not NotImplemented:
                res = create_fmpq_mpoly(self.ctx)
                fmpq_mpoly_sub_fmpq(res.val, (<fmpq_mpoly>self).val, (<fmpq>other).val, res.ctx.val)
                return res
        return NotImplemented

    def __rsub__(self, other):
        cdef fmpq_mpoly res
        other = any_as_fmpq(other)
        if other is not NotImplemented:
            res = create_fmpq_mpoly(self.ctx)
            fmpq_mpoly_sub_fmpq(res.val, (<fmpq_mpoly>self).val, (<fmpq>other).val, res.ctx.val)
            return -res
        return NotImplemented

    def __isub__(self, other):
        if typecheck(other, fmpq_mpoly):
            if (<fmpq_mpoly>self).ctx is not (<fmpq_mpoly>other).ctx:
                raise IncompatibleContextError(f"{(<fmpq_mpoly>self).ctx} is not {(<fmpq_mpoly>other).ctx}")
            fmpq_mpoly_sub((<fmpq_mpoly>self).val, (<fmpq_mpoly>self).val, (<fmpq_mpoly>other).val, self.ctx.val)
            return self
        else:
            other = any_as_fmpq(other)
            if other is not NotImplemented:
                fmpq_mpoly_sub_fmpq((<fmpq_mpoly>self).val, (<fmpq_mpoly>self).val, (<fmpq>other).val, self.ctx.val)
                return self
        return NotImplemented

    def __mul__(self, other):
        cdef fmpq_mpoly res
        if typecheck(other, fmpq_mpoly):
            if (<fmpq_mpoly>self).ctx is not (<fmpq_mpoly>other).ctx:
                raise IncompatibleContextError(f"{(<fmpq_mpoly>self).ctx} is not {(<fmpq_mpoly>other).ctx}")
            res = create_fmpq_mpoly(self.ctx)
            fmpq_mpoly_mul(res.val, (<fmpq_mpoly>self).val, (<fmpq_mpoly>other).val, res.ctx.val)
            return res
        else:
            other = any_as_fmpq(other)
            if other is not NotImplemented:
                res = create_fmpq_mpoly(self.ctx)
                fmpq_mpoly_scalar_mul_fmpq(res.val, (<fmpq_mpoly>self).val, (<fmpq>other).val, res.ctx.val)
                return res
        return NotImplemented

    def __rmul__(self, other):
        cdef fmpq_mpoly res
        other = any_as_fmpq(other)
        if other is not NotImplemented:
            res = create_fmpq_mpoly(self.ctx)
            fmpq_mpoly_scalar_mul_fmpq(res.val, (<fmpq_mpoly>self).val, (<fmpq>other).val, res.ctx.val)
            return res
        return NotImplemented

    def __imul__(self, other):
        if typecheck(other, fmpq_mpoly):
            if (<fmpq_mpoly>self).ctx is not (<fmpq_mpoly>other).ctx:
                raise IncompatibleContextError(f"{(<fmpq_mpoly>self).ctx} is not {(<fmpq_mpoly>other).ctx}")
            fmpq_mpoly_mul((<fmpq_mpoly>self).val, (<fmpq_mpoly>self).val, (<fmpq_mpoly>other).val, self.ctx.val)
            return self
        else:
            other = any_as_fmpq(other)
            if other is not NotImplemented:
                fmpq_mpoly_scalar_mul_fmpq(self.val, (<fmpq_mpoly>self).val, (<fmpq>other).val, self.ctx.val)
                return self
        return NotImplemented

    def __pow__(self, other, modulus):
        cdef fmpq_mpoly res
        if modulus is not None:
            raise NotImplementedError
        other = any_as_fmpz(other)
        if other is NotImplemented:
            return other
        if other < 0:
            raise ValueError("cannot raise fmpq_mpoly to negative power")
        res = create_fmpq_mpoly(self.ctx)
        if fmpq_mpoly_pow_fmpz(res.val, (<fmpq_mpoly>self).val, (<fmpz>other).val, res.ctx.val) == 0:
            raise ValueError("unreasonably large polynomial")  # pragma: no cover
        return res

    def __divmod__(self, other):
        cdef fmpq_mpoly res, res2
        if typecheck(other, fmpq_mpoly):
            if not other:
                raise ZeroDivisionError("fmpq_mpoly division by zero")
            elif (<fmpq_mpoly>self).ctx is not (<fmpq_mpoly>other).ctx:
                raise IncompatibleContextError(f"{(<fmpq_mpoly>self).ctx} is not {(<fmpq_mpoly>other).ctx}")
            res = create_fmpq_mpoly(self.ctx)
            res2 = create_fmpq_mpoly(self.ctx)
            fmpq_mpoly_divrem(res.val, res2.val, (<fmpq_mpoly>self).val, (<fmpq_mpoly>other).val, res.ctx.val)
            return (res, res2)
        else:
            other = any_as_fmpq(other)
            if other is not NotImplemented:
                other = fmpq_mpoly(other, self.ctx)
                if not other:
                    raise ZeroDivisionError("fmpq_mpoly division by zero")
                res = create_fmpq_mpoly(self.ctx)
                res2 = create_fmpq_mpoly(self.ctx)
                fmpq_mpoly_divrem(res.val, res2.val, (<fmpq_mpoly>self).val, (<fmpq_mpoly>other).val, res.ctx.val)
                return (res, res2)
        return NotImplemented

    def __rdivmod__(self, other):
        cdef fmpq_mpoly res, res2
        if not self:
            raise ZeroDivisionError("fmpq_mpoly division by zero")
        other = any_as_fmpq(other)
        if other is not NotImplemented:
            other = fmpq_mpoly(other, self.ctx)
            res = create_fmpq_mpoly(self.ctx)
            res2 = create_fmpq_mpoly(self.ctx)
            fmpq_mpoly_divrem(res.val, res2.val, (<fmpq_mpoly>other).val, (<fmpq_mpoly>self).val, res.ctx.val)
            return (res, res2)
        return NotImplemented

    def __floordiv__(self, other):
        cdef fmpq_mpoly res
        if typecheck(other, fmpq_mpoly):
            if not other:
                raise ZeroDivisionError("fmpq_mpoly division by zero")
            elif (<fmpq_mpoly>self).ctx is not (<fmpq_mpoly>other).ctx:
                raise IncompatibleContextError(f"{(<fmpq_mpoly>self).ctx} is not {(<fmpq_mpoly>other).ctx}")
            res = create_fmpq_mpoly(self.ctx)
            fmpq_mpoly_div(res.val, (<fmpq_mpoly>self).val, (<fmpq_mpoly>other).val, res.ctx.val)
            return res
        else:
            other = any_as_fmpq(other)
            if other is not NotImplemented:
                if not other:
                    raise ZeroDivisionError("fmpq_mpoly division by zero")
                other = fmpq_mpoly(other, self.ctx)
                res = create_fmpq_mpoly(self.ctx)
                fmpq_mpoly_div(res.val, (<fmpq_mpoly>self).val, (<fmpq_mpoly>other).val, res.ctx.val)
                return res
        return NotImplemented

    def __rfloordiv__(self, other):
        cdef fmpq_mpoly res
        if not self:
            raise ZeroDivisionError("fmpq_mpoly division by zero")
        other = any_as_fmpq(other)
        if other is not NotImplemented:
            other = fmpq_mpoly(other, self.ctx)
            res = create_fmpq_mpoly(self.ctx)
            fmpq_mpoly_div(res.val, (<fmpq_mpoly>other).val, self.val, res.ctx.val)
            return res
        return NotImplemented

    def __truediv__(self, other):
        cdef:
            fmpq_mpoly res

        if typecheck(other, fmpq_mpoly):
            if not other:
                raise ZeroDivisionError("fmpq_mpoly division by zero")
            elif self.ctx is not (<fmpq_mpoly>other).ctx:
                raise IncompatibleContextError(f"{self.ctx} is not {(<fmpq_mpoly>other).ctx}")

            res = create_fmpq_mpoly(self.ctx)
            if fmpq_mpoly_divides(res.val, self.val, (<fmpq_mpoly>other).val, self.ctx.val):
                return res
            else:
                raise DomainError("fmpq_mpoly division is not exact")
        else:
            o = any_as_fmpq(other)
            if o is NotImplemented:
                return NotImplemented
            elif not o:
                raise ZeroDivisionError("fmpq_mpoly division by zero")
            res = create_fmpq_mpoly(self.ctx)
            fmpq_mpoly_scalar_div_fmpq(res.val, self.val, (<fmpq>o).val, self.ctx.val)
            return res

    def __rtruediv__(self, other):
        cdef fmpq_mpoly res
        if not self:
            raise ZeroDivisionError("fmpq_mpoly division by zero")
        o = any_as_fmpq(other)
        if o is NotImplemented:
            return NotImplemented
        res = create_fmpq_mpoly(self.ctx)
        fmpq_mpoly_set_fmpq(res.val, (<fmpq>o).val, self.ctx.val)
        return res / self

    def __mod__(self, other):
        return divmod(self, other)[1]

    def __rmod__(self, other):
        return divmod(other, self)[1]

    def __call__(self, *args) -> fmpq:
        cdef:
            fmpq_vec V
            fmpq vres
            slong nvars = self.ctx.nvars(), nargs = len(args)

        if nargs < nvars:
            raise ValueError("Not enough arguments provided")
        elif nargs > nvars:
            raise ValueError("More arguments provided than variables")

        args_fmpq = tuple(any_as_fmpq(v) for v in args)
        for arg in args_fmpq:
            if arg is NotImplemented:
                raise TypeError(f"Cannot coerce argument ('{arg}') to fmpq")

        V = fmpq_vec(args_fmpq, double_indirect=True)
        vres = fmpq.__new__(fmpq)
        if fmpq_mpoly_evaluate_all_fmpq(vres.val, self.val, V.double_indirect, self.ctx.val) == 0:
            raise ValueError("Unreasonably large polynomial")  # pragma: no cover
        return vres

    def subs(self, dict_args) -> fmpq_mpoly:
        """
        Partial evaluate this polynomial.
        """
        cdef:
            fmpq_mpoly res
            fmpq_mpoly res2
            slong i, nargs

        partial_args = tuple((i, dict_args[x]) for i, x in enumerate(self.ctx.names()) if x in dict_args)
        nargs = len(partial_args)

        # If we've been provided with an invalid keyword arg then the length of our filter
        # args will be less than what we've been provided with.
        # If the length is equal to the number of variables then all arguments have been provided
        # otherwise we need to do partial application
        if nargs < len(dict_args):
            raise ValueError("Unknown keyword argument provided")

        args_fmpq = tuple(any_as_fmpq(v) for _, v in partial_args)
        for arg in args_fmpq:
            if arg is NotImplemented:
                raise TypeError(f"Cannot coerce argument ('{arg}') to fmpq")

        # Partial application with args in Z. We evaluate the polynomial one variable at a time
        res = create_fmpq_mpoly(self.ctx)
        res2 = create_fmpq_mpoly(self.ctx)

        fmpq_mpoly_set(res2.val, self.val, self.ctx.val)
        for (i, _), arg in zip(partial_args, args_fmpq):
            if fmpq_mpoly_evaluate_one_fmpq(res.val, res2.val, i, (<fmpq>arg).val, self.ctx.val) == 0:
                raise ValueError("Unreasonably large polynomial")  # pragma: no cover
            fmpq_mpoly_set(res2.val, res.val, self.ctx.val)
        return res

    def compose(self, *args) -> fmpq_mpoly:
        """
        Compose this polynomial with other fmpq_mpolys. Takes positional arguments.
        """
        cdef:
            fmpq_mpoly res
            fmpq_mpoly_ctx res_ctx
            fmpq_mpoly_vec C
            slong i, nvars = self.ctx.nvars(), nargs = len(args)

        if nargs < nvars:
            raise ValueError("Not enough arguments provided")
        elif nargs > nvars:
            raise ValueError("More arguments provided than variables")
        elif not all(typecheck(arg, fmpq_mpoly) for arg in args):
            raise TypeError("All arguments must be fmpq_mpolys")

        res_ctx = (<fmpq_mpoly> args[0]).ctx
        if not all((<fmpq_mpoly> args[i]).ctx is res_ctx for i in range(1, len(args))):
            raise IncompatibleContextError("All arguments must share the same context")

        C = fmpq_mpoly_vec(args, self.ctx, double_indirect=True)
        res = create_fmpq_mpoly(self.ctx)
        if fmpq_mpoly_compose_fmpq_mpoly(res.val, self.val, C.double_indirect, self.ctx.val, res_ctx.val) == 0:
            raise ValueError("Unreasonably large polynomial")  # pragma: no cover
        return res

    def context(self):
        """
        Return the context object for this polynomials.

            >>> from flint import Ordering
            >>> ctx = fmpq_mpoly_ctx.get_context(2, Ordering.lex, 'x')
            >>> p = ctx.from_dict({(0, 1): 2})
            >>> ctx is p.context()
            True
        """
        return self.ctx

    def is_one(self):
        return fmpq_mpoly_is_one(self.val, self.ctx.val)

    def coefficient(self, slong i):
        """
        Return the coefficient at index `i`.

            >>> from flint import Ordering
            >>> ctx = fmpq_mpoly_ctx.get_context(2, Ordering.lex, 'x')
            >>> p = ctx.from_dict({(0, 1): 2, (1, 1): 3})
            >>> p.coefficient(1)
            2
        """
        cdef fmpq v
        if not 0 <= i < fmpq_mpoly_length(self.val, self.ctx.val):
            raise IndexError("term index out of range")
        else:
            v = fmpq.__new__(fmpq)
            fmpq_mpoly_get_term_coeff_fmpq(v.val, self.val, i, self.ctx.val)
            return v

    def exponent_tuple(self, slong i):
        """
        Return the exponent vector at index `i` as a tuple.

            >>> from flint import Ordering
            >>> ctx = fmpq_mpoly_ctx.get_context(2, Ordering.lex, 'x')
            >>> p = ctx.from_dict({(0, 1): 2, (1, 1): 3})
            >>> p.exponent_tuple(1)
            (0, 1)
        """
        cdef:
            slong nvars = self.ctx.nvars()

        if not 0 <= i < fmpq_mpoly_length(self.val, self.ctx.val):
            raise IndexError("term index out of range")
        res = fmpz_vec(nvars, double_indirect=True)
        fmpq_mpoly_get_term_exp_fmpz(res.double_indirect, self.val, i, self.ctx.val)
        return res.to_tuple()

    def degrees(self):
        """
        Return a dictionary of variable name to degree.

            >>> from flint import Ordering
            >>> ctx = fmpq_mpoly_ctx.get_context(4, Ordering.lex, 'x')
            >>> p = ctx.from_dict({(1, 0, 0, 0): 1, (0, 2, 0, 0): 2, (0, 0, 3, 0): 3})
            >>> p.degrees()
            {'x0': 1, 'x1': 2, 'x2': 3, 'x3': 0}
        """
        cdef:
            slong nvars = self.ctx.nvars()

        res = fmpz_vec(nvars, double_indirect=True)
        fmpq_mpoly_degrees_fmpz(res.double_indirect, self.val, self.ctx.val)
        return dict(zip(self.ctx.names(), res.to_tuple()))

    def total_degree(self):
        """
        Return the total degree.

            >>> from flint import Ordering
            >>> ctx = fmpq_mpoly_ctx.get_context(4, Ordering.lex, 'x')
            >>> p = ctx.from_dict({(1, 0, 0, 0): 1, (0, 2, 0, 0): 2, (0, 0, 3, 0): 3})
            >>> p.total_degree()
            3
        """
        cdef fmpz res = fmpz()
        fmpq_mpoly_total_degree_fmpz((<fmpz> res).val, self.val, self.ctx.val)
        return res

    def repr(self):
        return f"{self.ctx}.from_dict({self.to_dict()})"

    def str(self):
        cdef bytes s = <bytes> fmpq_mpoly_get_str_pretty(self.val, self.ctx.c_names, self.ctx.val)
        # res = s.decode().replace("+", " + ").replace("-", " - ")
        res = s.decode()
        if res.startswith("- "):
            res = "-" + res[3:]
        return res

    def gcd(self, other):
        """
        Return the gcd of self and other.

            >>> from flint import Ordering
            >>> ctx = fmpq_mpoly_ctx.get_context(2, Ordering.lex, 'x')
            >>> f = ctx.from_dict({(1, 1): 4, (0, 0): 1})
            >>> g = ctx.from_dict({(0, 1): 2, (1, 0): 2})
            >>> (f * g).gcd(f)
            x0*x1 + 1/4

        """
        cdef fmpq_mpoly res
        if not typecheck(other, fmpq_mpoly):
            raise TypeError("argument must be a fmpq_mpoly")
        elif (<fmpq_mpoly>self).ctx is not (<fmpq_mpoly>other).ctx:
            raise IncompatibleContextError(f"{(<fmpq_mpoly>self).ctx} is not {(<fmpq_mpoly>other).ctx}")
        res = create_fmpq_mpoly(self.ctx)
        fmpq_mpoly_gcd(res.val, (<fmpq_mpoly>self).val, (<fmpq_mpoly>other).val, res.ctx.val)
        return res

    def sqrt(self):
        """
        Return the square root of self.

            >>> from flint import Ordering
            >>> ctx = fmpq_mpoly_ctx.get_context(2, Ordering.lex, 'x')
            >>> f = ctx.from_dict({(1, 1): 4, (0, 0): 1})
            >>> (f * f).sqrt()
            4*x0*x1 + 1
        """
        cdef fmpq_mpoly res
        res = create_fmpq_mpoly(self.ctx)

        if fmpq_mpoly_sqrt(res.val, self.val, self.ctx.val):
            return res
        else:
            raise ValueError("polynomial is not a perfect square")

    def factor(self):
        """
        Factors self into irreducible factors, returning a tuple
        (c, factors) where c is the content of the coefficients and
        factors is a list of (poly, exp) pairs.

            >>> from flint import Ordering
            >>> Zm = fmpq_mpoly
            >>> ctx = fmpq_mpoly_ctx.get_context(3, Ordering.lex, 'x,y,z')
            >>> p1 = Zm("2*x + 4", ctx)
            >>> p2 = Zm("3*x*z +  + 3*x + 3*z + 3", ctx)
            >>> (p1 * p2).factor()
            (6, [(z + 1, 1), (x + 2, 1), (x + 1, 1)])
            >>> (p2 * p1 * p2).factor()
            (18, [(z + 1, 2), (x + 2, 1), (x + 1, 2)])
        """
        cdef:
            fmpq_mpoly_factor_t fac
            fmpq_mpoly u

        fmpq_mpoly_factor_init(fac, self.ctx.val)
        fmpq_mpoly_factor(fac, self.val, self.ctx.val)
        res = [0] * fac.num

        for i in range(fac.num):
            u = create_fmpq_mpoly(self.ctx)
            fmpq_mpoly_set((<fmpq_mpoly>u).val, &fac.poly[i], self.ctx.val)

            c = fmpz.__new__(fmpz)
            fmpz_init_set((<fmpz>c).val, &fac.exp[i])

            res[i] = (u, c)

        c = fmpq.__new__(fmpq)
        fmpq_set((<fmpq>c).val, fac.constant)
        fmpq_mpoly_factor_clear(fac, self.ctx.val)
        return c, res

    def factor_squarefree(self):
        """
        Factors self into irreducible factors, returning a tuple
        (c, factors) where c is the content of the coefficients and
        factors is a list of (poly, exp) pairs.

            >>> from flint import Ordering
            >>> Zm = fmpq_mpoly
            >>> ctx = fmpq_mpoly_ctx.get_context(3, Ordering.lex, 'x,y,z')
            >>> p1 = Zm("2*x + 4", ctx)
            >>> p2 = Zm("3*x*y + 3*x + 3*y + 3", ctx)
            >>> (p1 * p2).factor_squarefree()
            (6, [(y + 1, 1), (x^2 + 3*x + 2, 1)])
            >>> (p1 * p2 * p1).factor_squarefree()
            (12, [(y + 1, 1), (x + 1, 1), (x + 2, 2)])
        """
        cdef:
            fmpq_mpoly_factor_t fac
            fmpq_mpoly u

        fmpq_mpoly_factor_init(fac, self.ctx.val)
        fmpq_mpoly_factor_squarefree(fac, self.val, self.ctx.val)
        res = [0] * fac.num

        for i in range(fac.num):
            u = create_fmpq_mpoly(self.ctx)
            fmpq_mpoly_init(u.val, u.ctx.val)
            fmpq_mpoly_set((<fmpq_mpoly>u).val, &fac.poly[i], self.ctx.val)

            c = fmpz.__new__(fmpz)
            fmpz_init_set((<fmpz>c).val, &fac.exp[i])

            res[i] = (u, c)

        c = fmpq.__new__(fmpq)
        fmpq_set((<fmpq>c).val, fac.constant)
        fmpq_mpoly_factor_clear(fac, self.ctx.val)
        return c, res

    def derivative(self, var):
        """
        Return the derivative of this polynomial with respect to the provided variable.
        The argument can either be the variable as a string, or the index of the
        variable in the context.

            >>> from flint import Ordering
            >>> ctx = fmpq_mpoly_ctx.get_context(2, Ordering.lex, 'x')
            >>> p = ctx.from_dict({(0, 3): 2, (2, 1): 3})
            >>> p
            3*x0^2*x1 + 2*x1^3
            >>> p.derivative("x0")
            6*x0*x1
            >>> p.derivative(1)
            3*x0^2 + 6*x1^2

        """
        cdef:
            fmpq_mpoly res
            slong i = self.ctx.variable_to_index(var)

        res = create_fmpq_mpoly(self.ctx)

        fmpq_mpoly_derivative(res.val, self.val, i, self.ctx.val)
        return res

    def integral(self, var):
        """
        Return the integral of this polynomial with respect to the provided variable The argument can either be the
        variable as a string, or the index of the variable in the context.

            >>> from flint import Ordering
            >>> ctx = fmpq_mpoly_ctx.get_context(2, Ordering.lex, 'x')
            >>> p = ctx.from_dict({(0, 3): 2, (2, 1): 3})
            >>> p
            3*x0^2*x1 + 2*x1^3
            >>> p.integral("x0")
            x0^3*x1 + 2*x0*x1^3
            >>> p.integral(1)
            3/2*x0^2*x1^2 + 1/2*x1^4

        """
        cdef:
            fmpq_mpoly res
            slong i = self.ctx.variable_to_index(var)

        res = create_fmpq_mpoly(self.ctx)

        fmpq_mpoly_integral(res.val, self.val, i, self.ctx.val)
        return res


cdef class fmpq_mpoly_vec:
    """
    A class representing a vector of fmpq_mpolys. Not present in FLINT.
    """

    def __cinit__(self, iterable_or_len, fmpq_mpoly_ctx ctx, bint double_indirect = False):
        if isinstance(iterable_or_len, int):
            self.length = iterable_or_len
        else:
            self.length = len(iterable_or_len)

        self.ctx = ctx

        self.val = <fmpq_mpoly_struct *> libc.stdlib.malloc(self.length * sizeof(fmpq_mpoly_struct))
        for i in range(self.length):
            fmpq_mpoly_init(&self.val[i], self.ctx.val)

        if double_indirect:
            self.double_indirect = <fmpq_mpoly_struct **> libc.stdlib.malloc(self.length * sizeof(fmpq_mpoly_struct *))
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
            fmpq_mpoly_clear(&self.val[i], self.ctx.val)
        libc.stdlib.free(self.val)

    def __getitem__(self, x):
        if not isinstance(x, int):
            raise TypeError("index is not integer")
        elif not 0 <= x < self.length:
            raise IndexError("index out of range")

        cdef fmpq_mpoly z = create_fmpq_mpoly(self.ctx)
        fmpq_mpoly_set(z.val, &self.val[x], self.ctx.val)
        return z

    def __setitem__(self, x, y):
        if not isinstance(x, int):
            raise TypeError("index is not integer")
        elif not 0 <= x < self.length:
            raise IndexError("index out of range")
        elif not typecheck(y, fmpq_mpoly):
            raise TypeError("argument is not fmpq_mpoly")
        elif (<fmpq_mpoly>y).ctx is not self.ctx:
            raise IncompatibleContextError(f"{(<fmpq_mpoly>y).ctx} is not {self.ctx}")

        fmpq_mpoly_set(&self.val[x], (<fmpq_mpoly>y).val, self.ctx.val)

    def __len__(self):
        return self.length

    def __str__(self):
        s = [None] * self.length
        for i in range(self.length):
            x = create_fmpq_mpoly(self.ctx)
            fmpq_mpoly_set((<fmpq_mpoly> x).val, &self.val[i], self.ctx.val)
            s[i] = str(x)
        return f"[{', '.join(s)}]"

    def __repr__(self):
        return f"fmpq_mpoly_vec({self}, ctx={self.ctx})"

    def to_tuple(self):
        return tuple(self[i] for i in range(self.length))
