from flint.flint_base.flint_base cimport (
    flint_mpoly,
    flint_mpoly_context,
    Ordering,
    ordering_py_to_c,
    ordering_c_to_py,
)

from flint.utils.typecheck cimport typecheck
from flint.utils.flint_exceptions import DomainError, IncompatibleContextError

from flint.types.fmpz cimport any_as_fmpz, fmpz
from flint.types.fmpz_vec cimport fmpz_vec

from flint.flintlib.fmpz cimport fmpz_set
from flint.flintlib.fmpz_mpoly cimport (
    fmpz_mpoly_add,
    fmpz_mpoly_add_fmpz,
    fmpz_mpoly_clear,
    fmpz_mpoly_compose_fmpz_mpoly,
    fmpz_mpoly_ctx_init,
    fmpz_mpoly_degrees_fmpz,
    fmpz_mpoly_derivative,
    fmpz_mpoly_div,
    fmpz_mpoly_divrem,
    fmpz_mpoly_equal,
    fmpz_mpoly_evaluate_all_fmpz,
    fmpz_mpoly_evaluate_one_fmpz,
    fmpz_mpoly_gcd,
    fmpz_mpoly_gen,
    fmpz_mpoly_get_coeff_fmpz_fmpz,
    fmpz_mpoly_get_str_pretty,
    fmpz_mpoly_get_term,
    fmpz_mpoly_get_term_coeff_fmpz,
    fmpz_mpoly_get_term_exp_fmpz,
    fmpz_mpoly_integral,
    fmpz_mpoly_is_one,
    fmpz_mpoly_is_zero,
    fmpz_mpoly_length,
    fmpz_mpoly_mul,
    fmpz_mpoly_neg,
    fmpz_mpoly_pow_fmpz,
    fmpz_mpoly_push_term_fmpz_ffmpz,
    fmpz_mpoly_scalar_divides_fmpz,
    fmpz_mpoly_scalar_mul_fmpz,
    fmpz_mpoly_set,
    fmpz_mpoly_set_coeff_fmpz_fmpz,
    fmpz_mpoly_set_fmpz,
    fmpz_mpoly_set_str_pretty,
    fmpz_mpoly_set_term_coeff_fmpz,
    fmpz_mpoly_sort_terms,
    fmpz_mpoly_sqrt_heap,
    fmpz_mpoly_sub,
    fmpz_mpoly_sub_fmpz,
    fmpz_mpoly_total_degree_fmpz,
    fmpz_mpoly_vec_clear,
    fmpz_mpoly_vec_entry,
    fmpz_mpoly_vec_init,
)
from flint.flintlib.fmpz_mpoly_factor cimport (
    fmpz_mpoly_factor,
    fmpz_mpoly_factor_clear,
    fmpz_mpoly_factor_init,
    fmpz_mpoly_factor_squarefree,
    fmpz_mpoly_factor_t,
)

from cpython.object cimport Py_EQ, Py_NE
cimport cython
cimport libc.stdlib

cdef dict _fmpz_mpoly_ctx_cache = {}


cdef class fmpz_mpoly_ctx(flint_mpoly_context):
    """
    A class for storing the polynomial context

    :param nvars: The number of variables in the ring
    :param ordering:  The term order for the ring
    :param names:  A tuple containing the names of the variables of the ring.

    Do not construct one of these directly, use `fmpz_mpoly_ctx.get_context`.
    """

    _ctx_cache = _fmpz_mpoly_ctx_cache

    def __init__(self, slong nvars, ordering, names):
        fmpz_mpoly_ctx_init(self.val, nvars, ordering_py_to_c(ordering))
        super().__init__(nvars, names)

    def nvars(self):
        """
        Return the number of variables in the context

            >>> from flint import Ordering
            >>> ctx = fmpz_mpoly_ctx.get_context(4, Ordering.lex, 'x')
            >>> ctx.nvars()
            4
        """
        return self.val.minfo.nvars

    def ordering(self):
        """
        Return the term order of the context object.

            >>> from flint import Ordering
            >>> ctx = fmpz_mpoly_ctx.get_context(4, Ordering.deglex, 'w')
            >>> ctx.ordering()
            <Ordering.deglex: 1>
        """
        return ordering_c_to_py(self.val.minfo.ord)

    def gen(self, slong i):
        """
        Return the `i`th generator of the polynomial ring

            >>> from flint import Ordering
            >>> ctx = fmpz_mpoly_ctx.get_context(3, Ordering.degrevlex, 'z')
            >>> ctx.gen(1)
            z1
        """
        cdef fmpz_mpoly res
        if not 0 <= i < self.val.minfo.nvars:
            raise IndexError("generator index out of range")
        res = create_fmpz_mpoly(self)
        fmpz_mpoly_gen(res.val, i, res.ctx.val)
        return res

    def constant(self, z):
        """
        Create a constant polynomial in this context
        """
        cdef fmpz_mpoly res
        z = any_as_fmpz(z)
        if z is NotImplemented:
            raise TypeError("cannot coerce argument to fmpz")
        res = create_fmpz_mpoly(self)
        fmpz_mpoly_set_fmpz(res.val, (<fmpz>z).val, res.ctx.val)
        return res

    def from_dict(self, d):
        """
        Create a fmpz_mpoly from a dictionary in this context.

        The dictionary's keys are tuples of ints (or anything that implicitly converts
        to fmpz) representing exponents, and corresponding values of fmpz.

            >>> from flint import Ordering
            >>> ctx = fmpz_mpoly_ctx.get_context(2, Ordering.lex, 'x,y')
            >>> ctx.from_dict({(1,0):2, (1,1):3, (0,1):1})
            3*x*y + 2*x + y
        """
        cdef:
            fmpz_vec exp_vec
            slong i, nvars = self.nvars()
            fmpz_mpoly res

        if not isinstance(d, dict):
            raise ValueError("expected a dictionary")

        res = create_fmpz_mpoly(self)

        for i, (k, v) in enumerate(d.items()):
            v = any_as_fmpz(v)
            if v is NotImplemented:
                raise TypeError(f"cannot coerce coefficient '{v}' to fmpz")
            elif len(k) != nvars:
                raise ValueError(f"expected {nvars} exponents, got {len(k)}")

            exp_vec = fmpz_vec(k)

            if v:
                fmpz_mpoly_push_term_fmpz_ffmpz(res.val, (<fmpz>v).val, exp_vec.val, self.val)

        fmpz_mpoly_sort_terms(res.val, self.val)
        return res


cdef class fmpz_mpoly(flint_mpoly):
    """
    The *fmpz_mpoly* type represents sparse multivariate polynomials over
    the integers.
    """

    def __cinit__(self):
        self._init = False

    def __dealloc__(self):
        if self._init:
            fmpz_mpoly_clear(self.val, self.ctx.val)
            self._init = False

    def __init__(self, val=0, ctx=None):
        if typecheck(val, fmpz_mpoly):
            if ctx is None or ctx == (<fmpz_mpoly>val).ctx:
                init_fmpz_mpoly(self, (<fmpz_mpoly>val).ctx)
                fmpz_mpoly_set(self.val, (<fmpz_mpoly>val).val, self.ctx.val)
            else:
                raise IncompatibleContextError(f"{ctx} is not {(<fmpz_mpoly>val).ctx}")
        elif isinstance(val, dict):
            if ctx is None:
                if len(val) == 0:
                    raise ValueError("Need context for zero polynomial")
                k = list(val.keys())[0]
                if not isinstance(k, tuple):
                    raise ValueError("Dict should be keyed with tuples of integers")
                ctx = fmpz_mpoly_ctx.get_context(len(k))
            x = ctx.from_dict(val)
            # XXX: this copy is silly, have a ctx function that assigns an fmpz_mpoly_t
            init_fmpz_mpoly(self, ctx)
            fmpz_mpoly_set(self.val, (<fmpz_mpoly>x).val, self.ctx.val)
        elif isinstance(val, str):
            if ctx is None:
                raise ValueError("Cannot parse a polynomial without context")
            val = bytes(val, 'utf-8')
            init_fmpz_mpoly(self, ctx)
            fmpz_mpoly_set_str_pretty(self.val, val, self.ctx.c_names, self.ctx.val)
            fmpz_mpoly_sort_terms(self.val, self.ctx.val)
        else:
            v = any_as_fmpz(val)
            if v is NotImplemented:
                raise TypeError("cannot create fmpz_mpoly from type %s" % type(val))
            if ctx is None:
                raise ValueError("Need context to convert  fmpz to fmpz_mpoly")
            init_fmpz_mpoly(self, ctx)
            fmpz_mpoly_set_fmpz(self.val, (<fmpz>v).val, self.ctx.val)

    def __bool__(self):
        return not fmpz_mpoly_is_zero(self.val, self.ctx.val)

    def __richcmp__(self, other, int op):
        if not (op == Py_EQ or op == Py_NE):
            return NotImplemented
        elif other is None:
            return op == Py_NE
        elif typecheck(self, fmpz_mpoly) and typecheck(other, fmpz_mpoly):
            if (<fmpz_mpoly>self).ctx is (<fmpz_mpoly>other).ctx:
                return (op == Py_NE) ^ bool(
                    fmpz_mpoly_equal((<fmpz_mpoly>self).val, (<fmpz_mpoly>other).val, (<fmpz_mpoly>self).ctx.val)
                )
            else:
                return op == Py_NE
        else:
            return NotImplemented

    def __len__(self):
        return fmpz_mpoly_length(self.val, self.ctx.val)

    def __getitem__(self, x):
        """
        Return the term at index `x` if `x` is an `int`, or the term with the exponent
        vector `x` if `x` is a tuple of `int`s or `fmpz`s.

            >>> from flint import Ordering
            >>> ctx = fmpz_mpoly_ctx.get_context(2, Ordering.lex, 'x')
            >>> p = ctx.from_dict({(0, 1): 2, (1, 1): 3})
            >>> p[1]
            2*x1
            >>> p[1, 1]
            3

        """
        cdef:
            slong nvars = self.ctx.nvars()

        if isinstance(x, int):
            if not 0 <= x < fmpz_mpoly_length(self.val, self.ctx.val):
                raise IndexError("term index out of range")
            res = create_fmpz_mpoly(self.ctx)
            fmpz_mpoly_get_term((<fmpz_mpoly> res).val, self.val, x, self.ctx.val)
            return res
        elif isinstance(x, tuple):
            if len(x) != nvars:
                raise ValueError("exponent vector provided does not match number of variables")
            res = fmpz()
            exp_vec = fmpz_vec(x, double_indirect=True)
            fmpz_mpoly_get_coeff_fmpz_fmpz((<fmpz>res).val, self.val, exp_vec.double_indirect, self.ctx.val)
            return res
        else:
            raise TypeError("index is not integer or tuple")

    def __setitem__(self, x, y):
        """
        Set the coefficient at index `x` to `y` if `x` is an `int`, or the term with
        the exponent vector `x` if `x` is a tuple of `int`s or `fmpz`s.

            >>> from flint import Ordering
            >>> ctx = fmpz_mpoly_ctx.get_context(2, Ordering.lex, 'x')
            >>> p = ctx.from_dict({(0, 1): 2, (1, 1): 3})
            >>> p[1] = 10
            >>> p[1, 1] = 20
            >>> p
            20*x0*x1 + 10*x1

        """
        cdef:
            slong nvars = self.ctx.nvars()

        coeff = any_as_fmpz(y)
        if coeff is NotImplemented:
            raise TypeError("provided coefficient not coercible to fmpz")

        if isinstance(x, int):
            if not 0 <= x < fmpz_mpoly_length(self.val, self.ctx.val):
                raise IndexError("term index out of range")
            fmpz_mpoly_set_term_coeff_fmpz(self.val, x, (<fmpz>coeff).val, self.ctx.val)
        elif isinstance(x, tuple):
            if len(x) != nvars:
                raise ValueError("exponent vector provided does not match number of variables")
            exp_vec = fmpz_vec(x, double_indirect=True)
            fmpz_mpoly_set_coeff_fmpz_fmpz(self.val, (<fmpz>coeff).val, exp_vec.double_indirect, self.ctx.val)
        else:
            raise TypeError("index is not integer or tuple")

    def __pos__(self):
        return self

    def __neg__(self):
        cdef fmpz_mpoly res
        res = create_fmpz_mpoly(self.ctx)
        fmpz_mpoly_neg(res.val, (<fmpz_mpoly>self).val, res.ctx.val)
        return res

    def __add__(self, other):
        cdef fmpz_mpoly res
        if typecheck(other, fmpz_mpoly):
            if (<fmpz_mpoly>self).ctx is not (<fmpz_mpoly>other).ctx:
                raise IncompatibleContextError(f"{(<fmpz_mpoly>self).ctx} is not {(<fmpz_mpoly>other).ctx}")
            res = create_fmpz_mpoly(self.ctx)
            fmpz_mpoly_add(res.val, (<fmpz_mpoly>self).val, (<fmpz_mpoly>other).val, res.ctx.val)
            return res
        else:
            other = any_as_fmpz(other)
            if other is not NotImplemented:
                res = create_fmpz_mpoly(self.ctx)
                fmpz_mpoly_add_fmpz(res.val, (<fmpz_mpoly>self).val, (<fmpz>other).val, self.ctx.val)
                return res
        return NotImplemented

    def __radd__(self, other):
        cdef fmpz_mpoly res
        other = any_as_fmpz(other)
        if other is not NotImplemented:
            res = create_fmpz_mpoly(self.ctx)
            fmpz_mpoly_add_fmpz(res.val, (<fmpz_mpoly>self).val, (<fmpz>other).val, self.ctx.val)
            return res
        return NotImplemented

    def __iadd__(self, other):
        if typecheck(other, fmpz_mpoly):
            if (<fmpz_mpoly>self).ctx is not (<fmpz_mpoly>other).ctx:
                raise IncompatibleContextError(f"{(<fmpz_mpoly>self).ctx} is not {(<fmpz_mpoly>other).ctx}")
            fmpz_mpoly_add((<fmpz_mpoly>self).val, (<fmpz_mpoly>self).val, (<fmpz_mpoly>other).val, self.ctx.val)
            return self
        else:
            zval = any_as_fmpz(other)
            if zval is not NotImplemented:
                fmpz_mpoly_add_fmpz((<fmpz_mpoly>self).val, (<fmpz_mpoly>self).val, (<fmpz>zval).val, self.ctx.val)
                return self
        return NotImplemented

    def __sub__(self, other):
        cdef fmpz_mpoly res
        if typecheck(other, fmpz_mpoly):
            if (<fmpz_mpoly>self).ctx is not (<fmpz_mpoly>other).ctx:
                raise IncompatibleContextError(f"{(<fmpz_mpoly>self).ctx} is not {(<fmpz_mpoly>other).ctx}")
            res = create_fmpz_mpoly(self.ctx)
            fmpz_mpoly_sub(res.val, (<fmpz_mpoly>self).val, (<fmpz_mpoly>other).val, res.ctx.val)
            return res
        else:
            other = any_as_fmpz(other)
            if other is not NotImplemented:
                res = create_fmpz_mpoly(self.ctx)
                fmpz_mpoly_sub_fmpz(res.val, (<fmpz_mpoly>self).val, (<fmpz>other).val, self.ctx.val)
                return res
        return NotImplemented

    def __rsub__(self, other):
        cdef fmpz_mpoly res
        other = any_as_fmpz(other)
        if other is not NotImplemented:
            res = create_fmpz_mpoly(self.ctx)
            fmpz_mpoly_sub_fmpz(res.val, (<fmpz_mpoly>self).val, (<fmpz>other).val, res.ctx.val)
            return -res
        return NotImplemented

    def __isub__(self, other):
        if typecheck(other, fmpz_mpoly):
            if (<fmpz_mpoly>self).ctx is not (<fmpz_mpoly>other).ctx:
                raise IncompatibleContextError(f"{(<fmpz_mpoly>self).ctx} is not {(<fmpz_mpoly>other).ctx}")
            fmpz_mpoly_sub((<fmpz_mpoly>self).val, (<fmpz_mpoly>self).val, (<fmpz_mpoly>other).val, self.ctx.val)
            return self
        else:
            other = any_as_fmpz(other)
            if other is not NotImplemented:
                fmpz_mpoly_sub_fmpz((<fmpz_mpoly>self).val, (<fmpz_mpoly>self).val, (<fmpz>other).val, self.ctx.val)
                return self
        return NotImplemented

    def __mul__(self, other):
        cdef fmpz_mpoly res
        if typecheck(other, fmpz_mpoly):
            if (<fmpz_mpoly>self).ctx is not (<fmpz_mpoly>other).ctx:
                raise IncompatibleContextError(f"{(<fmpz_mpoly>self).ctx} is not {(<fmpz_mpoly>other).ctx}")
            res = create_fmpz_mpoly(self.ctx)
            fmpz_mpoly_mul(res.val, (<fmpz_mpoly>self).val, (<fmpz_mpoly>other).val, res.ctx.val)
            return res
        else:
            other = any_as_fmpz(other)
            if other is not NotImplemented:
                res = create_fmpz_mpoly(self.ctx)
                fmpz_mpoly_scalar_mul_fmpz(res.val, (<fmpz_mpoly>self).val, (<fmpz>other).val, res.ctx.val)
                return res
        return NotImplemented

    def __rmul__(self, other):
        cdef fmpz_mpoly res
        other = any_as_fmpz(other)
        if other is not NotImplemented:
            res = create_fmpz_mpoly(self.ctx)
            fmpz_mpoly_scalar_mul_fmpz(res.val, (<fmpz_mpoly>self).val, (<fmpz>other).val, res.ctx.val)
            return res
        return NotImplemented

    def __imul__(self, other):
        if typecheck(other, fmpz_mpoly):
            if (<fmpz_mpoly>self).ctx is not (<fmpz_mpoly>other).ctx:
                raise IncompatibleContextError(f"{(<fmpz_mpoly>self).ctx} is not {(<fmpz_mpoly>other).ctx}")
            fmpz_mpoly_mul((<fmpz_mpoly>self).val, (<fmpz_mpoly>self).val, (<fmpz_mpoly>other).val, self.ctx.val)
            return self
        else:
            other = any_as_fmpz(other)
            if other is not NotImplemented:
                fmpz_mpoly_scalar_mul_fmpz(self.val, (<fmpz_mpoly>self).val, (<fmpz>other).val, self.ctx.val)
                return self
        return NotImplemented

    def __pow__(self, other, modulus):
        cdef fmpz_mpoly res
        if modulus is not None:
            raise NotImplementedError
        other = any_as_fmpz(other)
        if other is NotImplemented:
            return other
        if other < 0:
            raise ValueError("cannot raise fmpz_mpoly to negative power")
        res = create_fmpz_mpoly(self.ctx)
        if fmpz_mpoly_pow_fmpz(res.val, (<fmpz_mpoly>self).val, (<fmpz>other).val, res.ctx.val) == 0:
            raise ValueError("unreasonably large polynomial")  # pragma: no cover
        return res

    def __divmod__(self, other):
        cdef fmpz_mpoly res, res2
        if typecheck(other, fmpz_mpoly):
            if not other:
                raise ZeroDivisionError("fmpz_mpoly division by zero")
            elif (<fmpz_mpoly>self).ctx is not (<fmpz_mpoly>other).ctx:
                raise IncompatibleContextError(f"{(<fmpz_mpoly>self).ctx} is not {(<fmpz_mpoly>other).ctx}")
            res = create_fmpz_mpoly(self.ctx)
            res2 = create_fmpz_mpoly(self.ctx)
            fmpz_mpoly_divrem(res.val, res2.val, (<fmpz_mpoly>self).val, (<fmpz_mpoly>other).val, res.ctx.val)
            return (res, res2)
        else:
            other = any_as_fmpz(other)
            if other is not NotImplemented:
                other = fmpz_mpoly(other, self.ctx)
                if not other:
                    raise ZeroDivisionError("fmpz_mpoly division by zero")
                res = create_fmpz_mpoly(self.ctx)
                res2 = create_fmpz_mpoly(self.ctx)
                fmpz_mpoly_divrem(res.val, res2.val, (<fmpz_mpoly>self).val, (<fmpz_mpoly>other).val, res.ctx.val)
                return (res, res2)
        return NotImplemented

    def __rdivmod__(self, other):
        cdef fmpz_mpoly res, res2
        if not self:
            raise ZeroDivisionError("fmpz_mpoly division by zero")
        other = any_as_fmpz(other)
        if other is not NotImplemented:
            other = fmpz_mpoly(other, self.ctx)
            res = create_fmpz_mpoly(self.ctx)
            res2 = create_fmpz_mpoly(self.ctx)
            fmpz_mpoly_divrem(res.val, res2.val, (<fmpz_mpoly>other).val, (<fmpz_mpoly>self).val, res.ctx.val)
            return (res, res2)
        return NotImplemented

    def __floordiv__(self, other):
        cdef fmpz_mpoly res
        if typecheck(other, fmpz_mpoly):
            if not other:
                raise ZeroDivisionError("fmpz_mpoly division by zero")
            elif (<fmpz_mpoly>self).ctx is not (<fmpz_mpoly>other).ctx:
                raise IncompatibleContextError(f"{(<fmpz_mpoly>self).ctx} is not {(<fmpz_mpoly>other).ctx}")
            res = create_fmpz_mpoly(self.ctx)
            fmpz_mpoly_div(res.val, (<fmpz_mpoly>self).val, (<fmpz_mpoly>other).val, res.ctx.val)
            return res
        else:
            other = any_as_fmpz(other)
            if other is not NotImplemented:
                if not other:
                    raise ZeroDivisionError("fmpz_mpoly division by zero")
                other = fmpz_mpoly(other, self.ctx)
                res = create_fmpz_mpoly(self.ctx)
                fmpz_mpoly_div(res.val, (<fmpz_mpoly>self).val, (<fmpz_mpoly>other).val, res.ctx.val)
                return res
        return NotImplemented

    def __rfloordiv__(self, other):
        cdef fmpz_mpoly res
        if not self:
            raise ZeroDivisionError("fmpz_mpoly division by zero")
        other = any_as_fmpz(other)
        if other is not NotImplemented:
            other = fmpz_mpoly(other, self.ctx)
            res = create_fmpz_mpoly(self.ctx)
            fmpz_mpoly_div(res.val, (<fmpz_mpoly>other).val,  self.val, res.ctx.val)
            return res
        return NotImplemented

    def __truediv__(self, other):
        cdef:
            fmpz_mpoly res

        if typecheck(other, fmpz_mpoly):
            res, r = divmod(self, other)
            if not r:
                return res
            else:
                raise DomainError("fmpz_mpoly division is not exact")
        else:
            o = any_as_fmpz(other)
            if o is NotImplemented:
                return NotImplemented
            elif not o:
                raise ZeroDivisionError("fmpz_mpoly division by zero")
            res = create_fmpz_mpoly(self.ctx)
            if fmpz_mpoly_scalar_divides_fmpz(res.val, self.val, (<fmpz>o).val, self.ctx.val):
                return res
            else:
                raise DomainError("fmpz_mpoly division is not exact")

    def __rtruediv__(self, other):
        cdef fmpz_mpoly res
        if not self:
            raise ZeroDivisionError("fmpz_mpoly division by zero")
        o = any_as_fmpz(other)
        if o is NotImplemented:
            return NotImplemented
        res = create_fmpz_mpoly(self.ctx)
        fmpz_mpoly_set_fmpz(res.val, (<fmpz>o).val, self.ctx.val)
        return res / self

    def __mod__(self, other):
        return divmod(self, other)[1]

    def __rmod__(self, other):
        return divmod(other, self)[1]

    def __call__(self, *args, **kwargs):
        cdef:
            fmpz_mpoly res
            fmpz_mpoly res2
            fmpz_mpoly_ctx res_ctx
            fmpz_vec V
            fmpz vres
            fmpz_mpoly_vec C
            slong i, nvars = self.ctx.nvars(), nargs = len(args)

        if not args and not kwargs:
            return self
        elif args and kwargs:
            raise ValueError("only supply positional or keyword arguments")
        elif len(args) < nvars and not kwargs:
            raise ValueError("partial application requires keyword arguments")

        if kwargs:
            # Sort and filter the provided keyword args to only valid ones
            partial_args = tuple((i, kwargs[x]) for i, x in enumerate(self.ctx.names()) if x in kwargs)
            nargs = len(partial_args)

            # If we've been provided with an invalid keyword arg then the length of our filter
            # args will be less than what we've been provided with.
            # If the length is equal to the number of variables then all arguments have been provided
            # otherwise we need to do partial application
            if nargs < len(kwargs):
                raise ValueError("unknown keyword argument provided")
            elif nargs == nvars:
                args = tuple(arg for _, arg in partial_args)
                kwargs = None
        elif nargs > nvars:
            raise ValueError("more arguments provided than variables")

        if kwargs:
            args_fmpz = tuple(any_as_fmpz(v) for _, v in partial_args)
        else:
            args_fmpz = tuple(any_as_fmpz(v) for v in args)

        all_fmpz = NotImplemented not in args_fmpz

        if args and all_fmpz:
            # Normal evaluation
            V = fmpz_vec(args_fmpz, double_indirect=True)
            vres = fmpz.__new__(fmpz)
            if fmpz_mpoly_evaluate_all_fmpz(vres.val, self.val, V.double_indirect, self.ctx.val) == 0:
                raise ValueError("unreasonably large polynomial")  # pragma: no cover
            return vres
        elif kwargs and all_fmpz:
            # Partial application with args in Z. We evaluate the polynomial one variable at a time
            res = create_fmpz_mpoly(self.ctx)
            res2 = create_fmpz_mpoly(self.ctx)

            fmpz_mpoly_set(res2.val, self.val, self.ctx.val)
            for (i, _), arg in zip(partial_args, args_fmpz):
                if fmpz_mpoly_evaluate_one_fmpz(res.val, res2.val, i, (<fmpz>arg).val, self.ctx.val) == 0:
                    raise ValueError("unreasonably large polynomial")  # pragma: no cover
                fmpz_mpoly_set(res2.val, res.val, self.ctx.val)
            return res
        elif args and not all_fmpz:
            # Complete function composition
            if not all(typecheck(arg, fmpz_mpoly) for arg in args):
                raise TypeError("all arguments must be fmpz_mpolys")
            res_ctx = (<fmpz_mpoly> args[0]).ctx
            if not all((<fmpz_mpoly> args[i]).ctx is res_ctx for i in range(1, len(args))):
                raise IncompatibleContextError("all arguments must share the same context")

            C = fmpz_mpoly_vec(args, self.ctx, double_indirect=True)
            res = create_fmpz_mpoly(self.ctx)
            if fmpz_mpoly_compose_fmpz_mpoly(res.val, self.val, C.double_indirect, self.ctx.val, res_ctx.val) == 0:
                raise ValueError("unreasonably large polynomial")  # pragma: no cover
            return res
        else:
            # Partial function composition. We do this by composing with all arguments, however the ones
            # that have not be provided are set to the trivial monomial. This is why we require all the
            # polynomial itself, and all arguments exist in the same context, otherwise we have no way of
            # finding the correct monomial to use.
            if not all(typecheck(arg, fmpz_mpoly) for _, arg in partial_args):
                raise TypeError("all arguments must be fmpz_mpolys for partial composition")
            if not all((<fmpz_mpoly> arg).ctx is self.ctx for _, arg in partial_args):
                raise IncompatibleContextError("all arguments must share the same context")

            polys = [None] * nvars
            for i, poly in partial_args:
                polys[i] = poly

            for i in range(nvars):
                if polys[i] is None:
                    vec = [0] * nvars
                    vec[i] = 1
                    polys[i] = self.ctx.from_dict({tuple(vec): 1})

            C = fmpz_mpoly_vec(polys, self.ctx, double_indirect=True)
            res = create_fmpz_mpoly(self.ctx)
            if fmpz_mpoly_compose_fmpz_mpoly(res.val, self.val, C.double_indirect, self.ctx.val, self.ctx.val) == 0:
                raise ValueError("unreasonably large polynomial")  # pragma: no cover
            return res

    def context(self):
        """
        Return the context object for this polynomials.

            >>> from flint import Ordering
            >>> ctx = fmpz_mpoly_ctx.get_context(2, Ordering.lex, 'x')
            >>> p = ctx.from_dict({(0, 1): 2})
            >>> ctx is p.context()
            True
        """
        return self.ctx

    def is_one(self):
        return fmpz_mpoly_is_one(self.val, self.ctx.val)

    def coefficient(self, slong i):
        """
        Return the coefficient at index `i`.

            >>> from flint import Ordering
            >>> ctx = fmpz_mpoly_ctx.get_context(2, Ordering.lex, 'x')
            >>> p = ctx.from_dict({(0, 1): 2, (1, 1): 3})
            >>> p.coefficient(1)
            2
        """
        cdef fmpz v
        if not 0 <= i < fmpz_mpoly_length(self.val, self.ctx.val):
            raise IndexError("term index out of range")
        else:
            v = fmpz.__new__(fmpz)
            fmpz_mpoly_get_term_coeff_fmpz(v.val, self.val, i, self.ctx.val)
            return v

    def exponent_tuple(self, slong i):
        """
        Return the exponent vector at index `i` as a tuple.

            >>> from flint import Ordering
            >>> ctx = fmpz_mpoly_ctx.get_context(2, Ordering.lex, 'x')
            >>> p = ctx.from_dict({(0, 1): 2, (1, 1): 3})
            >>> p.exponent_tuple(1)
            (0, 1)
        """
        cdef:
            slong nvars = self.ctx.nvars()

        if not 0 <= i < fmpz_mpoly_length(self.val, self.ctx.val):
            raise IndexError("term index out of range")
        res = fmpz_vec(nvars, double_indirect=True)
        fmpz_mpoly_get_term_exp_fmpz(res.double_indirect, self.val, i, self.ctx.val)
        return res.to_tuple()

    def degrees(self):
        """
        Return a dictionary of variable name to degree.

            >>> from flint import Ordering
            >>> ctx = fmpz_mpoly_ctx.get_context(4, Ordering.lex, 'x')
            >>> p = ctx.from_dict({(1, 0, 0, 0): 1, (0, 2, 0, 0): 2, (0, 0, 3, 0): 3})
            >>> p.degrees()
            {'x0': 1, 'x1': 2, 'x2': 3, 'x3': 0}
        """
        cdef:
            slong nvars = self.ctx.nvars()

        res = fmpz_vec(nvars, double_indirect=True)
        fmpz_mpoly_degrees_fmpz(res.double_indirect, self.val, self.ctx.val)
        return dict(zip(self.ctx.names(), res.to_tuple()))

    def total_degree(self):
        """
        Return the total degree.

            >>> from flint import Ordering
            >>> ctx = fmpz_mpoly_ctx.get_context(4, Ordering.lex, 'x')
            >>> p = ctx.from_dict({(1, 0, 0, 0): 1, (0, 2, 0, 0): 2, (0, 0, 3, 0): 3})
            >>> p.total_degree()
            3
        """
        cdef fmpz res = fmpz()
        fmpz_mpoly_total_degree_fmpz((<fmpz> res).val, self.val, self.ctx.val)
        return res

    def repr(self):
        return f"{self.ctx}.from_dict({self.to_dict()})"

    def str(self):
        cdef bytes s = <bytes> fmpz_mpoly_get_str_pretty(self.val, self.ctx.c_names, self.ctx.val)
        res = s.decode().replace("+", " + ").replace("-", " - ")
        if res.startswith(" - "):
            res = "-" + res[3:]
        return res

    def gcd(self, other):
        """
        Return the gcd of self and other.

            >>> from flint import Ordering
            >>> ctx = fmpz_mpoly_ctx.get_context(2, Ordering.lex, 'x')
            >>> f = ctx.from_dict({(1, 1): 4, (0, 0): 1})
            >>> g = ctx.from_dict({(0, 1): 2, (1, 0): 2})
            >>> (f * g).gcd(f)
            4*x0*x1 + 1
        """
        cdef fmpz_mpoly res
        if not typecheck(other, fmpz_mpoly):
            raise TypeError("argument must be a fmpz_mpoly")
        elif (<fmpz_mpoly>self).ctx is not (<fmpz_mpoly>other).ctx:
            raise IncompatibleContextError(f"{(<fmpz_mpoly>self).ctx} is not {(<fmpz_mpoly>other).ctx}")
        res = create_fmpz_mpoly(self.ctx)
        fmpz_mpoly_gcd(res.val, (<fmpz_mpoly>self).val, (<fmpz_mpoly>other).val, res.ctx.val)
        return res

    def sqrt(self, assume_perfect_square: bool = False):
        """
        Return the square root of self.
        If self is known to be a perfect square provide `assume_perfect_square=True` for a more efficient
        result. If self is not a square root the result is not guaranteed to be correct.

            >>> from flint import Ordering
            >>> ctx = fmpz_mpoly_ctx.get_context(2, Ordering.lex, 'x')
            >>> f = ctx.from_dict({(1, 1): 4, (0, 0): 1})
            >>> (f * f).sqrt()
            4*x0*x1 + 1
        """
        cdef fmpz_mpoly res
        res = create_fmpz_mpoly(self.ctx)

        if fmpz_mpoly_sqrt_heap(res.val, self.val, self.ctx.val, not assume_perfect_square):
            return res
        else:
            raise ValueError("polynomial is not a perfect square")

    def factor(self):
        """
        Factors self into irreducible factors, returning a tuple
        (c, factors) where c is the content of the coefficients and
        factors is a list of (poly, exp) pairs.

            >>> from flint import Ordering
            >>> Zm = fmpz_mpoly
            >>> ctx = fmpz_mpoly_ctx.get_context(3, Ordering.lex, 'x,y,z')
            >>> p1 = Zm("2*x + 4", ctx)
            >>> p2 = Zm("3*x*z +  + 3*x + 3*z + 3", ctx)
            >>> (p1 * p2).factor()
            (6, [(z + 1, 1), (x + 2, 1), (x + 1, 1)])
            >>> (p2 * p1 * p2).factor()
            (18, [(z + 1, 2), (x + 2, 1), (x + 1, 2)])
        """
        cdef:
            fmpz_mpoly_factor_t fac
            fmpz c
            fmpz_mpoly u

        fmpz_mpoly_factor_init(fac, self.ctx.val)
        fmpz_mpoly_factor(fac, self.val, self.ctx.val)
        res = [0] * fac.num

        for i in range(fac.num):
            u = create_fmpz_mpoly(self.ctx)
            fmpz_mpoly_set((<fmpz_mpoly>u).val, &fac.poly[i], self.ctx.val)

            c = fmpz.__new__(fmpz)
            fmpz_set((<fmpz>c).val, &fac.exp[i])

            res[i] = (u, c)

        c = fmpz.__new__(fmpz)
        fmpz_set((<fmpz>c).val, fac.constant)
        fmpz_mpoly_factor_clear(fac, self.ctx.val)
        return c, res

    def factor_squarefree(self):
        """
        Factors self into irreducible factors, returning a tuple
        (c, factors) where c is the content of the coefficients and
        factors is a list of (poly, exp) pairs.

            >>> from flint import Ordering
            >>> Zm = fmpz_mpoly
            >>> ctx = fmpz_mpoly_ctx.get_context(3, Ordering.lex, 'x,y,z')
            >>> p1 = Zm("2*x + 4", ctx)
            >>> p2 = Zm("3*x*y + 3*x + 3*y + 3", ctx)
            >>> (p1 * p2).factor_squarefree()
            (6, [(y + 1, 1), (x^2 + 3*x + 2, 1)])
            >>> (p1 * p2 * p1).factor_squarefree()
            (12, [(y + 1, 1), (x + 1, 1), (x + 2, 2)])
        """
        cdef:
            fmpz_mpoly_factor_t fac
            fmpz c
            fmpz_mpoly u

        fmpz_mpoly_factor_init(fac, self.ctx.val)
        fmpz_mpoly_factor_squarefree(fac, self.val, self.ctx.val)
        res = [0] * fac.num

        for i in range(fac.num):
            u = create_fmpz_mpoly(self.ctx)
            fmpz_mpoly_init(u.val, u.ctx.val)
            fmpz_mpoly_set((<fmpz_mpoly>u).val, &fac.poly[i], self.ctx.val)

            c = fmpz.__new__(fmpz)
            fmpz_set((<fmpz>c).val, &fac.exp[i])

            res[i] = (u, c)

        c = fmpz.__new__(fmpz)
        fmpz_set((<fmpz>c).val, fac.constant)
        fmpz_mpoly_factor_clear(fac, self.ctx.val)
        return c, res

    # TODO: Rethink context conversions, particularly the proposed methods in #132
    # def coerce_to_context(self, ctx):
    #     cdef:
    #         fmpz_mpoly res
    #         slong *C
    #         slong i

    #     if not typecheck(ctx, fmpz_mpoly_ctx):
    #         raise ValueError("provided context is not a fmpz_mpoly_ctx")

    #     if self.ctx is ctx:
    #         return self

    #     C = <slong *> libc.stdlib.malloc(self.ctx.val.minfo.nvars * sizeof(slong *))
    #     if C is NULL:
    #         raise MemoryError("malloc returned a null pointer")
    #     res = create_fmpz_mpoly(self.ctx)

    #     vars = {x: i for i, x in enumerate(ctx.py_names)}
    #     for i, var in enumerate(self.ctx.py_names):
    #         C[i] = <slong>vars[var]

    #     fmpz_mpoly_compose_fmpz_mpoly_gen(res.val, self.val, C, self.ctx.val, (<fmpz_mpoly_ctx>ctx).val)

    #     libc.stdlib.free(C)
    #     return res

    def derivative(self, var):
        """
        Return the derivative of this polynomial with respect to the provided variable.
        The argument and either be the variable as a string, or the index of the
        variable in the context.

            >>> from flint import Ordering
            >>> ctx = fmpz_mpoly_ctx.get_context(2, Ordering.lex, 'x')
            >>> p = ctx.from_dict({(0, 3): 2, (2, 1): 3})
            >>> p
            3*x0^2*x1 + 2*x1^3
            >>> p.derivative("x0")
            6*x0*x1
            >>> p.derivative(1)
            3*x0^2 + 6*x1^2

        """
        cdef:
            fmpz_mpoly res
            slong i = 0

        if isinstance(var, str):
            vars = {x: i for i, x in enumerate(self.ctx.names())}
            if var not in vars:
                raise ValueError("variable not in context")
            else:
                i = vars[var]
        elif isinstance(var, int):
            if not 0 <= var < self.ctx.nvars():
                raise IndexError("generator index out of range")
            i = <slong>var
        else:
            raise TypeError("invalid variable type")

        res = create_fmpz_mpoly(self.ctx)

        fmpz_mpoly_derivative(res.val, self.val, i, self.ctx.val)
        return res

    def integral(self, var):
        """
        Return the integral of this polynomial*B with respect to the provided variable
        where B is minimal. The argument and either be the variable as a string, or
        the index of the variable in the context.

            >>> from flint import Ordering
            >>> ctx = fmpz_mpoly_ctx.get_context(2, Ordering.lex, 'x')
            >>> p = ctx.from_dict({(0, 3): 2, (2, 1): 3})
            >>> p
            3*x0^2*x1 + 2*x1^3
            >>> p.integral("x0")
            (1, x0^3*x1 + 2*x0*x1^3)
            >>> p.integral(1)
            (2, 3*x0^2*x1^2 + x1^4)

        """
        cdef:
            fmpz_mpoly res
            fmpz scale
            slong i = 0

        if isinstance(var, str):
            vars = {x: i for i, x in enumerate(self.ctx.names())}
            if var not in vars:
                raise ValueError("variable not in context")
            else:
                i = vars[var]
        elif isinstance(var, int):
            if not 0 <= var < self.ctx.nvars():
                raise IndexError("generator index out of range")
            i = <slong>var
        else:
            raise TypeError("invalid variable type")

        res = create_fmpz_mpoly(self.ctx)

        scale = fmpz()

        fmpz_mpoly_integral(res.val, scale.val, self.val, i, self.ctx.val)
        return scale, res


cdef class fmpz_mpoly_vec:
    """
    A class representing a vector of fmpz_mpolys.
    """

    def __cinit__(self, iterable_or_len, fmpz_mpoly_ctx ctx, bint double_indirect = False):
        if isinstance(iterable_or_len, int):
            length = iterable_or_len
        else:
            length = len(iterable_or_len)

        self.ctx = ctx
        fmpz_mpoly_vec_init(self.val, length, self.ctx.val)

        if double_indirect:
            self.double_indirect = <fmpz_mpoly_struct **> libc.stdlib.malloc(length * sizeof(fmpz_mpoly_struct *))
            if self.double_indirect is NULL:
                raise MemoryError("malloc returned a null pointer")  # pragma: no cover

            for i in range(length):
                self.double_indirect[i] = fmpz_mpoly_vec_entry(self.val, i)
        else:
            self.double_indirect = NULL

    def __init__(self, iterable_or_len, _, double_indirect: bool = False):
        if not isinstance(iterable_or_len, int):
            for i, x in enumerate(iterable_or_len):
                self[i] = x

    def __dealloc__(self):
        libc.stdlib.free(self.double_indirect)
        fmpz_mpoly_vec_clear(self.val, self.ctx.val)

    def __getitem__(self, x):
        if not isinstance(x, int):
            raise TypeError("index is not integer")
        elif not 0 <= x < self.val.length:
            raise IndexError("index out of range")

        cdef fmpz_mpoly z = create_fmpz_mpoly(self.ctx)
        fmpz_mpoly_set(z.val, fmpz_mpoly_vec_entry(self.val, x), self.ctx.val)
        return z

    def __setitem__(self, x, y):
        if not isinstance(x, int):
            raise TypeError("index is not integer")
        elif not 0 <= x < self.val.length:
            raise IndexError("index out of range")
        elif not typecheck(y, fmpz_mpoly):
            raise TypeError("argument is not fmpz_mpoly")
        elif (<fmpz_mpoly>y).ctx is not self.ctx:
            raise IncompatibleContextError(f"{(<fmpz_mpoly>y).ctx} is not {self.ctx}")

        fmpz_mpoly_set(fmpz_mpoly_vec_entry(self.val, x), (<fmpz_mpoly>y).val, self.ctx.val)

    def __len__(self):
        return self.val.length

    def __str__(self):
        s = [None] * self.val.length
        for i in range(self.val.length):
            x = create_fmpz_mpoly(self.ctx)
            fmpz_mpoly_set(x.val, fmpz_mpoly_vec_entry(self.val, i), self.ctx.val)
            s[i] = str(x)
        return f"[{', '.join(s)}]"

    def __repr__(self):
        return f"fmpz_mpoly_vec({self}, ctx={self.ctx})"

    def to_tuple(self):
        return tuple(self[i] for i in range(self.val.length))
