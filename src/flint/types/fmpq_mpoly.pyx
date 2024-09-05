from flint.flint_base.flint_base cimport (
    flint_mpoly,
    flint_mpoly_context,
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

from flint.flintlib.functions.fmpq cimport fmpq_set, fmpq_one
from flint.flintlib.functions.fmpq_mpoly cimport (
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
    fmpq_mpoly_equal_fmpq,
    fmpq_mpoly_equal_fmpz,
    fmpq_mpoly_evaluate_all_fmpq,
    fmpq_mpoly_evaluate_one_fmpq,
    fmpq_mpoly_gcd,
    fmpq_mpoly_gen,
    fmpq_mpoly_get_coeff_fmpq_fmpz,
    fmpq_mpoly_get_str_pretty,
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
    fmpq_mpoly_scalar_mul_fmpq,
    fmpq_mpoly_set,
    fmpq_mpoly_set_coeff_fmpq_fmpz,
    fmpq_mpoly_set_fmpq,
    fmpq_mpoly_set_str_pretty,
    fmpq_mpoly_sort_terms,
    fmpq_mpoly_sqrt,
    fmpq_mpoly_sub,
    fmpq_mpoly_sub_fmpq,
    fmpq_mpoly_total_degree_fmpz,
)
from flint.flintlib.functions.fmpq_mpoly_factor cimport (
    fmpq_mpoly_factor,
    fmpq_mpoly_factor_clear,
    fmpq_mpoly_factor_init,
    fmpq_mpoly_factor_squarefree,
    fmpq_mpoly_factor_t,
)

from flint.flintlib.functions.fmpz cimport fmpz_init_set
from flint.flintlib.functions.fmpz_mpoly cimport fmpz_mpoly_set

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

    def any_as_scalar(self, other):
        if isinstance(other, int):
            return any_as_fmpq(other)
        elif typecheck(other, fmpz):
            return any_as_fmpq(other)
        elif typecheck(other, fmpq):
            res = fmpq.__new__(fmpq)
            fmpq_set((<fmpq>res).val, (<fmpq>other).val)
            return res
        else:
            return NotImplemented

    def scalar_as_mpoly(self, other: fmpq):
        # non-fmpq scalars should first be converted via self.any_as_scalar
        return self.constant(<fmpq>other)

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
            raise ValueError("expected a dictionary")

        res = create_fmpq_mpoly(self)

        for k, v in d.items():
            o = any_as_fmpq(v)
            if o is NotImplemented:
                raise TypeError(f"cannot coerce coefficient '{v}' to fmpq")
            elif len(k) != nvars:
                raise ValueError(f"expected {nvars} exponents, got {len(k)}")

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
                    f"Provided context ('{ctx}') and provided fmpz_mpoly ('{val}') don't share the same number of "
                    "variables"
                )
            init_fmpq_mpoly(self, ctx)
            fmpz_mpoly_set(self.val.zpoly, (<fmpz_mpoly>val).val, (<fmpz_mpoly> val).ctx.val)
            fmpq_one(self.val.content)
            fmpq_mpoly_reduce(self.val, self.ctx.val)
        elif isinstance(val, dict):
            if ctx is None:
                raise ValueError("a context is required to create a fmpq_mpoly from a dict")
            x = ctx.from_dict(val)
            # XXX: this copy is silly, have a ctx function that assigns an fmpz_mpoly_t
            init_fmpq_mpoly(self, ctx)
            fmpq_mpoly_set(self.val, (<fmpq_mpoly>x).val, self.ctx.val)
        elif isinstance(val, str):
            if ctx is None:
                raise ValueError("cannot parse a polynomial without context")
            val = val.encode("ascii")
            init_fmpq_mpoly(self, ctx)
            if fmpq_mpoly_set_str_pretty(self.val, val, self.ctx.c_names, self.ctx.val) == -1:
                raise ValueError("unable to parse fmpq_mpoly from string")
            fmpq_mpoly_sort_terms(self.val, self.ctx.val)
        else:
            v = any_as_fmpq(val)
            if v is NotImplemented:
                raise TypeError("cannot create fmpz_mpoly from type %s" % type(val))
            if ctx is None:
                raise ValueError("need context to convert  fmpz to fmpq_mpoly")
            init_fmpq_mpoly(self, ctx)
            fmpq_mpoly_set_fmpq(self.val, (<fmpq>v).val, self.ctx.val)

    def __bool__(self):
        return not fmpq_mpoly_is_zero(self.val, self.ctx.val)

    def is_zero(self):
        return <bint>fmpq_mpoly_is_zero(self.val, self.ctx.val)

    def is_one(self):
        return <bint>fmpq_mpoly_is_one(self.val, self.ctx.val)

    def __richcmp__(self, other, int op):
        if not (op == Py_EQ or op == Py_NE):
            return NotImplemented
        elif other is None:
            return op == Py_NE
        elif typecheck(self, fmpq_mpoly) and typecheck(other, fmpq_mpoly):
            if (<fmpq_mpoly>self).ctx is (<fmpq_mpoly>other).ctx:
                return (op == Py_NE) ^ <bint>fmpq_mpoly_equal(self.val, (<fmpq_mpoly>other).val, self.ctx.val)
            else:
                return op == Py_NE
        elif typecheck(other, fmpq):
            return (op == Py_NE) ^ <bint>fmpq_mpoly_equal_fmpq(self.val, (<fmpq>other).val, self.ctx.val)
        elif typecheck(other, fmpz):
            return (op == Py_NE) ^ <bint>fmpq_mpoly_equal_fmpz(self.val, (<fmpz>other).val, self.ctx.val)
        elif isinstance(other, int):
            other = any_as_fmpz(other)
            if other is NotImplemented:
                return NotImplemented
            return (op == Py_NE) ^ <bint>fmpq_mpoly_equal_fmpz(self.val, (<fmpz>other).val, self.ctx.val)
        else:
            return NotImplemented

    def __len__(self):
        return fmpq_mpoly_length(self.val, self.ctx.val)

    def __getitem__(self, x):
        """
        Return the coefficient of the term with the exponent vector `x`.
        Always returns a value, missing keys will return `0`.
        Negative exponents are made positive.

            >>> from flint import Ordering
            >>> ctx = fmpq_mpoly_ctx.get_context(2, Ordering.lex, 'x')
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

        res = fmpq()
        exp_vec = fmpz_vec(x, double_indirect=True)
        fmpq_mpoly_get_coeff_fmpq_fmpz((<fmpq>res).val, self.val, exp_vec.double_indirect, self.ctx.val)
        return res

    def __setitem__(self, x, y):
        """
        Set the coefficient of the term with the exponent vector `x` to `y`.
        Will always set a value, missing keys will create a new term.
        Negative exponents are made positive.

            >>> from flint import Ordering
            >>> ctx = fmpq_mpoly_ctx.get_context(2, Ordering.lex, 'x')
            >>> p = ctx.from_dict({(0, 1): 2, (1, 1): 3})
            >>> p[1, 1] = 20
            >>> p
            20*x0*x1 + 2*x1

        """
        cdef:
            slong nvars = self.ctx.nvars()

        coeff = any_as_fmpq(y)
        if coeff is NotImplemented:
            raise TypeError("provided coefficient not coercible to fmpq")
        elif not isinstance(x, tuple):
            raise TypeError("exponent vector index is not a tuple")
        elif len(x) != nvars:
            raise ValueError("exponent vector provided does not match number of variables")

        exp_vec = fmpz_vec(x, double_indirect=True)
        fmpq_mpoly_set_coeff_fmpq_fmpz(self.val, (<fmpq>coeff).val, exp_vec.double_indirect, self.ctx.val)

    def __neg__(self):
        cdef fmpq_mpoly res
        res = create_fmpq_mpoly(self.ctx)
        fmpq_mpoly_neg(res.val, (<fmpq_mpoly>self).val, res.ctx.val)
        return res

    cdef _add_scalar_(self, arg):
        cdef fmpq_mpoly res
        cdef fmpq other = <fmpq>arg
        res = create_fmpq_mpoly(self.ctx)
        fmpq_mpoly_add_fmpq(res.val, self.val, other.val, self.ctx.val)
        return res

    cdef _add_mpoly_(self, arg):
        cdef fmpq_mpoly res, other = <fmpq_mpoly>arg
        res = create_fmpq_mpoly(self.ctx)
        fmpq_mpoly_add(res.val, self.val, other.val, res.ctx.val)
        return res

    cdef _sub_scalar_(self, arg):
        cdef fmpq_mpoly res
        cdef fmpq other = <fmpq>arg
        res = create_fmpq_mpoly(self.ctx)
        fmpq_mpoly_sub_fmpq(res.val, self.val, other.val, self.ctx.val)
        return res

    cdef _sub_mpoly_(self, arg):
        cdef fmpq_mpoly res, other = <fmpq_mpoly>arg
        res = create_fmpq_mpoly(self.ctx)
        fmpq_mpoly_sub(res.val, self.val, other.val, res.ctx.val)
        return res

    cdef _mul_scalar_(self, arg):
        cdef fmpq_mpoly res
        cdef fmpq other = <fmpq>arg
        res = create_fmpq_mpoly(self.ctx)
        fmpq_mpoly_scalar_mul_fmpq(res.val, self.val, other.val, self.ctx.val)
        return res

    cdef _mul_mpoly_(self, arg):
        cdef fmpq_mpoly res, other = <fmpq_mpoly>arg
        res = create_fmpq_mpoly(self.ctx)
        fmpq_mpoly_mul(res.val, self.val, other.val, res.ctx.val)
        return res

    cdef _pow_(self, arg):
        cdef fmpq_mpoly res
        cdef fmpz other = <fmpz>arg
        res = create_fmpq_mpoly(self.ctx)
        if fmpq_mpoly_pow_fmpz(res.val, self.val, other.val, res.ctx.val) == 0:
            raise ValueError("unreasonably large polynomial")  # pragma: no cover
        return res

    cdef _divmod_mpoly_(self, arg):
        cdef fmpq_mpoly quotient, remainder, other = <fmpq_mpoly>arg
        quotient = create_fmpq_mpoly(self.ctx)
        remainder = create_fmpq_mpoly(self.ctx)
        fmpq_mpoly_divrem(quotient.val, remainder.val, self.val, other.val, self.ctx.val)
        return (quotient, remainder)

    cdef _floordiv_mpoly_(self, arg):
        cdef fmpq_mpoly quotient, other = <fmpq_mpoly>arg
        quotient = create_fmpq_mpoly(self.ctx)
        fmpq_mpoly_div(quotient.val, self.val, other.val, self.ctx.val)
        return quotient

    cdef _truediv_mpoly_(self, arg):
        cdef fmpq_mpoly quotient, other = <fmpq_mpoly>arg
        quotient = create_fmpq_mpoly(self.ctx)
        if fmpq_mpoly_divides(quotient.val, self.val, other.val, self.ctx.val):
            return quotient
        else:
            raise DomainError("fmpq_mpoly division is not exact")

    cdef _mod_mpoly_(self, arg):
        cdef fmpq_mpoly quotient, remainder, other = <fmpq_mpoly>arg
        quotient = create_fmpq_mpoly(self.ctx)
        remainder = create_fmpq_mpoly(self.ctx)
        fmpq_mpoly_divrem(quotient.val, remainder.val, self.val, other.val, self.ctx.val)
        return remainder

    cdef _rsub_scalar_(self, arg):
        cdef fmpq_mpoly res
        cdef fmpq other = <fmpq>arg
        res = create_fmpq_mpoly(self.ctx)
        fmpq_mpoly_sub_fmpq(res.val, self.val, other.val, self.ctx.val)
        fmpq_mpoly_neg(res.val, res.val, res.ctx.val)
        return res

    cdef _rsub_mpoly_(self, arg):
        return (<fmpq_mpoly>arg)._sub_mpoly_(self)

    cdef _rdivmod_mpoly_(self, arg):
        return (<fmpq_mpoly>arg)._divmod_mpoly_(self)

    cdef _rfloordiv_mpoly_(self, arg):
        return (<fmpq_mpoly>arg)._floordiv_mpoly_(self)

    cdef _rtruediv_mpoly_(self, arg):
        return (<fmpq_mpoly>arg)._truediv_mpoly_(self)

    cdef _rmod_mpoly_(self, arg):
        return (<fmpq_mpoly>arg)._mod_mpoly_(self)

    cdef _iadd_scalar_(self, arg):
        cdef fmpq other = <fmpq>arg
        fmpq_mpoly_add_fmpq(self.val, self.val, other.val, self.ctx.val)

    cdef _isub_scalar_(self, arg):
        cdef fmpq other = <fmpq>arg
        fmpq_mpoly_sub_fmpq(self.val, self.val, other.val, self.ctx.val)

    cdef _imul_scalar_(self, arg):
        cdef fmpq other = <fmpq>arg
        fmpq_mpoly_scalar_mul_fmpq(self.val, self.val, other.val, self.ctx.val)

    cdef _iadd_mpoly_(self, arg):
        cdef fmpq_mpoly other = <fmpq_mpoly>arg
        fmpq_mpoly_add(self.val, self.val, other.val, self.ctx.val)

    cdef _isub_mpoly_(self, arg):
        cdef fmpq_mpoly other = <fmpq_mpoly>arg
        fmpq_mpoly_sub(self.val, self.val, other.val, self.ctx.val)

    cdef _imul_mpoly_(self, arg):
        cdef fmpq_mpoly other = <fmpq_mpoly>arg
        fmpq_mpoly_mul(self.val, self.val, other.val, self.ctx.val)

    def __call__(self, *args) -> fmpq:
        cdef:
            fmpq_vec V
            fmpq vres
            slong nvars = self.ctx.nvars(), nargs = len(args)

        if nargs != nvars:
            raise ValueError("number of generators does not match number of arguments")

        V = fmpq_vec(args, double_indirect=True)
        vres = fmpq.__new__(fmpq)
        if fmpq_mpoly_evaluate_all_fmpq(vres.val, self.val, V.double_indirect, self.ctx.val) == 0:
            raise ValueError("unreasonably large polynomial")  # pragma: no cover
        return vres

    def monoms(self):
        """
        Return the exponent vectors of each term as a tuple of fmpz.

            >>> from flint import Ordering
            >>> ctx = fmpq_mpoly_ctx.get_context(2, Ordering.lex, 'x')
            >>> f = ctx.from_dict({(0, 0): 1, (1, 0): 2, (0, 1): 3, (1, 1): 4})
            >>> f.monoms()
            [(1, 1), (1, 0), (0, 1), (0, 0)]

        """
        cdef:
            slong i, nvars = self.ctx.nvars()
            fmpz_vec vec = fmpz_vec(nvars, double_indirect=True)

        res = []
        for i in range(len(self)):
            fmpq_mpoly_get_term_exp_fmpz(vec.double_indirect, self.val, i, self.ctx.val)
            res.append(vec.to_tuple())

        return res

    def coeffs(self):
        """
        Return the coefficients of each term as a fmpq.

            >>> from flint import Ordering
            >>> ctx = fmpq_mpoly_ctx.get_context(2, Ordering.lex, 'x')
            >>> f = ctx.from_dict({(0, 0): 1, (1, 0): 2, (0, 1): 3, (1, 1): 4})
            >>> f.coeffs()
            [4, 2, 3, 1]

        """
        cdef:
            fmpq coeff
            slong i

        res = []
        for i in range(len(self)):
            coeff = fmpq.__new__(fmpq)
            fmpq_mpoly_get_term_coeff_fmpq(coeff.val, self.val, i, self.ctx.val)
            res.append(coeff)

        return res

    # def terms(self):
    #     """
    #     Return the terms of this polynomial as a list of fmpq_mpolys.

    #         >>> from flint import Ordering
    #         >>> ctx = fmpq_mpoly_ctx.get_context(2, Ordering.lex, 'x')
    #         >>> f = ctx.from_dict({(0, 0): 1, (1, 0): 2, (0, 1): 3, (1, 1): 4})
    #         >>> f.terms()
    #         [4*x0*x1, 2*x0, 3*x1, 1]

    #     """
    #     cdef:
    #         fmpq_mpoly term
    #         slong i

    #     res = []
    #     for i in range(len(self)):
    #         term = create_fmpq_mpoly(self.ctx)
    #         fmpq_mpoly_get_term(term.val, self.val, i, self.ctx.val)
    #         res.append(term)

    #     return res

    def subs(self, dict_args) -> fmpq_mpoly:
        """
        Partial evaluate this polynomial with select constants. Keys must be generator names or generator indices,
        all values must be fmpq.

            >>> from flint import Ordering
            >>> ctx = fmpq_mpoly_ctx.get_context(2, Ordering.lex, 'x')
            >>> f = ctx.from_dict({(0, 0): 1, (1, 0): 2, (0, 1): 3, (1, 1): 4})
            >>> f.subs({"x1": 0})
            2*x0 + 1

        """
        cdef:
            fmpq_mpoly res
            slong i

        args = tuple((self.ctx.variable_to_index(k), any_as_fmpq(v)) for k, v in dict_args.items())
        for _, v in args:
            if v is NotImplemented:
                raise TypeError(f"cannot coerce argument ('{v}') to fmpq")

        # Partial application with args in Q. We evaluate the polynomial one variable at a time
        res = create_fmpq_mpoly(self.ctx)

        fmpq_mpoly_set(res.val, self.val, self.ctx.val)
        for i, arg in args:
            if fmpq_mpoly_evaluate_one_fmpq(res.val, res.val, i, (<fmpq>arg).val, self.ctx.val) == 0:
                raise ValueError("unreasonably large polynomial")  # pragma: no cover
        return res

    def compose(self, *args, ctx=None) -> fmpq_mpoly:
        """
        Compose this polynomial with other fmpq_mpolys. All arguments must share the same context, it may different
        from this polynomials context.

            >>> from flint import Ordering
            >>> ctx = fmpq_mpoly_ctx.get_context(1, Ordering.lex, 'x')
            >>> ctx1 = fmpq_mpoly_ctx.get_context(2, Ordering.lex, 'y')
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
            fmpq_mpoly res
            fmpq_mpoly_ctx res_ctx
            fmpq_mpoly_vec C
            slong nvars = self.ctx.nvars(), nargs = len(args)

        if nargs < nvars:
            raise ValueError("not enough arguments provided")
        elif nargs > nvars:
            raise ValueError("more arguments provided than variables")
        elif self.ctx.nvars() == 0 and ctx is None:
            raise ValueError("a context must be provided when composing a polynomial with no generators")
        elif not all(typecheck(arg, fmpq_mpoly) for arg in args):
            raise TypeError("all arguments must be fmpq_mpolys")

        if ctx is None:
            res_ctx = (<fmpq_mpoly> args[0]).ctx
        elif typecheck(ctx, fmpq_mpoly_ctx):
            res_ctx = <fmpq_mpoly_ctx>ctx
        else:
            raise TypeError(f"provided context ({ctx}) is not a fmpq_mpoly_ctx")

        if not all((<fmpq_mpoly> arg).ctx is res_ctx for arg in args):
            raise IncompatibleContextError(
                "all arguments must share the " + ("same" if ctx is not None else "provided") + " context"
            )

        C = fmpq_mpoly_vec(args, res_ctx, double_indirect=True)
        res = create_fmpq_mpoly(res_ctx)
        if fmpq_mpoly_compose_fmpq_mpoly(res.val, self.val, C.double_indirect, self.ctx.val, res_ctx.val) == 0:
            raise ValueError("unreasonably large polynomial")  # pragma: no cover
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

    def monomial(self, slong i):
        """
        Return the exponent vector at index `i` as a tuple.

            >>> from flint import Ordering
            >>> ctx = fmpq_mpoly_ctx.get_context(2, Ordering.lex, 'x')
            >>> p = ctx.from_dict({(0, 1): 2, (1, 1): 3})
            >>> p.monomial(1)
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
            (1, 2, 3, 0)
        """
        cdef:
            slong nvars = self.ctx.nvars()

        res = fmpz_vec(nvars, double_indirect=True)
        fmpq_mpoly_degrees_fmpz(res.double_indirect, self.val, self.ctx.val)
        return res.to_tuple()

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

    def leading_coefficient(self):
        """
        Leading coefficient in the monomial ordering.

            >>> from flint import Ordering
            >>> ctx = fmpq_mpoly_ctx(2, Ordering.lex, ['x', 'y'])
            >>> x, y = ctx.gens()
            >>> p = 2*x*y + 3*x + 4*y**2 + 5
            >>> p
            2*x*y + 3*x + 4*y^2 + 5
            >>> p.leading_coefficient()
            2

        """
        if fmpq_mpoly_is_zero(self.val, self.ctx.val):
            return fmpq(0)
        else:
            return self.coefficient(0)

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
            raise DomainError("polynomial is not a perfect square")

    def factor(self):
        """
        Factors self into irreducible factors, returning a tuple
        (c, factors) where c is the content of the coefficients and
        factors is a list of (poly, exp) pairs.

            >>> from flint import Ordering
            >>> Zm = fmpq_mpoly
            >>> ctx = fmpq_mpoly_ctx.get_context(3, Ordering.lex, 'x,y,z')
            >>> p1 = Zm("2*x + 4", ctx)
            >>> p2 = Zm("3*x*z + 3*x + 3*z + 3", ctx)
            >>> (p1 * p2).factor()
            (6, [(z + 1, 1), (x + 2, 1), (x + 1, 1)])
            >>> (p2 * p1 * p2).factor()
            (18, [(z + 1, 2), (x + 2, 1), (x + 1, 2)])
        """
        cdef:
            fmpq_mpoly_factor_t fac
            fmpq_mpoly u

        fmpq_mpoly_factor_init(fac, self.ctx.val)
        if not fmpq_mpoly_factor(fac, self.val, self.ctx.val):
            raise RuntimeError("factorisation failed")
        res = [0] * fac.num

        for i in range(fac.num):
            u = create_fmpq_mpoly(self.ctx)
            fmpq_mpoly_set((<fmpq_mpoly>u).val, &fac.poly[i], self.ctx.val)

            c = fmpz.__new__(fmpz)
            fmpz_init_set((<fmpz>c).val, &fac.exp[i])

            res[i] = (u, int(c))

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

            res[i] = (u, int(c))

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

    def __richcmp__(self, other, int op):
        if not (op == Py_EQ or op == Py_NE):
            return NotImplemented
        elif typecheck(other, fmpq_mpoly_vec):
            if (<fmpq_mpoly_vec>self).ctx is (<fmpq_mpoly_vec>other).ctx and len(self) == len(other):
                return (op == Py_NE) ^ all(x == y for x, y in zip(self, other))
            else:
                return op == Py_NE
        else:
            return NotImplemented

    def to_tuple(self):
        return tuple(self[i] for i in range(self.length))
