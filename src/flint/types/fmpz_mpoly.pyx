from flint.flint_base.flint_base cimport (
    flint_mpoly,
    flint_mpoly_context,
    ordering_py_to_c,
    ordering_c_to_py,
)

from flint.utils.typecheck cimport typecheck
from flint.utils.flint_exceptions import DomainError, IncompatibleContextError

from flint.types.fmpz cimport any_as_fmpz, fmpz
from flint.types.fmpz_vec cimport fmpz_vec

from flint.flintlib.types.flint cimport slong
from flint.flintlib.types.fmpz cimport fmpz_mpoly_vec_entry
from flint.flintlib.functions.fmpz cimport fmpz_set

from flint.flintlib.functions.fmpz_mpoly cimport (
    fmpz_mpoly_add,
    fmpz_mpoly_add_fmpz,
    fmpz_mpoly_buchberger_naive,
    fmpz_mpoly_buchberger_naive_with_limits,
    fmpz_mpoly_clear,
    fmpz_mpoly_compose_fmpz_mpoly,
    fmpz_mpoly_compose_fmpz_mpoly_gen,
    fmpz_mpoly_ctx_init,
    fmpz_mpoly_deflate,
    fmpz_mpoly_deflation,
    fmpz_mpoly_degrees_fmpz,
    fmpz_mpoly_derivative,
    fmpz_mpoly_discriminant,
    fmpz_mpoly_div,
    fmpz_mpoly_divides,
    fmpz_mpoly_divrem,
    fmpz_mpoly_equal,
    fmpz_mpoly_equal_fmpz,
    fmpz_mpoly_evaluate_all_fmpz,
    fmpz_mpoly_evaluate_one_fmpz,
    fmpz_mpoly_gcd,
    fmpz_mpoly_gen,
    fmpz_mpoly_get_coeff_fmpz_fmpz,
    fmpz_mpoly_get_str_pretty,
    fmpz_mpoly_get_term_coeff_fmpz,
    fmpz_mpoly_get_term_exp_fmpz,
    fmpz_mpoly_inflate,
    fmpz_mpoly_integral,
    fmpz_mpoly_is_one,
    fmpz_mpoly_is_zero,
    fmpz_mpoly_length,
    fmpz_mpoly_mul,
    fmpz_mpoly_neg,
    fmpz_mpoly_pow_fmpz,
    fmpz_mpoly_push_term_fmpz_ffmpz,
    fmpz_mpoly_push_term_ui_ffmpz,
    fmpz_mpoly_reduction_primitive_part,
    fmpz_mpoly_resultant,
    fmpz_mpoly_scalar_divexact_fmpz,
    fmpz_mpoly_scalar_divides_fmpz,
    fmpz_mpoly_scalar_mul_fmpz,
    fmpz_mpoly_set,
    fmpz_mpoly_set_coeff_fmpz_fmpz,
    fmpz_mpoly_set_fmpz,
    fmpz_mpoly_set_str_pretty,
    fmpz_mpoly_sort_terms,
    fmpz_mpoly_spoly,
    fmpz_mpoly_sqrt_heap,
    fmpz_mpoly_sub,
    fmpz_mpoly_sub_fmpz,
    fmpz_mpoly_term_content,
    fmpz_mpoly_total_degree_fmpz,
    fmpz_mpoly_vec_autoreduction,
    fmpz_mpoly_vec_autoreduction_groebner,
    fmpz_mpoly_vec_clear,
    fmpz_mpoly_vec_init,
    fmpz_mpoly_vec_is_autoreduced,
    fmpz_mpoly_vec_is_groebner,
)
from flint.flintlib.functions.fmpz_mpoly_factor cimport (
    fmpz_mpoly_factor,
    fmpz_mpoly_factor_clear,
    fmpz_mpoly_factor_init,
    fmpz_mpoly_factor_squarefree,
    fmpz_mpoly_factor_t,
)
from flint.flintlib.functions.fmpz_vec cimport _fmpz_vec_content

from cpython.object cimport Py_EQ, Py_NE
cimport libc.stdlib

cdef dict _fmpz_mpoly_ctx_cache = {}


cdef class fmpz_mpoly_ctx(flint_mpoly_context):
    """
    A class for storing the polynomial context

    :param names:  A tuple containing the names of the variables of the ring.
    :param ordering:  The term order for the ring.

    Do not construct one of these directly, use ``fmpz_mpoly_ctx.get``.
    """

    _ctx_cache = _fmpz_mpoly_ctx_cache

    @classmethod
    def _new_(cls, names, ordering):
        cdef fmpz_mpoly_ctx self = cls.__new__(cls)
        super()._new_(self, names)
        fmpz_mpoly_ctx_init(self.val, len(names), ordering_py_to_c(ordering))
        return self

    def _any_as_scalar(self, other):
        if isinstance(other, int):
            return any_as_fmpz(other)
        elif typecheck(other, fmpz):
            res = fmpz.__new__(fmpz)
            fmpz_set((<fmpz>res).val, (<fmpz>other).val)
            return res
        else:
            return NotImplemented

    def _scalar_as_mpoly(self, other: fmpz):
        # non-fmpz scalars should first be converted via self._any_as_scalar
        return self.constant(<fmpz>other)

    def nvars(self):
        """
        Return the number of variables in the context

            >>> ctx = fmpz_mpoly_ctx.get(('x', 4), 'lex')
            >>> ctx.nvars()
            4
        """
        return self.val.minfo.nvars

    def ordering(self):
        """
        Return the term order of the context object.

            >>> ctx = fmpz_mpoly_ctx.get(('w', 4), 'deglex')
            >>> ctx.ordering()
            <Ordering.deglex: 'deglex'>
        """
        return ordering_c_to_py(self.val.minfo.ord)

    def gen(self, slong i):
        """
        Return the ``i`` th generator of the polynomial ring

            >>> ctx = fmpz_mpoly_ctx.get(('z', 3), 'degrevlex')
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

            >>> ctx = fmpz_mpoly_ctx.get(('x', 'y'), 'lex')
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
            o = any_as_fmpz(v)
            if o is NotImplemented:
                raise TypeError(f"cannot coerce coefficient '{v}' to fmpz")
            elif len(k) != nvars:
                raise ValueError(f"expected {nvars} exponents, got {len(k)}")

            exp_vec = fmpz_vec(k)

            if o:
                fmpz_mpoly_push_term_fmpz_ffmpz(res.val, (<fmpz>o).val, exp_vec.val, self.val)

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
                raise ValueError("a context is required to create a fmpz_mpoly from a dict")
            x = ctx.from_dict(val)
            # XXX: this copy is silly, have a ctx function that assigns an fmpz_mpoly_t
            init_fmpz_mpoly(self, ctx)
            fmpz_mpoly_set(self.val, (<fmpz_mpoly>x).val, self.ctx.val)
        elif isinstance(val, str):
            if ctx is None:
                raise ValueError("cannot parse a polynomial without context")
            val = bytes(val, 'utf-8')
            init_fmpz_mpoly(self, ctx)
            if fmpz_mpoly_set_str_pretty(self.val, val, self.ctx.c_names, self.ctx.val) == -1:
                raise ValueError("unable to parse fmpz_mpoly from string")
            fmpz_mpoly_sort_terms(self.val, self.ctx.val)
        else:
            v = any_as_fmpz(val)
            if v is NotImplemented:
                raise TypeError("cannot create fmpz_mpoly from type %s" % type(val))
            if ctx is None:
                raise ValueError("need context to convert  fmpz to fmpz_mpoly")
            init_fmpz_mpoly(self, ctx)
            fmpz_mpoly_set_fmpz(self.val, (<fmpz>v).val, self.ctx.val)

    def __bool__(self):
        return not fmpz_mpoly_is_zero(self.val, self.ctx.val)

    def is_zero(self):
        return <bint>fmpz_mpoly_is_zero(self.val, self.ctx.val)

    def is_one(self):
        return <bint>fmpz_mpoly_is_one(self.val, self.ctx.val)

    def is_constant(self):
        """
        Returns True if this is a constant polynomial.

        >>> ctx = fmpz_mpoly_ctx.get(['x', 'y'])
        >>> x, y = ctx.gens()
        >>> p = x**2 + y
        >>> p.is_constant()
        False
        >>> (0*p + 1).is_constant()
        True
        """
        return self.total_degree() <= 0

    def __richcmp__(self, other, int op):
        if not (op == Py_EQ or op == Py_NE):
            return NotImplemented
        elif other is None:
            return op == Py_NE
        elif typecheck(other, fmpz_mpoly):
            if (<fmpz_mpoly>self).ctx is (<fmpz_mpoly>other).ctx:
                return (op == Py_NE) ^ <bint>fmpz_mpoly_equal(self.val, (<fmpz_mpoly>other).val, self.ctx.val)
            else:
                return op == Py_NE
        elif typecheck(other, fmpz):
            return (op == Py_NE) ^ <bint>fmpz_mpoly_equal_fmpz(self.val, (<fmpz>other).val, self.ctx.val)
        elif isinstance(other, int):
            other = any_as_fmpz(other)
            if other is NotImplemented:
                return NotImplemented
            return (op == Py_NE) ^ <bint>fmpz_mpoly_equal_fmpz(self.val, (<fmpz>other).val, self.ctx.val)
        else:
            return NotImplemented

    def __len__(self):
        return fmpz_mpoly_length(self.val, self.ctx.val)

    def __getitem__(self, x):
        """
        Return the coefficient of the term with the exponent vector ``x``.
        Always returns a value, missing keys will return ``0``.
        Negative exponents are made positive.

            >>> ctx = fmpz_mpoly_ctx.get(('x', 2), 'lex')
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
        fmpz_mpoly_get_coeff_fmpz_fmpz((<fmpz>res).val, self.val, exp_vec.double_indirect, self.ctx.val)
        return res

    def __setitem__(self, x, y):
        """
        Set the coefficient of the term with the exponent vector ``x`` to ``y``.
        Will always set a value, missing keys will create a new term.
        Negative exponents are made positive.

            >>> ctx = fmpz_mpoly_ctx.get(('x', 2), 'lex')
            >>> p = ctx.from_dict({(0, 1): 2, (1, 1): 3})
            >>> p[1, 1] = 20
            >>> p
            20*x0*x1 + 2*x1

        """
        cdef:
            slong nvars = self.ctx.nvars()

        coeff = any_as_fmpz(y)
        if coeff is NotImplemented:
            raise TypeError("provided coefficient not coercible to fmpz")
        elif not isinstance(x, tuple):
            raise TypeError("exponent vector index is not a tuple")
        elif len(x) != nvars:
            raise ValueError("exponent vector provided does not match number of variables")

        exp_vec = fmpz_vec(x, double_indirect=True)
        fmpz_mpoly_set_coeff_fmpz_fmpz(self.val, (<fmpz>coeff).val, exp_vec.double_indirect, self.ctx.val)

    def __neg__(self):
        cdef fmpz_mpoly res
        res = create_fmpz_mpoly(self.ctx)
        fmpz_mpoly_neg(res.val, (<fmpz_mpoly>self).val, res.ctx.val)
        return res

    cdef _add_scalar_(self, arg):
        cdef fmpz_mpoly res
        cdef fmpz other = <fmpz>arg
        res = create_fmpz_mpoly(self.ctx)
        fmpz_mpoly_add_fmpz(res.val, self.val, other.val, self.ctx.val)
        return res

    cdef _sub_scalar_(self, arg):
        cdef fmpz_mpoly res
        cdef fmpz other = <fmpz>arg
        res = create_fmpz_mpoly(self.ctx)
        fmpz_mpoly_sub_fmpz(res.val, self.val, other.val, self.ctx.val)
        return res

    cdef _mul_scalar_(self, arg):
        cdef fmpz_mpoly res
        cdef fmpz other = <fmpz>arg
        res = create_fmpz_mpoly(self.ctx)
        fmpz_mpoly_scalar_mul_fmpz(res.val, self.val, other.val, self.ctx.val)
        return res

    cdef _pow_(self, arg):
        cdef fmpz_mpoly res
        cdef fmpz other = <fmpz>arg
        res = create_fmpz_mpoly(self.ctx)
        if fmpz_mpoly_pow_fmpz(res.val, self.val, other.val, res.ctx.val) == 0:
            raise ValueError("unreasonably large polynomial")  # pragma: no cover
        return res

    cdef _add_mpoly_(self, arg):
        cdef fmpz_mpoly res, other = <fmpz_mpoly>arg
        res = create_fmpz_mpoly(self.ctx)
        fmpz_mpoly_add(res.val, self.val, other.val, res.ctx.val)
        return res

    cdef _sub_mpoly_(self, arg):
        cdef fmpz_mpoly res, other = <fmpz_mpoly>arg
        res = create_fmpz_mpoly(self.ctx)
        fmpz_mpoly_sub(res.val, self.val, other.val, res.ctx.val)
        return res

    cdef _mul_mpoly_(self, arg):
        cdef fmpz_mpoly res, other = <fmpz_mpoly>arg
        res = create_fmpz_mpoly(self.ctx)
        fmpz_mpoly_mul(res.val, self.val, other.val, res.ctx.val)
        return res

    cdef _divmod_mpoly_(self, arg):
        cdef fmpz_mpoly quotient, remainder, other = <fmpz_mpoly>arg
        quotient = create_fmpz_mpoly(self.ctx)
        remainder = create_fmpz_mpoly(self.ctx)
        fmpz_mpoly_divrem(quotient.val, remainder.val, self.val, other.val, self.ctx.val)
        return (quotient, remainder)

    cdef _floordiv_mpoly_(self, arg):
        cdef fmpz_mpoly quotient, other = <fmpz_mpoly>arg
        quotient = create_fmpz_mpoly(self.ctx)
        fmpz_mpoly_div(quotient.val, self.val, other.val, self.ctx.val)
        return quotient

    cdef _truediv_scalar_(self, arg):
        cdef fmpz_mpoly quotient,
        cdef fmpz other = <fmpz>arg
        quotient = create_fmpz_mpoly(self.ctx)
        if fmpz_mpoly_scalar_divides_fmpz(quotient.val, self.val, other.val, self.ctx.val):
            return quotient
        else:
            raise DomainError("fmpz_mpoly division is not exact")

    cdef _divexact_scalar_(self, arg):
        cdef fmpz_mpoly quotient,
        cdef fmpz other = <fmpz>arg
        quotient = create_fmpz_mpoly(self.ctx)
        fmpz_mpoly_scalar_divexact_fmpz(quotient.val, self.val, other.val, self.ctx.val)
        return quotient

    cdef _truediv_mpoly_(self, arg):
        cdef fmpz_mpoly quotient, other = <fmpz_mpoly>arg
        quotient = create_fmpz_mpoly(self.ctx)
        if fmpz_mpoly_divides(quotient.val, self.val, other.val, self.ctx.val):
            return quotient
        else:
            raise DomainError("fmpz_mpoly division is not exact")

    cdef _mod_mpoly_(self, arg):
        cdef fmpz_mpoly quotient, remainder, other = <fmpz_mpoly>arg
        quotient = create_fmpz_mpoly(self.ctx)
        remainder = create_fmpz_mpoly(self.ctx)
        fmpz_mpoly_divrem(quotient.val, remainder.val, self.val, other.val, self.ctx.val)
        return remainder

    cdef _rsub_scalar_(self, arg):
        cdef fmpz_mpoly res
        cdef fmpz other = <fmpz>arg
        res = create_fmpz_mpoly(self.ctx)
        fmpz_mpoly_sub_fmpz(res.val, self.val, other.val, self.ctx.val)
        fmpz_mpoly_neg(res.val, res.val, res.ctx.val)
        return res

    cdef _rsub_mpoly_(self, arg):
        return (<fmpz_mpoly>arg)._sub_mpoly_(self)

    cdef _rdivmod_mpoly_(self, arg):
        return (<fmpz_mpoly>arg)._divmod_mpoly_(self)

    cdef _rfloordiv_mpoly_(self, arg):
        return (<fmpz_mpoly>arg)._floordiv_mpoly_(self)

    cdef _rtruediv_mpoly_(self, arg):
        return (<fmpz_mpoly>arg)._truediv_mpoly_(self)

    cdef _rmod_mpoly_(self, arg):
        return (<fmpz_mpoly>arg)._mod_mpoly_(self)

    cdef _iadd_scalar_(self, arg):
        cdef fmpz other = <fmpz>arg
        fmpz_mpoly_add_fmpz(self.val, self.val, other.val, self.ctx.val)

    cdef _isub_scalar_(self, arg):
        cdef fmpz other = <fmpz>arg
        fmpz_mpoly_sub_fmpz(self.val, self.val, other.val, self.ctx.val)

    cdef _imul_scalar_(self, arg):
        cdef fmpz other = <fmpz>arg
        fmpz_mpoly_scalar_mul_fmpz(self.val, self.val, other.val, self.ctx.val)

    cdef _iadd_mpoly_(self, arg):
        cdef fmpz_mpoly other = <fmpz_mpoly>arg
        fmpz_mpoly_add(self.val, self.val, other.val, self.ctx.val)

    cdef _isub_mpoly_(self, arg):
        cdef fmpz_mpoly other = <fmpz_mpoly>arg
        fmpz_mpoly_sub(self.val, self.val, other.val, self.ctx.val)

    cdef _imul_mpoly_(self, arg):
        cdef fmpz_mpoly other = <fmpz_mpoly>arg
        fmpz_mpoly_mul(self.val, self.val, other.val, self.ctx.val)

    def __call__(self, *args) -> fmpz:
        cdef:
            fmpz_vec V
            fmpz vres
            slong nvars = self.ctx.nvars(), nargs = len(args)

        if nargs != nvars:
            raise ValueError("number of generators does not match number of arguments")

        V = fmpz_vec(args, double_indirect=True)
        vres = fmpz.__new__(fmpz)
        if fmpz_mpoly_evaluate_all_fmpz(vres.val, self.val, V.double_indirect, self.ctx.val) == 0:
            raise ValueError("unreasonably large polynomial")  # pragma: no cover
        return vres

    def monoms(self):
        """
        Return the exponent vectors of each term as a tuple of fmpz.

            >>> ctx = fmpz_mpoly_ctx.get(('x', 2), 'lex')
            >>> f = ctx.from_dict({(0, 0): 1, (1, 0): 2, (0, 1): 3, (1, 1): 4})
            >>> f.monoms()
            [(1, 1), (1, 0), (0, 1), (0, 0)]

        """
        cdef:
            slong i, nvars = self.ctx.nvars()
            fmpz_vec vec = fmpz_vec(nvars, double_indirect=True)

        res = []
        for i in range(len(self)):
            fmpz_mpoly_get_term_exp_fmpz(vec.double_indirect, self.val, i, self.ctx.val)
            res.append(tuple(vec))

        return res

    def coeffs(self):
        """
        Return the coefficients of each term as a fmpz

            >>> ctx = fmpz_mpoly_ctx.get(('x', 2), 'lex')
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
            fmpz_mpoly_get_term_coeff_fmpz(coeff.val, self.val, i, self.ctx.val)
            res.append(coeff)

        return res

    def subs(self, dict_args) -> fmpz_mpoly:
        """
        Partial evaluate this polynomial with select constants. Keys must be generator names or generator indices,
        all values must be fmpz.

            >>> ctx = fmpz_mpoly_ctx.get(('x', 2), 'lex')
            >>> f = ctx.from_dict({(0, 0): 1, (1, 0): 2, (0, 1): 3, (1, 1): 4})
            >>> f.subs({"x1": 0})
            2*x0 + 1

        """
        cdef:
            fmpz_mpoly res
            slong i

        args = tuple((self.ctx.variable_to_index(k), any_as_fmpz(v)) for k, v in dict_args.items())
        for _, v in args:
            if v is NotImplemented:
                raise TypeError(f"cannot coerce argument ('{v}') to fmpz")

        # Partial application with args in Z. We evaluate the polynomial one variable at a time
        res = create_fmpz_mpoly(self.ctx)

        fmpz_mpoly_set(res.val, self.val, self.ctx.val)
        for i, arg in args:
            if fmpz_mpoly_evaluate_one_fmpz(res.val, res.val, i, (<fmpz>arg).val, self.ctx.val) == 0:
                raise ValueError("unreasonably large polynomial")  # pragma: no cover
        return res

    def compose(self, *args, ctx=None) -> fmpz_mpoly:
        """
        Compose this polynomial with other fmpz_mpolys. All arguments must share the same context, it may different
        from this polynomials context.

            >>> ctx = fmpz_mpoly_ctx.get(('x',), 'lex')
            >>> ctx1 = fmpz_mpoly_ctx.get(('y', 2), 'lex')
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
            fmpz_mpoly res
            fmpz_mpoly_ctx res_ctx
            fmpz_mpoly_vec C
            slong nvars = self.ctx.nvars(), nargs = len(args)

        if nargs < nvars:
            raise ValueError("not enough arguments provided")
        elif nargs > nvars:
            raise ValueError("more arguments provided than variables")
        elif self.ctx.nvars() == 0 and ctx is None:
            raise ValueError("a context must be provided when composing a polynomial with no generators")
        elif not all(typecheck(arg, fmpz_mpoly) for arg in args):
            raise TypeError("all arguments must be fmpz_mpolys")

        if ctx is None:
            res_ctx = (<fmpz_mpoly> args[0]).ctx
        elif typecheck(ctx, fmpz_mpoly_ctx):
            res_ctx = <fmpz_mpoly_ctx>ctx
        else:
            raise TypeError(f"provided context ({ctx}) is not a fmpz_mpoly_ctx")

        if not all((<fmpz_mpoly> arg).ctx is res_ctx for arg in args):
            raise IncompatibleContextError(
                "all arguments must share the " + ("same" if ctx is not None else "provided") + " context"
            )

        C = fmpz_mpoly_vec(args, res_ctx, double_indirect=True)
        res = create_fmpz_mpoly(res_ctx)
        if fmpz_mpoly_compose_fmpz_mpoly(res.val, self.val, C.double_indirect, self.ctx.val, res_ctx.val) == 0:
            raise ValueError("unreasonably large polynomial")  # pragma: no cover
        return res

    def context(self):
        """
        Return the context object for this polynomials.

            >>> ctx = fmpz_mpoly_ctx.get(('x', 2), 'lex')
            >>> p = ctx.from_dict({(0, 1): 2})
            >>> ctx is p.context()
            True
        """
        return self.ctx

    def coefficient(self, slong i):
        """
        Return the coefficient at index ``i``.

            >>> ctx = fmpz_mpoly_ctx.get(('x', 2), 'lex')
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

    def monomial(self, slong i):
        """
        Return the exponent vector at index ``i`` as a tuple.

            >>> ctx = fmpz_mpoly_ctx.get(('x', 2), 'lex')
            >>> p = ctx.from_dict({(0, 1): 2, (1, 1): 3})
            >>> p.monomial(1)
            (0, 1)
        """
        cdef:
            slong nvars = self.ctx.nvars()

        if not 0 <= i < fmpz_mpoly_length(self.val, self.ctx.val):
            raise IndexError("term index out of range")
        res = fmpz_vec(nvars, double_indirect=True)
        fmpz_mpoly_get_term_exp_fmpz(res.double_indirect, self.val, i, self.ctx.val)
        return tuple(res)

    def degrees(self):
        """
        Return a tuple of degrees.

            >>> ctx = fmpz_mpoly_ctx.get(('x', 4), 'lex')
            >>> p = sum(x**i for i, x in enumerate(ctx.gens()))
            >>> p
            x1 + x2^2 + x3^3 + 1
            >>> p.degrees()
            (0, 1, 2, 3)
        """
        cdef:
            slong nvars = self.ctx.nvars()

        res = fmpz_vec(nvars, double_indirect=True)
        fmpz_mpoly_degrees_fmpz(res.double_indirect, self.val, self.ctx.val)
        return tuple(res)

    def total_degree(self):
        """
        Return the total degree.

            >>> ctx = fmpz_mpoly_ctx.get(('x', 4), 'lex')
            >>> p = ctx.from_dict({(1, 0, 0, 0): 1, (0, 2, 0, 0): 2, (0, 0, 3, 0): 3})
            >>> p.total_degree()
            3
        """
        cdef fmpz res = fmpz()
        fmpz_mpoly_total_degree_fmpz((<fmpz> res).val, self.val, self.ctx.val)
        return res

    def leading_coefficient(self):
        """
        Leading coefficient in the monomial ordering.

            >>> ctx = fmpz_mpoly_ctx.get(('x', 'y'), 'lex')
            >>> x, y = ctx.gens()
            >>> p = 2*x*y + 3*x + 4*y**2 + 5
            >>> p
            2*x*y + 3*x + 4*y^2 + 5
            >>> p.leading_coefficient()
            2

        """
        if fmpz_mpoly_is_zero(self.val, self.ctx.val):
            return fmpz(0)
        else:
            return self.coefficient(0)

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

            >>> ctx = fmpz_mpoly_ctx.get(('x', 2), 'lex')
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

    def content(self):
        """
        Return the GCD of the coefficients of ``self``.

            >>> ctx = fmpz_mpoly_ctx.get(('x', 2), 'lex')
            >>> x0, x1 = ctx.gens()
            >>> f = 3 * x0**2 * x1 + 6 * x0 * x1
            >>> f.content()
            3
        """
        cdef fmpz res = fmpz()
        _fmpz_vec_content(res.val, self.val.coeffs, self.val.length)
        return res

    def term_content(self):
        """
        Return the GCD of the terms of ``self``. If ``self`` is zero, then the result will
        be zero, otherwise it will be a monomial with positive coefficient.

            >>> ctx = fmpz_mpoly_ctx.get(('x', 2), 'lex')
            >>> x0, x1 = ctx.gens()
            >>> f = 3 * x0**2 * x1 + 6 * x0 * x1
            >>> f.term_content()
            3*x0*x1
        """
        cdef fmpz_mpoly res = create_fmpz_mpoly(self.ctx)
        fmpz_mpoly_term_content(res.val, self.val, self.ctx.val)
        return res

    def resultant(self, other, var):
        """
        Return the resultant of ``self`` and ``other`` with respect to variable ``var``.

            >>> ctx = fmpz_mpoly_ctx.get(('x', 2), 'lex')
            >>> x0, x1 = ctx.gens()
            >>> f = x0**2 * x1 + x0 * x1
            >>> g = x0 + x1
            >>> f.resultant(g, 'x1')
            x0^3 + x0^2
        """
        cdef:
            fmpz_mpoly res
            slong i

        if not typecheck(other, fmpz_mpoly):
            raise TypeError("argument must be a fmpz_mpoly")
        elif (<fmpz_mpoly>self).ctx is not (<fmpz_mpoly>other).ctx:
            raise IncompatibleContextError(f"{(<fmpz_mpoly>self).ctx} is not {(<fmpz_mpoly>other).ctx}")

        i = self.ctx.variable_to_index(var)
        res = create_fmpz_mpoly(self.ctx)
        if not fmpz_mpoly_resultant(res.val, self.val, (<fmpz_mpoly>other).val, i, self.ctx.val):
            raise RuntimeError(f"failed to compute resultant with respect to {var}")
        return res

    def discriminant(self, var):
        """
        Return the discriminant of ``self`` with respect to variable ``var``.

            >>> ctx = fmpz_mpoly_ctx.get(('x', 2), 'lex')
            >>> x0, x1 = ctx.gens()
            >>> f = (x0 + x1)**2 + 1
            >>> f.discriminant('x1')
            -4

        """
        cdef:
            fmpz_mpoly res
            slong i

        i = self.ctx.variable_to_index(var)
        res = create_fmpz_mpoly(self.ctx)
        if not fmpz_mpoly_discriminant(res.val, self.val, i, self.ctx.val):
            raise RuntimeError(f"failed to compute discriminant with respect to {var}")
        return res

    def primitive(self):
        """
        Return the content and primitive of ``self``.

            >>> ctx = fmpz_mpoly_ctx.get(('x', 2), 'lex')
            >>> x, y = ctx.gens()
            >>> f = 4*x + 2*x*y
            >>> f.primitive()
            (2, x0*x1 + 2*x0)
        """
        cdef fmpz content = self.content()
        return content, self._divexact_scalar_(content)

    def sqrt(self, assume_perfect_square: bool = False):
        """
        Return the square root of self.
        If self is known to be a perfect square provide ``assume_perfect_square=True`` for a more efficient
        result. If self is not a square root the result is not guaranteed to be correct.

            >>> ctx = fmpz_mpoly_ctx.get(('x', 2), 'lex')
            >>> f = ctx.from_dict({(1, 1): 4, (0, 0): 1})
            >>> (f * f).sqrt()
            4*x0*x1 + 1
        """
        cdef fmpz_mpoly res
        res = create_fmpz_mpoly(self.ctx)

        if fmpz_mpoly_sqrt_heap(res.val, self.val, self.ctx.val, not assume_perfect_square):
            return res
        else:
            raise DomainError("polynomial is not a perfect square")

    def factor(self):
        """
        Factors self into irreducible factors, returning a tuple
        (c, factors) where c is the content of the coefficients and
        factors is a list of (poly, exp) pairs.

            >>> Zm = fmpz_mpoly
            >>> ctx = fmpz_mpoly_ctx.get(('x', 'y', 'z'), 'lex')
            >>> p1 = Zm("2*x + 4", ctx)
            >>> p2 = Zm("3*x*z + 3*x + 3*z + 3", ctx)
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
        if not fmpz_mpoly_factor(fac, self.val, self.ctx.val):
            raise RuntimeError("factorisation failed")
        res = [0] * fac.num

        for i in range(fac.num):
            u = create_fmpz_mpoly(self.ctx)
            fmpz_mpoly_set((<fmpz_mpoly>u).val, &fac.poly[i], self.ctx.val)

            c = fmpz.__new__(fmpz)
            fmpz_set((<fmpz>c).val, &fac.exp[i])

            res[i] = (u, int(c))

        c = fmpz.__new__(fmpz)
        fmpz_set((<fmpz>c).val, fac.constant)
        fmpz_mpoly_factor_clear(fac, self.ctx.val)
        return c, res

    def factor_squarefree(self):
        """
        Factors self into square-free factors, returning a tuple
        (c, factors) where c is the content of the coefficients and
        factors is a list of (poly, exp) pairs.

            >>> Zm = fmpz_mpoly
            >>> ctx = fmpz_mpoly_ctx.get(('x', 'y', 'z'), 'lex')
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

            res[i] = (u, int(c))

        c = fmpz.__new__(fmpz)
        fmpz_set((<fmpz>c).val, fac.constant)
        fmpz_mpoly_factor_clear(fac, self.ctx.val)
        return c, res

    cdef _compose_gens_(self, ctx, slong *mapping):
        cdef fmpz_mpoly res = create_fmpz_mpoly(ctx)
        fmpz_mpoly_compose_fmpz_mpoly_gen(
            res.val,
            self.val,
            mapping,
            self.ctx.val,
            (<fmpz_mpoly_ctx>ctx).val
        )

        return res

    def derivative(self, var):
        """
        Return the derivative of this polynomial with respect to the provided variable.
        The argument can either be the variable as a string, or the index of the
        variable in the context.

            >>> ctx = fmpz_mpoly_ctx.get(('x', 2), 'lex')
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
            slong i = self.ctx.variable_to_index(var)

        res = create_fmpz_mpoly(self.ctx)

        fmpz_mpoly_derivative(res.val, self.val, i, self.ctx.val)
        return res

    def integral(self, var):
        """
        Return the integral of this polynomial*B with respect to the provided variable where B is minimal. That is,
        given p this method returns (B, P) such that P' = B*p or P/B is the formal integral of p.

        The argument can either be the variable as a string, or the index of the variable in the context.

            >>> ctx = fmpz_mpoly_ctx.get(('x', 2), 'lex')
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
            slong i = self.ctx.variable_to_index(var)

        res = create_fmpz_mpoly(self.ctx)

        scale = fmpz()

        fmpz_mpoly_integral(res.val, scale.val, self.val, i, self.ctx.val)
        return scale, res

    def spoly(self, g):
        """
        Compute the S-polynomial of ``self`` and ``g``, scaled to an integer polynomial
        by computing the LCM of the leading coefficients.

            >>> ctx = fmpz_mpoly_ctx.get(('x', 'y'), 'lex')
            >>> f = ctx.from_dict({(2, 0): 1, (0, 1): -1})
            >>> g = ctx.from_dict({(3, 0): 1, (1, 0): -1})
            >>> f.spoly(g)
            -x*y + x
        """
        cdef fmpz_mpoly res = create_fmpz_mpoly(self.ctx)

        if not typecheck(g, fmpz_mpoly):
            raise TypeError(f"expected fmpz_mpoly, got {type(g)}")

        self.ctx.compatible_context_check((<fmpz_mpoly>g).ctx)
        fmpz_mpoly_spoly(res.val, self.val, (<fmpz_mpoly>g).val, self.ctx.val)
        return res

    def reduction_primitive_part(self, vec):
        """
        Compute the the primitive part of the reduction (remainder of multivariate
        quasi-division with remainder) with respect to the polynomials ``vec``.

            >>> ctx = fmpz_mpoly_ctx.get(('x', 'y'), 'lex')
            >>> x, y = ctx.gens()
            >>> f = 2 * x**3 -x**2 * y + y**3 + 3 * y
            >>> g1 = x**2 + y**2 + 1
            >>> g2 = x * y - 2
            >>> vec = fmpz_mpoly_vec([g1, g2], ctx)
            >>> vec
            fmpz_mpoly_vec([x^2 + y^2 + 1, x*y - 2], ctx=fmpz_mpoly_ctx(2, '<Ordering.lex: 'lex'>', ('x', 'y')))
            >>> f.reduction_primitive_part(vec)
            x - y^3
        """
        cdef fmpz_mpoly res = create_fmpz_mpoly(self.ctx)
        if not typecheck(vec, fmpz_mpoly_vec):
            raise TypeError(f"expected fmpz_mpoly, got {type(vec)}")

        self.ctx.compatible_context_check((<fmpz_mpoly_vec>vec).ctx)
        fmpz_mpoly_reduction_primitive_part(res.val, self.val, (<fmpz_mpoly_vec>vec).val, self.ctx.val)
        return res

    def inflate(self, N: list[int]) -> fmpz_mpoly:
        """
        Compute the inflation of ``self`` for a provided ``N``, that is return ``q``
        such that ``q(X) = p(X^N)``.

            >>> ctx = fmpz_mpoly_ctx.get(('x', 'y'), 'lex')
            >>> x, y = ctx.gens()
            >>> f = x + y + 1
            >>> f.inflate([2, 3])
            x^2 + y^3 + 1
        """

        cdef nvars = self.ctx.nvars()

        if nvars != len(N):
            raise ValueError(f"expected list of length {nvars}, got {len(N)}")
        elif any(n < 0 for n in N):
            raise ValueError("all inflate strides must be non-negative")

        cdef:
            fmpz_vec shift = fmpz_vec(nvars)
            fmpz_vec stride = fmpz_vec(N)
            fmpz_mpoly res = create_fmpz_mpoly(self.ctx)

        fmpz_mpoly_inflate(res.val, self.val, shift.val, stride.val, self.ctx.val)
        return res

    def deflate(self, N: list[int]) -> fmpz_mpoly:
        """
        Compute the deflation of ``self`` for a provided ``N``, that is return ``q``
        such that ``q(X) = p(X^(1/N))``.

            >>> ctx = fmpz_mpoly_ctx.get(('x', 'y'), 'lex')
            >>> x, y = ctx.gens()
            >>> f = x**3 * y + x * y**4 + x * y
            >>> f.deflate([2, 3])
            x + y + 1
        """
        cdef slong nvars = self.ctx.nvars()

        if nvars != len(N):
            raise ValueError(f"expected list of length {nvars}, got {len(N)}")

        cdef:
            fmpz_vec shift = fmpz_vec(nvars)
            fmpz_vec stride = fmpz_vec(N)
            fmpz_mpoly res = create_fmpz_mpoly(self.ctx)

        fmpz_mpoly_deflate(res.val, self.val, shift.val, stride.val, self.ctx.val)
        return res

    def deflation(self) -> tuple[fmpz_mpoly, list[int]]:
        """
        Compute the deflation of ``self``, that is ``p(X^(1/N))`` for maximal N.

            >>> ctx = fmpz_mpoly_ctx.get(('x', 'y'), 'lex')
            >>> x, y = ctx.gens()
            >>> f = x**2 * y**2 + x * y**2
            >>> q, N = f.deflation()
            >>> q, N
            (x^2*y + x*y, [1, 2])
            >>> q.inflate(N) == f
            True
        """
        cdef:
            slong nvars = self.ctx.nvars()
            fmpz_vec shift = fmpz_vec(nvars)
            fmpz_vec stride = fmpz_vec(nvars)
            fmpz_mpoly res = create_fmpz_mpoly(self.ctx)

        fmpz_mpoly_deflation(shift.val, stride.val, self.val, self.ctx.val)

        for i in range(nvars):
            stride[i] = shift[i].gcd(stride[i])
            shift[i] = 0

        fmpz_mpoly_deflate(res.val, self.val, shift.val, stride.val, self.ctx.val)

        return res, list(stride)

    def deflation_monom(self) -> tuple[fmpz_mpoly, list[int], fmpz_mpoly]:
        """
        Compute the exponent vector ``N`` and monomial ``m`` such that ``p(X^(1/N))
        = m * q(X^N)`` for maximal N. The returned monomial allows the undo-ing of the
        deflation.

            >>> ctx = fmpz_mpoly_ctx.get(('x', 'y'), 'lex')
            >>> x, y = ctx.gens()
            >>> f = x**3 * y + x * y**4 + x * y
            >>> fd, N, m = f.deflation_monom()
            >>> fd, N, m
            (x + y + 1, [2, 3], x*y)
            >>> m * fd.inflate(N)
            x^3*y + x*y^4 + x*y
        """
        cdef:
            slong nvars = self.ctx.nvars()
            fmpz_mpoly res = create_fmpz_mpoly(self.ctx)
            fmpz_mpoly monom = create_fmpz_mpoly(self.ctx)
            fmpz_vec shift = fmpz_vec(nvars)
            fmpz_vec stride = fmpz_vec(nvars)

        fmpz_mpoly_deflation(shift.val, stride.val, self.val, self.ctx.val)
        fmpz_mpoly_push_term_ui_ffmpz(monom.val, 1, fmpz_vec(shift).val, self.ctx.val)
        fmpz_mpoly_deflate(res.val, self.val, shift.val, stride.val, self.ctx.val)

        return res, list(stride), monom

    def deflation_index(self) -> tuple[list[int], list[int]]:
        """
        Compute the exponent vectors ``N`` and ``I`` such that ``p(X^(1/N)) = X^I *
        q(X^N)`` for maximal N. Importantly the deflation itself is not computed
        here. The returned exponent vector ``I`` is the shift that was applied to the
        exponents. It is the exponent vector of the monomial returned by
        ``deflation_monom``.

            >>> ctx = fmpz_mpoly_ctx.get(('x', 'y'), 'lex')
            >>> x, y = ctx.gens()
            >>> f = x**3 * y + x * y**4 + x * y
            >>> N, I = f.deflation_index()
            >>> N, I
            ([2, 3], [1, 1])
            >>> f_deflated = f.deflate(N)
            >>> f_deflated
            x + y + 1
            >>> m = ctx.term(exp_vec=I)
            >>> m
            x*y
            >>> m * f_deflated.inflate(N)
            x^3*y + x*y^4 + x*y
        """
        cdef:
            slong nvars = self.ctx.nvars()
            fmpz_vec shift = fmpz_vec(nvars)
            fmpz_vec stride = fmpz_vec(nvars)

        fmpz_mpoly_deflation(shift.val, stride.val, self.val, self.ctx.val)
        return list(stride), list(shift)


cdef class fmpz_mpoly_vec:
    """
    A class representing a vector of fmpz_mpolys.
    """

    def __cinit__(
            self, iterable_or_len, fmpz_mpoly_ctx ctx, bint double_indirect = False
    ):
        if isinstance(iterable_or_len, int):
            length = iterable_or_len
        else:
            length = len(iterable_or_len)

        self.ctx = ctx
        fmpz_mpoly_vec_init(self.val, length, self.ctx.val)

        if double_indirect:
            self.double_indirect = <fmpz_mpoly_struct **> libc.stdlib.malloc(
                length * sizeof(fmpz_mpoly_struct *)
            )
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

        fmpz_mpoly_set(
            fmpz_mpoly_vec_entry(self.val, x), (<fmpz_mpoly>y).val, self.ctx.val
        )

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

    def __richcmp__(self, other, int op):
        if not (op == Py_EQ or op == Py_NE):
            return NotImplemented
        elif typecheck(other, fmpz_mpoly_vec):
            if (
                    (<fmpz_mpoly_vec>self).ctx is (<fmpz_mpoly_vec>other).ctx
                    and len(self) == len(other)
            ):
                return (op == Py_NE) ^ all(x == y for x, y in zip(self, other))
            else:
                return op == Py_NE
        else:
            return NotImplemented

    def __iter__(self):
        cdef fmpz_mpoly z
        for i in range(self.val.length):
            z = create_fmpz_mpoly(self.ctx)
            fmpz_mpoly_set(z.val, fmpz_mpoly_vec_entry(self.val, i), self.ctx.val)
            yield z

    def is_groebner(self, other=None) -> bool:
        """
        Check if self is a Gröbner basis. If ``other`` is not None then check if self
        is a Gröbner basis for ``other``.

            >>> ctx = fmpz_mpoly_ctx.get(('x', 'y'), 'lex')
            >>> x, y = ctx.gens()
            >>> f = x**2 - y
            >>> g = x**3 - x
            >>> k = x*y - x
            >>> h = y**2 - y
            >>> vec = fmpz_mpoly_vec([f, k, h], ctx)
            >>> vec.is_groebner()
            True
            >>> vec.is_groebner(fmpz_mpoly_vec([f, g], ctx))
            True
            >>> vec.is_groebner(fmpz_mpoly_vec([f, x**3], ctx))
            False
        """
        if other is None:
            return <bint>fmpz_mpoly_vec_is_groebner(self.val, NULL, self.ctx.val)
        elif typecheck(other, fmpz_mpoly_vec):
            self.ctx.compatible_context_check((<fmpz_mpoly_vec>other).ctx)
            return <bint>fmpz_mpoly_vec_is_groebner(
                self.val, (<fmpz_mpoly_vec>other).val, self.ctx.val
            )
        else:
            raise TypeError(
                f"expected either None or a fmpz_mpoly_vec, got {type(other)}"
            )

    def is_autoreduced(self) -> bool:
        """
        Check if self is auto-reduced (or inter-reduced).

            >>> ctx = fmpz_mpoly_ctx.get(('x', 'y'), 'lex')
            >>> x, y = ctx.gens()
            >>> f2 = 3*x**2 - y
            >>> k = x*y - x
            >>> k2 = 3*x*y - 3*x
            >>> h = y**2 - y
            >>> vec = fmpz_mpoly_vec([f2, k, h], ctx)
            >>> vec.is_autoreduced()
            True
            >>> vec = fmpz_mpoly_vec([f2, k2, h], ctx)
            >>> vec.is_autoreduced()
            False

        """
        return <bint>fmpz_mpoly_vec_is_autoreduced(self.val, self.ctx.val)

    def autoreduction(self, groebner=False) -> fmpz_mpoly_vec:
        """
        Compute the autoreduction of ``self``. If ``groebner`` is True and ``self`` is a
        Gröbner basis, compute the reduced reduced Gröbner basis of ``self``, throws an
        ``RuntimeError`` otherwise.

            >>> ctx = fmpz_mpoly_ctx.get(('x', 'y'), 'lex')
            >>> x, y = ctx.gens()
            >>> f2 = 3*x**2 - y
            >>> k2 = 3*x*y - 3*x
            >>> h = y**2 - y
            >>> vec = fmpz_mpoly_vec([f2, k2, h], ctx)
            >>> vec.is_autoreduced()
            False
            >>> vec2 = vec.autoreduction()
            >>> vec2.is_autoreduced()
            True
            >>> vec2
            fmpz_mpoly_vec([3*x^2 - y, x*y - x, y^2 - y], ctx=fmpz_mpoly_ctx(2, '<Ordering.lex: 'lex'>', ('x', 'y')))
        """

        cdef fmpz_mpoly_vec h = fmpz_mpoly_vec(0, self.ctx)

        if groebner:
            if not self.is_groebner():
                raise RuntimeError(
                    "reduced Gröbner basis construction requires that ``self`` is a "
                    "Gröbner basis."
                )
            fmpz_mpoly_vec_autoreduction_groebner(h.val, self.val, self.ctx.val)
        else:
            fmpz_mpoly_vec_autoreduction(h.val, self.val, self.ctx.val)

        return h

    def buchberger_naive(self, limits=None):
        """
        Compute the Gröbner basis of ``self`` using a naive implementation of
        Buchberger’s algorithm.

        Provide ``limits`` in the form of a tuple of ``(ideal_len_limit, poly_len_limit,
        poly_bits_limit)`` to halt execution if the length of the ideal basis set exceeds
        ``ideal_len_limit``, the length of any polynomial exceeds ``poly_len_limit``, or the
        size of the coefficients of any polynomial exceeds ``poly_bits_limit``.

        If limits is provided return a tuple of ``(result, success)``. If ``success`` is
        False then ``result`` is a valid basis for ``self``, but it may not be a Gröbner
        basis.

        NOTE: This function is exposed only for convenience, it is a naive
        implementation and does not compute a reduced basis. To construct a reduced
        basis use ``autoreduce``.

            >>> ctx = fmpz_mpoly_ctx.get(('x', 'y'), 'lex')
            >>> x, y = ctx.gens()
            >>> f = x**2 - y
            >>> g = x**3*y - x
            >>> vec = fmpz_mpoly_vec([f, g], ctx)
            >>> vec.is_groebner()
            False
            >>> vec2 = vec.buchberger_naive()
            >>> vec2
            fmpz_mpoly_vec([x^2 - y, x^3*y - x, x*y^2 - x, y^3 - y], ctx=fmpz_mpoly_ctx(2, '<Ordering.lex: 'lex'>', ('x', 'y')))
            >>> vec.buchberger_naive(limits=(2, 2, 512))
            (fmpz_mpoly_vec([x^2 - y, x^3*y - x], ctx=fmpz_mpoly_ctx(2, '<Ordering.lex: 'lex'>', ('x', 'y'))), False)
            >>> vec2.is_autoreduced()
            False
            >>> vec2.autoreduction()
            fmpz_mpoly_vec([x^2 - y, y^3 - y, x*y^2 - x], ctx=fmpz_mpoly_ctx(2, '<Ordering.lex: 'lex'>', ('x', 'y')))
        """

        cdef:
            fmpz_mpoly_vec g = fmpz_mpoly_vec(0, self.ctx)
            slong ideal_len_limit, poly_len_limit, poly_bits_limit

        if limits is not None:
            ideal_len_limit, poly_len_limit, poly_bits_limit = limits
            if fmpz_mpoly_buchberger_naive_with_limits(
                    g.val,
                    self.val,
                    ideal_len_limit,
                    poly_len_limit,
                    poly_bits_limit,
                    self.ctx.val
            ):
                return g, True
            else:
                return g, False
        else:
            fmpz_mpoly_buchberger_naive(g.val, self.val, self.ctx.val)
            return g
