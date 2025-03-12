from flint.flint_base.flint_base cimport (
    flint_mpoly,
    flint_mod_mpoly_context,
    ordering_py_to_c,
    ordering_c_to_py,
)

from flint.flint_base.flint_base import FLINT_RELEASE

from flint.utils.typecheck cimport typecheck
from flint.utils.flint_exceptions import DomainError, IncompatibleContextError

from flint.types.fmpz cimport any_as_fmpz, fmpz
from flint.types.fmpz_vec cimport fmpz_vec
from flint.types.fmpz_mod cimport fmpz_mod
from flint.types.nmod cimport nmod

from flint.flintlib.functions.fmpz cimport fmpz_set
from flint.flintlib.functions.fmpz_mod_mpoly cimport (
    fmpz_mod_mpoly_add,
    fmpz_mod_mpoly_add_fmpz,
    fmpz_mod_mpoly_clear,
    fmpz_mod_mpoly_compose_fmpz_mod_mpoly,
    fmpz_mod_mpoly_ctx_get_modulus,
    fmpz_mod_mpoly_ctx_init,
    fmpz_mod_mpoly_deflate,
    fmpz_mod_mpoly_deflation,
    fmpz_mod_mpoly_degrees_fmpz,
    fmpz_mod_mpoly_derivative,
    fmpz_mod_mpoly_discriminant,
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
    fmpz_mod_mpoly_get_term_coeff_fmpz,
    fmpz_mod_mpoly_get_term_exp_fmpz,
    fmpz_mod_mpoly_inflate,
    fmpz_mod_mpoly_is_one,
    fmpz_mod_mpoly_is_zero,
    fmpz_mod_mpoly_length,
    fmpz_mod_mpoly_mul,
    fmpz_mod_mpoly_neg,
    fmpz_mod_mpoly_pow_fmpz,
    fmpz_mod_mpoly_push_term_fmpz_ffmpz,
    fmpz_mod_mpoly_push_term_ui_ffmpz,
    fmpz_mod_mpoly_resultant,
    fmpz_mod_mpoly_scalar_mul_fmpz,
    fmpz_mod_mpoly_set,
    fmpz_mod_mpoly_set_coeff_fmpz_fmpz,
    fmpz_mod_mpoly_set_fmpz,
    fmpz_mod_mpoly_set_str_pretty,
    fmpz_mod_mpoly_sort_terms,
    fmpz_mod_mpoly_sqrt,
    fmpz_mod_mpoly_sub,
    fmpz_mod_mpoly_sub_fmpz,
    fmpz_mod_mpoly_term_content,
    fmpz_mod_mpoly_total_degree_fmpz,
)
from flint.flintlib.functions.fmpz_mod_mpoly_factor cimport (
    fmpz_mod_mpoly_factor,
    fmpz_mod_mpoly_factor_clear,
    fmpz_mod_mpoly_factor_init,
    fmpz_mod_mpoly_factor_squarefree,
    fmpz_mod_mpoly_factor_t,
)
from flint.flintlib.functions.compat cimport compat_fmpz_mod_mpoly_compose_fmpz_mod_mpoly_gen

from flint.types.fmpz_mpoly cimport fmpz_mpoly_ctx, fmpz_mpoly


from cpython.object cimport Py_EQ, Py_NE
cimport libc.stdlib

cdef dict _fmpz_mod_mpoly_ctx_cache = {}


cdef class fmpz_mod_mpoly_ctx(flint_mod_mpoly_context):
    """
    A class for storing the polynomial context

    :param names:  A tuple containing the names of the variables of the ring.
    :param ordering:  The term order for the ring.
    :param modulus:  The modulus for the ring.

    Do not construct one of these directly, use ``fmpz_mod_mpoly_ctx.get``.
    """

    _ctx_cache = _fmpz_mod_mpoly_ctx_cache

    @classmethod
    def _new_(cls, names, ordering, modulus):
        cdef fmpz_mod_mpoly_ctx self = cls.__new__(cls)
        cdef fmpz m
        if not typecheck(modulus, fmpz):
            m = any_as_fmpz(modulus)
            if m is NotImplemented:
                raise TypeError(f"modulus ({modulus}) is not coercible to fmpz")
        else:
            m = modulus

        super()._new_(self, names, m.is_prime())
        fmpz_mod_mpoly_ctx_init(self.val, len(names), ordering_py_to_c(ordering), m.val)
        return self

    def _any_as_scalar(self, other):
        if isinstance(other, int):
            return any_as_fmpz(other)
        elif typecheck(other, nmod):
            return any_as_fmpz((<nmod>other).val)
        elif typecheck(other, fmpz):
            res = fmpz.__new__(fmpz)
            fmpz_set((<fmpz>res).val, (<fmpz>other).val)
            return res
        elif typecheck(other, fmpz_mod):
            res = fmpz.__new__(fmpz)
            fmpz_set((<fmpz>res).val, (<fmpz_mod>other).val)
            return res
        else:
            return NotImplemented

    def _scalar_as_mpoly(self, other: fmpz):
        # non-fmpz scalars should first be converted via self._any_as_scalar
        return self.constant(<fmpz>other)

    def nvars(self):
        """
        Return the number of variables in the context

            >>> ctx = fmpz_mod_mpoly_ctx.get(('x', 4), 11, 'lex')
            >>> ctx.nvars()
            4
        """
        return self.val.minfo.nvars

    def ordering(self):
        """
        Return the term order of the context object.

            >>> ctx = fmpz_mod_mpoly_ctx.get(('w', 4), 11, 'deglex')
            >>> ctx.ordering()
            <Ordering.deglex: 'deglex'>
        """
        return ordering_c_to_py(self.val.minfo.ord)

    def modulus(self):
        """
        Return the modulus of the context object.

            >>> ctx = fmpz_mod_mpoly_ctx.get(('w', 4), 2, 'deglex')
            >>> ctx.modulus()
            2

        """
        cdef fmpz m = fmpz.__new__(fmpz)
        fmpz_mod_mpoly_ctx_get_modulus(m.val, self.val)
        return m

    def gen(self, slong i):
        """
        Return the ``i`` th generator of the polynomial ring

            >>> ctx = fmpz_mod_mpoly_ctx.get(('z', 3), 11, 'degrevlex')
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

            >>> ctx = fmpz_mod_mpoly_ctx.get(('x', 'y'), 11, 'lex')
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
            coeff_scalar = self._any_as_scalar(coeff)
            if coeff_scalar is NotImplemented:
                raise TypeError(f"cannot coerce {repr(coeff)} to nmod_mpoly coefficient")

            fmpz_mod_mpoly_push_term_fmpz_ffmpz(
                res.val,
                (<fmpz>coeff_scalar).val,
                exp_vec.val,
                self.val
            )

        fmpz_mod_mpoly_sort_terms(res.val, self.val)
        return res


cdef class fmpz_mod_mpoly(flint_mpoly):
    """
    The *fmpz_mod_mpoly* type represents sparse multivariate polynomials over
    the integers modulo ``n``, for large ``n``.
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

    def _division_check(self, other):
        super()._division_check(other)
        if not self.ctx.is_prime():
            raise DomainError("division with non-prime modulus is not supported")

    def __bool__(self):
        return not fmpz_mod_mpoly_is_zero(self.val, self.ctx.val)

    def is_zero(self):
        return <bint>fmpz_mod_mpoly_is_zero(self.val, self.ctx.val)

    def is_one(self):
        return <bint>fmpz_mod_mpoly_is_one(self.val, self.ctx.val)

    def is_constant(self):
        """
        Returns True if this is a constant polynomial.

        >>> R = fmpz_mod_mpoly_ctx.get(['x', 'y'], modulus=11)
        >>> x, y = R.gens()
        >>> x.is_constant()
        False
        >>> (0*x + 1).is_constant()
        True
        """
        return self.total_degree() <= 0

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
        Return the coefficient of the term with the exponent vector ``x``.
        Always returns a value, missing keys will return ``0``.
        Negative exponents are made positive.

            >>> ctx = fmpz_mod_mpoly_ctx.get(('x', 2), 11, 'lex')
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
        Set the coefficient of the term with the exponent vector ``x`` to ``y``.
        Will always set a value, missing keys will create a new term.
        Negative exponents are made positive.

            >>> ctx = fmpz_mod_mpoly_ctx.get(('x', 2), 11, 'lex')
            >>> p = ctx.from_dict({(0, 1): 2, (1, 1): 3})
            >>> p[1, 1] = 20
            >>> p
            9*x0*x1 + 2*x1

        """
        cdef:
            slong nvars = self.ctx.nvars()

        if not isinstance(x, tuple):
            raise TypeError("exponent vector index is not a tuple")
        elif len(x) != nvars:
            raise ValueError("exponent vector provided does not match number of variables")
        exp_vec = fmpz_vec(x, double_indirect=True)

        coeff = self.ctx._any_as_scalar(y)
        if coeff is NotImplemented:
            raise TypeError("provided coefficient not coercible to fmpz")
        fmpz_mod_mpoly_set_coeff_fmpz_fmpz(self.val, (<fmpz>coeff).val, exp_vec.double_indirect, self.ctx.val)

    def __neg__(self):
        cdef fmpz_mod_mpoly res
        res = create_fmpz_mod_mpoly(self.ctx)
        fmpz_mod_mpoly_neg(res.val, (<fmpz_mod_mpoly>self).val, res.ctx.val)
        return res

    cdef _add_scalar_(self, arg):
        cdef fmpz_mod_mpoly res
        cdef fmpz other = <fmpz>arg
        res = create_fmpz_mod_mpoly(self.ctx)
        fmpz_mod_mpoly_add_fmpz(res.val, self.val, other.val, self.ctx.val)
        return res

    cdef _sub_scalar_(self, arg):
        cdef fmpz_mod_mpoly res
        cdef fmpz other = <fmpz>arg
        res = create_fmpz_mod_mpoly(self.ctx)
        fmpz_mod_mpoly_sub_fmpz(res.val, self.val, other.val, self.ctx.val)
        return res

    cdef _mul_scalar_(self, arg):
        cdef fmpz_mod_mpoly res
        cdef fmpz other = <fmpz>arg
        res = create_fmpz_mod_mpoly(self.ctx)
        fmpz_mod_mpoly_scalar_mul_fmpz(res.val, self.val, other.val, self.ctx.val)
        return res

    cdef _pow_(self, arg):
        cdef fmpz_mod_mpoly res
        cdef fmpz other = <fmpz>arg
        res = create_fmpz_mod_mpoly(self.ctx)
        if fmpz_mod_mpoly_pow_fmpz(res.val, self.val, other.val, res.ctx.val) == 0:
            raise ValueError("unreasonably large polynomial")  # pragma: no cover
        return res

    cdef _add_mpoly_(self, arg):
        cdef fmpz_mod_mpoly res, other = <fmpz_mod_mpoly>arg
        res = create_fmpz_mod_mpoly(self.ctx)
        fmpz_mod_mpoly_add(res.val, self.val, other.val, res.ctx.val)
        return res

    cdef _sub_mpoly_(self, arg):
        cdef fmpz_mod_mpoly res, other = <fmpz_mod_mpoly>arg
        res = create_fmpz_mod_mpoly(self.ctx)
        fmpz_mod_mpoly_sub(res.val, self.val, other.val, res.ctx.val)
        return res

    cdef _mul_mpoly_(self, arg):
        cdef fmpz_mod_mpoly res, other = <fmpz_mod_mpoly>arg
        res = create_fmpz_mod_mpoly(self.ctx)
        fmpz_mod_mpoly_mul(res.val, self.val, other.val, res.ctx.val)
        return res

    cdef _divmod_mpoly_(self, arg):
        cdef fmpz_mod_mpoly quotient, remainder, other = <fmpz_mod_mpoly>arg
        quotient = create_fmpz_mod_mpoly(self.ctx)
        remainder = create_fmpz_mod_mpoly(self.ctx)
        fmpz_mod_mpoly_divrem(quotient.val, remainder.val, self.val, other.val, self.ctx.val)
        return (quotient, remainder)

    cdef _floordiv_mpoly_(self, arg):
        cdef fmpz_mod_mpoly quotient, other = <fmpz_mod_mpoly>arg
        quotient = create_fmpz_mod_mpoly(self.ctx)
        fmpz_mod_mpoly_div(quotient.val, self.val, other.val, self.ctx.val)
        return quotient

    cdef _truediv_mpoly_(self, arg):
        cdef fmpz_mod_mpoly quotient, other = <fmpz_mod_mpoly>arg
        quotient = create_fmpz_mod_mpoly(self.ctx)
        if fmpz_mod_mpoly_divides(quotient.val, self.val, other.val, self.ctx.val):
            return quotient
        else:
            raise DomainError("fmpz_mod_mpoly division is not exact")

    cdef _mod_mpoly_(self, arg):
        cdef fmpz_mod_mpoly quotient, remainder, other = <fmpz_mod_mpoly>arg
        quotient = create_fmpz_mod_mpoly(self.ctx)
        remainder = create_fmpz_mod_mpoly(self.ctx)
        fmpz_mod_mpoly_divrem(quotient.val, remainder.val, self.val, other.val, self.ctx.val)
        return remainder

    cdef _rsub_scalar_(self, arg):
        cdef fmpz_mod_mpoly res
        cdef fmpz other = <fmpz>arg
        res = create_fmpz_mod_mpoly(self.ctx)
        fmpz_mod_mpoly_sub_fmpz(res.val, self.val, other.val, self.ctx.val)
        fmpz_mod_mpoly_neg(res.val, res.val, res.ctx.val)
        return res

    cdef _rsub_mpoly_(self, arg):
        return (<fmpz_mod_mpoly>arg)._sub_mpoly_(self)

    cdef _rdivmod_mpoly_(self, arg):
        return (<fmpz_mod_mpoly>arg)._divmod_mpoly_(self)

    cdef _rfloordiv_mpoly_(self, arg):
        return (<fmpz_mod_mpoly>arg)._floordiv_mpoly_(self)

    cdef _rtruediv_mpoly_(self, arg):
        return (<fmpz_mod_mpoly>arg)._truediv_mpoly_(self)

    cdef _rmod_mpoly_(self, arg):
        return (<fmpz_mod_mpoly>arg)._mod_mpoly_(self)

    cdef _iadd_scalar_(self, arg):
        cdef fmpz other = <fmpz>arg
        fmpz_mod_mpoly_add_fmpz(self.val, self.val, other.val, self.ctx.val)

    cdef _isub_scalar_(self, arg):
        cdef fmpz other = <fmpz>arg
        fmpz_mod_mpoly_sub_fmpz(self.val, self.val, other.val, self.ctx.val)

    cdef _imul_scalar_(self, arg):
        cdef fmpz other = <fmpz>arg
        fmpz_mod_mpoly_scalar_mul_fmpz(self.val, self.val, other.val, self.ctx.val)

    cdef _iadd_mpoly_(self, arg):
        cdef fmpz_mod_mpoly other = <fmpz_mod_mpoly>arg
        fmpz_mod_mpoly_add(self.val, self.val, other.val, self.ctx.val)

    cdef _isub_mpoly_(self, arg):
        cdef fmpz_mod_mpoly other = <fmpz_mod_mpoly>arg
        fmpz_mod_mpoly_sub(self.val, self.val, other.val, self.ctx.val)

    cdef _imul_mpoly_(self, arg):
        cdef fmpz_mod_mpoly other = <fmpz_mod_mpoly>arg
        fmpz_mod_mpoly_mul(self.val, self.val, other.val, self.ctx.val)

    def __call__(self, *args) -> fmpz:
        cdef:
            fmpz_vec V
            fmpz vres
            slong nvars = self.ctx.nvars(), nargs = len(args)

        if nargs != nvars:
            raise ValueError("number of generators does not match number of arguments")
        args = [self.ctx._any_as_scalar(x) for x in args]
        V = fmpz_vec(args, double_indirect=True)
        vres = fmpz.__new__(fmpz)
        fmpz_mod_mpoly_evaluate_all_fmpz(vres.val, self.val, V.double_indirect, self.ctx.val)
        return vres

    def monoms(self):
        """
        Return the exponent vectors of each term as a tuple of fmpz.

            >>> ctx = fmpz_mod_mpoly_ctx.get(('x', 2), 11, 'lex')
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
            res.append(tuple(vec))

        return res

    def coeffs(self):
        """
        Return the coefficients of each term as a fmpz

            >>> ctx = fmpz_mod_mpoly_ctx.get(('x', 2), 11, 'lex')
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

    def subs(self, dict_args) -> fmpz_mod_mpoly:
        """
        Partial evaluate this polynomial with select constants. Keys must be generator names or generator indices,
        all values must be fmpz.

            >>> ctx = fmpz_mod_mpoly_ctx.get(('x', 2), 11, 'lex')
            >>> f = ctx.from_dict({(0, 0): 1, (1, 0): 2, (0, 1): 3, (1, 1): 4})
            >>> f.subs({"x1": 0})
            2*x0 + 1

        """
        cdef:
            fmpz_mod_mpoly res
            slong i

        args = tuple((self.ctx.variable_to_index(k), self.ctx._any_as_scalar(v)) for k, v in dict_args.items())
        for (_, v), old in zip(args, dict_args.values()):
            if v is NotImplemented:
                raise TypeError(f"cannot coerce {type(old)} to fmpz")

        # Partial application with args in Z. We evaluate the polynomial one variable at a time
        res = create_fmpz_mod_mpoly(self.ctx)

        fmpz_mod_mpoly_set(res.val, self.val, self.ctx.val)
        for i, arg in args:
            fmpz_mod_mpoly_evaluate_one_fmpz(res.val, res.val, i, (<fmpz>arg).val, self.ctx.val)
        return res

    def compose(self, *args, ctx=None) -> fmpz_mod_mpoly:
        """
        Compose this polynomial with other fmpz_mod_mpolys. All arguments must share the same context, it may different
        from this polynomials context.

            >>> ctx = fmpz_mod_mpoly_ctx.get(('x',), 11, 'lex')
            >>> ctx1 = fmpz_mod_mpoly_ctx.get(('y', 2), 11, 'lex')
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
            slong nvars = self.ctx.nvars(), nargs = len(args)

        if nargs != nvars:
            raise ValueError("number of generators does not match number of arguments")
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

            >>> ctx = fmpz_mod_mpoly_ctx.get(('x', 2), 11, 'lex')
            >>> p = ctx.from_dict({(0, 1): 2})
            >>> ctx is p.context()
            True
        """
        return self.ctx

    def coefficient(self, slong i):
        """
        Return the coefficient at index ``i``.

            >>> ctx = fmpz_mod_mpoly_ctx.get(('x', 2), 11, 'lex')
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
        Return the exponent vector at index ``i`` as a tuple.

            >>> ctx = fmpz_mod_mpoly_ctx.get(('x', 2), 11, 'lex')
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
        return tuple(res)

    def degrees(self):
        """
        Return a dictionary of variable name to degree.

            >>> ctx = fmpz_mod_mpoly_ctx.get(('x', 4), 11, 'lex')
            >>> p = sum(x**i for i, x in enumerate(ctx.gens()))
            >>> p
            x1 + x2^2 + x3^3 + 1
            >>> p.degrees()
            (0, 1, 2, 3)
        """
        cdef:
            slong nvars = self.ctx.nvars()

        res = fmpz_vec(nvars, double_indirect=True)
        fmpz_mod_mpoly_degrees_fmpz(res.double_indirect, self.val, self.ctx.val)
        return tuple(res)

    def total_degree(self):
        """
        Return the total degree.

            >>> ctx = fmpz_mod_mpoly_ctx.get(('x', 4), 11, 'lex')
            >>> p = ctx.from_dict({(1, 0, 0, 0): 1, (0, 2, 0, 0): 2, (0, 0, 3, 0): 3})
            >>> p.total_degree()
            3
        """
        cdef fmpz res = fmpz()
        fmpz_mod_mpoly_total_degree_fmpz((<fmpz> res).val, self.val, self.ctx.val)
        return res

    def leading_coefficient(self):
        """
        Leading coefficient in the monomial ordering.

            >>> ctx = fmpz_mod_mpoly_ctx.get(('x', 'y'), 11, 'lex')
            >>> x, y = ctx.gens()
            >>> p = 2*x*y + 3*x + 4*y**2 + 5
            >>> p
            2*x*y + 3*x + 4*y^2 + 5
            >>> p.leading_coefficient()
            2

        """
        if fmpz_mod_mpoly_is_zero(self.val, self.ctx.val):
            return fmpz(0)
        else:
            return self.coefficient(0)

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

            >>> ctx = fmpz_mod_mpoly_ctx.get(('x', 2), 11, 'lex')
            >>> f = ctx.from_dict({(1, 1): 4, (0, 0): 1})
            >>> g = ctx.from_dict({(0, 1): 2, (1, 0): 2})
            >>> (f * g).gcd(f)
            x0*x1 + 3
        """
        cdef fmpz_mod_mpoly res
        if not typecheck(other, fmpz_mod_mpoly):
            raise TypeError("argument must be a fmpz_mod_mpoly")
        elif (<fmpz_mod_mpoly>self).ctx is not (<fmpz_mod_mpoly>other).ctx:
            raise IncompatibleContextError(f"{(<fmpz_mod_mpoly>self).ctx} is not {(<fmpz_mod_mpoly>other).ctx}")
        elif not self.ctx.is_prime():
            raise DomainError("gcd with non-prime modulus is not supported")
        res = create_fmpz_mod_mpoly(self.ctx)
        fmpz_mod_mpoly_gcd(res.val, (<fmpz_mod_mpoly>self).val, (<fmpz_mod_mpoly>other).val, res.ctx.val)
        return res

    def term_content(self):
        """
        Return the GCD of the terms of ``self``. If ``self`` is zero, then the result will
        be zero, otherwise it will be a monomial with positive coefficient.

            >>> ctx = fmpz_mod_mpoly_ctx.get(('x', 2), 11, 'lex')
            >>> x0, x1 = ctx.gens()
            >>> f = 3 * x0**2 * x1 + 6 * x0 * x1
            >>> f.term_content()
            x0*x1
        """

        cdef fmpz_mod_mpoly res = create_fmpz_mod_mpoly(self.ctx)
        fmpz_mod_mpoly_term_content(res.val, self.val, self.ctx.val)
        return res

    def resultant(self, other, var):
        """
        Return the resultant of ``self`` and ``other`` with respect to variable ``var``.

            >>> ctx = fmpz_mod_mpoly_ctx.get(('x', 2), 11, 'lex')
            >>> x0, x1 = ctx.gens()
            >>> f = x0**2 * x1 + x0 * x1
            >>> g = x0 + x1
            >>> f.resultant(g, 'x1')
            x0^3 + x0^2
        """
        cdef:
            fmpz_mod_mpoly res
            slong i

        if not typecheck(other, fmpz_mod_mpoly):
            raise TypeError("argument must be a fmpz_mod_mpoly")
        elif (<fmpz_mod_mpoly>self).ctx is not (<fmpz_mod_mpoly>other).ctx:
            raise IncompatibleContextError(f"{(<fmpz_mod_mpoly>self).ctx} is not {(<fmpz_mod_mpoly>other).ctx}")

        i = self.ctx.variable_to_index(var)
        res = create_fmpz_mod_mpoly(self.ctx)
        if not fmpz_mod_mpoly_resultant(res.val, self.val, (<fmpz_mod_mpoly>other).val, i, self.ctx.val):
            raise RuntimeError(f"failed to compute resultant with respect to {var}")
        return res

    def discriminant(self, var):
        """
        Return the discriminant of ``self`` with respect to variable ``var``.

            >>> ctx = fmpz_mod_mpoly_ctx.get(('x', 2), 11, 'lex')
            >>> x0, x1 = ctx.gens()
            >>> f = (x0 + x1)**2 + 1
            >>> f.discriminant('x1')
            7

        """
        cdef:
            fmpz_mod_mpoly res
            slong i

        i = self.ctx.variable_to_index(var)
        res = create_fmpz_mod_mpoly(self.ctx)
        if not fmpz_mod_mpoly_discriminant(res.val, self.val, i, self.ctx.val):
            raise RuntimeError(f"failed to compute discriminant with respect to {var}")
        return res

    def sqrt(self):
        """
        Return the square root of self.

            >>> ctx = fmpz_mod_mpoly_ctx.get(('x', 2), 11, 'lex')
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
            raise DomainError("polynomial is not a perfect square")

    def factor(self):
        """
        Factors self into irreducible factors, returning a tuple
        (c, factors) where c is the content of the coefficients and
        factors is a list of (poly, exp) pairs.

            >>> Zm = fmpz_mod_mpoly
            >>> ctx = fmpz_mod_mpoly_ctx.get(('x', 'y', 'z'), 11, 'lex')
            >>> p1 = Zm("2*x + 4", ctx)
            >>> p2 = Zm("3*x*z + 3*x + 3*z + 3", ctx)
            >>> (p1 * p2).factor()
            (6, [(z + 1, 1), (x + 1, 1), (x + 2, 1)])
            >>> (p2 * p1 * p2).factor()
            (7, [(z + 1, 2), (x + 2, 1), (x + 1, 2)])
        """
        cdef:
            fmpz_mod_mpoly_factor_t fac
            fmpz c
            fmpz_mod_mpoly u

        if not self.ctx.is_prime():
            raise DomainError("factorisation with non-prime modulus is not supported")

        fmpz_mod_mpoly_factor_init(fac, self.ctx.val)
        if not fmpz_mod_mpoly_factor(fac, self.val, self.ctx.val):
            raise RuntimeError("factorisation failed")
        res = [0] * fac.num

        for i in range(fac.num):
            u = create_fmpz_mod_mpoly(self.ctx)
            fmpz_mod_mpoly_set((<fmpz_mod_mpoly>u).val, &fac.poly[i], self.ctx.val)

            c = fmpz.__new__(fmpz)
            fmpz_set((<fmpz>c).val, &fac.exp[i])

            res[i] = (u, int(c))

        c = fmpz.__new__(fmpz)
        fmpz_set((<fmpz>c).val, fac.constant)
        fmpz_mod_mpoly_factor_clear(fac, self.ctx.val)
        return c, res

    def factor_squarefree(self):
        """
        Factors self into irreducible factors, returning a tuple
        (c, factors) where c is the content of the coefficients and
        factors is a list of (poly, exp) pairs.

            >>> Zm = fmpz_mod_mpoly
            >>> ctx = fmpz_mod_mpoly_ctx.get(('x', 'y', 'z'), 11, 'lex')
            >>> p1 = Zm("2*x + 4", ctx)
            >>> p2 = Zm("3*x*y + 3*x + 3*y + 3", ctx)
            >>> (p1 * p2).factor_squarefree()
            (6, [(y + 1, 1), (x^2 + 3*x + 2, 1)])
            >>> (p1 * p2 * p1).factor_squarefree()
            (1, [(y + 1, 1), (x + 1, 1), (x + 2, 2)])
        """
        cdef:
            fmpz_mod_mpoly_factor_t fac
            fmpz c
            fmpz_mod_mpoly u

        if not self.ctx.is_prime():
            raise DomainError("factorisation with non-prime modulus is not supported")

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

            res[i] = (u, int(c))

        c = fmpz.__new__(fmpz)
        fmpz_set((<fmpz>c).val, fac.constant)
        fmpz_mod_mpoly_factor_clear(fac, self.ctx.val)
        return c, res

    def derivative(self, var):
        """
        Return the derivative of this polynomial with respect to the provided variable.
        The argument can either be the variable as a string, or the index of the
        variable in the context.

            >>> ctx = fmpz_mod_mpoly_ctx.get(('x', 2), 11, 'lex')
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

    def inflate(self, N: list[int]) -> fmpz_mod_mpoly:
        """
        Compute the inflation of ``self`` for a provided ``N``, that is return ``q``
        such that ``q(X) = p(X^N)``.

            >>> ctx = fmpz_mod_mpoly_ctx.get(('x', 'y'), 11, 'lex')
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
            fmpz_mod_mpoly res = create_fmpz_mod_mpoly(self.ctx)

        fmpz_mod_mpoly_inflate(res.val, self.val, shift.val, stride.val, self.ctx.val)
        return res

    def deflate(self, N: list[int]) -> fmpz_mod_mpoly:
        """
        Compute the deflation of ``self`` for a provided ``N``, that is return ``q``
        such that ``q(X) = p(X^(1/N))``.

            >>> ctx = fmpz_mod_mpoly_ctx.get(('x', 'y'), 11, 'lex')
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
            fmpz_mod_mpoly res = create_fmpz_mod_mpoly(self.ctx)

        fmpz_mod_mpoly_deflate(res.val, self.val, shift.val, stride.val, self.ctx.val)
        return res

    def deflation(self) -> tuple[fmpz_mod_mpoly, list[int]]:
        """
        Compute the deflation of ``self``, that is ``p(X^(1/N))`` for maximal
        N. Returns ``q, N`` such that ``self == q.inflate(N)``.

            >>> ctx = fmpz_mod_mpoly_ctx.get(('x', 'y'), 11, 'lex')
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
            fmpz_mod_mpoly res = create_fmpz_mod_mpoly(self.ctx)

        fmpz_mod_mpoly_deflation(shift.val, stride.val, self.val, self.ctx.val)

        for i in range(nvars):
            stride[i] = shift[i].gcd(stride[i])
            shift[i] = 0

        fmpz_mod_mpoly_deflate(res.val, self.val, shift.val, stride.val, self.ctx.val)

        return res, list(stride)

    def deflation_monom(self) -> tuple[fmpz_mod_mpoly, list[int], fmpz_mod_mpoly]:
        """
        Compute the exponent vector ``N`` and monomial ``m`` such that ``p(X^(1/N))
        = m * q(X^N)`` for maximal N. The returned monomial allows the undo-ing of the
        deflation.

            >>> ctx = fmpz_mod_mpoly_ctx.get(('x', 'y'), 11, 'lex')
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
            fmpz_mod_mpoly res = create_fmpz_mod_mpoly(self.ctx)
            fmpz_mod_mpoly monom = create_fmpz_mod_mpoly(self.ctx)
            fmpz_vec shift = fmpz_vec(nvars)
            fmpz_vec stride = fmpz_vec(nvars)

        fmpz_mod_mpoly_deflation(shift.val, stride.val, self.val, self.ctx.val)
        fmpz_mod_mpoly_push_term_ui_ffmpz(monom.val, 1, fmpz_vec(shift).val, self.ctx.val)
        fmpz_mod_mpoly_deflate(res.val, self.val, shift.val, stride.val, self.ctx.val)

        return res, list(stride), monom

    def deflation_index(self) -> tuple[list[int], list[int]]:
        """
        Compute the exponent vectors ``N`` and ``I`` such that ``p(X^(1/N)) = X^I *
        q(X^N)`` for maximal N. Importantly the deflation itself is not computed
        here. The returned exponent vector ``I`` is the shift that was applied to the
        exponents. It is the exponent vector of the monomial returned by
        ``deflation_monom``.

            >>> ctx = fmpz_mod_mpoly_ctx.get(('x', 'y'), 11, 'lex')
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

        fmpz_mod_mpoly_deflation(shift.val, stride.val, self.val, self.ctx.val)
        return list(stride), list(shift)

    cdef _compose_gens_(self, ctx, slong *mapping):
        # FIXME: Remove this when FLINT < 3.2 is dropped
        cdef fmpz_mod_mpoly res
        if FLINT_RELEASE >= 30200:
            res = create_fmpz_mod_mpoly(ctx)
            compat_fmpz_mod_mpoly_compose_fmpz_mod_mpoly_gen(
                res.val,
                self.val,
                mapping,
                self.ctx.val,
                (<fmpz_mod_mpoly_ctx>ctx).val
            )
            return res

        cdef:
            fmpz_mpoly_ctx mpoly_ctx = fmpz_mpoly_ctx.from_context(self.context())
            fmpz_mpoly_ctx res_ctx = fmpz_mpoly_ctx.from_context(ctx)

            fmpz_mpoly poly = mpoly_ctx.from_dict(self.to_dict())
            fmpz_mpoly res1 = poly._compose_gens_(res_ctx, mapping)

        return ctx.from_dict(res1.to_dict())


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

    def __iter__(self):
        cdef fmpz_mod_mpoly z
        for i in range(self.val.length):
            z = create_fmpz_mod_mpoly(self.ctx)
            fmpz_mod_mpoly_set(z.val, &self.val[i], self.ctx.val)
            yield z
