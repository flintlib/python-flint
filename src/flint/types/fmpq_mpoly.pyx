from cpython.dict cimport PyDict_Size, PyDict_Check, PyDict_Next
from cpython.tuple cimport PyTuple_Check, PyTuple_GET_SIZE

from flint.types.fmpq cimport any_as_fmpq, fmpq, fmpq_set_any_ref
from flint.types.fmpz cimport fmpz, fmpz_set_any_ref, any_as_fmpz
from flint.types.fmpz_mpoly cimport fmpz_mpoly, fmpz_mpoly_ctx

from flint.utils.typecheck cimport typecheck
from flint.utils.conversion cimport str_from_chars

from flint.flint_base.flint_base cimport flint_mpoly_context

from flint.flintlib.flint cimport FMPZ_UNKNOWN
from flint.flintlib.fmpq cimport fmpq_init, fmpq_clear, fmpq_is_zero, fmpq_set, fmpq_one
from flint.flintlib.fmpz cimport fmpz_clear, fmpz_init, fmpz_set
from flint.flintlib.fmpz_mpoly cimport fmpz_mpoly_set
from flint.flintlib.fmpq_mpoly_factor cimport fmpq_mpoly_factor_t, \
    fmpq_mpoly_factor_init, fmpq_mpoly_factor, fmpq_mpoly_factor_clear, \
    fmpq_mpoly_factor_squarefree

cdef extern from *:
    """
    /* An ugly hack to get around the ugly hack of renaming fmpq to avoid a c/python name collision */
    typedef fmpq fmpq_struct;
    """


cimport cython
cimport libc.stdlib

cdef dict _fmpq_mpoly_ctx_cache = {}


@cython.auto_pickle(False)
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
        if ordering == "lex":
            fmpq_mpoly_ctx_init(self.val, nvars, ordering_t.ORD_LEX)
        elif ordering == "deglex":
            fmpq_mpoly_ctx_init(self.val, nvars, ordering_t.ORD_DEGLEX)
        elif ordering == "degrevlex":
            fmpq_mpoly_ctx_init(self.val, nvars, ordering_t.ORD_DEGREVLEX)
        else:
            raise ValueError("Unimplemented term order %s" % ordering)

        super().__init__(nvars, names)

    cpdef slong nvars(self):
        """
        Return the number of variables in the context

            >>> ctx = fmpq_mpoly_ctx.get_context(4, "lex", 'x')
            >>> ctx.nvars()
            4
        """
        return self.val.zctx.minfo.nvars

    cpdef ordering(self):
        """
        Return the term order of the context object.

            >>> ctx = fmpq_mpoly_ctx.get_context(4, "deglex", 'w')
            >>> ctx.ordering()
            'deglex'
        """
        if self.val.zctx.minfo.ord == ordering_t.ORD_LEX:
            return "lex"
        if self.val.zctx.minfo.ord == ordering_t.ORD_DEGLEX:
            return "deglex"
        if self.val.zctx.minfo.ord == ordering_t.ORD_DEGREVLEX:
            return "degrevlex"

    def gen(self, slong i):
        """
        Return the `i`th generator of the polynomial ring

            >>> ctx = fmpq_mpoly_ctx.get_context(3, 'degrevlex', 'z')
            >>> ctx.gen(1)
            z1
        """
        cdef fmpq_mpoly res
        assert i >= 0 and i < self.val.zctx.minfo.nvars
        res = fmpq_mpoly.__new__(fmpq_mpoly)
        res.ctx = self
        fmpq_mpoly_init(res.val, res.ctx.val)
        res._init = True
        fmpq_mpoly_gen(res.val, i, res.ctx.val)
        return res

    def constant(self, z):
        """
        Create a constant polynomial in this context
        """
        cdef fmpq_mpoly res
        z = any_as_fmpq(z)
        if z is NotImplemented:
            raise ValueError("A constant fmpq_mpoly is a fmpq")
        res = fmpq_mpoly.__new__(fmpq_mpoly)
        res.ctx = self
        fmpq_mpoly_init(res.val, res.ctx.val)
        res._init = True
        fmpq_mpoly_set_fmpq(res.val, (<fmpq>z).val, res.ctx.val)
        return res

    def fmpq_mpoly_from_dict(self, d):
        """
        Create a fmpz_mpoly from a dictionary.

        The dictionary's keys are tuples of ints (or anything that implicitly converts
        to fmpz) representing exponents, and corresponding coefficient values of fmpq.

            >>> ctx = fmpq_mpoly_ctx.get_context(2,'lex','x,y')
            >>> ctx.fmpq_mpoly_from_dict({(1,0):2, (1,1):3, (0,1):1})
            3*x*y + 2*x + y
        """
        cdef long n
        cdef fmpq_t coefficient
        cdef int xtype
        cdef fmpz_struct *exponents
        cdef fmpz_struct **exp_ptr
        cdef int nvars = self.nvars()
        cdef int i,j
        cdef int count
        cdef fmpq_mpoly res

        if not PyDict_Check(d):
            raise ValueError("expected a dictionary")
        n = PyDict_Size(d)
        exponents = <fmpz_struct *> libc.stdlib.calloc(nvars, sizeof(fmpz_struct))
        if exponents == NULL:
            raise MemoryError()
        exp_ptr = <fmpz_struct **> libc.stdlib.calloc(nvars, sizeof(fmpz_struct *))
        if exp_ptr == NULL:
            libc.stdlib.free(exponents)
            raise MemoryError()
        for i in range(nvars):
            fmpz_init(exponents + i)
            exp_ptr[i] = exponents + i
        # fmpq_init(coefficient)
        res = fmpq_mpoly.__new__(fmpq_mpoly)
        res.ctx = self
        fmpq_mpoly_init(res.val, res.ctx.val)
        res._init = True
        count = 0
        for k, v in d.items():
            xtype = fmpq_set_any_ref(coefficient, v)
            if xtype == FMPZ_UNKNOWN:
                for i in range(nvars):
                    fmpz_clear(exponents + i)
                libc.stdlib.free(exponents)
                libc.stdlib.free(exp_ptr)
                raise TypeError("invalid coefficient type %s" % type(v))
            if not PyTuple_Check(k):
                for i in range(nvars):
                    fmpz_clear(exponents + i)
                libc.stdlib.free(exponents)
                raise TypeError("Expected tuple of ints as key not %s" % type(k))
            if PyTuple_GET_SIZE(k) != nvars:
                for i in range(nvars):
                    fmpz_clear(exponents + i)
                libc.stdlib.free(exponents)
                raise TypeError("Expected exponent tuple of length %d" % nvars)
            for i, tup in enumerate(k):
                xtype = fmpz_set_any_ref(exponents + i, tup)
                if xtype == FMPZ_UNKNOWN:
                    for i in range(nvars):
                        fmpz_clear(exponents + i)
                    libc.stdlib.free(exponents)
                    raise TypeError("Invalid exponent type %s" % type(tup))
            #Todo lobby for fmpz_mpoly_push_term_fmpz_ffmpz
            if not fmpq_is_zero(coefficient):
                fmpq_mpoly_push_term_fmpq_fmpz(res.val, coefficient, exp_ptr, self.val)
                # _fmpq_mpoly_push_exp_ffmpz(res.val, exponents, self.val)
                # fmpq_mpoly_set_term_coeff_fmpz(res.val, count, coefficient, self.val)
                count += 1
        for i in range(nvars):
            fmpz_clear(exponents + i)
            fmpq_mpoly_sort_terms(res.val, self.val)
            fmpq_mpoly_reduce(res.val, self.val)
        return res


cdef class fmpq_mpoly(flint_mpoly):
    """
    The *fmpz_poly* type represents sparse multivariate polynomials over
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
                raise ValueError("Cannot automatically coerce contexts")
        if typecheck(val, fmpz_mpoly):
            if ctx is None:
                ctx = fmpq_mpoly_ctx.from_context((<fmpz_mpoly>val).ctx)
            elif ctx != fmpq_mpoly_ctx.from_context((<fmpz_mpoly>val).ctx):
                raise ValueError("Cannot automatically coerce contexts")
            init_fmpq_mpoly(self, ctx)
            fmpz_mpoly_set(self.val.zpoly, (<fmpz_mpoly>val).val, (<fmpz_mpoly> val).ctx.val)
            fmpq_one(self.val.content)
            fmpq_mpoly_reduce(self.val, self.ctx.val)
        elif isinstance(val, dict):
            if ctx is None:
                if len(val) == 0:
                    raise ValueError("Need context for zero polynomial")
                k = list(val.keys())[0]
                if not isinstance(k, tuple):
                    raise ValueError("Dict should be keyed with tuples of integers")
                ctx = fmpq_mpoly_ctx.get_context(len(k))
            x = ctx.fmpq_mpoly_from_dict(val)
            #XXX this copy is silly, have a ctx function that assigns an fmpz_mpoly_t
            init_fmpq_mpoly(self, ctx)
            fmpq_mpoly_set(self.val, (<fmpq_mpoly>x).val, self.ctx.val)
        elif isinstance(val, str):
            if ctx is None:
                raise ValueError("Cannot parse a polynomial without context")
            val = bytes(val, 'utf-8')
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

    def __nonzero__(self):
        return not fmpq_mpoly_is_zero(self.val, self.ctx.val)

    def __bool__(self):
        return not fmpq_mpoly_is_zero(self.val, self.ctx.val)

    def is_one(self):
        return fmpq_mpoly_is_one(self.val, self.ctx.val)

    def __richcmp__(self, other, int op):
        if op != 2 and op != 3:
            return NotImplemented
        if typecheck(self, fmpq_mpoly) and typecheck(other, fmpq_mpoly):
            if (<fmpq_mpoly>self).ctx is (<fmpq_mpoly>other).ctx:
                if op == 2:
                    return bool(fmpq_mpoly_equal((<fmpq_mpoly>self).val, (<fmpq_mpoly>other).val, (<fmpq_mpoly>self).ctx.val))
                else:
                    return not bool(fmpq_mpoly_equal((<fmpq_mpoly>self).val, (<fmpq_mpoly>other).val, (<fmpq_mpoly>self).ctx.val))
            else:
                if op == 2:
                    return False
                else:
                    return True
        if op == 2:
            return not bool(self - other)
        else:
            return bool(self - other)

    def context(self):
        return self.ctx

    def __len__(self):
        return fmpq_mpoly_length(self.val, self.ctx.val)

    def coefficient(self, slong i):
        cdef fmpq v
        if i < 0 or i >= fmpq_mpoly_length(self.val, self.ctx.val):
            return fmpq(0)
        else:
            v = fmpq.__new__(fmpz)
            fmpq_mpoly_get_term_coeff_fmpq(v.val, self.val, i, self.ctx.val)
            return v

    def exponent_tuple(self, slong i):
        cdef slong j, nvars
        cdef fmpz_struct ** tmp
        if i < 0 or i >= fmpq_mpoly_length(self.val, self.ctx.val):
            raise ValueError
        nvars = self.ctx.nvars()
        res = tuple(fmpz() for j in range(nvars))
        tmp = <fmpz_struct **> libc.stdlib.malloc(nvars * sizeof(fmpz_struct **))
        try:
            for j in range(nvars):
                tmp[j] = &((<fmpz> (res[j])).val[0])
            fmpq_mpoly_get_term_exp_fmpz(tmp, self.val, i, self.ctx.val)
        finally:
            libc.stdlib.free(tmp)
        return res

    def repr(self):
        return self.str() + "  (nvars=%s, ordering=%s names=%s)" % (self.ctx.nvars(), self.ctx.ordering(), self.ctx.py_names)

    def str(self):
        cdef char * s = fmpq_mpoly_get_str_pretty(self.val, self.ctx.c_names, self.ctx.val)
        try:
            res = str_from_chars(s)
        finally:
            libc.stdlib.free(s)
        if res.startswith(" - "):
            res = "-" + res[3:]
        return res

    def __neg__(self):
        cdef fmpq_mpoly res
        res = create_fmpq_mpoly(self.ctx)
        fmpq_mpoly_neg(res.val, (<fmpq_mpoly>self).val, res.ctx.val)
        return res

    def __add__(self, other):
        cdef fmpq_mpoly res
        if typecheck(other, fmpq_mpoly):
            if (<fmpq_mpoly>self).ctx is not (<fmpq_mpoly>other).ctx:
                return NotImplemented
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
                return NotImplemented
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
                return NotImplemented
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
                return NotImplemented
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
                return NotImplemented
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
                return NotImplemented
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
            raise ValueError("cannot raise fmpz_mpoly to negative power")
        res = create_fmpq_mpoly(self.ctx)
        if fmpq_mpoly_pow_fmpz(res.val, (<fmpq_mpoly>self).val, (<fmpz>other).val, res.ctx.val) == 0:
            raise ValueError("unreasonably large polynomial")
        return res

    def __divmod__(self, other):
        cdef fmpq_mpoly res, res2
        if typecheck(other, fmpq_mpoly):
            if (<fmpq_mpoly>self).ctx is not (<fmpq_mpoly>other).ctx:
                return NotImplemented
            res = create_fmpq_mpoly(self.ctx)
            res2 = create_fmpq_mpoly(self.ctx)
            fmpq_mpoly_divrem(res.val, res2.val, (<fmpq_mpoly>self).val, (<fmpq_mpoly>other).val, res.ctx.val)
            return (res, res2)
        else:
            other = any_as_fmpq(other)
            if other is not NotImplemented:
                other= fmpq_mpoly(other, self.ctx)
                res = create_fmpq_mpoly(self.ctx)
                res2 = create_fmpq_mpoly(self.ctx)
                fmpq_mpoly_divrem(res.val, res2.val, (<fmpq_mpoly>self).val, (<fmpq_mpoly>other).val, res.ctx.val)
                return (res, res2)
        return NotImplemented

    def __rdivmod__(self, other):
        cdef fmpq_mpoly res, res2
        other = any_as_fmpq(other)
        if other is not NotImplemented:
            other = fmpq_mpoly(other, self.ctx)
            res = create_fmpq_mpoly(self.ctx)
            res2 = create_fmpq_mpoly(self.ctx)
            fmpq_mpoly_divrem(res.val, res2.val, (<fmpq_mpoly>other).val, (<fmpq_mpoly>self).val, res.ctx.val)
            return res
        return NotImplemented

    def __floordiv__(self, other):
        cdef fmpq_mpoly res
        if typecheck(other, fmpq_mpoly):
            if (<fmpq_mpoly>self).ctx is not (<fmpq_mpoly>other).ctx:
                return NotImplemented
            res = create_fmpq_mpoly(self.ctx)
            fmpq_mpoly_div(res.val, (<fmpq_mpoly>self).val, (<fmpq_mpoly>other).val, res.ctx.val)
            return res
        else:
            other = any_as_fmpq(other)
            if other is not NotImplemented:
                other = fmpq_mpoly(other, self.ctx)
                res = create_fmpq_mpoly(self.ctx)
                fmpq_mpoly_div(res.val, (<fmpq_mpoly>self).val, (<fmpq_mpoly>other).val, res.ctx.val)
                return res
        return NotImplemented

    def __rfloordiv__(self, other):
        cdef fmpq_mpoly res
        other = any_as_fmpq(other)
        if other is not NotImplemented:
            other = fmpq_mpoly(other, self.ctx)
            res = create_fmpq_mpoly(self.ctx)
            fmpq_mpoly_div(res.val, (<fmpq_mpoly>other).val, self.val, res.ctx.val)
            return res
        return NotImplemented

    def __mod__(self, other):
        return divmod(self, other)[1]

    def gcd(self, other):
        cdef fmpq_mpoly res
        assert isinstance(other, fmpq_mpoly)
        if (<fmpq_mpoly>self).ctx is not (<fmpq_mpoly>other).ctx:
            return NotImplemented
        res = create_fmpq_mpoly(self.ctx)
        fmpq_mpoly_gcd(res.val, (<fmpq_mpoly>self).val, (<fmpq_mpoly>other).val, res.ctx.val)
        return res

    def factor(self):
        """
        Factors self into irreducible factors, returning a tuple
        (c, factors) where c is the content of the coefficients and
        factors is a list of (poly, exp) pairs.

            >>> Zm = fmpq_mpoly
            >>> ctx = fmpq_mpoly_ctx.get_context(3, 'lex', 'x,y,z')
            >>> p1 = Zm("2*x + 4", ctx)
            >>> p2 = Zm("3*x*z +  + 3*x + 3*z + 3", ctx)
            >>> (p1 * p2).factor()
            (6, [(z + 1, 1), (x + 2, 1), (x + 1, 1)])
            >>> (p2 * p1 * p2).factor()
            (18, [(z + 1, 2), (x + 2, 1), (x + 1, 2)])
        """
        cdef fmpq_mpoly_factor_t fac
        cdef int i
        cdef fmpq c
        cdef fmpz exp
        cdef fmpq_mpoly u
        fmpq_mpoly_factor_init(fac, self.ctx.val)
        fmpq_mpoly_factor(fac, self.val, self.ctx.val)
        res = [0] * fac.num
        for i in range(fac.num):
            u = fmpq_mpoly.__new__(fmpq_mpoly)
            u.ctx = self.ctx
            fmpq_mpoly_init(u.val, u.ctx.val)
            u._init = True
            fmpq_mpoly_set((<fmpq_mpoly>u).val, &fac.poly[i], self.ctx.val)
            exp = fmpz.__new__(fmpz)
            fmpz_set((<fmpz>exp).val, fac.exp + i)
            res[i] = (u, exp)
        c = fmpq.__new__(fmpq)
        fmpq_set((<fmpq>c).val, fac.constant)
        fmpq_mpoly_factor_clear(fac, self.ctx.val)
        return c, res

    def factor_squarefree(self):
        """
        Factors self into irreducible factors, returning a tuple
        (c, factors) where c is the content of the coefficients and
        factors is a list of (poly, exp) pairs.

            >>> Zm = fmpq_mpoly
            >>> ctx = fmpq_mpoly_ctx.get_context(3, 'lex', 'x,y,z')
            >>> p1 = Zm("2*x + 4", ctx)
            >>> p2 = Zm("3*x*y + 3*x + 3*y + 3", ctx)
            >>> (p1 * p2).factor_squarefree()
            (6, [(y + 1, 1), (x^2 + 3*x + 2, 1)])
            >>> (p1 * p2 * p1).factor_squarefree()
            (12, [(y + 1, 1), (x + 1, 1), (x + 2, 2)])
        """
        cdef fmpq_mpoly_factor_t fac
        cdef int i
        cdef fmpq c
        cdef fmpz exp
        cdef fmpq_mpoly u
        fmpq_mpoly_factor_init(fac, self.ctx.val)
        fmpq_mpoly_factor_squarefree(fac, self.val, self.ctx.val)
        res = [0] * fac.num
        for 0 <= i < fac.num:
            u = fmpq_mpoly.__new__(fmpq_mpoly)
            u.ctx = self.ctx
            fmpq_mpoly_init(u.val, u.ctx.val)
            u._init = True
            fmpq_mpoly_set((<fmpq_mpoly>u).val, &fac.poly[i], self.ctx.val)
            exp = fmpz.__new__(fmpz)
            fmpz_set((<fmpz>exp).val, fac.exp + i)
            res[i] = (u, exp)
        c = fmpq.__new__(fmpq)
        fmpq_set((<fmpq>c).val, fac.constant)
        fmpq_mpoly_factor_clear(fac, self.ctx.val)
        return c, res
