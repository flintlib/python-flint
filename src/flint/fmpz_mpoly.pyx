cdef any_as_fmpz_mpoly(x):
    cdef fmpz_mpoly res
    """
    if typecheck(x, fmpz_poly):
        return x
    elif typecheck(x, fmpz):
        res = fmpz_poly.__new__(fmpz_poly)
        fmpz_poly_set_fmpz(res.val, (<fmpz>x).val)
        return res
    elif PY_MAJOR_VERSION < 3 and PyInt_Check(<PyObject*>x):
        res = fmpz_poly.__new__(fmpz_poly)
        fmpz_poly_set_si(res.val, PyInt_AS_LONG(<PyObject*>x))
        return res
    elif PyLong_Check(<PyObject*>x):
        res = fmpz_poly.__new__(fmpz_poly)
        t = fmpz(x)
        fmpz_poly_set_fmpz(res.val, (<fmpz>t).val)
        return res
    """
    return NotImplemented

"""
cdef fmpz_poly_set_list(fmpz_poly_t poly, list val):
    cdef long i, n
    cdef fmpz_t x
    n = PyList_GET_SIZE(<PyObject*>val)
    fmpz_poly_fit_length(poly, n)
    fmpz_init(x)
    for i from 0 <= i < n:
        if typecheck(val[i], fmpz):
            fmpz_poly_set_coeff_fmpz(poly, i, (<fmpz>(val[i])).val)
        elif fmpz_set_python(x, val[i]):
            fmpz_poly_set_coeff_fmpz(poly, i, x)
        else:
            fmpz_clear(x)
            raise TypeError("unsupported coefficient in list")
    fmpz_clear(x)
"""

cdef dict _fmpz_mpoly_ctx_cache = {}

@cython.auto_pickle(False)
cdef class fmpz_mpoly_ctx:
    cdef fmpz_mpoly_ctx_t val

    def __init__(self, slong nvars, ordering="lex"):
        assert nvars >= 1
        fmpz_mpoly_ctx_init(self.val, nvars, ordering_t.ORD_LEX)

    cpdef slong nvars(self):
        return self.val.minfo.nvars

    cpdef ordering(self):
        return "lex"

cdef get_fmpz_mpoly_context(slong nvars=1, ordering=None):
    if nvars <= 0:
        nvars = 1
    if ordering is None:
        ordering = "lex"
    key = (nvars, ordering)
    ctx = _fmpz_mpoly_ctx_cache.get(key)
    if ctx is None:
        ctx = fmpz_mpoly_ctx(nvars, ordering)
        _fmpz_mpoly_ctx_cache[key] = ctx
    return ctx

cdef _fmpz_mpoly_set2(fmpz_mpoly_t out, fmpz_mpoly_ctx_t outctx, fmpz_mpoly_t inp, fmpz_mpoly_ctx_t inpctx):
    cdef slong * C
    cdef slong i
    cdef slong inpvars, outvars
    if outctx == inpctx:
        fmpz_mpoly_set(out, inp, inpctx)
    else:
        inpvars = inpctx.minfo.nvars
        outvars = inpctx.minfo.nvars
        C = <slong *> libc.stdlib.malloc(inpvars * sizeof(slong *))
        for i in range(min(outvars, inpvars)):
            C[i] = i
        for i in range(outvars, inpvars):
            C[i] = -1
        fmpz_mpoly_compose_fmpz_mpoly_gen(out, inp, C, inpctx, outctx)
        libc.stdlib.free(C)

def coerce_fmpz_mpolys(*args):
    cdef fmpz_mpoly_ctx ctx
    ctx = get_fmpz_mpoly_context()
    if not args:
        return ctx, []
    args = list(args)
    if typecheck(args[0], fmpz_mpoly):
        ctx = (<fmpz_mpoly> args[0]).ctx
        if all(typecheck(args[i], fmpz_mpoly) and (<fmpz_mpoly> args[i]).ctx is ctx for i in range(1, len(args))):
            return ctx, args
    for i in range(len(args)):
        if not typecheck(args[i], fmpz_mpoly):
            args[i] = fmpz_mpoly(args[i])
    nvars = max((<fmpz_mpoly>pol).ctx.nvars() for pol in args)
    ctx = get_fmpz_mpoly_context(nvars)
    args2 = [fmpz_mpoly() for i in range(len(args))]
    for i in range(len(args)):
        (<fmpz_mpoly> args2[i]).ctx = ctx
        _fmpz_mpoly_set2((<fmpz_mpoly> args2[i]).val, ctx.val, (<fmpz_mpoly> args[i]).val, (<fmpz_mpoly> args[i]).ctx.val)
    return ctx, args2


# todo: store cached context objects externally
cdef class fmpz_mpoly(flint_mpoly):
    """
    The *fmpz_poly* type represents sparse multivariate polynomials over
    the integers.
    """

    cdef fmpz_mpoly_t val
    cdef fmpz_mpoly_ctx ctx
    cdef bint _init

    def __cinit__(self):
        self._init = False

    def __dealloc__(self):
        if self._init:
            fmpz_mpoly_clear(self.val, self.ctx.val)
            self._init = False

    def __init__(self, val=0, slong nvars=-1, ordering=None):
        if typecheck(val, fmpz_mpoly):
            if nvars == -1 and ordering is None:
                self.ctx = (<fmpz_mpoly>val).ctx
                fmpz_mpoly_init(self.val, self.ctx.val)
                self._init = True
                fmpz_mpoly_set(self.val, (<fmpz_mpoly>val).val, self.ctx.val)
            else:
                self.ctx = get_fmpz_mpoly_context(nvars, ordering)
                fmpz_mpoly_init(self.val, self.ctx.val)
                self._init = True
                _fmpz_mpoly_set2(self.val, self.ctx.val, (<fmpz_mpoly>val).val, (<fmpz_mpoly>val).ctx.val)
        else:
            v = any_as_fmpz(val)
            if v is NotImplemented:
                raise TypeError("cannot create fmpz_mpoly from type %s" % type(val))
            self.ctx = get_fmpz_mpoly_context(nvars, ordering)
            fmpz_mpoly_init(self.val, self.ctx.val)
            self._init = True
            fmpz_mpoly_set_fmpz(self.val, (<fmpz>v).val, self.ctx.val)

    def __nonzero__(self):
        return not fmpz_mpoly_is_zero(self.val, self.ctx.val)

    def is_one(self):
        return fmpz_mpoly_is_one(self.val, self.ctx.val)

    def __richcmp__(self, other, int op):
        if op != 2 and op != 3:
            return NotImplemented
        if typecheck(self, fmpz_mpoly) and typecheck(other, fmpz_mpoly):
            if (<fmpz_mpoly>self).ctx is (<fmpz_mpoly>other).ctx:
                if op == 2:
                    return bool(fmpz_mpoly_equal((<fmpz_mpoly>self).val, (<fmpz_mpoly>other).val, (<fmpz_mpoly>self).ctx.val))
                else:
                    return not bool(fmpz_mpoly_equal((<fmpz_mpoly>self).val, (<fmpz_mpoly>other).val, (<fmpz_mpoly>self).ctx.val))
        if op == 2:
            return not bool(self - other)
        else:
            return bool(self - other)

    def __len__(self):
        return fmpz_mpoly_length(self.val, self.ctx.val)

    def __hash__(self):
        s = str(self)
        i = s.index("(nvars")
        s = s[:i]
        return hash(s)

    def coefficient(self, slong i):
        cdef fmpz v
        if i < 0 or i >= fmpz_mpoly_length(self.val, self.ctx.val):
            return fmpz(0)
        else:
            v = fmpz.__new__(fmpz)
            fmpz_mpoly_get_term_coeff_fmpz(v.val, self.val, i, self.ctx.val)
            return v

    def leading_coefficient(self):
        return self.coefficient(0)

    def exponent_tuple(self, slong i):
        cdef slong j, nvars
        cdef fmpz_struct ** tmp
        if i < 0 or i >= fmpz_mpoly_length(self.val, self.ctx.val):
            raise ValueError
        nvars = self.ctx.nvars()
        res = tuple(fmpz() for j in range(nvars))
        tmp = <fmpz_struct **> libc.stdlib.malloc(nvars * sizeof(fmpz_struct **))
        try:
            for j in range(nvars):
                tmp[j] = &((<fmpz> (res[j])).val[0])
            fmpz_mpoly_get_term_exp_fmpz(tmp, self.val, i, self.ctx.val)
        finally:
            libc.stdlib.free(tmp)
        return res

    @staticmethod
    def gen(slong i, slong nvars=-1, ordering=None):
        cdef fmpz_mpoly res
        assert i >= 0
        if nvars <= 0:
            nvars = i + 1
        assert i < nvars
        res = fmpz_mpoly.__new__(fmpz_mpoly)
        res.ctx = get_fmpz_mpoly_context(nvars, ordering)
        fmpz_mpoly_init(res.val, res.ctx.val)
        res._init = True
        fmpz_mpoly_gen(res.val, i, res.ctx.val)
        return res

    @staticmethod
    def gens(slong n, ordering=None):
        # todo: (i, n)? or just (i)?
        return tuple(fmpz_mpoly.gen(i, n) for i in range(n))

    def repr(self):
        cdef char * s = fmpz_mpoly_get_str_pretty(self.val, NULL, self.ctx.val)
        try:
            res = str_from_chars(s)
        finally:
            libc.stdlib.free(s)
        res = res.replace("+", " + ")
        res = res.replace("-", " - ")
        if res.startswith(" - "):
            res = "-" + res[3:]
        return res + "  (nvars=%s, ordering=%s)" % (self.ctx.nvars(), self.ctx.ordering())

    def str(self):
        return self.repr()

    def __neg__(self):
        cdef fmpz_mpoly res
        res = fmpz_mpoly.__new__(fmpz_mpoly)
        res.ctx = (<fmpz_mpoly>self).ctx
        fmpz_mpoly_init(res.val, res.ctx.val)
        res._init = True
        fmpz_mpoly_neg(res.val, (<fmpz_mpoly>self).val, res.ctx.val)
        return res

    def __add__(self, other):
        cdef fmpz_mpoly res
        if typecheck(self, fmpz_mpoly):
            if typecheck(other, fmpz_mpoly):
                if (<fmpz_mpoly>self).ctx is not (<fmpz_mpoly>other).ctx:
                    ctx, (self, other) = coerce_fmpz_mpolys(self, other)
                res = fmpz_mpoly.__new__(fmpz_mpoly)
                res.ctx = (<fmpz_mpoly>self).ctx
                fmpz_mpoly_init(res.val, res.ctx.val)
                res._init = True
                fmpz_mpoly_add(res.val, (<fmpz_mpoly>self).val, (<fmpz_mpoly>other).val, res.ctx.val)
                return res
            else:
                other = any_as_fmpz(other)
                if other is not NotImplemented:
                    res = fmpz_mpoly.__new__(fmpz_mpoly)
                    res.ctx = (<fmpz_mpoly>self).ctx
                    fmpz_mpoly_init(res.val, res.ctx.val)
                    res._init = True
                    fmpz_mpoly_add_fmpz(res.val, (<fmpz_mpoly>self).val, (<fmpz>other).val, res.ctx.val)
                    return res
        else:
            self = any_as_fmpz(self)
            if self is not NotImplemented:
                res = fmpz_mpoly.__new__(fmpz_mpoly)
                res.ctx = (<fmpz_mpoly>other).ctx
                fmpz_mpoly_init(res.val, res.ctx.val)
                res._init = True
                fmpz_mpoly_add_fmpz(res.val, (<fmpz_mpoly>other).val, (<fmpz>self).val, res.ctx.val)
                return res
        return NotImplemented

    def __sub__(self, other):
        cdef fmpz_mpoly res
        if typecheck(self, fmpz_mpoly):
            if typecheck(other, fmpz_mpoly):
                if (<fmpz_mpoly>self).ctx is not (<fmpz_mpoly>other).ctx:
                    ctx, (self, other) = coerce_fmpz_mpolys(self, other)
                res = fmpz_mpoly.__new__(fmpz_mpoly)
                res.ctx = (<fmpz_mpoly>self).ctx
                fmpz_mpoly_init(res.val, res.ctx.val)
                res._init = True
                fmpz_mpoly_sub(res.val, (<fmpz_mpoly>self).val, (<fmpz_mpoly>other).val, res.ctx.val)
                return res
            else:
                other = any_as_fmpz(other)
                if other is not NotImplemented:
                    res = fmpz_mpoly.__new__(fmpz_mpoly)
                    res.ctx = (<fmpz_mpoly>self).ctx
                    fmpz_mpoly_init(res.val, res.ctx.val)
                    res._init = True
                    fmpz_mpoly_sub_fmpz(res.val, (<fmpz_mpoly>self).val, (<fmpz>other).val, res.ctx.val)
                    return res
        else:
            self = any_as_fmpz(self)
            if self is not NotImplemented:
                res = fmpz_mpoly.__new__(fmpz_mpoly)
                res.ctx = (<fmpz_mpoly>other).ctx
                fmpz_mpoly_init(res.val, res.ctx.val)
                res._init = True
                fmpz_mpoly_sub_fmpz(res.val, (<fmpz_mpoly>other).val, (<fmpz>self).val, res.ctx.val)
                fmpz_mpoly_neg(res.val, res.val, res.ctx.val)
                return res
        return NotImplemented

    def __mul__(self, other):
        cdef fmpz_mpoly res
        if typecheck(self, fmpz_mpoly):
            if typecheck(other, fmpz_mpoly):
                if (<fmpz_mpoly>self).ctx is not (<fmpz_mpoly>other).ctx:
                    ctx, (self, other) = coerce_fmpz_mpolys(self, other)
                res = fmpz_mpoly.__new__(fmpz_mpoly)
                res.ctx = (<fmpz_mpoly>self).ctx
                fmpz_mpoly_init(res.val, res.ctx.val)
                res._init = True
                fmpz_mpoly_mul(res.val, (<fmpz_mpoly>self).val, (<fmpz_mpoly>other).val, res.ctx.val)
                return res
            else:
                other = any_as_fmpz(other)
                if other is not NotImplemented:
                    res = fmpz_mpoly.__new__(fmpz_mpoly)
                    res.ctx = (<fmpz_mpoly>self).ctx
                    fmpz_mpoly_init(res.val, res.ctx.val)
                    res._init = True
                    fmpz_mpoly_scalar_mul_fmpz(res.val, (<fmpz_mpoly>self).val, (<fmpz>other).val, res.ctx.val)
                    return res
        else:
            self = any_as_fmpz(self)
            if self is not NotImplemented:
                res = fmpz_mpoly.__new__(fmpz_mpoly)
                res.ctx = (<fmpz_mpoly>other).ctx
                fmpz_mpoly_init(res.val, res.ctx.val)
                res._init = True
                fmpz_mpoly_scalar_mul_fmpz(res.val, (<fmpz_mpoly>other).val, (<fmpz>self).val, res.ctx.val)
                return res
        return NotImplemented

    def __pow__(self, other, modulus):
        cdef fmpz_mpoly res
        if modulus is not None:
            raise NotImplementedError
        if typecheck(self, fmpz_mpoly):
            other = any_as_fmpz(other)
            if other is NotImplemented:
                return other
            if other < 0:
                raise ValueError("cannot raise fmpz_mpoly to negative power")
            res = fmpz_mpoly.__new__(fmpz_mpoly)
            res.ctx = (<fmpz_mpoly>self).ctx
            fmpz_mpoly_init(res.val, res.ctx.val)
            res._init = True
            if fmpz_mpoly_pow_fmpz(res.val, (<fmpz_mpoly>self).val, (<fmpz>other).val, res.ctx.val) == 0:
                raise ValueError("unreasonably large polynomial")
            return res
        return NotImplemented

    def __divmod__(self, other):
        cdef fmpz_mpoly res, res2
        if typecheck(self, fmpz_mpoly):
            if typecheck(other, fmpz_mpoly):
                if (<fmpz_mpoly>self).ctx is not (<fmpz_mpoly>other).ctx:
                    ctx, (self, other) = coerce_fmpz_mpolys(self, other)
                res = fmpz_mpoly.__new__(fmpz_mpoly)
                res.ctx = (<fmpz_mpoly>self).ctx
                fmpz_mpoly_init(res.val, res.ctx.val)
                res._init = True
                res2 = fmpz_mpoly.__new__(fmpz_mpoly)
                res2.ctx = (<fmpz_mpoly>self).ctx
                fmpz_mpoly_init(res2.val, res2.ctx.val)
                res2._init = True
                fmpz_mpoly_divrem(res.val, res2.val, (<fmpz_mpoly>self).val, (<fmpz_mpoly>other).val, res.ctx.val)
                return (res, res2)
        return NotImplemented

    def __floordiv__(self, other):
        cdef fmpz_mpoly res
        if typecheck(self, fmpz_mpoly):
            if typecheck(other, fmpz_mpoly):
                if (<fmpz_mpoly>self).ctx is not (<fmpz_mpoly>other).ctx:
                    ctx, (self, other) = coerce_fmpz_mpolys(self, other)
                res = fmpz_mpoly.__new__(fmpz_mpoly)
                res.ctx = (<fmpz_mpoly>self).ctx
                fmpz_mpoly_init(res.val, res.ctx.val)
                res._init = True
                fmpz_mpoly_div(res.val, (<fmpz_mpoly>self).val, (<fmpz_mpoly>other).val, res.ctx.val)
                return res
        return NotImplemented

    def __mod__(self, other):
        return divmod(self, other)[1]

    def gcd(self, other):
        cdef fmpz_mpoly res
        assert isinstance(other, fmpz_mpoly)
        if (<fmpz_mpoly>self).ctx is not (<fmpz_mpoly>other).ctx:
            ctx, (self, other) = coerce_fmpz_mpolys(self, other)
        res = fmpz_mpoly.__new__(fmpz_mpoly)
        res.ctx = (<fmpz_mpoly>self).ctx
        fmpz_mpoly_init(res.val, res.ctx.val)
        res._init = True
        fmpz_mpoly_gcd(res.val, (<fmpz_mpoly>self).val, (<fmpz_mpoly>other).val, res.ctx.val)
        return res

    def __call__(self, *args):
        cdef fmpz_mpoly res
        cdef fmpz_mpoly_ctx res_ctx
        cdef fmpz_struct ** V
        cdef fmpz vres
        cdef fmpz_mpoly_struct ** C
        cdef slong i, nvars, nargs
        other = tuple(args)
        nargs = len(args)
        nvars = self.ctx.nvars()
        # todo: should extend with generators instead?
        if nargs < nvars:
            args = args + (0,) * (nvars - nargs)
        if nargs > nvars:
            args = args[:nvars]
        args_fmpz = [any_as_fmpz(v) for v in args]
        # todo: for combination, compose
        # todo: if fewer than number of variables, evaluate partially?
        if NotImplemented not in args_fmpz:
            V = <fmpz_struct **> libc.stdlib.malloc(nvars * sizeof(fmpz_struct *))
            try:
                for i in range(nvars):
                    V[i] = &((<fmpz> args_fmpz[i]).val[0])
                vres = fmpz.__new__(fmpz)
                if fmpz_mpoly_evaluate_all_fmpz(vres.val, self.val, V, self.ctx.val) == 0:
                    raise ValueError("unreasonably large polynomial")
                return vres
            finally:
                libc.stdlib.free(V)
        else:
            res_ctx, args = coerce_fmpz_mpolys(*args)
            C = <fmpz_mpoly_struct **> libc.stdlib.malloc(nvars * sizeof(fmpz_mpoly_struct *))
            try:
                for i in range(nvars):
                    C[i] = &((<fmpz_mpoly> args[i]).val[0])
                res = fmpz_mpoly.__new__(fmpz_mpoly)
                res.ctx = res_ctx
                fmpz_mpoly_init(res.val, res.ctx.val)
                res._init = True
                if fmpz_mpoly_compose_fmpz_mpoly(res.val, self.val, C, self.ctx.val, res_ctx.val) == 0:
                    raise ValueError("unreasonably large polynomial")
                return res
            finally:
                libc.stdlib.free(C)

    '''
    def factor(self):
        """
        Factors self into irreducible factors, returning a tuple
        (c, factors) where c is the content of the coefficients and
        factors is a list of (poly, exp) pairs.

        """
        cdef fmpz_mpoly_factor_t fac
        cdef int i
        cdef fmpz c
        cdef fmpz_mpoly u
        fmpz_mpoly_factor_init(fac, self.ctx.val)
        fmpz_mpoly_factor(fac, self.val, 1, self.ctx.val)
        res = [0] * fac.length
        for 0 <= i < fac.length:
            u = fmpz_mpoly.__new__(fmpz_mpoly)
            u.ctx = self.ctx
            fmpz_mpoly_init(u.val, u.ctx.val)
            u._init = True
            fmpz_mpoly_set((<fmpz_mpoly>u).val, &fac.poly[i], self.ctx.val)
            c = fmpz.__new__(fmpz)
            fmpz_set((<fmpz>c).val, &fac.exp[i])
            res[i] = (u, c)
        c = fmpz.__new__(fmpz)
        fmpz_set((<fmpz>c).val, fac.content)   # should be & with ...
        fmpz_mpoly_factor_clear(fac, self.ctx.val)
        return c, res
    '''

