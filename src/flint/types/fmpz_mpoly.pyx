from cpython.version cimport PY_MAJOR_VERSION
from cpython.dict cimport PyDict_Size, PyDict_Check, PyDict_Next
from cpython.tuple cimport PyTuple_Check, PyTuple_GET_SIZE
from flint.flintlib.fmpz cimport fmpz_init, fmpz_clear, fmpz_is_zero
from flint.flintlib.fmpq_mpoly cimport fmpq_mpoly_t, fmpq_mpoly_add_fmpq, fmpq_mpoly_sub_fmpq, fmpq_mpoly_scalar_mul_fmpq
from flint.flintlib.fmpq_mpoly cimport fmpq_mpoly_scalar_div_fmpq
from flint.flintlib.flint cimport *
from flint.flintlib.fmpq cimport fmpq_numref, fmpq_denref
from flint.flintlib.fmpz_mpoly_q cimport  fmpz_mpoly_q_div_fmpz

from flint.utils.conversion cimport str_from_chars
from flint.utils.typecheck cimport typecheck
from flint.flint_base.flint_base cimport flint_mpoly
from flint.flint_base.flint_base cimport flint_mpoly_context
from flint.types.fmpz cimport any_as_fmpz
from flint.types.fmpz cimport fmpz, fmpz_set_any_ref
from flint.types.fmpq cimport any_as_fmpq, fmpq_set_any_ref
from flint.types.fmpq_mpoly cimport fmpq_mpoly
from flint.types.fmpz_mpoly_q cimport fmpz_mpoly_q


cimport cython
cimport libc.stdlib
from flint.flintlib.fmpz cimport fmpz_set
from flint.flintlib.fmpz_mpoly cimport *
from flint.flintlib.fmpz_mpoly_factor cimport *

cdef extern from *:
    """
    /* An ugly hack to get around the ugly hack of renaming fmpq to avoid a c/python name collision */
    typedef fmpq fmpq_struct;
    """


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
        if ordering == "lex":
            fmpz_mpoly_ctx_init(self.val, nvars, ordering_t.ORD_LEX)
        elif ordering == "deglex":
            fmpz_mpoly_ctx_init(self.val, nvars, ordering_t.ORD_DEGLEX)
        elif ordering == "degrevlex":
            fmpz_mpoly_ctx_init(self.val, nvars, ordering_t.ORD_DEGREVLEX)
        else:
            raise ValueError("Unimplemented term order %s" % ordering)
        super().__init__(nvars, names)

    cpdef slong nvars(self):
        """
        Return the number of variables in the context

            >>> ctx = fmpz_mpoly_ctx.get_context(4, "lex", 'x')
            >>> ctx.nvars()
            4
        """
        return self.val.minfo.nvars

    cpdef ordering(self):
        """
        Return the term order of the context object.

            >>> ctx = fmpz_mpoly_ctx.get_context(4, "deglex", 'w')
            >>> ctx.ordering()
            'deglex'
        """
        if self.val.minfo.ord == ordering_t.ORD_LEX:
            return "lex"
        if self.val.minfo.ord == ordering_t.ORD_DEGLEX:
            return "deglex"
        if self.val.minfo.ord == ordering_t.ORD_DEGREVLEX:
            return "degrevlex"

    def gen(self, slong i):
        """
        Return the `i`th generator of the polynomial ring

            >>> ctx = fmpz_mpoly_ctx.get_context(3, 'degrevlex', 'z')
            >>> ctx.gen(1)
            z1
        """
        cdef fmpz_mpoly res
        assert i >= 0 and i < self.val.minfo.nvars
        res = fmpz_mpoly.__new__(fmpz_mpoly)
        res.ctx = self
        fmpz_mpoly_init(res.val, res.ctx.val)
        res._init = True
        fmpz_mpoly_gen(res.val, i, res.ctx.val)
        return res

    def constant(self, z):
        """
        Create a constant polynomial in this context
        """
        cdef fmpz_mpoly res
        z = any_as_fmpz(z)
        if z is NotImplemented:
            raise ValueError("A constant fmpz_mpoly is a fmpz")
        res = fmpz_mpoly.__new__(fmpz_mpoly)
        res.ctx = self
        fmpz_mpoly_init(res.val, res.ctx.val)
        res._init = True
        fmpz_mpoly_set_fmpz(res.val, (<fmpz>z).val, res.ctx.val)
        return res

    def fmpz_mpoly_from_dict(self, d):
        """
        Create a fmpz_mpoly from a dictionary.

        The dictionary's keys are tuples of ints (or anything that implicitly converts
        to fmpz) representing exponents, and corresponding values of fmpz.

            >>> ctx = fmpz_mpoly_ctx.get_context(2,'lex','x,y')
            >>> ctx.fmpz_mpoly_from_dict({(1,0):2, (1,1):3, (0,1):1})
            3*x*y + 2*x + y
        """
        cdef long n
        cdef fmpz_t coefficient
        cdef fmpz_struct *exponents
        cdef int xtype
        cdef int nvars = self.nvars()
        cdef int i,j
        cdef int count
        cdef fmpz_mpoly res

        if not PyDict_Check(d):
            raise ValueError("expected a dictionary")
        n = PyDict_Size(d)
        fmpz_init(coefficient)
        exponents = <fmpz_struct *> libc.stdlib.calloc(nvars, sizeof(fmpz_struct))
        if exponents == NULL:
            raise MemoryError()
        for i in range(nvars):
            fmpz_init(exponents + i)
        fmpz_init(coefficient)
        res = fmpz_mpoly.__new__(fmpz_mpoly)
        res.ctx = self
        fmpz_mpoly_init(res.val, res.ctx.val)
        res._init = True
        count = 0
        for k, v in d.items():
            xtype = fmpz_set_any_ref(coefficient, v)
            if xtype == FMPZ_UNKNOWN:
                for i in range(nvars):
                    fmpz_clear(exponents + i)
                libc.stdlib.free(exponents)
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
            if not fmpz_is_zero(coefficient):
                _fmpz_mpoly_push_exp_ffmpz(res.val, exponents, self.val)
                fmpz_mpoly_set_term_coeff_fmpz(res.val, count, coefficient, self.val)
                count += 1
        for i in range(nvars):
            fmpz_clear(exponents + i)
        fmpz_clear(coefficient)
        fmpz_mpoly_sort_terms(res.val, self.val)
        return res


cdef class fmpz_mpoly(flint_mpoly):
    """
    The *fmpz_poly* type represents sparse multivariate polynomials over
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
                raise ValueError("Cannot automatically coerce contexts")
        elif isinstance(val, dict):
            if ctx is None:
                if len(val) == 0:
                    raise ValueError("Need context for zero polynomial")
                k = list(val.keys())[0]
                if not isinstance(k, tuple):
                    raise ValueError("Dict should be keyed with tuples of integers")
                ctx = fmpz_mpoly_ctx.get_context(len(k))
            x = ctx.fmpz_mpoly_from_dict(val)
            #XXX this copy is silly, have a ctx function that assigns an fmpz_mpoly_t
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

    def context(self):
        return self.ctx

    def __nonzero__(self):
        return not fmpz_mpoly_is_zero(self.val, self.ctx.val)

    def __bool__(self):
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
            else:
                if op == 2:
                    return False
                else:
                    return True
        if op == 2:
            return not bool(self - other)
        else:
            return bool(self - other)

    def __len__(self):
        return fmpz_mpoly_length(self.val, self.ctx.val)

    def coefficient(self, slong i):
        cdef fmpz v
        if i < 0 or i >= fmpz_mpoly_length(self.val, self.ctx.val):
            return fmpz(0)
        else:
            v = fmpz.__new__(fmpz)
            fmpz_mpoly_get_term_coeff_fmpz(v.val, self.val, i, self.ctx.val)
            return v

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

    def repr(self):
        return self.str() + "  (nvars=%s, ordering=%s names=%s)" % (self.ctx.nvars(), self.ctx.ordering(), self.ctx.py_names)

    def str(self):
        cdef char * s = fmpz_mpoly_get_str_pretty(self.val, self.ctx.c_names, self.ctx.val)
        try:
            res = str_from_chars(s)
        finally:
            libc.stdlib.free(s)
        res = res.replace("+", " + ")
        res = res.replace("-", " - ")
        if res.startswith(" - "):
            res = "-" + res[3:]
        return res

    def __neg__(self):
        cdef fmpz_mpoly res
        res = create_fmpz_mpoly(self.ctx)
        fmpz_mpoly_neg(res.val, (<fmpz_mpoly>self).val, res.ctx.val)
        return res

    def __add__(self, other):
        cdef fmpz_mpoly res
        cdef fmpq_mpoly qres
        cdef fmpz_t z_other
        cdef fmpq_t q_other
        cdef int xtype
        if typecheck(other, fmpz_mpoly):
            if (<fmpz_mpoly>self).ctx is not (<fmpz_mpoly>other).ctx:
                return NotImplemented
            res = create_fmpz_mpoly(self.ctx)
            fmpz_mpoly_add(res.val, (<fmpz_mpoly>self).val, (<fmpz_mpoly>other).val, res.ctx.val)
            return res
        else:
            xtype = fmpz_set_any_ref(z_other, other)
            if xtype != FMPZ_UNKNOWN:
                res = create_fmpz_mpoly(self.ctx)
                fmpz_mpoly_add_fmpz(res.val, (<fmpz_mpoly>self).val, z_other, res.ctx.val)
                return res
            xtype = fmpq_set_any_ref(q_other, other)
            if xtype != FMPZ_UNKNOWN:
                qres = fmpq_mpoly(self)
                fmpq_mpoly_add_fmpq(qres.val, qres.val, q_other, qres.ctx.val)
                return qres
        return NotImplemented

    def __radd__(self, other):
        cdef fmpz_mpoly res
        cdef fmpq_mpoly qres
        cdef fmpz_t z_other
        cdef fmpq_t q_other
        cdef int xtype

        xtype = fmpz_set_any_ref(z_other, other)
        if xtype != FMPZ_UNKNOWN:
            res = create_fmpz_mpoly(self.ctx)
            fmpz_mpoly_add_fmpz(res.val, (<fmpz_mpoly>self).val, z_other, res.ctx.val)
            return res
        xtype = fmpq_set_any_ref(q_other, other)
        if xtype != FMPZ_UNKNOWN:
            qres = fmpq_mpoly(self)
            fmpq_mpoly_add_fmpq(qres.val, qres.val, q_other, qres.ctx.val)
            return qres
        return NotImplemented

    def __iadd__(self, other):
        if typecheck(other, fmpz_mpoly):
            if (<fmpz_mpoly>self).ctx is not (<fmpz_mpoly>other).ctx:
                return NotImplemented
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
        cdef fmpq_mpoly qres
        cdef fmpz_t z_other
        cdef fmpq_t q_other
        cdef int xtype
        if typecheck(other, fmpz_mpoly):
            if (<fmpz_mpoly>self).ctx is not (<fmpz_mpoly>other).ctx:
                return NotImplemented
            res = create_fmpz_mpoly(self.ctx)
            fmpz_mpoly_sub(res.val, (<fmpz_mpoly>self).val, (<fmpz_mpoly>other).val, res.ctx.val)
            return res
        else:
            xtype = fmpz_set_any_ref(z_other, other)
            if xtype != FMPZ_UNKNOWN:
                res = create_fmpz_mpoly(self.ctx)
                fmpz_mpoly_sub_fmpz(res.val, (<fmpz_mpoly>self).val, z_other, res.ctx.val)
                return res
            xtype = fmpq_set_any_ref(q_other, other)
            if xtype != FMPZ_UNKNOWN:
                qres = fmpq_mpoly(self)
                fmpq_mpoly_sub_fmpq(qres.val, qres.val, q_other, qres.ctx.val)
                return qres
        return NotImplemented

    def __rsub__(self, other):
        cdef fmpz_mpoly res
        cdef fmpq_mpoly qres
        cdef fmpz_t z_other
        cdef fmpq_t q_other
        cdef int xtype

        xtype = fmpz_set_any_ref(z_other, other)
        if xtype != FMPZ_UNKNOWN:
            res = create_fmpz_mpoly(self.ctx)
            fmpz_mpoly_sub_fmpz(res.val, (<fmpz_mpoly>self).val, z_other, res.ctx.val)
            return -res
        xtype = fmpq_set_any_ref(q_other, other)
        if xtype != FMPZ_UNKNOWN:
            qres = fmpq_mpoly(self)
            fmpq_mpoly_sub_fmpq(qres.val, qres.val, q_other, qres.ctx.val)
            return -qres
        return NotImplemented

    def __isub__(self, other):
        if typecheck(other, fmpz_mpoly):
            if (<fmpz_mpoly>self).ctx is not (<fmpz_mpoly>other).ctx:
                return NotImplemented
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
        cdef fmpq_mpoly qres
        cdef fmpz_t z_other
        cdef fmpq_t q_other
        cdef int xtype

        if typecheck(other, fmpz_mpoly):
            if (<fmpz_mpoly>self).ctx is not (<fmpz_mpoly>other).ctx:
                return NotImplemented
            res = create_fmpz_mpoly(self.ctx)
            fmpz_mpoly_mul(res.val, (<fmpz_mpoly>self).val, (<fmpz_mpoly>other).val, res.ctx.val)
            return res
        else:
            xtype = fmpz_set_any_ref(z_other, other)
            if xtype != FMPZ_UNKNOWN:
                res = create_fmpz_mpoly(self.ctx)
                fmpz_mpoly_scalar_mul_fmpz(res.val, (<fmpz_mpoly>self).val, z_other, res.ctx.val)
                return res
            xtype = fmpq_set_any_ref(q_other, other)
            if xtype != FMPZ_UNKNOWN:
                qres = fmpq_mpoly(self)
                fmpq_mpoly_scalar_mul_fmpq(qres.val, qres.val, q_other, qres.ctx.val)
                return qres
        return NotImplemented

    def __rmul__(self, other):
        cdef fmpz_mpoly res
        cdef fmpq_mpoly qres
        cdef fmpz_t z_other
        cdef fmpq_t q_other
        cdef int xtype

        xtype = fmpz_set_any_ref(z_other, other)
        if xtype != FMPZ_UNKNOWN:
            res = create_fmpz_mpoly(self.ctx)
            fmpz_mpoly_scalar_mul_fmpz(res.val, (<fmpz_mpoly>self).val, z_other, res.ctx.val)
            return res
        xtype = fmpq_set_any_ref(q_other, other)
        if xtype != FMPZ_UNKNOWN:
            qres = fmpq_mpoly(self)
            fmpq_mpoly_scalar_mul_fmpq(qres.val, qres.val, q_other, qres.ctx.val)
            return qres
        return NotImplemented

    def __imul__(self, other):
        if typecheck(other, fmpz_mpoly):
            if (<fmpz_mpoly>self).ctx is not (<fmpz_mpoly>other).ctx:
                return NotImplemented
            fmpz_mpoly_mul((<fmpz_mpoly>self).val, (<fmpz_mpoly>self).val, (<fmpz_mpoly>other).val, self.ctx.val)
            return self
        else:
            other = any_as_fmpz(other)
            if other is not NotImplemented:
                fmpz_mpoly_scalar_mul_fmpz(self.val, (<fmpz_mpoly>self).val, (<fmpz>other).val, self.ctx.val)
                return self
            other = any_as_fmpq(other)
            if other is not NotImplemented:
                return fmpq_mpoly(self) * other
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
            raise ValueError("unreasonably large polynomial")
        return res

    def __divmod__(self, other):
        cdef fmpz_mpoly res, res2
        if typecheck(other, fmpz_mpoly):
            if not other:
                raise ZeroDivisionError("fmpz_mpoly_divison by zero")
            if (<fmpz_mpoly>self).ctx is not (<fmpz_mpoly>other).ctx:
                return NotImplemented
            res = create_fmpz_mpoly(self.ctx)
            res2 = create_fmpz_mpoly(self.ctx)
            fmpz_mpoly_divrem(res.val, res2.val, (<fmpz_mpoly>self).val, (<fmpz_mpoly>other).val, res.ctx.val)
            return (res, res2)
        else:
            other = any_as_fmpz(other)
            if other is not NotImplemented:
                other= fmpz_mpoly(other, self.ctx)
                if not other:
                    raise ZeroDivisionError("fmpz_mpoly divison by zero")
                res = create_fmpz_mpoly(self.ctx)
                res2 = create_fmpz_mpoly(self.ctx)
                fmpz_mpoly_divrem(res.val, res2.val, (<fmpz_mpoly>self).val, (<fmpz_mpoly>other).val, res.ctx.val)
                return (res, res2)
        return NotImplemented

    def __rdivmod__(self, other):
        cdef fmpz_mpoly res, res2
        if not self:
            raise ZeroDivisionError("fmpz_mpoly divison by zero")
        other = any_as_fmpz(other)
        if other is not NotImplemented:
            other = fmpz_mpoly(other, self.ctx)
            res = create_fmpz_mpoly(self.ctx)
            res2 = create_fmpz_mpoly(self.ctx)
            fmpz_mpoly_divrem(res.val, res2.val, (<fmpz_mpoly>other).val, (<fmpz_mpoly>self).val, res.ctx.val)
            return res
        return NotImplemented

    def __floordiv__(self, other):
        cdef fmpz_mpoly res
        if typecheck(other, fmpz_mpoly):
            if not other:
                raise ZeroDivisionError("fmpz_mpoly division by zero")
            if (<fmpz_mpoly>self).ctx is not (<fmpz_mpoly>other).ctx:
                return NotImplemented
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

    def __rfloordiv__(self,other):
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
        cdef fmpq_mpoly qres
        cdef fmpq_t q_other
        cdef int xtype

        if typecheck(other, fmpz_mpoly):
            if not other:
                raise ZeroDivisionError("fmpz_mpoly division by zero")
            if self.ctx is not (<fmpz_mpoly> other).ctx:
                return NotImplemented
            return fmpz_mpoly_q(self, other)
        xtype = fmpq_set_any_ref(q_other, other)
        if xtype != FMPZ_UNKNOWN:
            if not other:
                raise ZeroDivisionError("fmpz_mpoly division by zero")
            qres = fmpq_mpoly(self)
            fmpq_mpoly_scalar_div_fmpq(qres.val, qres.val, q_other, qres.ctx.val)
            return qres

    def __rtruediv__(self, other):
        cdef fmpz_mpoly num
        cdef fmpz_mpoly_q ret
        cdef fmpz_t z_other
        cdef fmpq_t q_other
        cdef int xtype

        if not self:
            raise ZeroDivisionError("fmpz_mpoly division by zero")
        if typecheck(other, fmpz_mpoly):
            if self.ctx is not (<fmpz_mpoly> other).ctx:
                return NotImplemented
            return fmpz_mpoly_q(other, self)
        xtype = fmpz_set_any_ref(z_other, other)
        if xtype != FMPZ_UNKNOWN:
            num = create_fmpz_mpoly(self.ctx)
            fmpz_mpoly_set_fmpz(num.val, z_other, self.ctx.val)
            return fmpz_mpoly_q(num, self)
        xtype = fmpq_set_any_ref(q_other, other)
        if xtype != FMPZ_UNKNOWN:
            num = create_fmpz_mpoly(self.ctx)
            fmpz_mpoly_set_fmpz(num.val, fmpq_numref(q_other),self.ctx.val)
            ret = fmpz_mpoly_q(num, self)
            fmpz_mpoly_q_div_fmpz(ret.fraction, ret.fraction, fmpq_denref(q_other), self.ctx.val)
            return ret
        return NotImplemented

    def __mod__(self, other):
        return divmod(self, other)[1]

    def gcd(self, other):
        cdef fmpz_mpoly res
        assert isinstance(other, fmpz_mpoly)
        if (<fmpz_mpoly>self).ctx is not (<fmpz_mpoly>other).ctx:
            return NotImplemented
        res = create_fmpz_mpoly(self.ctx)
        fmpz_mpoly_gcd(res.val, (<fmpz_mpoly>self).val, (<fmpz_mpoly>other).val, res.ctx.val)
        return res

    def __call__(self, *args, **kwargs):
        cdef:
            fmpz_mpoly res
            fmpz_mpoly res2
            fmpz_mpoly_ctx res_ctx
            fmpz_struct ** V
            fmpz vres
            fmpz_mpoly_struct ** C
            slong i, nvars = self.ctx.nvars(), nargs = len(args)

        if args and kwargs:
            raise ValueError("only supply positional or keyword arguments")

        if kwargs:
            # Sort and filter the provided keyword args to only valid ones
            args = tuple((i, kwargs[x]) for i, x in enumerate(self.ctx.names()) if x in kwargs)
            nargs = len(args)

            # If we've been provided with an invalid keyword arg then the length of our filter
            # args will be less than what we've been provided with.
            # If the length is equal to the number of variables then all arguments have been provided
            # otherwise we need to do partial application
            if nargs < len(kwargs):
                raise ValueError("unknown keyword argument provided")
            elif nargs == nvars:
                args = tuple(x for _, x in args)
                partial = False
            else:
                partial = True
        elif nargs > nvars:
            raise ValueError("more arguments provided than variables")
        elif nargs < nvars:
            # Positional partial application
            args = tuple((i, x) for i, x in enumerate(args))
            partial = True

        if partial:
            args_fmpz = tuple(any_as_fmpz(v) for _, v in args)
        else:
            args_fmpz = tuple(any_as_fmpz(v) for v in args)

        all_fmpz = NotImplemented not in args_fmpz
        # if NotImplemented in args_fmpz:
        #     raise NotImplementedError("an argument was unable to be coerced to fmpz")

        if not partial and all_fmpz:
            try:
                V = <fmpz_struct **> libc.stdlib.malloc(nvars * sizeof(fmpz_struct *))
                for i in range(nvars):
                    V[i] = &((<fmpz> args_fmpz[i]).val[0])
                vres = fmpz.__new__(fmpz)
                if fmpz_mpoly_evaluate_all_fmpz(vres.val, self.val, V, self.ctx.val) == 0:
                    raise ValueError("unreasonably large polynomial")
                return vres
            finally:
                libc.stdlib.free(V)
        elif partial and all_fmpz:
            res = fmpz_mpoly.__new__(fmpz_mpoly)
            res2 = fmpz_mpoly.__new__(fmpz_mpoly)
            res.ctx = self.ctx
            res2.ctx = self.ctx

            fmpz_mpoly_set(res2.val, self.val, self.ctx.val)
            for (i, _), arg in zip(args, args_fmpz):
                if fmpz_mpoly_evaluate_one_fmpz(res.val, res2.val, i, (<fmpz>arg).val, self.ctx.val) == 0:
                    raise ValueError("unreasonably large polynomial")
                fmpz_mpoly_set(res2.val, res.val, self.ctx.val)
            return res
        elif not partial and not all_fmpz:
            res_ctx, args = coerce_fmpz_mpolys(args)
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
        else:
            raise NotImplemented("partial composition not implemented")
            # res_ctx, new_args = coerce_fmpz_mpolys(tuple(x for _, x in args))

            # polys = [None] * nvars
            # for (i, _), poly in zip(args, new_args):
            #     polys[i] = poly

            # for i in range(nvars):
            #     if polys[i] is None:
            #         polys[i] = fmpz_mpoly({ : })

    def factor(self):
        """
        Factors self into irreducible factors, returning a tuple
        (c, factors) where c is the content of the coefficients and
        factors is a list of (poly, exp) pairs.

            >>> Zm = fmpz_mpoly
            >>> ctx = fmpz_mpoly_ctx.get_context(3, 'lex', 'x,y,z')
            >>> p1 = Zm("2*x + 4", ctx)
            >>> p2 = Zm("3*x*z +  + 3*x + 3*z + 3", ctx)
            >>> (p1 * p2).factor()
            (6, [(z + 1, 1), (x + 2, 1), (x + 1, 1)])
            >>> (p2 * p1 * p2).factor()
            (18, [(z + 1, 2), (x + 2, 1), (x + 1, 2)])
        """
        cdef fmpz_mpoly_factor_t fac
        cdef int i
        cdef fmpz c
        cdef fmpz_mpoly u
        fmpz_mpoly_factor_init(fac, self.ctx.val)
        fmpz_mpoly_factor(fac, self.val, self.ctx.val)
        res = [0] * fac.num
        for 0 <= i < fac.num:
            u = fmpz_mpoly.__new__(fmpz_mpoly)
            u.ctx = self.ctx
            fmpz_mpoly_init(u.val, u.ctx.val)
            u._init = True
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

            >>> Zm = fmpz_mpoly
            >>> ctx = fmpz_mpoly_ctx.get_context(3, 'lex', 'x,y,z')
            >>> p1 = Zm("2*x + 4", ctx)
            >>> p2 = Zm("3*x*y + 3*x + 3*y + 3", ctx)
            >>> (p1 * p2).factor_squarefree()
            (6, [(y + 1, 1), (x^2 + 3*x + 2, 1)])
            >>> (p1 * p2 * p1).factor_squarefree()
            (12, [(y + 1, 1), (x + 1, 1), (x + 2, 2)])
        """
        cdef fmpz_mpoly_factor_t fac
        cdef int i
        cdef fmpz c
        cdef fmpz_mpoly u
        fmpz_mpoly_factor_init(fac, self.ctx.val)
        fmpz_mpoly_factor_squarefree(fac, self.val, self.ctx.val)
        res = [0] * fac.num
        for 0 <= i < fac.num:
            u = fmpz_mpoly.__new__(fmpz_mpoly)
            u.ctx = self.ctx
            fmpz_mpoly_init(u.val, u.ctx.val)
            u._init = True
            fmpz_mpoly_set((<fmpz_mpoly>u).val, &fac.poly[i], self.ctx.val)
            c = fmpz.__new__(fmpz)
            fmpz_set((<fmpz>c).val, &fac.exp[i])
            res[i] = (u, c)
        c = fmpz.__new__(fmpz)
        fmpz_set((<fmpz>c).val, fac.constant)
        fmpz_mpoly_factor_clear(fac, self.ctx.val)
        return c, res

    def coerce_to_context(self, ctx):
        cdef:
            fmpz_mpoly outpoly
            slong *C
            slong i

        if not typecheck(ctx, fmpz_mpoly_ctx):
            raise ValueError("provided context is not a fmpz_mpoly_ctx")

        if self.ctx is ctx:
            return self

        C = <slong *> libc.stdlib.malloc(self.ctx.val.minfo.nvars * sizeof(slong *))
        outpoly = fmpz_mpoly.__new__(fmpz_mpoly)
        init_fmpz_mpoly(outpoly, ctx)

        vars = {x: i for i, x in enumerate(ctx.py_names)}
        for i, var in enumerate(self.ctx.py_names):
            C[i] = <slong>vars[var]

        fmpz_mpoly_compose_fmpz_mpoly_gen(outpoly.val, self.val, C, self.ctx.val, (<fmpz_mpoly_ctx>ctx).val)

        libc.stdlib.free(C)
        return outpoly
