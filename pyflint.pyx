"""
Python wrapper for FLINT - Fast Library for Number Theory
http://www.flintlib.org/
"""

cimport flint
cimport stdlib

cdef flint_rand_t global_random_state
flint_randinit(global_random_state)

cdef extern from "Python.h":
    ctypedef void PyObject
    ctypedef void PyTypeObject
    ctypedef long Py_ssize_t
    int PyObject_TypeCheck(object, PyTypeObject*)
    int PyInt_Check(PyObject *o)
    PyObject* PyInt_FromLong(long ival)
    int PyLong_Check(PyObject *o)
    long PyInt_AS_LONG(PyObject *io)
    Py_ssize_t PyList_GET_SIZE(PyObject *list)

cdef inline bint typecheck(object ob, object tp):
    return PyObject_TypeCheck(ob, <PyTypeObject*>tp)

cdef inline int fmpz_set_python(fmpz_t x, obj):
    if PyInt_Check(<PyObject*>obj):
        fmpz_set_si(x, PyInt_AS_LONG(<PyObject*>obj))
        return 1
    if PyLong_Check(<PyObject*>obj):
        s = "%x" % obj
        fmpz_set_str(x, s, 16)      # XXX: slow
        return 1
    return 0

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
            raise ValueError("unsupported coefficient in list")
    fmpz_clear(x)

cdef inline any_as_fmpz_poly(x):
    """
    Tries to convert x to an fmpz_poly
    """
    cdef fmpz_poly res
    if typecheck(x, fmpz_poly):
        return x
    elif typecheck(x, fmpz):
        res = fmpz_poly.__new__(fmpz_poly)
        fmpz_poly_set_fmpz(res.val, (<fmpz>x).val)
        return res
    elif PyInt_Check(<PyObject*>x):
        res = fmpz_poly.__new__(fmpz_poly)
        fmpz_poly_set_si(res.val, PyInt_AS_LONG(<PyObject*>x))
        return res
    elif PyLong_Check(<PyObject*>x):
        res = fmpz_poly.__new__(fmpz_poly)
        t = fmpz(x)   # XXX: slow
        fmpz_poly_set_fmpz(res.val, (<fmpz>t).val)
        return res
    return NotImplemented

DEF FMPZ_UNKNOWN = 0
DEF FMPZ_REF = 1
DEF FMPZ_TMP = 2

cdef inline int fmpz_set_any_ref(fmpz_struct *x, obj):
    if typecheck(obj, fmpz):
        x[0] = (<fmpz>obj).val[0]
        return FMPZ_REF
    if PyInt_Check(<PyObject*>obj):
        fmpz_init(x)
        fmpz_set_si(x, PyInt_AS_LONG(<PyObject*>obj))
        return FMPZ_TMP
    if PyLong_Check(<PyObject*>obj):
        fmpz_init(x)
        s = "%x" % obj             # XXX: slow
        fmpz_set_str(x, s, 16)
        return FMPZ_TMP
    return FMPZ_UNKNOWN

cdef inline any_as_fmpz(obj):
    cdef fmpz_struct x[1]
    cdef bint xtype
    cdef fmpz v
    xtype = fmpz_set_any_ref(x, obj)
    if xtype == FMPZ_REF:
        v = fmpz.__new__(fmpz)
        fmpz_set(v.val, x)
        return v
    elif xtype == FMPZ_TMP:
        v = fmpz.__new__(fmpz)
        fmpz_clear(v.val)
        v.val[0] = x[0]
        return v
    else:
        return NotImplemented

cdef any_as_fmpz_mat(obj):
    if typecheck(obj, fmpz_mat):
        return obj
    return NotImplemented


cdef fmpz_get_intlong(fmpz_t x):
    """
    Convert fmpz_t to a Python int or long.
    """
    cdef char * s
    if COEFF_IS_MPZ(x[0]):
        s = fmpz_get_str(NULL, 16, x)   # XXX: slow
        v = int(s, 16)
        stdlib.free(s)
        return v
    else:
        return <long>x[0]

cdef class fmpz:

    cdef fmpz_t val

    def __cinit__(self):
        fmpz_init(self.val)

    def __dealloc__(self):
        fmpz_clear(self.val)

    def __init__(self, val=None):
        cdef long x
        if val is not None:
            if typecheck(val, fmpz):
                fmpz_set(self.val, (<fmpz>val).val)
            else:
                fmpz_set_any_ref(self.val, val)   # XXX


    # XXX: improve!
    def __int__(self):
        return fmpz_get_intlong(self.val)

    def __long__(self):
        return long(fmpz_get_intlong(self.val))

    def __index__(self):
        return fmpz_get_intlong(self.val)

    def __richcmp__(s, t, int op):
        cdef bint res = 0
        cdef long tl
        cdef fmpz_struct tval[1]
        cdef fmpz_struct *sval = (<fmpz>s).val
        cdef int ttype
        if PyInt_Check(<PyObject*>t):
            tl = PyInt_AS_LONG(<PyObject*>t)
            if   op == 2: res = fmpz_cmp_si(sval, tl) == 0
            elif op == 3: res = fmpz_cmp_si(sval, tl) != 0
            elif op == 0: res = fmpz_cmp_si(sval, tl) < 0
            elif op == 1: res = fmpz_cmp_si(sval, tl) <= 0
            elif op == 4: res = fmpz_cmp_si(sval, tl) > 0
            elif op == 5: res = fmpz_cmp_si(sval, tl) >= 0
        else:
            ttype = fmpz_set_any_ref(tval, t)
            if ttype != FMPZ_UNKNOWN:
                if   op == 2: res = fmpz_equal(sval, tval)
                elif op == 3: res = not fmpz_equal(sval, tval)
                elif op == 0: res = fmpz_cmp(sval, tval) < 0
                elif op == 1: res = fmpz_cmp(sval, tval) <= 0
                elif op == 4: res = fmpz_cmp(sval, tval) > 0
                elif op == 5: res = fmpz_cmp(sval, tval) >= 0
            if ttype == FMPZ_TMP:
                fmpz_clear(tval)
            if ttype == FMPZ_UNKNOWN:
                return NotImplemented
        return res

    def __str__(self):
        cdef char * s = fmpz_get_str(NULL, 10, self.val)
        try:
            res = s
        finally:
            stdlib.free(s)
        return res

    def __repr__(self):
        return "fmpz(%s)" % self.__str__()

    def __pos__(self):
        return self

    def __neg__(self):
        cdef fmpz res = fmpz.__new__(fmpz)
        fmpz_neg(res.val, self.val)
        return res

    def __abs__(self):
        if fmpz_sgn(self.val) >= 0:
            return self
        return -self

    def __add__(s, t):
        cdef fmpz_struct sval[1]
        cdef fmpz_struct tval[1]
        cdef int stype = FMPZ_UNKNOWN
        cdef int ttype = FMPZ_UNKNOWN
        u = NotImplemented
        stype = fmpz_set_any_ref(sval, s)
        if stype != FMPZ_UNKNOWN:
            ttype = fmpz_set_any_ref(tval, t)
            if ttype != FMPZ_UNKNOWN:
                u = fmpz.__new__(fmpz)
                fmpz_add((<fmpz>u).val, sval, tval)
        if stype == FMPZ_TMP: fmpz_clear(sval)
        if ttype == FMPZ_TMP: fmpz_clear(tval)
        return u

    def __sub__(s, t):
        cdef fmpz_struct sval[1]
        cdef fmpz_struct tval[1]
        cdef int stype = FMPZ_UNKNOWN
        cdef int ttype = FMPZ_UNKNOWN
        u = NotImplemented
        stype = fmpz_set_any_ref(sval, s)
        if stype != FMPZ_UNKNOWN:
            ttype = fmpz_set_any_ref(tval, t)
            if ttype != FMPZ_UNKNOWN:
                u = fmpz.__new__(fmpz)
                fmpz_sub((<fmpz>u).val, sval, tval)
        if stype == FMPZ_TMP: fmpz_clear(sval)
        if ttype == FMPZ_TMP: fmpz_clear(tval)
        return u

    def __mul__(s, t):
        cdef fmpz_struct sval[1]
        cdef fmpz_struct tval[1]
        cdef int stype = FMPZ_UNKNOWN
        cdef int ttype = FMPZ_UNKNOWN
        u = NotImplemented
        stype = fmpz_set_any_ref(sval, s)
        if stype != FMPZ_UNKNOWN:
            ttype = fmpz_set_any_ref(tval, t)
            if ttype != FMPZ_UNKNOWN:
                u = fmpz.__new__(fmpz)
                fmpz_mul((<fmpz>u).val, sval, tval)
        if stype == FMPZ_TMP: fmpz_clear(sval)
        if ttype == FMPZ_TMP: fmpz_clear(tval)
        return u

    def __floordiv__(s, t):
        cdef fmpz_struct sval[1]
        cdef fmpz_struct tval[1]
        cdef int stype = FMPZ_UNKNOWN
        cdef int ttype = FMPZ_UNKNOWN
        u = NotImplemented
        ttype = fmpz_set_any_ref(tval, t)
        if ttype != FMPZ_UNKNOWN:
            if fmpz_is_zero(tval):
                if ttype == FMPZ_TMP:
                    fmpz_clear(tval)
                raise ZeroDivisionError("fmpz division by zero")
            stype = fmpz_set_any_ref(sval, s)
            if stype != FMPZ_UNKNOWN:
                u = fmpz.__new__(fmpz)
                fmpz_fdiv_q((<fmpz>u).val, sval, tval)
        if stype == FMPZ_TMP: fmpz_clear(sval)
        if ttype == FMPZ_TMP: fmpz_clear(tval)
        return u

    def __mod__(s, t):
        cdef fmpz_struct sval[1]
        cdef fmpz_struct tval[1]
        cdef int stype = FMPZ_UNKNOWN
        cdef int ttype = FMPZ_UNKNOWN
        u = NotImplemented
        ttype = fmpz_set_any_ref(tval, t)
        if ttype != FMPZ_UNKNOWN:
            if fmpz_is_zero(tval):
                if ttype == FMPZ_TMP:
                    fmpz_clear(tval)
                raise ZeroDivisionError("fmpz division by zero")
            stype = fmpz_set_any_ref(sval, s)
            if stype != FMPZ_UNKNOWN:
                u = fmpz.__new__(fmpz)
                fmpz_fdiv_r((<fmpz>u).val, sval, tval)
        if stype == FMPZ_TMP: fmpz_clear(sval)
        if ttype == FMPZ_TMP: fmpz_clear(tval)
        return u

    def __divmod__(s, t):
        cdef fmpz_struct sval[1]
        cdef fmpz_struct tval[1]
        cdef int stype = FMPZ_UNKNOWN
        cdef int ttype = FMPZ_UNKNOWN
        u = NotImplemented
        ttype = fmpz_set_any_ref(tval, t)
        if ttype != FMPZ_UNKNOWN:
            if fmpz_is_zero(tval):
                if ttype == FMPZ_TMP:
                    fmpz_clear(tval)
                raise ZeroDivisionError("fmpz division by zero")
            stype = fmpz_set_any_ref(sval, s)
            if stype != FMPZ_UNKNOWN:
                u1 = fmpz.__new__(fmpz)
                u2 = fmpz.__new__(fmpz)
                fmpz_fdiv_qr((<fmpz>u1).val, (<fmpz>u2).val, sval, tval)
                u = u1, u2
        if stype == FMPZ_TMP: fmpz_clear(sval)
        if ttype == FMPZ_TMP: fmpz_clear(tval)
        return u

    def __pow__(s, t, m):
        cdef fmpz_struct sval[1]
        cdef int stype = FMPZ_UNKNOWN
        cdef ulong exp
        u = NotImplemented
        if m is not None:
            raise NotImplementedError("modular exponentiation")
        stype = fmpz_set_any_ref(sval, s)
        if stype != FMPZ_UNKNOWN:
            c = t
            u = fmpz.__new__(fmpz)
            fmpz_pow_ui((<fmpz>u).val, sval, c)
        if stype == FMPZ_TMP: fmpz_clear(sval)
        return u



cdef class fmpz_poly:

    cdef fmpz_poly_t val

    def __cinit__(self):
        fmpz_poly_init(self.val)

    def __dealloc__(self):
        fmpz_poly_clear(self.val)

    def __init__(self, val=None):
        if val is not None:
            if typecheck(val, fmpz_poly):
                fmpz_poly_set(self.val, (<fmpz_poly>val).val)
            elif isinstance(val, list):
                fmpz_poly_set_list(self.val, val)
            else:
                raise TypeError("cannot create fmpz_poly from input of type %s", type(val))

    cpdef long length(self):
        return fmpz_poly_length(self.val)

    cpdef long degree(self):
        return fmpz_poly_degree(self.val)

    def __richcmp__(self, other, int op):
        cdef bint r
        if op != 2 and op != 3:
            raise TypeError("polynomials cannot be ordered")
        self = any_as_fmpz_poly(self)
        if self is NotImplemented:
            return self
        other = any_as_fmpz_poly(other)
        if other is NotImplemented:
            return other
        r = fmpz_poly_equal((<fmpz_poly>self).val, (<fmpz_poly>other).val)
        if op == 3:
            r = not r
        return r

    def coeffs(self):
        cdef long i, n
        n = self.length()
        L = [fmpz() for i in range(n)]
        for i from 0 <= i < n:
            fmpz_poly_get_coeff_fmpz((<fmpz>(L[i])).val, self.val, i)
        return L

    def __getitem__(self, long i):
        cdef fmpz x
        x = fmpz()
        if i < 0:
            return x
        fmpz_poly_get_coeff_fmpz(x.val, self.val, i)
        return x

    def __setitem__(self, long i, x):
        if i < 0:
            raise ValueError("cannot assign to index < 0 of polynomial")
        v = fmpz(x)  # XXX
        fmpz_poly_set_coeff_fmpz(self.val, i, (<fmpz>v).val)

    def __str__(self):
        cdef char * s = fmpz_poly_get_str_pretty(self.val, "x")
        try:
            res = s
        finally:
            stdlib.free(s)
        return res

    def __repr__(self):
        return "fmpz_poly(%s)" % map(int, self.coeffs())

    def __pos__(self):
        return self

    def __neg__(self):
        cdef fmpz_poly res = fmpz_poly.__new__(fmpz_poly)
        fmpz_poly_neg(res.val, self.val)
        return res

    def __add__(self, other):
        cdef fmpz_poly res
        self = any_as_fmpz_poly(self)
        if self is NotImplemented:
            return self
        other = any_as_fmpz_poly(other)
        if other is NotImplemented:
            return other
        res = fmpz_poly.__new__(fmpz_poly)
        fmpz_poly_add(res.val, (<fmpz_poly>self).val, (<fmpz_poly>other).val)
        return res

    def __sub__(self, other):
        cdef fmpz_poly res
        self = any_as_fmpz_poly(self)
        if self is NotImplemented:
            return self
        other = any_as_fmpz_poly(other)
        if other is NotImplemented:
            return other
        res = fmpz_poly.__new__(fmpz_poly)
        fmpz_poly_sub(res.val, (<fmpz_poly>self).val, (<fmpz_poly>other).val)
        return res

    def __mul__(self, other):
        cdef fmpz_poly res
        self = any_as_fmpz_poly(self)
        if self is NotImplemented:
            return self
        other = any_as_fmpz_poly(other)
        if other is NotImplemented:
            return other
        res = fmpz_poly.__new__(fmpz_poly)
        fmpz_poly_mul(res.val, (<fmpz_poly>self).val, (<fmpz_poly>other).val)
        return res

    def __floordiv__(self, other):
        cdef fmpz_poly res
        self = any_as_fmpz_poly(self)
        if self is NotImplemented:
            return self
        other = any_as_fmpz_poly(other)
        if other is NotImplemented:
            return other
        if fmpz_poly_is_zero((<fmpz_poly>other).val):
            raise ZeroDivisionError("fmpz_poly division by 0")
        res = fmpz_poly.__new__(fmpz_poly)
        fmpz_poly_div(res.val, (<fmpz_poly>self).val, (<fmpz_poly>other).val)
        return res

    def __mod__(self, other):
        cdef fmpz_poly res
        self = any_as_fmpz_poly(self)
        if self is NotImplemented:
            return self
        other = any_as_fmpz_poly(other)
        if other is NotImplemented:
            return other
        if fmpz_poly_is_zero((<fmpz_poly>other).val):
            raise ZeroDivisionError("fmpz_poly division by 0")
        res = fmpz_poly.__new__(fmpz_poly)
        fmpz_poly_rem(res.val, (<fmpz_poly>self).val, (<fmpz_poly>other).val)
        return res

    def __divmod__(self, other):
        cdef fmpz_poly P, Q
        self = any_as_fmpz_poly(self)
        if self is NotImplemented:
            return self
        other = any_as_fmpz_poly(other)
        if other is NotImplemented:
            return other
        if fmpz_poly_is_zero((<fmpz_poly>other).val):
            raise ZeroDivisionError("fmpz_poly divmod by 0")
        P = fmpz_poly.__new__(fmpz_poly)
        Q = fmpz_poly.__new__(fmpz_poly)
        fmpz_poly_divrem(P.val, Q.val, (<fmpz_poly>self).val, (<fmpz_poly>other).val)
        return P, Q

    def __pow__(fmpz_poly self, ulong exp, mod):
        cdef fmpz_poly res
        if mod is not None:
            raise NotImplementedError("fmpz_poly modular exponentiation")
        res = fmpz_poly.__new__(fmpz_poly)
        fmpz_poly_pow(res.val, self.val, exp)
        return res


cdef class fmpz_mat:

    cdef fmpz_mat_t val

    def __cinit__(self):
        fmpz_mat_init(self.val, 0, 0)

    def __dealloc__(self):
        fmpz_mat_clear(self.val)

    def __init__(self, *args):
        cdef long m, n, i, j
        if len(args) == 1:
            val = args[0]
            if typecheck(val, fmpz_mat):
                fmpz_mat_init_set(self.val, (<fmpz_mat>val).val)
            else:
                raise TypeError("cannot create fmpz_mat from input of type %s" % type(val))
        elif len(args) == 2:
            m, n = args
            fmpz_mat_init(self.val, m, n)
        elif len(args) == 3:
            m, n, entries = args
            fmpz_mat_init(self.val, m, n)
            entries = list(entries)
            if len(entries) != m*n:
                raise ValueError("list of entries has the wrong length")
            for i from 0 <= i < m:
                for j from 0 <= j < n:
                    # XXX: slow
                    x = fmpz(entries[i*n + j])
                    fmpz_set(fmpz_mat_entry(self.val, i, j), (<fmpz>x).val)
        else:
            raise ValueError("fmpz_mat: expected 1-3 arguments")

    def __richcmp__(s, t, int op):
        cdef bint r
        if op != 2 and op != 3:
            raise TypeError("matrices cannot be ordered")
        s = any_as_fmpz_mat(s)
        if t is NotImplemented:
            return s
        t = any_as_fmpz_mat(t)
        if t is NotImplemented:
            return t
        r = fmpz_mat_equal((<fmpz_mat>s).val, (<fmpz_mat>t).val)
        if op == 3:
            r = not r
        return r

    cpdef long nrows(self):
        return fmpz_mat_nrows(self.val)

    cpdef long ncols(self):
        return fmpz_mat_ncols(self.val)

    def __repr__(self):
        return "fmpz_mat(%i, %i, [%s])" % (self.nrows(), self.ncols(),
            (", ".join(map(str, self.entries()))))

    def __str__(self):
        tab = self.table()
        if len(tab) == 0 or len(tab[0]) == 0:
            return "[]"
        tab = [map(str, row) for row in tab]
        widths = []
        for i in xrange(len(tab[0])):
            w = max([len(row[i]) for row in tab])
            widths.append(w)
        for i in xrange(len(tab)):
            tab[i] = [s.rjust(widths[j]) for j, s in enumerate(tab[i])]
            tab[i] = "[" + (", ".join(tab[i])) + "]"
        return "\n".join(tab)

    def __getitem__(self, index):
        cdef long i, j
        cdef fmpz x
        i, j = index
        if i < 0 or i >= self.nrows() or j < 0 or j >= self.ncols():
            raise ValueError("index %i,%i exceeds matrix dimensions" % (i, j))
        x = fmpz.__new__(fmpz)
        fmpz_set(x.val, fmpz_mat_entry(self.val, i, j))
        return x

    def __setitem__(self, index, value):
        cdef long i, j
        i, j = index
        if i < 0 or i >= self.nrows() or j < 0 or j >= self.ncols():
            raise ValueError("index %i,%i exceeds matrix dimensions" % (i, j))
        c = fmpz(value)  # XXX
        fmpz_set(fmpz_mat_entry(self.val, i, j), (<fmpz>c).val)

    def entries(self):
        cdef long i, j, m, n
        cdef fmpz t
        m = self.nrows()
        n = self.ncols()
        L = [None] * (m * n)
        for i from 0 <= i < m:
            for j from 0 <= j < n:
                t = fmpz.__new__(fmpz)
                fmpz_set(t.val, fmpz_mat_entry(self.val, i, j))
                L[i*n + j] = t
        return L

    def table(self):
        cdef long i, m, n
        m = self.nrows()
        n = self.ncols()
        L = self.entries()
        return [L[i*n:(i+1)*n] for i in range(m)]

    def det(self):
        cdef fmpz d
        if not fmpz_mat_is_square(self.val):
            raise ValueError("matrix must be square")
        d = fmpz.__new__(fmpz)
        fmpz_mat_det(d.val, self.val)
        return d

    def __add__(s, t):
        cdef fmpz_mat u
        cdef fmpz_mat_struct *sval, *tval
        sm = any_as_fmpz_mat(s)
        if sm is NotImplemented:
            return sm
        tm = any_as_fmpz_mat(t)
        if tm is NotImplemented:
            return tm
        sval = (<fmpz_mat>sm).val
        tval = (<fmpz_mat>tm).val
        if (fmpz_mat_nrows(sval) != fmpz_mat_nrows(tval) or
           fmpz_mat_ncols(sval) != fmpz_mat_ncols(tval)):
            raise ValueError("incompatible shapes for matrix addition")
        u = fmpz_mat.__new__(fmpz_mat)
        fmpz_mat_init(u.val, fmpz_mat_nrows(sval), fmpz_mat_ncols(sval))
        fmpz_mat_add(u.val, sval, tval)
        return u

    def __sub__(s, t):
        cdef fmpz_mat u
        cdef fmpz_mat_struct *sval, *tval
        sm = any_as_fmpz_mat(s)
        if sm is NotImplemented:
            return sm
        tm = any_as_fmpz_mat(t)
        if tm is NotImplemented:
            return tm
        sval = (<fmpz_mat>sm).val
        tval = (<fmpz_mat>tm).val
        if (fmpz_mat_nrows(sval) != fmpz_mat_nrows(tval) or
           fmpz_mat_ncols(sval) != fmpz_mat_ncols(tval)):
            raise ValueError("incompatible shapes for matrix subtraction")
        u = fmpz_mat.__new__(fmpz_mat)
        fmpz_mat_init(u.val, fmpz_mat_nrows(sval), fmpz_mat_ncols(sval))
        fmpz_mat_sub(u.val, sval, tval)
        return u

    cdef __mul_fmpz(self, fmpz c):
        cdef fmpz_mat u
        u = fmpz_mat.__new__(fmpz_mat)
        fmpz_mat_init(u.val, fmpz_mat_nrows(self.val), fmpz_mat_ncols(self.val))
        fmpz_mat_scalar_mul_fmpz(u.val, self.val, c.val)
        return u

    def __mul__(s, t):
        cdef fmpz_mat u
        cdef fmpz_mat_struct *sval, *tval
        cdef int ttype
        if typecheck(s, fmpz_mat) and typecheck(t, fmpz_mat):
            sval = (<fmpz_mat>s).val
            tval = (<fmpz_mat>t).val
            if fmpz_mat_ncols(sval) != fmpz_mat_nrows(tval):
                raise ValueError("incompatible shapes for matrix multiplication")
            u = fmpz_mat.__new__(fmpz_mat)
            fmpz_mat_init(u.val, fmpz_mat_nrows(sval), fmpz_mat_ncols(tval))
            fmpz_mat_mul(u.val, sval, tval)
            return u
        else:
            if typecheck(t, fmpz_mat):
                s, t = t, s
            c = any_as_fmpz(t)
            if c is not NotImplemented:
                return (<fmpz_mat>s).__mul_fmpz(c)
        return NotImplemented

    @classmethod
    def randtest(cls, m, n, bits):
        cdef fmpz_mat mat = fmpz_mat(m, n)
        fmpz_mat_randtest(mat.val, global_random_state, bits)
        return mat

    @classmethod
    def randbits(cls, m, n, bits):
        cdef fmpz_mat mat = fmpz_mat(m, n)
        fmpz_mat_randbits(mat.val, global_random_state, bits)
        return mat

    @classmethod
    def randrank(cls, m, n, rank, bits):
        cdef fmpz_mat mat = fmpz_mat(m, n)
        fmpz_mat_randrank(mat.val, global_random_state, rank, bits)
        return mat

    def rank(self):
        return fmpz_mat_rank(self.val)

