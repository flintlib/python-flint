"""
Python wrapper for FLINT - Fast Library for Number Theory
http://www.flintlib.org/
"""

cimport flint
cimport libc.stdlib
cimport cython

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


#----------------------------------------------------------------------------#
#                                                                            #
#                        Various utilities                                   #
#                                                                            #
#----------------------------------------------------------------------------#

def matrix_to_str(tab):
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

cdef fmpz_get_intlong(fmpz_t x):
    """
    Convert fmpz_t to a Python int or long.
    """
    cdef char * s
    if COEFF_IS_MPZ(x[0]):
        s = fmpz_get_str(NULL, 16, x)   # XXX: slow
        v = int(s, 16)
        libc.stdlib.free(s)
        return v
    else:
        return <long>x[0]

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

cdef inline any_as_fmpz_poly(x):
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
            raise TypeError("unsupported coefficient in list")
    fmpz_clear(x)

cdef inline any_as_fmpz_mat(obj):
    if typecheck(obj, fmpz_mat):
        return obj
    return NotImplemented

cdef inline any_as_fmpq(obj):
    if typecheck(obj, fmpq):
        return obj
    z = any_as_fmpz(obj)
    if z is NotImplemented:
        return z
    q = fmpq.__new__(fmpq)
    fmpz_set(fmpq_numref((<fmpq>q).val), (<fmpz>z).val)
    fmpz_set_ui(fmpq_denref((<fmpq>q).val), 1)
    return q

cdef inline any_as_fmpq_poly(obj):
    if typecheck(obj, fmpq_poly):
        return obj
    x = any_as_fmpz(obj)
    if x is not NotImplemented:
        r = fmpq_poly.__new__(fmpq_poly)
        fmpq_poly_set_fmpz((<fmpq_poly>r).val, (<fmpz>x).val)
        return r
    x = any_as_fmpz_poly(obj)
    if x is not NotImplemented:
        r = fmpq_poly.__new__(fmpq_poly)
        fmpq_poly_set_fmpz_poly((<fmpq_poly>r).val, (<fmpz_poly>x).val)
        return r
    x = any_as_fmpq(obj)
    if x is not NotImplemented:
        r = fmpq_poly.__new__(fmpq_poly)
        fmpq_poly_set_fmpq((<fmpq_poly>r).val, (<fmpq>x).val)
        return r
    return NotImplemented

cdef fmpq_poly_set_list(fmpq_poly_t poly, list val):
    cdef long i, n
    n = PyList_GET_SIZE(<PyObject*>val)
    fmpq_poly_fit_length(poly, n)
    for i from 0 <= i < n:
        c = val[i]
        x = any_as_fmpz(c)
        if x is not NotImplemented:
            fmpq_poly_set_coeff_fmpz(poly, i, (<fmpz>x).val)
            continue
        x = any_as_fmpq(c)
        if x is not NotImplemented:
            fmpq_poly_set_coeff_fmpq(poly, i, (<fmpq>x).val)
            continue
        raise TypeError("unsupported coefficient in list")

cdef inline any_as_fmpq_mat(obj):
    if typecheck(obj, fmpq_mat):
        return obj
    if typecheck(obj, fmpz_mat):
        return fmpq_mat(obj)
    return NotImplemented

cdef int any_as_nmod(mp_limb_t * val, obj, nmod_t mod) except -1:
    cdef int success
    cdef fmpz_t t
    if typecheck(obj, nmod):
        if (<nmod>obj).mod.n != mod.n:
            raise ValueError("cannot coerce integers mod n with different n")
        val[0] = (<nmod>obj).val
        return 1
    z = any_as_fmpz(obj)
    if z is not NotImplemented:
        val[0] = fmpz_fdiv_ui((<fmpz>z).val, mod.n)
        return 1
    q = any_as_fmpq(obj)
    if q is not NotImplemented:
        # XXX: write an fmpq_mod_ui or fmpq_get_nmod in flint
        fmpz_init(t)
        fmpz_set_ui(t, mod.n)
        success = fmpq_mod_fmpz(t, (<fmpq>q).val, t)
        val[0] = fmpz_get_ui(t)
        fmpz_clear(t)
        if not success:
            raise ZeroDivisionError("%s does not exist mod %i!" % (q, mod.n))
        return 1
    return 0

cdef any_as_nmod_poly(obj, nmod_t mod):
    cdef nmod_poly r
    cdef mp_limb_t v
    # XXX: should check that modulus is the same here, and not all over the place
    if typecheck(obj, nmod_poly):
        return obj
    if any_as_nmod(&v, obj, mod):
        r = nmod_poly.__new__(nmod_poly)
        nmod_poly_init(r.val, mod.n)
        nmod_poly_set_coeff_ui(r.val, 0, v)
        return r
    x = any_as_fmpz_poly(obj)
    if x is not NotImplemented:
        r = nmod_poly.__new__(nmod_poly)
        nmod_poly_init(r.val, mod.n)   # XXX: create flint _nmod_poly_set_modulus for this?
        fmpz_poly_get_nmod_poly(r.val, (<fmpz_poly>x).val)
        return r
    return NotImplemented

cdef nmod_poly_set_list(nmod_poly_t poly, list val):
    cdef long i, n
    cdef nmod_t mod
    cdef mp_limb_t v
    nmod_init(&mod, nmod_poly_modulus(poly)) # XXX
    n = PyList_GET_SIZE(<PyObject*>val)
    nmod_poly_fit_length(poly, n)
    for i from 0 <= i < n:
        c = val[i]
        if any_as_nmod(&v, val[i], mod):
            nmod_poly_set_coeff_ui(poly, i, v)
        else:
            raise TypeError("unsupported coefficient in list")

cdef any_as_nmod_mat(obj, nmod_t mod):
    cdef nmod_mat r
    cdef mp_limb_t v
    if typecheck(obj, nmod_mat):
        return obj
    x = any_as_fmpz_mat(obj)
    if x is not NotImplemented:
        r = nmod_mat.__new__(nmod_mat)
        nmod_mat_init(r.val, fmpz_mat_nrows((<fmpz_mat>x).val),
                             fmpz_mat_ncols((<fmpz_mat>x).val), mod.n)
        fmpz_mat_get_nmod_mat(r.val, (<fmpz_mat>x).val)
        return r
    return NotImplemented

#----------------------------------------------------------------------------#
#                                                                            #
#                                  fmpz                                      #
#                                                                            #
#----------------------------------------------------------------------------#


cdef class fmpz:
    """
    The fmpz type represents multiprecision integers.

        >>> fmpz(3) ** 25
        fmpz(847288609443)

    """

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
                if fmpz_set_any_ref(self.val, val) == FMPZ_UNKNOWN: # XXX
                    raise TypeError("cannot create fmpz from type", type(val))

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
        cdef fmpz_struct * sval
        cdef int ttype
        sval = &((<fmpz>s).val[0])
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
            res = str(s)
        finally:
            libc.stdlib.free(s)
        return res

    def __repr__(self):
        return "fmpz(%s)" % self.__str__()

    def __nonzero__(self):
        return not fmpz_is_zero(self.val)

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

    def gcd(self, other):
        """
        Returns the greatest common divisor of self and other.

            >>> fmpz(30).gcd(45)
            fmpz(15)
        """
        cdef fmpz_struct tval[1]
        cdef int ttype = FMPZ_UNKNOWN
        ttype = fmpz_set_any_ref(tval, other)
        if ttype == FMPZ_UNKNOWN:
            raise TypeError("input must be an integer")
        u = fmpz.__new__(fmpz)
        fmpz_gcd((<fmpz>u).val, self.val, tval)
        if ttype == FMPZ_TMP:
            fmpz_clear(tval)
        return u

    def factor(self):
        """
        Factors self into pseudoprimes, returning a list of
        (prime, exp) pairs. The sign is ignored.

            >>> fmpz(5040).factor()
            [(fmpz(2), 4), (fmpz(3), 2), (fmpz(5), 1), (fmpz(7), 1)]

        Warning: this is a preliminary implementation and should only
        be used for integers that fit in a single word, or which are
        larger but only have very small prime factors. It will
        crash if trying to factor a large integer with large
        prime factors.
        """
        cdef fmpz_factor_t fac
        cdef int i
        fmpz_factor_init(fac)
        fmpz_factor(fac, self.val)
        res = [0] * fac.num
        for 0 <= i < fac.num:
            u = fmpz.__new__(fmpz)
            fmpz_set((<fmpz>u).val, &fac.p[i])
            exp = <long> fac.exp[i]
            res[i] = (u, exp)
        fmpz_factor_clear(fac)
        return res


#----------------------------------------------------------------------------#
#                                                                            #
#                               fmpz_poly                                    #
#                                                                            #
#----------------------------------------------------------------------------#

cdef class fmpz_poly:
    """
    The fmpz_poly type represents dense univariate polynomials over
    the integers.

        >>> fmpz_poly([1,2,3]) ** 3
        fmpz_poly([1, 6, 21, 44, 63, 54, 27])
        >>> divmod(fmpz_poly([2,0,1,1,6]), fmpz_poly([3,5,7]))
        (fmpz_poly([]), fmpz_poly([2, 0, 1, 1, 6]))
    """

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

    def __len__(self):
        return fmpz_poly_length(self.val)

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

    def __iter__(self):
        cdef long i, n
        n = self.length()
        for i from 0 <= i < n:
            yield self[i]

    def coeffs(self):
        cdef long i, n
        cdef list L
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
            res = str(s)
        finally:
            libc.stdlib.free(s)
        return res

    def __repr__(self):
        return "fmpz_poly(%s)" % map(int, self.coeffs())

    def __nonzero__(self):
        return not fmpz_poly_is_zero(self.val)

    def __call__(self, other):
        t = any_as_fmpz(other)
        if t is not NotImplemented:
            v = fmpz.__new__(fmpz)
            fmpz_poly_evaluate_fmpz((<fmpz>v).val, self.val, (<fmpz>t).val)
            return v
        t = any_as_fmpz_poly(other)
        if t is not NotImplemented:
            v = fmpz_poly.__new__(fmpz_poly)
            fmpz_poly_compose((<fmpz_poly>v).val, self.val, (<fmpz_poly>t).val)
            return v
        t = any_as_fmpq(other)
        if t is not NotImplemented:
            return fmpq_poly(self)(t)
        t = any_as_fmpq_poly(other)
        if t is not NotImplemented:
            return fmpq_poly(self)(t)
        raise TypeError("cannot call fmpz_poly with input of type %s", type(other))

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

    def gcd(self, other):
        """
        Returns the greatest common divisor of self and other.

            >>> A = fmpz_poly([2,0,1,0,5]); B = fmpz_poly([2,3,4])
            >>> (A*B).gcd(B)
            fmpz_poly([2, 3, 4])

        """
        cdef fmpz_poly res
        other = any_as_fmpz_poly(other)
        if other is NotImplemented:
            raise TypeError("cannot convert input to fmpz_poly")
        res = fmpz_poly.__new__(fmpz_poly)
        fmpz_poly_gcd(res.val, self.val, (<fmpz_poly>other).val)
        return res

    def factor(self):
        """
        Factors self into irreducible factors, returning a tuple
        (c, factors) where c is the content of the coefficients and
        factors is a list of (poly, exp) pairs.

            >>> (-73 * fmpz_poly([1,2,3]) ** 3 * fmpz_poly([5,6,7,8,9]) ** 8).factor()
            (fmpz(-73), [(fmpz_poly([1, 2, 3]), 3), (fmpz_poly([5, 6, 7, 8, 9]), 8)])
            >>> chebyshev_t_polynomial(6).factor()
            (fmpz(1), [(fmpz_poly([-1, 0, 2]), 1), (fmpz_poly([1, 0, -16, 0, 16]), 1)])
            >>> (fmpz_poly([-1,1])**100).factor()
            (fmpz(1), [(fmpz_poly([-1, 1]), 100)])
            >>> fmpz_poly([1,2,3,4,5,6]).factor()
            (fmpz(1), [(fmpz_poly([1, 2, 3, 4, 5, 6]), 1)])

        """
        cdef fmpz_poly_factor_t fac
        cdef int i
        fmpz_poly_factor_init(fac)
        fmpz_poly_factor_zassenhaus(fac, self.val)
        res = [0] * fac.num
        for 0 <= i < fac.num:
            u = fmpz_poly.__new__(fmpz_poly)
            fmpz_poly_set((<fmpz_poly>u).val, &fac.p[i])
            exp = fac.exp[i]
            res[i] = (u, exp)
        c = fmpz.__new__(fmpz)
        fmpz_set((<fmpz>c).val, &fac.c)
        fmpz_poly_factor_clear(fac)
        return c, res

#----------------------------------------------------------------------------#
#                                                                            #
#                               fmpz_mat                                     #
#                                                                            #
#----------------------------------------------------------------------------#

cdef class fmpz_mat:
    """
    The fmpz_mat type represents dense matrices over the integers.

    An fmpz_mat can be constructed empty, from a list of entries
    in row-major order, or from an existing matrix::

        >>> fmpz_mat(2, 4)
        fmpz_mat(2, 4, [0, 0, 0, 0, 0, 0, 0, 0])
        >>> A = fmpz_mat(2, 4, range(8))
        >>> print(A)
        [0, 1, 2, 3]
        [4, 5, 6, 7]
        >>> fmpz_mat(A) == A
        True

    Entries can be accessed and set::

        >>> A[0,1]
        fmpz(1)
        >>> A[1,3] = 8
        >>> A
        fmpz_mat(2, 4, [0, 1, 2, 3, 4, 5, 6, 8])

    Arithmetic operations are supported::

        >>> A * 3
        fmpz_mat(2, 4, [0, 3, 6, 9, 12, 15, 18, 24])
        >>> A + (-A)
        fmpz_mat(2, 4, [0, 0, 0, 0, 0, 0, 0, 0])
        >>> A * fmpz_mat(4, 3, range(12))
        fmpz_mat(2, 3, [42, 48, 54, 123, 146, 169])
        >>> A * fmpz_mat(3, 3, range(9))
        Traceback (most recent call last):
          ...
        ValueError: incompatible shapes for matrix multiplication

    The ~ operator returns the inverse of the matrix, which will
    be of type fmpq_mat::

        >>> print(~fmpz_mat(3,3,[1,2,4,0,1,1,2,-1,0]))
        [-1/3,  4/3,  2/3]
        [-2/3,  8/3,  1/3]
        [ 2/3, -5/3, -1/3]

    """

    cdef fmpz_mat_t val

    def __cinit__(self):
        fmpz_mat_init(self.val, 0, 0)

    def __dealloc__(self):
        fmpz_mat_clear(self.val)

    @cython.embedsignature(False)
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

    def __nonzero__(self):
        return not fmpz_mat_is_zero(self.val)

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
        """
        Returns the number of rows of self.
        """
        return fmpz_mat_nrows(self.val)

    cpdef long ncols(self):
        """
        Returns the number of columns of self.
        """
        return fmpz_mat_ncols(self.val)

    def __repr__(self):
        return "fmpz_mat(%i, %i, [%s])" % (self.nrows(), self.ncols(),
            (", ".join(map(str, self.entries()))))

    def __str__(self):
        return matrix_to_str(self.table())

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
        """
        Returns self as a flat list of fmpz entries,
        output in row-major order.

            >>> A = fmpz_mat(2, 4, range(8))
            >>> A.entries()
            [fmpz(0), fmpz(1), fmpz(2), fmpz(3), fmpz(4), fmpz(5), fmpz(6), fmpz(7)]

        """
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
        """
        Returns self as a nested list of fmpz entries,
        output in row-major order.

            >>> A = fmpz_mat(2, 4, range(8))
            >>> A.table()
            [[fmpz(0), fmpz(1), fmpz(2), fmpz(3)], [fmpz(4), fmpz(5), fmpz(6), fmpz(7)]]
        """
        cdef long i, m, n
        m = self.nrows()
        n = self.ncols()
        L = self.entries()
        return [L[i*n:(i+1)*n] for i in range(m)]

    def det(self):
        """
        Returns the determinant of self as an fmpz.

            >>> A = fmpz_mat(3, 3, range(9))
            >>> A.det()
            fmpz(0)
            >>> A[2,2] = 10
            >>> A.det()
            fmpz(-6)
            >>> (A * A).det()
            fmpz(36)
            >>> fmpz_mat(0, 0).det()
            fmpz(1)

        """
        cdef fmpz d
        if not fmpz_mat_is_square(self.val):
            raise ValueError("matrix must be square")
        d = fmpz.__new__(fmpz)
        fmpz_mat_det(d.val, self.val)
        return d

    def __pos__(self):
        return self

    def __neg__(self):
        cdef fmpz_mat t = fmpz_mat(self)
        fmpz_mat_neg(t.val, t.val)   # XXX
        return t

    def __add__(s, t):
        cdef fmpz_mat u
        cdef fmpz_mat_struct *sval
        cdef fmpz_mat_struct *tval
        sm = any_as_fmpz_mat(s)
        if sm is NotImplemented:
            return sm
        tm = any_as_fmpz_mat(t)
        if tm is NotImplemented:
            return tm
        sval = &(<fmpz_mat>s).val[0]
        tval = &(<fmpz_mat>t).val[0]
        if (fmpz_mat_nrows(sval) != fmpz_mat_nrows(tval) or
           fmpz_mat_ncols(sval) != fmpz_mat_ncols(tval)):
            raise ValueError("incompatible shapes for matrix addition")
        u = fmpz_mat.__new__(fmpz_mat)
        fmpz_mat_init(u.val, fmpz_mat_nrows(sval), fmpz_mat_ncols(sval))
        fmpz_mat_add(u.val, sval, tval)
        return u

    def __sub__(s, t):
        cdef fmpz_mat u
        cdef fmpz_mat_struct *sval
        cdef fmpz_mat_struct *tval
        sm = any_as_fmpz_mat(s)
        if sm is NotImplemented:
            return sm
        tm = any_as_fmpz_mat(t)
        if tm is NotImplemented:
            return tm
        sval = &(<fmpz_mat>s).val[0]
        tval = &(<fmpz_mat>t).val[0]
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
        cdef fmpz_mat_struct *sval
        cdef fmpz_mat_struct *tval
        cdef int ttype
        if typecheck(s, fmpz_mat) and typecheck(t, fmpz_mat):
            sval = &(<fmpz_mat>s).val[0]
            tval = &(<fmpz_mat>t).val[0]
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
            c = any_as_fmpq(t)
            if c is not NotImplemented:
                # XXX: improve this
                return fmpq_mat(s) * t
        return NotImplemented

    def __div__(fmpz_mat s, t):
        return s * (1 / fmpq(t))

    def __truediv__(fmpz_mat s, t):
        return s.__div__(t)

    @classmethod
    def randtest(cls, ulong m, ulong n, ulong bits):
        """
        Returns a random (m, n) matrix with non-uniformly chosen
        entries up to the specified number of bits in size. Small and
        zero entries are generated with increased probability.

            >>> print(fmpz_mat.randtest(3, 2, 100))   # doctest: +SKIP
            [        5442103460568216, 1906839377153448]
            [-37778931862922801979391,                0]
            [                       0,                1]

        """
        cdef fmpz_mat mat = fmpz_mat(m, n)
        fmpz_mat_randtest(mat.val, global_random_state, bits)
        return mat

    @classmethod
    def randbits(cls, ulong m, ulong n, ulong bits):
        """
        Returns a random (m, n) matrix with uniformly chosen entries up
        to the specified number of bits in size.

            >>> print(fmpz_mat.randbits(3, 2, 100))   # doctest: +SKIP
            [  502804798116524380422349115480, -93136769489619409388141424916]
            [-1201505897735399116254292047234, 145439343004100224514654363320]
            [ 1183889243483733739229662952032, 632771515833349927306121868518]

        """
        cdef fmpz_mat mat = fmpz_mat(m, n)
        fmpz_mat_randbits(mat.val, global_random_state, bits)
        return mat

    @classmethod
    def randrank(cls, ulong m, ulong n, ulong rank, ulong bits):
        """
        Returns a random sparse (m, n) matrix of the specified rank
        with entries up to the specified number of bits in size.

            >>> print(flint.fmpz_mat.randrank(3,6,2,20))   # doctest: +SKIP
            [0, 484749, 0, 0, 0,     0]
            [0,      0, 0, 0, 0,     0]
            [0,      0, 0, 0, 0, -2048]

        """
        cdef fmpz_mat mat
        if rank > m or rank > n:
            raise ValueError("impossible rank")
        mat = fmpz_mat(m, n)
        fmpz_mat_randrank(mat.val, global_random_state, rank, bits)
        return mat

    def rank(self):
        """
        Returns the rank of self.

            >>> A = fmpz_mat(3, 3, range(9))
            >>> A.rank()
            2
            >>> A[2,2] = 10
            >>> A.rank()
            3
        """
        return fmpz_mat_rank(self.val)

    def __invert__(self):
        # XXX: write flint function for this
        cdef fmpz_mat_t tmp
        cdef fmpq_mat u
        cdef fmpz_t den
        if not fmpz_mat_is_square(self.val):
            raise ValueError("matrix must be square")
        fmpz_mat_init_set(tmp, self.val)
        fmpz_init(den)
        try:
            fmpz_mat_inv(tmp, den, self.val)
            if fmpz_is_zero(den):
                raise ZeroDivisionError("matrix is singular")
            u = fmpq_mat.__new__(fmpq_mat)
            fmpq_mat_init(u.val, fmpz_mat_nrows(self.val), fmpz_mat_ncols(self.val))
            fmpq_mat_set_fmpz_mat_div_fmpz(u.val, tmp, den)
            return u
        finally:
            fmpz_clear(den)
            fmpz_mat_clear(tmp)

    def transpose(self):
        """
        Returns the transpose of self.

            >>> fmpz_mat(2,3,range(6)).transpose()
            fmpz_mat(3, 2, [0, 3, 1, 4, 2, 5])
        """
        cdef fmpz_mat u
        u = fmpz_mat.__new__(fmpz_mat)
        fmpz_mat_init(u.val, fmpz_mat_ncols(self.val), fmpz_mat_nrows(self.val))
        fmpz_mat_transpose(u.val, self.val)
        return u

    def solve(self, other):
        """
        Given self = A and other = B, returns a matrix X and
        denominator den such that A*(X/den) = B, assuming that
        self is square and invertible.

        TODO: support Dixon solver

            >>> A = fmpz_mat(2, 2, [1,4,8,3])
            >>> B = fmpz_mat(2, 3, range(6))
            >>> X, d = A.solve(B)
            >>> A*X == B*d
            True
            >>> fmpz_mat(2, 2, [1,0,2,0]).solve(B)
            Traceback (most recent call last):
              ...
            ZeroDivisionError: singular matrix in solve()
            >>> A.solve(fmpz_mat(1, 2, [2,3]))
            Traceback (most recent call last):
              ...
            ValueError: need a square system and compatible right hand side

        """
        cdef fmpz_mat u
        cdef fmpz d
        cdef int result
        t = any_as_fmpz_mat(other)
        if t is NotImplemented:
            raise TypeError("cannot convert input to fmpz_mat")
        if (fmpz_mat_nrows(self.val) != fmpz_mat_ncols(self.val) or
            fmpz_mat_nrows(self.val) != fmpz_mat_nrows((<fmpz_mat>t).val)):
            raise ValueError("need a square system and compatible right hand side")
        u = fmpz_mat.__new__(fmpz_mat)
        fmpz_mat_init(u.val, fmpz_mat_nrows((<fmpz_mat>t).val),
            fmpz_mat_ncols((<fmpz_mat>t).val))
        d = fmpz.__new__(fmpz)
        result = fmpz_mat_solve(u.val, d.val, self.val, (<fmpz_mat>t).val)
        if not result:
            raise ZeroDivisionError("singular matrix in solve()")
        return u, d

    def rref(self, inplace=False):
        """
        Computes the reduced row echelon form (rref) of self,
        either returning a new copy or modifying self in-place.
        Returns (rref, denominator, rank).

            >>> A = fmpz_mat(3,3,range(9))
            >>> A.rref()
            (fmpz_mat(3, 3, [3, 0, -3, 0, 3, 6, 0, 0, 0]), fmpz(3), 2)
            >>> A.rref(inplace=True)
            (fmpz_mat(3, 3, [3, 0, -3, 0, 3, 6, 0, 0, 0]), fmpz(3), 2)
            >>> A
            fmpz_mat(3, 3, [3, 0, -3, 0, 3, 6, 0, 0, 0])

        """
        cdef fmpz d
        if inplace:
            res = self
        else:
            res = fmpz_mat.__new__(fmpz_mat)
            fmpz_mat_init((<fmpz_mat>res).val, fmpz_mat_nrows(self.val), fmpz_mat_ncols(self.val))
        d = fmpz.__new__(fmpz)
        rank = fmpz_mat_rref((<fmpz_mat>res).val, d.val, self.val)
        return res, d, rank

    def nullspace(self):
        """
        Computes a basis for the nullspace of self. Returns (X, nullity)
        where nullity is the rank of the nullspace of self and X is a
        matrix whose first (nullity) columns are linearly independent,
        and such that self * X = 0.

            >>> A = fmpz_mat(3,5,range(1,16))
            >>> X, nullity = A.nullspace()
            >>> A.rank(), nullity, X.rank()
            (2, 3, 3)
            >>> A * X
            fmpz_mat(3, 5, [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
            >>> print(X)
            [  5,  10,  15, 0, 0]
            [-10, -15, -20, 0, 0]
            [  5,   0,   0, 0, 0]
            [  0,   5,   0, 0, 0]
            [  0,   0,   5, 0, 0]

        """
        cdef fmpz_mat res
        res = fmpz_mat.__new__(fmpz_mat)
        fmpz_mat_init(res.val, fmpz_mat_ncols(self.val), fmpz_mat_ncols(self.val))
        nullity = fmpz_mat_nullspace(res.val, self.val)
        return res, nullity


#----------------------------------------------------------------------------#
#                                                                            #
#                                  fmpq                                      #
#                                                                            #
#----------------------------------------------------------------------------#

cdef class fmpq:
    """
    The fmpq type represents multiprecision rational numbers.

        >>> fmpq(1,7) + fmpq(50,51)
        fmpq(401,357)

    """

    cdef fmpq_t val

    def __cinit__(self):
        fmpq_init(self.val)

    def __dealloc__(self):
        fmpq_clear(self.val)

    def __init__(self, p=None, q=None):
        cdef long x
        if q is None:
            if p is None:
                return # zero
            else:
                p = any_as_fmpq(p)
                if p is NotImplemented:
                    raise ValueError("cannot create fmpq from object of type %s" % type(p))
                fmpq_set(self.val, (<fmpq>p).val)
                return
        p = any_as_fmpz(p)
        if p is NotImplemented:
            raise ValueError("cannot create fmpq from object of type %s" % type(p))
        q = any_as_fmpz(q)
        if q is NotImplemented:
            raise ValueError("cannot create fmpq from object of type %s" % type(q))
        if fmpz_is_zero((<fmpz>q).val):
            raise ZeroDivisionError("cannot create rational number with zero denominator")
        fmpz_set(fmpq_numref(self.val), (<fmpz>p).val)
        fmpz_set(fmpq_denref(self.val), (<fmpz>q).val)
        fmpq_canonicalise(self.val)

    def __richcmp__(s, t, int op):
        cdef bint res
        s = any_as_fmpq(s)
        if s is NotImplemented:
            return s
        t = any_as_fmpq(t)
        if t is NotImplemented:
            return t
        if op == 2 or op == 3:
            res = fmpq_equal((<fmpq>s).val, (<fmpq>t).val)
            if op == 3:
                res = not res
            return res
        else:
            raise NotImplementedError("fmpq comparisons")

    def numer(self):
        cdef fmpz x = fmpz.__new__(fmpz)
        fmpz_set(x.val, fmpq_numref(self.val))
        return x

    def denom(self):
        cdef fmpz x = fmpz.__new__(fmpz)
        fmpz_set(x.val, fmpq_denref(self.val))
        return x

    p = property(numer)
    q = property(denom)

    def __repr__(self):
        return "fmpq(%s,%s)" % (self.p, self.q)

    def __str__(self):
        return "%s/%s" % (self.p, self.q)

    def __nonzero__(self):
        return not fmpq_is_zero(self.val)

    def __pos__(self):
        return self

    def __neg__(self):
        cdef fmpq r = fmpq.__new__(fmpq)
        fmpq_neg(r.val, self.val)
        return r

    def __abs__(self):
        cdef fmpq r
        if fmpz_sgn(fmpq_numref(self.val)) >= 0:
            return self
        r = fmpq.__new__(fmpq)
        fmpq_neg(r.val, self.val)
        return r

    def __add__(s, t):
        cdef fmpq r
        s = any_as_fmpq(s)
        if s is NotImplemented:
            return s
        t = any_as_fmpq(t)
        if t is NotImplemented:
            return t
        r = fmpq.__new__(fmpq)
        fmpq_add(r.val, (<fmpq>s).val, (<fmpq>t).val)
        return r

    def __sub__(s, t):
        cdef fmpq r
        s = any_as_fmpq(s)
        if s is NotImplemented:
            return s
        t = any_as_fmpq(t)
        if t is NotImplemented:
            return t
        r = fmpq.__new__(fmpq)
        fmpq_sub(r.val, (<fmpq>s).val, (<fmpq>t).val)
        return r

    def __mul__(s, t):
        cdef fmpq r
        s = any_as_fmpq(s)
        if s is NotImplemented:
            return s
        t = any_as_fmpq(t)
        if t is NotImplemented:
            return t
        r = fmpq.__new__(fmpq)
        fmpq_mul(r.val, (<fmpq>s).val, (<fmpq>t).val)
        return r

    def __div__(s, t):
        cdef fmpq r
        s = any_as_fmpq(s)
        if s is NotImplemented:
            return s
        t = any_as_fmpq(t)
        if t is NotImplemented:
            return t
        if fmpq_is_zero((<fmpq>t).val):
            raise ZeroDivisionError("fmpq division by zero")
        r = fmpq.__new__(fmpq)
        fmpq_div(r.val, (<fmpq>s).val, (<fmpq>t).val)
        return r

    # __truediv__ = __div__ doesn't seem to work?
    def __truediv__(s, t):
        cdef fmpq r
        s = any_as_fmpq(s)
        if s is NotImplemented:
            return s
        t = any_as_fmpq(t)
        if t is NotImplemented:
            return t
        if fmpq_is_zero((<fmpq>t).val):
            raise ZeroDivisionError("fmpq division by zero")
        r = fmpq.__new__(fmpq)
        fmpq_div(r.val, (<fmpq>s).val, (<fmpq>t).val)
        return r


#----------------------------------------------------------------------------#
#                                                                            #
#                               fmpq_poly                                    #
#                                                                            #
#----------------------------------------------------------------------------#

cdef class fmpq_poly:
    """
    The fmpq_poly type represents dense univariate polynomials
    over the rational numbers. For efficiency reasons, an fmpq_poly is
    structurally an integer polynomial with a single common denominator.

        >>> fmpq_poly([1,2,3],5) ** 3
        fmpq_poly([1, 6, 21, 44, 63, 54, 27], 125)
        >>> print _
        27/125*x^6 + 54/125*x^5 + 63/125*x^4 + 44/125*x^3 + 21/125*x^2 + 6/125*x + 1/125
        >>> divmod(fmpq_poly([2,0,1,1,6]), fmpq_poly([3,5,7]))
        (fmpq_poly([38, -161, 294], 343), fmpq_poly([572, 293], 343))

    """

    cdef fmpq_poly_t val

    def __cinit__(self):
        fmpq_poly_init(self.val)

    def __dealloc__(self):
        fmpq_poly_clear(self.val)

    def __init__(self, p=None, q=None):
        if p is not None:
            if typecheck(p, fmpq_poly):
                fmpq_poly_set(self.val, (<fmpq_poly>p).val)
            elif typecheck(p, fmpz_poly):
                fmpq_poly_set_fmpz_poly(self.val, (<fmpz_poly>p).val)
            elif isinstance(p, list):
                fmpq_poly_set_list(self.val, p)
            else:
                raise TypeError("cannot create fmpq_poly from input of type %s", type(p))
        if q is not None:
            q = any_as_fmpz(q)
            if q is NotImplemented:
                raise TypeError("denominator must be an integer, got %s", type(q))
            if fmpz_is_zero((<fmpz>q).val):
                raise ZeroDivisionError("cannot create fmpq_poly with zero denominator")
            fmpq_poly_scalar_div_fmpz(self.val, self.val, (<fmpz>q).val)

    def __len__(self):
        return fmpq_poly_length(self.val)

    cpdef long length(self):
        return fmpq_poly_length(self.val)

    cpdef long degree(self):
        return fmpq_poly_degree(self.val)

    def __richcmp__(self, other, int op):
        cdef bint r
        if op != 2 and op != 3:
            raise TypeError("polynomials cannot be ordered")
        self = any_as_fmpq_poly(self)
        if self is NotImplemented:
            return self
        other = any_as_fmpq_poly(other)
        if other is NotImplemented:
            return other
        r = fmpq_poly_equal((<fmpq_poly>self).val, (<fmpq_poly>other).val)
        if op == 3:
            r = not r
        return r

    def numer(self):
        cdef fmpz_poly x = fmpz_poly.__new__(fmpz_poly)
        fmpq_poly_get_numerator(x.val, self.val)
        return x

    def denom(self):
        cdef fmpz x = fmpz.__new__(fmpz)
        fmpz_set(x.val, fmpq_poly_denref(self.val))
        return x

    p = property(numer)
    q = property(denom)

    def __iter__(self):
        cdef long i, n
        n = self.length()
        for i from 0 <= i < n:
            yield self[i]

    def coeffs(self):
        cdef long i, n
        cdef list L
        n = self.length()
        L = [fmpq() for i in range(n)]
        for i from 0 <= i < n:
            fmpq_poly_get_coeff_fmpq((<fmpq>(L[i])).val, self.val, i)
        return L

    def __getitem__(self, long i):
        cdef fmpq x
        x = fmpq()
        if i < 0:
            return x
        fmpq_poly_get_coeff_fmpq(x.val, self.val, i)
        return x

    def __setitem__(self, long i, x):
        if i < 0:
            raise ValueError("cannot assign to index < 0 of polynomial")
        v = fmpq(x)  # XXX
        fmpq_poly_set_coeff_fmpq(self.val, i, (<fmpq>v).val)

    def __repr__(self):
        #return "fmpq_poly(%r)" % self.coeffs()
        d = self.denom()
        n = self.numer()
        if d == 1:
            return "fmpq_poly(%s)" % map(int, n.coeffs())
        else:
            return "fmpq_poly(%s, %s)" % (map(int, n.coeffs()), d)

    def __str__(self):
        cdef char * s = fmpq_poly_get_str_pretty(self.val, "x")
        try:
            res = str(s)
        finally:
            libc.stdlib.free(s)
        return res

    def __nonzero__(self):
        return not fmpq_poly_is_zero(self.val)

    def __call__(self, other):
        t = any_as_fmpz(other)
        if t is not NotImplemented:
            v = fmpq.__new__(fmpq)
            fmpq_poly_evaluate_fmpz((<fmpq>v).val, self.val, (<fmpz>t).val)
            return v
        t = any_as_fmpq(other)
        if t is not NotImplemented:
            v = fmpq.__new__(fmpq)
            fmpq_poly_evaluate_fmpq((<fmpq>v).val, self.val, (<fmpq>t).val)
            return v
        t = any_as_fmpq_poly(other)
        if t is not NotImplemented:
            v = fmpq_poly.__new__(fmpq_poly)
            fmpq_poly_compose((<fmpq_poly>v).val, self.val, (<fmpq_poly>t).val)
            return v
        raise TypeError("cannot call fmpq_poly with input of type %s", type(other))

    def __pos__(self):
        return self

    def __neg__(self):
        cdef fmpq_poly res = fmpq_poly.__new__(fmpq_poly)
        fmpq_poly_neg(res.val, self.val)
        return res

    def __add__(s, t):
        cdef fmpq_poly r
        s = any_as_fmpq_poly(s)
        if s is NotImplemented:
            return s
        t = any_as_fmpq_poly(t)
        if t is NotImplemented:
            return t
        r = fmpq_poly.__new__(fmpq_poly)
        fmpq_poly_add(r.val, (<fmpq_poly>s).val, (<fmpq_poly>t).val)
        return r

    def __sub__(s, t):
        cdef fmpq_poly r
        s = any_as_fmpq_poly(s)
        if s is NotImplemented:
            return s
        t = any_as_fmpq_poly(t)
        if t is NotImplemented:
            return t
        r = fmpq_poly.__new__(fmpq_poly)
        fmpq_poly_sub(r.val, (<fmpq_poly>s).val, (<fmpq_poly>t).val)
        return r

    def __mul__(s, t):
        cdef fmpq_poly r
        s = any_as_fmpq_poly(s)
        if s is NotImplemented:
            return s
        t = any_as_fmpq_poly(t)
        if t is NotImplemented:
            return t
        r = fmpq_poly.__new__(fmpq_poly)
        fmpq_poly_mul(r.val, (<fmpq_poly>s).val, (<fmpq_poly>t).val)
        return r

    def __floordiv__(s, t):
        cdef fmpq_poly r
        s = any_as_fmpq_poly(s)
        if s is NotImplemented:
            return s
        t = any_as_fmpq_poly(t)
        if t is NotImplemented:
            return t
        if fmpq_poly_is_zero((<fmpq_poly>t).val):
            raise ZeroDivisionError("fmpq_poly division by 0")
        r = fmpq_poly.__new__(fmpq_poly)
        fmpq_poly_div(r.val, (<fmpq_poly>s).val, (<fmpq_poly>t).val)
        return r

    def __mod__(s, t):
        cdef fmpq_poly r
        s = any_as_fmpq_poly(s)
        if s is NotImplemented:
            return s
        t = any_as_fmpq_poly(t)
        if t is NotImplemented:
            return t
        if fmpq_poly_is_zero((<fmpq_poly>t).val):
            raise ZeroDivisionError("fmpq_poly division by 0")
        r = fmpq_poly.__new__(fmpq_poly)
        fmpq_poly_rem(r.val, (<fmpq_poly>s).val, (<fmpq_poly>t).val)
        return r

    def __div__(fmpq_poly s, t):
        cdef fmpq_poly r
        t = any_as_fmpq(t)
        if t is NotImplemented:
            return t
        if fmpq_is_zero((<fmpq>t).val):
            raise ZeroDivisionError("fmpq_poly scalar division by 0")
        r = fmpq_poly.__new__(fmpq_poly)
        # XXX: implement function in flint
        fmpq_poly_scalar_mul_fmpz(r.val, (<fmpq_poly>s).val, fmpq_denref((<fmpq>t).val))
        fmpq_poly_scalar_div_fmpz(r.val, (<fmpq_poly>r).val, fmpq_numref((<fmpq>t).val))
        return r

    # __truediv__ = __div__ doesn't seem to work?
    def __truediv__(fmpq_poly s, t):
        return fmpq_poly.__div__(s, t)

    def __divmod__(s, t):
        cdef fmpq_poly P, Q
        s = any_as_fmpq_poly(s)
        if s is NotImplemented:
            return s
        t = any_as_fmpq_poly(t)
        if t is NotImplemented:
            return t
        if fmpq_poly_is_zero((<fmpq_poly>t).val):
            raise ZeroDivisionError("fmpq_poly divmod by 0")
        P = fmpq_poly.__new__(fmpq_poly)
        Q = fmpq_poly.__new__(fmpq_poly)
        fmpq_poly_divrem(P.val, Q.val, (<fmpq_poly>s).val, (<fmpq_poly>t).val)
        return P, Q

    def __pow__(fmpq_poly self, ulong exp, mod):
        cdef fmpq_poly res
        if mod is not None:
            raise NotImplementedError("fmpz_poly modular exponentiation")
        res = fmpq_poly.__new__(fmpq_poly)
        fmpq_poly_pow(res.val, self.val, exp)
        return res

    def gcd(self, other):
        """
        Returns the greatest common divisor of self and other.

            >>> A = fmpq_poly([1,2,6],6); B = fmpq_poly([4,2,1],12)
            >>> (A * B).gcd(B)
            fmpq_poly([4, 2, 1])

        """
        cdef fmpq_poly res
        other = any_as_fmpq_poly(other)
        if other is NotImplemented:
            raise TypeError("cannot convert input to fmpq_poly")
        res = fmpq_poly.__new__(fmpq_poly)
        fmpq_poly_gcd(res.val, self.val, (<fmpq_poly>other).val)
        return res

    def factor(self):
        """
        Factors self into irreducible polynomials. Returns (c, factors)
        where c is the leading coefficient and factors is a list of
        (poly, exp) pairs with all poly monic.

            >>> legendre_polynomial(5).factor()
            (fmpq(63,8), [(fmpq_poly([0, 1]), 1), (fmpq_poly([15, 0, -70, 0, 63], 63), 1)])
            >>> (fmpq_poly([1,-1],10) ** 5 * fmpq_poly([1,2,3],7)).factor()
            (fmpq(-3,700000), [(fmpq_poly([1, 2, 3], 3), 1), (fmpq_poly([-1, 1]), 5)])

        """
        c, fac = self.numer().factor()
        c = fmpq(c)
        for i in range(len(fac)):
            base, exp = fac[i]
            lead = base[base.degree()]
            base = fmpq_poly(base, lead)
            c *= lead ** exp
            fac[i] = (base, exp)
        return c / self.denom(), fac

#----------------------------------------------------------------------------#
#                                                                            #
#                               fmpq_mat                                     #
#                                                                            #
#----------------------------------------------------------------------------#

cdef class fmpq_mat:
    """
    The fmpq_mat type represents dense matrices over the rational
    numbers.

        >>> A = fmpq_mat(3,3,[1,3,5,2,4,6,fmpq(2,3),2,4])
        >>> print (~A)
        [-3/1,  3/2,  3/2]
        [ 3/1, -1/2, -3/1]
        [-1/1,  0/1,  3/2]
        >>> print (~A) * A
        [1/1, 0/1, 0/1]
        [0/1, 1/1, 0/1]
        [0/1, 0/1, 1/1]

    """

    cdef fmpq_mat_t val

    def __cinit__(self):
        fmpq_mat_init(self.val, 0, 0)

    def __dealloc__(self):
        fmpq_mat_clear(self.val)

    @cython.embedsignature(False)
    def __init__(self, *args):
        cdef long m, n, i, j
        if len(args) == 1:
            val = args[0]
            if typecheck(val, fmpq_mat):
                # XXX: need fmpq_mat_init_set(self.val, (<fmpq_mat>val).val)
                fmpq_mat_init(self.val, fmpq_mat_nrows((<fmpq_mat>val).val),
                                        fmpq_mat_ncols((<fmpq_mat>val).val))
                fmpq_mat_set(self.val, (<fmpq_mat>val).val)
            elif typecheck(val, fmpz_mat):
                fmpq_mat_init(self.val, fmpz_mat_nrows((<fmpz_mat>val).val),
                                        fmpz_mat_ncols((<fmpz_mat>val).val))
                fmpq_mat_set_fmpz_mat(self.val, (<fmpz_mat>val).val)
            else:
                raise TypeError("cannot create fmpq_mat from input of type %s" % type(val))
        elif len(args) == 2:
            m, n = args
            fmpq_mat_init(self.val, m, n)
        elif len(args) == 3:
            m, n, entries = args
            fmpq_mat_init(self.val, m, n)
            entries = list(entries)
            if len(entries) != m*n:
                raise ValueError("list of entries has the wrong length")
            for i from 0 <= i < m:
                for j from 0 <= j < n:
                    # XXX: slow
                    x = fmpq(entries[i*n + j])
                    fmpq_set(fmpq_mat_entry(self.val, i, j), (<fmpq>x).val)
        else:
            raise ValueError("fmpq_mat: expected 1-3 arguments")

    def __nonzero__(self):
        return not fmpq_mat_is_zero(self.val)

    def __richcmp__(s, t, int op):
        cdef bint r
        if op != 2 and op != 3:
            raise TypeError("matrices cannot be ordered")
        s = any_as_fmpq_mat(s)
        if t is NotImplemented:
            return s
        t = any_as_fmpq_mat(t)
        if t is NotImplemented:
            return t
        r = fmpq_mat_equal((<fmpq_mat>s).val, (<fmpq_mat>t).val)
        if op == 3:
            r = not r
        return r

    cpdef long nrows(self):
        return fmpq_mat_nrows(self.val)

    cpdef long ncols(self):
        return fmpq_mat_ncols(self.val)

    def __repr__(self):
        return "fmpq_mat(%i, %i, %s)" % (self.nrows(), self.ncols(), self.entries())

    def __str__(self):
        return matrix_to_str(self.table())

    def __getitem__(self, index):
        cdef long i, j
        cdef fmpq x
        i, j = index
        if i < 0 or i >= self.nrows() or j < 0 or j >= self.ncols():
            raise ValueError("index %i,%i exceeds matrix dimensions" % (i, j))
        x = fmpq.__new__(fmpq)
        fmpq_set(x.val, fmpq_mat_entry(self.val, i, j))
        return x

    def __setitem__(self, index, value):
        cdef long i, j
        i, j = index
        if i < 0 or i >= self.nrows() or j < 0 or j >= self.ncols():
            raise ValueError("index %i,%i exceeds matrix dimensions" % (i, j))
        c = fmpq(value)  # XXX
        fmpq_set(fmpq_mat_entry(self.val, i, j), (<fmpq>c).val)

    def entries(self):
        cdef long i, j, m, n
        cdef fmpq t
        m = self.nrows()
        n = self.ncols()
        L = [None] * (m * n)
        for i from 0 <= i < m:
            for j from 0 <= j < n:
                t = fmpq.__new__(fmpq)
                fmpq_set(t.val, fmpq_mat_entry(self.val, i, j))
                L[i*n + j] = t
        return L

    def table(self):
        cdef long i, m, n
        m = self.nrows()
        n = self.ncols()
        L = self.entries()
        return [L[i*n:(i+1)*n] for i in range(m)]

    def det(self):
        """
        Returns the determinant of self as an fmpq.

            >>> (fmpq_mat(2,2,[1,2,3,4]) / 5).det()
            fmpq(-2,25)
        """
        cdef fmpq d
        if not fmpq_mat_is_square(self.val):
            raise ValueError("matrix must be square")
        d = fmpq.__new__(fmpq)
        fmpq_mat_det(d.val, self.val)
        return d

    def __pos__(self):
        return self

    def __neg__(self):
        cdef fmpq_mat t = fmpq_mat(self)
        fmpq_mat_neg(t.val, t.val)     # XXX
        return t

    def __add__(s, t):
        cdef fmpq_mat u
        cdef fmpq_mat_struct *sval
        cdef fmpq_mat_struct *tval
        s = any_as_fmpq_mat(s)
        if s is NotImplemented:
            return s
        t = any_as_fmpq_mat(t)
        if t is NotImplemented:
            return t
        sval = &(<fmpq_mat>s).val[0]
        tval = &(<fmpq_mat>t).val[0]
        if (fmpq_mat_nrows(sval) != fmpq_mat_nrows(tval) or
           fmpq_mat_ncols(sval) != fmpq_mat_ncols(tval)):
            raise ValueError("incompatible shapes for matrix addition")
        u = fmpq_mat.__new__(fmpq_mat)
        fmpq_mat_init(u.val, fmpq_mat_nrows(sval), fmpq_mat_ncols(sval))
        fmpq_mat_add(u.val, sval, tval)
        return u

    def __sub__(s, t):
        cdef fmpq_mat u
        cdef fmpq_mat_struct *sval
        cdef fmpq_mat_struct *tval
        s = any_as_fmpq_mat(s)
        if s is NotImplemented:
            return s
        t = any_as_fmpq_mat(t)
        if t is NotImplemented:
            return t
        sval = &(<fmpq_mat>s).val[0]
        tval = &(<fmpq_mat>t).val[0]
        if (fmpq_mat_nrows(sval) != fmpq_mat_nrows(tval) or
           fmpq_mat_ncols(sval) != fmpq_mat_ncols(tval)):
            raise ValueError("incompatible shapes for matrix subtraction")
        u = fmpq_mat.__new__(fmpq_mat)
        fmpq_mat_init(u.val, fmpq_mat_nrows(sval), fmpq_mat_ncols(sval))
        fmpq_mat_sub(u.val, sval, tval)
        return u

    cdef __mul_fmpz(self, fmpz c):
        cdef fmpq_mat u
        u = fmpq_mat.__new__(fmpq_mat)
        fmpq_mat_init(u.val, fmpq_mat_nrows(self.val), fmpq_mat_ncols(self.val))
        fmpq_mat_scalar_mul_fmpz(u.val, self.val, c.val)
        return u

    cdef __mul_fmpq(self, fmpq c):
        cdef fmpq_mat u
        u = fmpq_mat.__new__(fmpq_mat)
        fmpq_mat_init(u.val, fmpq_mat_nrows(self.val), fmpq_mat_ncols(self.val))
        fmpq_mat_scalar_mul_fmpz(u.val, self.val, fmpq_numref(c.val))
        fmpq_mat_scalar_div_fmpz(u.val, u.val, fmpq_denref(c.val))
        return u

    cdef __mul_fmpq_mat(self, fmpq_mat other):
        cdef fmpq_mat u
        if self.ncols() != other.nrows():
            raise ValueError("incompatible shapes for matrix multiplication")
        u = fmpq_mat.__new__(fmpq_mat)
        fmpq_mat_init(u.val, fmpq_mat_nrows(self.val), fmpq_mat_ncols(other.val))
        fmpq_mat_mul(u.val, self.val, other.val)
        return u

    cdef __mul_fmpz_mat(self, fmpz_mat other):
        cdef fmpq_mat u
        if self.ncols() != other.nrows():
            raise ValueError("incompatible shapes for matrix multiplication")
        u = fmpq_mat.__new__(fmpq_mat)
        fmpq_mat_init(u.val, fmpq_mat_nrows(self.val), fmpz_mat_ncols(other.val))
        fmpq_mat_mul_fmpz_mat(u.val, self.val, other.val)
        return u

    cdef __mul_r_fmpz_mat(self, fmpz_mat other):
        cdef fmpq_mat u
        if self.nrows() != other.ncols():
            raise ValueError("incompatible shapes for matrix multiplication")
        u = fmpq_mat.__new__(fmpq_mat)
        fmpq_mat_init(u.val, fmpz_mat_nrows(other.val), fmpq_mat_ncols(self.val))
        fmpq_mat_mul_r_fmpz_mat(u.val, other.val, self.val)
        return u

    def __mul__(s, t):
        cdef fmpz_mat u
        if typecheck(s, fmpq_mat):
            if typecheck(t, fmpq_mat):
                return (<fmpq_mat>s).__mul_fmpq_mat(t)
            elif typecheck(t, fmpz_mat):
                return (<fmpq_mat>s).__mul_fmpz_mat(t)
            else:
                c = any_as_fmpz(t)
                if c is not NotImplemented:
                    return (<fmpq_mat>s).__mul_fmpz(c)
                c = any_as_fmpq(t)
                if c is not NotImplemented:
                    return (<fmpq_mat>s).__mul_fmpq(c)
                return NotImplemented
        else:
            if typecheck(s, fmpz_mat):
                return (<fmpq_mat>t).__mul_r_fmpz_mat(s)
            else:
                c = any_as_fmpz(s)
                if c is not NotImplemented:
                    return (<fmpq_mat>t).__mul_fmpz(c)
                c = any_as_fmpq(s)
                if c is not NotImplemented:
                    return (<fmpq_mat>t).__mul_fmpq(c)
                return NotImplemented

    def __div__(fmpq_mat s, t):
        t = any_as_fmpq(t)
        if t is NotImplemented:
            return t
        return s * (1 / t)

    # __truediv__ = __div__ doesn't seem to work?
    def __truediv__(fmpq_mat s, t):
        return fmpq_mat.__div__(s, t)

    def __invert__(self):
        cdef fmpq_mat u
        if not fmpq_mat_is_square(self.val):
            raise ValueError("matrix must be square")
        u = fmpq_mat.__new__(fmpq_mat)
        fmpq_mat_init(u.val, fmpq_mat_nrows(self.val), fmpq_mat_ncols(self.val))
        if not fmpq_mat_inv(u.val, self.val):
            raise ZeroDivisionError("matrix is singular")
        return u

    def transpose(self):
        """
        Returns the transpose of self.
        
            >>> fmpq_mat(2,3,range(6)).transpose()
            fmpq_mat(3, 2, [fmpq(0,1), fmpq(3,1), fmpq(1,1), fmpq(4,1), fmpq(2,1), fmpq(5,1)])
        """
        cdef fmpq_mat u
        u = fmpq_mat.__new__(fmpq_mat)
        fmpq_mat_init(u.val, fmpq_mat_ncols(self.val), fmpq_mat_nrows(self.val))
        fmpq_mat_transpose(u.val, self.val)
        return u

    def solve(self, other, algorithm="fraction-free"):
        """
        Given self = A and other = B, returns a matrix X such
        that A*X = B, assuming that self is square and invertible.

        Algorithm can be "fraction-free" or "dixon"
        (faster for large matrices).

            >>> A = fmpq_mat(2, 2, [1,4,8,3])
            >>> B = fmpq_mat(2, 3, range(6))
            >>> X = A.solve(B)
            >>> print X
            [12/29, 13/29, 14/29]
            [-3/29,  4/29, 11/29]
            >>> A*X == B
            True
            >>> A.solve(B, algorithm='dixon') == X
            True
            >>> fmpq_mat(2, 2, [1,0,2,0]).solve(B)
            Traceback (most recent call last):
              ...
            ZeroDivisionError: singular matrix in solve()
            >>> A.solve(fmpq_mat(1, 2, [2,3]))
            Traceback (most recent call last):
              ...
            ValueError: need a square system and compatible right hand side

        """
        cdef fmpq_mat u
        cdef int result
        t = any_as_fmpq_mat(other)
        if t is NotImplemented:
            raise TypeError("cannot convert input to fmpq_mat")
        if (fmpq_mat_nrows(self.val) != fmpq_mat_ncols(self.val) or
            fmpq_mat_nrows(self.val) != fmpq_mat_nrows((<fmpq_mat>t).val)):
            raise ValueError("need a square system and compatible right hand side")
        u = fmpq_mat.__new__(fmpq_mat)
        fmpq_mat_init(u.val, fmpq_mat_nrows((<fmpq_mat>t).val),
            fmpq_mat_ncols((<fmpq_mat>t).val))
        if algorithm == "fraction-free":
            result = fmpq_mat_solve_fraction_free(u.val, self.val, (<fmpq_mat>t).val)
        else:
            result = fmpq_mat_solve_dixon(u.val, self.val, (<fmpq_mat>t).val)
        if not result:
            raise ZeroDivisionError("singular matrix in solve()")
        return u

    def rref(self, inplace=False):
        """
        Computes the reduced row echelon form (rref) of self,
        either returning a new copy or modifying self in-place.
        Returns (rref, rank).

            >>> A = fmpq_mat(3,3,range(9))
            >>> A.rref()
            (fmpq_mat(3, 3, [fmpq(1,1), fmpq(0,1), fmpq(-1,1), fmpq(0,1), fmpq(1,1), fmpq(2,1), fmpq(0,1), fmpq(0,1), fmpq(0,1)]), 2)
            >>> A.rref(inplace=True)
            (fmpq_mat(3, 3, [fmpq(1,1), fmpq(0,1), fmpq(-1,1), fmpq(0,1), fmpq(1,1), fmpq(2,1), fmpq(0,1), fmpq(0,1), fmpq(0,1)]), 2)
            >>> A
            fmpq_mat(3, 3, [fmpq(1,1), fmpq(0,1), fmpq(-1,1), fmpq(0,1), fmpq(1,1), fmpq(2,1), fmpq(0,1), fmpq(0,1), fmpq(0,1)])

        """
        if inplace:
            res = self
        else:
            res = fmpq_mat.__new__(fmpq_mat)
            fmpq_mat_init((<fmpq_mat>res).val, fmpq_mat_nrows(self.val), fmpq_mat_ncols(self.val))
        rank = fmpq_mat_rref((<fmpq_mat>res).val, self.val)
        return res, rank


#----------------------------------------------------------------------------#
#                                                                            #
#                               nmod                                         #
#                                                                            #
#----------------------------------------------------------------------------#

cdef class nmod:
    """
    The nmod type represents elements of Z/nZ for word-size n.

        >>> nmod(10,17) * 2
        nmod(3, 17)

    """

    cdef mp_limb_t val
    cdef nmod_t mod

    def __init__(self, val, mod):
        cdef mp_limb_t m
        m = mod
        nmod_init(&self.mod, m)
        if not any_as_nmod(&self.val, val, self.mod):
            raise TypeError("cannot create nmod from object of type %s" % type(val))

    def __repr__(self):
        return "nmod(%s, %s)" % (self.val, self.mod.n)

    def __str__(self):
        return str(self.val)

    def __int__(self):
        return self.val

    def __long__(self):
        return self.val

    def modulus(self):
        return self.mod.n

    def __richcmp__(s, t, int op):
        cdef mp_limb_t v
        cdef bint res
        if op != 2 and op != 3:
            raise TypeError("nmods cannot be ordered")
        if typecheck(s, nmod) and typecheck(t, nmod):
            res = ((<nmod>s).val == (<nmod>t).val) and \
                  ((<nmod>s).mod.n == (<nmod>t).mod.n)
            if op == 2:
                return res
            else:
                return not res
        return NotImplemented

    def __nonzero__(self):
        return self.val != 0

    def __pos__(self):
        return self

    def __neg__(self):
        cdef nmod r = nmod.__new__(nmod)
        r.mod = self.mod
        r.val = nmod_neg(self.val, self.mod)
        return r

    def __add__(s, t):
        cdef nmod r
        cdef mp_limb_t val
        if not typecheck(s, nmod):
            s, t = t, s
        if any_as_nmod(&val, t, (<nmod>s).mod):
            r = nmod.__new__(nmod)
            r.mod = (<nmod>s).mod
            r.val = nmod_add(val, (<nmod>s).val, r.mod)
            return r
        return NotImplemented

    def __sub__(s, t):
        cdef nmod r
        cdef mp_limb_t val
        if typecheck(s, nmod):
            if any_as_nmod(&val, t, (<nmod>s).mod):
                r = nmod.__new__(nmod)
                r.mod = (<nmod>s).mod
                r.val = nmod_sub((<nmod>s).val, val, r.mod)
                return r
        else:
            if any_as_nmod(&val, s, (<nmod>t).mod):
                r = nmod.__new__(nmod)
                r.mod = (<nmod>t).mod
                r.val = nmod_sub(val, (<nmod>t).val, r.mod)
                return r
        return NotImplemented

    def __mul__(s, t):
        cdef nmod r
        cdef mp_limb_t val
        if not typecheck(s, nmod):
            s, t = t, s
        if any_as_nmod(&val, t, (<nmod>s).mod):
            r = nmod.__new__(nmod)
            r.mod = (<nmod>s).mod
            r.val = nmod_mul(val, (<nmod>s).val, r.mod)
            return r
        return NotImplemented

    def __div__(s, t):
        cdef nmod r
        cdef mp_limb_t sval, tval, x
        cdef nmod_t mod
        if typecheck(s, nmod):
            mod = (<nmod>s).mod
            sval = (<nmod>s).val
            if not any_as_nmod(&tval, t, mod):
                return NotImplemented
        else:
            mod = (<nmod>t).mod
            tval = (<nmod>t).val
            if not any_as_nmod(&sval, s, mod):
                return NotImplemented
        if tval == 0:
            raise ZeroDivisionError("%s is not invertible mod %s" % (tval, mod.n))
        # XXX: check invertibility?
        x = nmod_div(sval, tval, mod)
        if x == 0:
            raise ZeroDivisionError("%s is not invertible mod %s" % (tval, mod.n))
        r = nmod.__new__(nmod)
        r.mod = mod
        r.val = x
        return r

    def __invert__(self):
        return (1 / self)   # XXX: speed up


#----------------------------------------------------------------------------#
#                                                                            #
#                               nmod_poly                                    #
#                                                                            #
#----------------------------------------------------------------------------#

cdef class nmod_poly:
    """
    The nmod_poly type represents dense univariate polynomials
    over Z/nZ for word-size n.

        >>> a = nmod_poly([5,1,10,14,8], 7)
        >>> print a
        x^4+3*x^2+x+5
        >>> print -a
        6*x^4+4*x^2+6*x+2
        >>> list(nmod_poly(range(3), 2))
        [nmod(0, 2), nmod(1, 2)]
        >>> nmod_poly([1, 2, 3], 23) ** 3
        nmod_poly([1L, 6L, 21L, 21L, 17L, 8L, 4L], 23)
        >>> divmod(nmod_poly([2,0,1,1,6],23), nmod_poly([3,5,7],23))
        (nmod_poly([4L, 0L, 14L], 23), nmod_poly([13L, 3L], 23))

    """

    cdef nmod_poly_t val

    #def __cinit__(self):

    def __dealloc__(self):
        nmod_poly_clear(self.val)

    def __init__(self, val=None, ulong mod=0):
        cdef ulong m2
        if typecheck(val, nmod_poly):
            m2 = nmod_poly_modulus((<nmod_poly>val).val)
            if m2 != mod:
                raise ValueError("different moduli!")
            nmod_poly_init(self.val, m2)
            nmod_poly_set(self.val, (<nmod_poly>val).val)
        else:
            if mod == 0:
                raise ValueError("a nonzero modulus is required")
            nmod_poly_init(self.val, mod)
            if typecheck(val, fmpz_poly):
                fmpz_poly_get_nmod_poly(self.val, (<fmpz_poly>val).val)
            elif typecheck(val, list):
                nmod_poly_set_list(self.val, val)
            else:
                raise TypeError("cannot create nmod_poly from input of type %s", type(val))

    def __len__(self):
        return nmod_poly_length(self.val)

    cpdef long length(self):
        return nmod_poly_length(self.val)

    cpdef long degree(self):
        return nmod_poly_degree(self.val)

    cpdef mp_limb_t modulus(self):
        return nmod_poly_modulus(self.val)

    def __richcmp__(s, t, int op):
        cdef mp_limb_t v
        cdef bint res
        if op != 2 and op != 3:
            raise TypeError("nmod_polyss cannot be ordered")
        if typecheck(s, nmod_poly) and typecheck(t, nmod_poly):
            if (<nmod_poly>s).val.mod.n != (<nmod_poly>t).val.mod.n:
                res = False
            else:
                res = nmod_poly_equal((<nmod_poly>s).val, (<nmod_poly>t).val)
            if op == 2:
                return res
            if op == 3:
                return not res
        return NotImplemented

    def __iter__(self):
        cdef long i, n
        n = self.length()
        for i from 0 <= i < n:
            yield self[i]

    def coeffs(self):
        cdef long i, n
        cdef list L
        cdef mp_limb_t m
        n = self.length()
        m = self.modulus()
        L = [nmod(0,m) for i in range(n)]   # XXX: speed up
        for i from 0 <= i < n:
            (<nmod>(L[i])).val = nmod_poly_get_coeff_ui(self.val, i)
        return L

    def __repr__(self):
        return "nmod_poly(%s, %s)" % (map(int, self.coeffs()), self.modulus())

    def __str__(self):
        return str(fmpz_poly(map(int, self.coeffs())))

    def __getitem__(self, long i):
        cdef nmod x
        x = nmod(0, self.modulus())
        if i < 0:
            return x
        x.val = nmod_poly_get_coeff_ui(self.val, i)
        return x

    def __setitem__(self, long i, x):
        cdef mp_limb_t v
        if i < 0:
            raise ValueError("cannot assign to index < 0 of polynomial")
        if any_as_nmod(&v, x, self.val.mod):
            nmod_poly_set_coeff_ui(self.val, i, v)
        else:
            raise TypeError("cannot set element of type %s" % type(x))

    def __nonzero__(self):
        return not nmod_poly_is_zero(self.val)

    def __call__(self, other):
        cdef mp_limb_t c
        if any_as_nmod(&c, other, self.val.mod):
            v = nmod(0, self.modulus())
            (<nmod>v).val = nmod_poly_evaluate_nmod(self.val, c)
            return v
        t = any_as_nmod_poly(other, self.val.mod)
        if t is not NotImplemented:
            r = nmod_poly.__new__(nmod_poly)
            nmod_poly_init_preinv((<nmod_poly>r).val, self.val.mod.n, self.val.mod.ninv)
            nmod_poly_compose((<nmod_poly>r).val, self.val, (<nmod_poly>t).val)
            return r
        raise TypeError("cannot call nmod_poly with input of type %s", type(other))

    def __pos__(self):
        return self

    def __neg__(self):
        cdef nmod_poly r = nmod_poly.__new__(nmod_poly)
        nmod_poly_init_preinv(r.val, self.val.mod.n, self.val.mod.ninv)
        nmod_poly_neg(r.val, self.val)
        return r

    def __add__(s, t):
        cdef nmod_poly r
        if typecheck(s, nmod_poly):
            t = any_as_nmod_poly(t, (<nmod_poly>s).val.mod)
            if t is NotImplemented:
                return t
        else:
            s = any_as_nmod_poly(s, (<nmod_poly>t).val.mod)
            if s is NotImplemented:
                return s
        if (<nmod_poly>s).val.mod.n != (<nmod_poly>t).val.mod.n:
            raise ValueError("cannot add nmod_polys with different moduli")
        r = nmod_poly.__new__(nmod_poly)
        nmod_poly_init_preinv(r.val, (<nmod_poly>t).val.mod.n, (<nmod_poly>t).val.mod.ninv)
        nmod_poly_add(r.val, (<nmod_poly>s).val, (<nmod_poly>t).val)
        return r

    def __sub__(s, t):
        cdef nmod_poly r
        if typecheck(s, nmod_poly):
            t = any_as_nmod_poly(t, (<nmod_poly>s).val.mod)
            if t is NotImplemented:
                return t
        else:
            s = any_as_nmod_poly(s, (<nmod_poly>t).val.mod)
            if s is NotImplemented:
                return s
        if (<nmod_poly>s).val.mod.n != (<nmod_poly>t).val.mod.n:
            raise ValueError("cannot subtract nmod_polys with different moduli")
        r = nmod_poly.__new__(nmod_poly)
        nmod_poly_init_preinv(r.val, (<nmod_poly>t).val.mod.n, (<nmod_poly>t).val.mod.ninv)
        nmod_poly_sub(r.val, (<nmod_poly>s).val, (<nmod_poly>t).val)
        return r

    def __mul__(s, t):
        cdef nmod_poly r
        if typecheck(s, nmod_poly):
            t = any_as_nmod_poly(t, (<nmod_poly>s).val.mod)
            if t is NotImplemented:
                return t
        else:
            s = any_as_nmod_poly(s, (<nmod_poly>t).val.mod)
            if s is NotImplemented:
                return s
        if (<nmod_poly>s).val.mod.n != (<nmod_poly>t).val.mod.n:
            raise ValueError("cannot multiply nmod_polys with different moduli")
        r = nmod_poly.__new__(nmod_poly)
        nmod_poly_init_preinv(r.val, (<nmod_poly>t).val.mod.n, (<nmod_poly>t).val.mod.ninv)
        nmod_poly_mul(r.val, (<nmod_poly>s).val, (<nmod_poly>t).val)
        return r

    # TODO: __div__, __truediv__

    def __floordiv__(s, t):
        cdef nmod_poly r
        if typecheck(s, nmod_poly):
            t = any_as_nmod_poly(t, (<nmod_poly>s).val.mod)
            if t is NotImplemented:
                return t
        else:
            s = any_as_nmod_poly(s, (<nmod_poly>t).val.mod)
            if s is NotImplemented:
                return s
        if (<nmod_poly>s).val.mod.n != (<nmod_poly>t).val.mod.n:
            raise ValueError("cannot divide nmod_polys with different moduli")
        if nmod_poly_is_zero((<nmod_poly>t).val):
            raise ZeroDivisionError("polynomial division by zero")
        r = nmod_poly.__new__(nmod_poly)
        nmod_poly_init_preinv(r.val, (<nmod_poly>t).val.mod.n, (<nmod_poly>t).val.mod.ninv)
        nmod_poly_div(r.val, (<nmod_poly>s).val, (<nmod_poly>t).val)
        return r

    def __divmod__(s, t):
        cdef nmod_poly P, Q
        if typecheck(s, nmod_poly):
            t = any_as_nmod_poly(t, (<nmod_poly>s).val.mod)
            if t is NotImplemented:
                return t
        else:
            s = any_as_nmod_poly(s, (<nmod_poly>t).val.mod)
            if s is NotImplemented:
                return s
        if (<nmod_poly>s).val.mod.n != (<nmod_poly>t).val.mod.n:
            raise ValueError("cannot divide nmod_polys with different moduli")
        if nmod_poly_is_zero((<nmod_poly>t).val):
            raise ZeroDivisionError("polynomial division by zero")
        P = nmod_poly.__new__(nmod_poly)
        Q = nmod_poly.__new__(nmod_poly)
        nmod_poly_init_preinv(P.val, (<nmod_poly>t).val.mod.n, (<nmod_poly>t).val.mod.ninv)
        nmod_poly_init_preinv(Q.val, (<nmod_poly>t).val.mod.n, (<nmod_poly>t).val.mod.ninv)
        nmod_poly_divrem(P.val, Q.val, (<nmod_poly>s).val, (<nmod_poly>t).val)
        return P, Q

    def __mod__(s, t):
        return divmod(s, t)[1]      # XXX

    def __pow__(nmod_poly self, ulong exp, mod):
        cdef nmod_poly res
        if mod is not None:
            raise NotImplementedError("nmod_poly modular exponentiation")
        res = nmod_poly.__new__(nmod_poly)
        nmod_poly_init_preinv(res.val, (<nmod_poly>self).val.mod.n, (<nmod_poly>self).val.mod.ninv)
        nmod_poly_pow(res.val, self.val, exp)
        return res

    def gcd(self, other):
        """
        Returns the monic greatest common divisor of self and other.

            >>> A = nmod_poly([1,2,3,4], 7); B = nmod_poly([4,1,5], 7)
            >>> (A * B).gcd(B) * 5
            nmod_poly([4L, 1L, 5L], 7)

        """
        cdef nmod_poly res
        other = any_as_nmod_poly(other, (<nmod_poly>self).val.mod)
        if other is NotImplemented:
            raise TypeError("cannot convert input to nmod_poly")
        if self.val.mod.n != (<nmod_poly>other).val.mod.n:
            raise ValueError("moduli must be the same")
        res = nmod_poly.__new__(nmod_poly)
        nmod_poly_init_preinv(res.val, self.val.mod.n, self.val.mod.ninv)
        nmod_poly_gcd(res.val, self.val, (<nmod_poly>other).val)
        return res

    def factor(self, algorithm=None):
        """
        Factors self into irreducible factors, returning a tuple
        (c, factors) where c is the leading coefficient and
        factors is a list of (poly, exp) pairs with all polynomials
        monic.

            >>> nmod_poly(range(10), 3).factor()
            (nmod(2, 3), [(nmod_poly([0L, 1L], 3), 1), (nmod_poly([2L, 1L], 3), 7)])
            >>> nmod_poly(range(10), 19).factor()
            (nmod(9, 19), [(nmod_poly([0L, 1L], 19), 1), (nmod_poly([3L, 7L, 2L, 15L, 1L], 19), 1), (nmod_poly([12L, 15L, 12L, 7L, 1L], 19), 1)])
            >>> nmod_poly(range(10), 53).factor()
            (nmod(9, 53), [(nmod_poly([0L, 1L], 53), 1), (nmod_poly([6L, 12L, 18L, 24L, 30L, 36L, 42L, 48L, 1L], 53), 1)])

        Algorithm can be None (default), 'berlekamp', or
        'cantor-zassenhaus'.

            >>> nmod_poly([3,2,1,2,3], 7).factor(algorithm='berlekamp')
            (nmod(3, 7), [(nmod_poly([2L, 1L], 7), 1), (nmod_poly([4L, 1L], 7), 1), (nmod_poly([1L, 4L, 1L], 7), 1)])
            >>> nmod_poly([3,2,1,2,3], 7).factor(algorithm='cantor-zassenhaus')
            (nmod(3, 7), [(nmod_poly([4L, 1L], 7), 1), (nmod_poly([2L, 1L], 7), 1), (nmod_poly([1L, 4L, 1L], 7), 1)])

        """
        cdef nmod_poly_factor_t fac
        cdef mp_limb_t lead
        cdef int i
        nmod_poly_factor_init(fac)
        if algorithm == 'berlekamp':
            lead = nmod_poly_factor_with_berlekamp(fac, self.val)
        elif algorithm == 'cantor-zassenhaus':
            lead = nmod_poly_factor_with_cantor_zassenhaus(fac, self.val)
        else:
            lead = nmod_poly_factor(fac, self.val)
        res = [None] * fac.num
        for 0 <= i < fac.num:
            u = nmod_poly.__new__(nmod_poly)
            nmod_poly_init_preinv((<nmod_poly>u).val,
                (<nmod_poly>self).val.mod.n, (<nmod_poly>self).val.mod.ninv)
            nmod_poly_set((<nmod_poly>u).val, &fac.p[i])
            exp = fac.exp[i]
            res[i] = (u, exp)
        c = nmod.__new__(nmod)
        (<nmod>c).mod = self.val.mod
        (<nmod>c).val = lead
        nmod_poly_factor_clear(fac)
        return c, res


#----------------------------------------------------------------------------#
#                                                                            #
#                               nmod_mat                                     #
#                                                                            #
#----------------------------------------------------------------------------#

cdef class nmod_mat:
    """
    The nmod_poly type represents dense matrices over Z/nZ for
    word-size n. Some operations may assume that n is a prime.
    """

    cdef nmod_mat_t val

    def __dealloc__(self):
        nmod_mat_clear(self.val)

    @cython.embedsignature(False)
    def __init__(self, *args):
        cdef long m, n, i, j
        cdef mp_limb_t mod
        if len(args) == 1:
            val = args[0]
            if typecheck(val, nmod_mat):
                nmod_mat_init_set(self.val, (<nmod_mat>val).val)
                return
        mod = args[-1]
        args = args[:-1]
        if mod == 0:
            raise ValueError("modulus must be nonzero")
        if len(args) == 1:
            val = args[0]
            if typecheck(val, fmpz_mat):
                nmod_mat_init(self.val, fmpz_mat_nrows((<fmpz_mat>val).val),
                    fmpz_mat_ncols((<fmpz_mat>val).val), mod)
                fmpz_mat_get_nmod_mat(self.val, (<fmpz_mat>val).val)
            else:
                raise TypeError("cannot create nmod_mat from input of type %s" % type(val))
        elif len(args) == 2:
            m, n = args
            nmod_mat_init(self.val, m, n, mod)
        elif len(args) == 3:
            m, n, entries = args
            nmod_mat_init(self.val, m, n, mod)
            entries = list(entries)
            if len(entries) != m*n:
                raise ValueError("list of entries has the wrong length")
            for i from 0 <= i < m:
                for j from 0 <= j < n:
                    x = nmod(entries[i*n + j], mod)         # XXX: slow
                    self.val.rows[i][j] = (<nmod>x).val
        else:
            raise ValueError("nmod_mat: expected 1-3 arguments plus modulus")

    def __nonzero__(self):
        return not nmod_mat_is_zero(self.val)

    def __richcmp__(s, t, int op):
        cdef mp_limb_t v
        cdef bint res
        if op != 2 and op != 3:
            raise TypeError("matrices cannot be ordered")
        if typecheck(s, nmod_mat) and typecheck(t, nmod_mat):
            if (<nmod_mat>s).val.mod.n != (<nmod_mat>t).val.mod.n:
                res = False
            else:
                res = nmod_mat_equal((<nmod_mat>s).val, (<nmod_mat>t).val)
            if op == 2:
                return res
            if op == 3:
                return not res
        return NotImplemented

    cpdef long nrows(self):
        return nmod_mat_nrows(self.val)

    cpdef long ncols(self):
        return nmod_mat_ncols(self.val)

    cpdef mp_limb_t modulus(self):
        return self.val.mod.n

    @classmethod
    def randtest(cls, ulong m, ulong n, ulong mod):
        """
        Returns a random (m, n) matrix with non-uniformly chosen
        entries. Zero entries are generated with increased probability.
        """
        cdef nmod_mat mat = nmod_mat(m, n, mod)
        nmod_mat_randtest(mat.val, global_random_state)
        return mat

    def __repr__(self):
        return "nmod_mat(%i, %i, %s, %i)" % (self.nrows(), self.ncols(),
            map(int, self.entries()), self.modulus())

    def __str__(self):
        return matrix_to_str(self.table())

    def entries(self):
        cdef long i, j, m, n
        cdef nmod t
        m = self.nrows()
        n = self.ncols()
        L = [None] * (m * n)
        for i from 0 <= i < m:
            for j from 0 <= j < n:
                # XXX
                t = nmod(nmod_mat_entry(self.val, i, j), self.val.mod.n)
                L[i*n + j] = t
        return L

    def table(self):
        cdef long i, m, n
        m = self.nrows()
        n = self.ncols()
        L = self.entries()
        return [L[i*n:(i+1)*n] for i in range(m)]

    def __getitem__(self, index):
        cdef long i, j
        cdef nmod x
        i, j = index
        if i < 0 or i >= self.nrows() or j < 0 or j >= self.ncols():
            raise ValueError("index %i,%i exceeds matrix dimensions" % (i, j))
        x = nmod(nmod_mat_entry(self.val, i, j), self.modulus()) # XXX: slow
        return x

    def __setitem__(self, index, value):
        cdef long i, j
        cdef mp_limb_t v
        i, j = index
        if i < 0 or i >= self.nrows() or j < 0 or j >= self.ncols():
            raise ValueError("index %i,%i exceeds matrix dimensions" % (i, j))
        if any_as_nmod(&v, value, self.val.mod):
            self.val.rows[i][j] = v
        else:
            raise ValueError("cannot set item of type %s" % type(value))

    def det(self):
        """
        Returns the determinant of self as an nmod.

            >>> nmod_mat(2,2,[1,2,3,4],17).det()
            nmod(15, 17)

        """
        if not nmod_mat_is_square(self.val):
            raise ValueError("matrix must be square")
        return nmod(nmod_mat_det(self.val), self.modulus())

    def rank(self):
        return nmod_mat_rank(self.val)

    def __pos__(self):
        return self

    def __neg__(self):
        cdef nmod_mat r = nmod_mat(self)   # XXX
        nmod_mat_neg(r.val, r.val)
        return r

    def __add__(s, t):
        cdef nmod_mat r
        cdef nmod_mat_struct *sv
        cdef nmod_mat_struct *tv
        if typecheck(s, nmod_mat):
            sv = &(<nmod_mat>s).val[0]
            t = any_as_nmod_mat(t, sv.mod)
            if t is NotImplemented:
                return t
            tv = &(<nmod_mat>t).val[0]
        else:
            tv = &(<nmod_mat>t).val[0]
            s = any_as_nmod_mat(s, tv.mod)
            if s is NotImplemented:
                return s
            sv = &(<nmod_mat>s).val[0]
        if sv.mod.n != tv.mod.n:
            raise ValueError("cannot add nmod_mats with different moduli")
        if sv.r != tv.r or sv.c != tv.c:
            raise ValueError("incompatible shapes for matrix addition")
        r = nmod_mat.__new__(nmod_mat)
        nmod_mat_init(r.val, sv.r, sv.c, sv.mod.n)
        nmod_mat_add(r.val, sv, tv)
        return r

    def __sub__(s, t):
        cdef nmod_mat r
        cdef nmod_mat_struct *sv
        cdef nmod_mat_struct *tv
        if typecheck(s, nmod_mat):
            sv = &(<nmod_mat>s).val[0]
            t = any_as_nmod_mat(t, sv.mod)
            if t is NotImplemented:
                return t
            tv = &(<nmod_mat>t).val[0]
        else:
            tv = &(<nmod_mat>t).val[0]
            s = any_as_nmod_mat(s, tv.mod)
            if s is NotImplemented:
                return s
            sv = &(<nmod_mat>s).val[0]
        if sv.mod.n != tv.mod.n:
            raise ValueError("cannot subtract nmod_mats with different moduli")
        if sv.r != tv.r or sv.c != tv.c:
            raise ValueError("incompatible shapes for matrix subtraction")
        r = nmod_mat.__new__(nmod_mat)
        nmod_mat_init(r.val, sv.r, sv.c, sv.mod.n)
        nmod_mat_sub(r.val, sv, tv)
        return r

    cdef __mul_nmod(self, mp_limb_t c):
        cdef nmod_mat r = nmod_mat.__new__(nmod_mat)
        nmod_mat_init(r.val, self.val.r, self.val.c, self.val.mod.n)
        nmod_mat_scalar_mul(r.val, self.val, c)
        return r

    def __mul__(s, t):
        cdef nmod_mat r
        cdef nmod_mat_struct *sv
        cdef nmod_mat_struct *tv
        cdef mp_limb_t c
        if typecheck(s, nmod_mat):
            sv = &(<nmod_mat>s).val[0]
            u = any_as_nmod_mat(t, sv.mod)
            if u is NotImplemented:
                if any_as_nmod(&c, t, sv.mod):
                    return (<nmod_mat>s).__mul_nmod(c)
                return NotImplemented
            tv = &(<nmod_mat>u).val[0]
        else:
            tv = &(<nmod_mat>t).val[0]
            u = any_as_nmod_mat(s, tv.mod)
            if u is NotImplemented:
                if any_as_nmod(&c, s, tv.mod):
                    return (<nmod_mat>t).__mul_nmod(c)
                return NotImplemented
            sv = &(<nmod_mat>u).val[0]
        if sv.mod.n != tv.mod.n:
            raise ValueError("cannot multiply nmod_mats with different moduli")
        if sv.c != tv.r:
            raise ValueError("incompatible shapes for matrix multiplication")
        r = nmod_mat.__new__(nmod_mat)
        nmod_mat_init(r.val, sv.r, tv.c, sv.mod.n)
        nmod_mat_mul(r.val, sv, tv)
        return r

    def __div__(nmod_mat s, t):
        cdef mp_limb_t v
        if not any_as_nmod(&v, t, s.val.mod):
            return NotImplemented
        t = nmod(v, s.val.mod.n)
        return s * (~t)

    # __truediv__ = __div__ doesn't seem to work?
    def __truediv__(nmod_mat s, t):
        return nmod_mat.__div__(s, t)

    def __invert__(self):
        cdef nmod_mat u
        if not nmod_mat_is_square(self.val):
            raise ValueError("matrix must be square")
        u = nmod_mat.__new__(nmod_mat)
        nmod_mat_init(u.val, nmod_mat_nrows(self.val),
            nmod_mat_ncols(self.val), self.val.mod.n)
        if not nmod_mat_inv(u.val, self.val):
            raise ZeroDivisionError("matrix is singular")
        return u

    def transpose(self):
        """
        Returns the transpose of self.

            >>> nmod_mat(2,3,range(6),7).transpose()
            nmod_mat(3, 2, [0L, 3L, 1L, 4L, 2L, 5L], 7)
        """
        cdef nmod_mat u
        u = nmod_mat.__new__(nmod_mat)
        nmod_mat_init(u.val, nmod_mat_ncols(self.val),
            nmod_mat_nrows(self.val), self.val.mod.n)
        nmod_mat_transpose(u.val, self.val)
        return u

    def solve(self, other):
        """
        Given self = A and other = B, returns a matrix X such
        that A*X = B, assuming that self is square and invertible.

            >>> A = nmod_mat(2, 2, [1,4,8,3], 11)
            >>> B = nmod_mat(2, 3, range(6), 11)
            >>> X = A.solve(B)
            >>> print X
            [8,  5, 2]
            [9, 10, 0]
            >>> A*X == B
            True
            >>> nmod_mat(2, 2, [1,0,2,0], 11).solve(B)
            Traceback (most recent call last):
              ...
            ZeroDivisionError: singular matrix in solve()
            >>> A.solve(nmod_mat(1, 2, [2,3], 11))
            Traceback (most recent call last):
              ...
            ValueError: need a square system and compatible right hand side

        """
        cdef nmod_mat u
        cdef int result
        t = any_as_nmod_mat(other, self.val.mod)
        if t is NotImplemented:
            raise TypeError("cannot convert input to nmod_mat")
        if (nmod_mat_nrows(self.val) != nmod_mat_ncols(self.val) or
            nmod_mat_nrows(self.val) != nmod_mat_nrows((<nmod_mat>t).val)):
            raise ValueError("need a square system and compatible right hand side")
        u = nmod_mat.__new__(nmod_mat)
        nmod_mat_init(u.val, nmod_mat_nrows((<nmod_mat>t).val),
            nmod_mat_ncols((<nmod_mat>t).val), self.val.mod.n)
        result = nmod_mat_solve(u.val, self.val, (<nmod_mat>t).val)
        if not result:
            raise ZeroDivisionError("singular matrix in solve()")
        return u

    def rref(self, inplace=False):
        """
        Computes the reduced row echelon form (rref) of self,
        either returning a new copy or modifying self in-place.
        Returns (rref, rank).

            >>> A = nmod_mat(3,3,range(9),23)
            >>> A.rref()
            (nmod_mat(3, 3, [1L, 0L, 22L, 0L, 1L, 2L, 0L, 0L, 0L], 23), 2)
            >>> A.rref(inplace=True)
            (nmod_mat(3, 3, [1L, 0L, 22L, 0L, 1L, 2L, 0L, 0L, 0L], 23), 2)
            >>> A
            nmod_mat(3, 3, [1L, 0L, 22L, 0L, 1L, 2L, 0L, 0L, 0L], 23)

        """
        if inplace:
            res = self
        else:
            res = nmod_mat.__new__(nmod_mat)
            nmod_mat_init_set((<nmod_mat>res).val, self.val)
        rank = nmod_mat_rref((<nmod_mat>res).val)
        return res, rank

    def nullspace(self):
        """
        Computes a basis for the nullspace of self. Returns (X, nullity)
        where nullity is the rank of the nullspace of self and X is a
        matrix whose first (nullity) columns are linearly independent,
        and such that self * X = 0.

            >>> A = nmod_mat(3,5,range(1,16),23)
            >>> X, nullity = A.nullspace()
            >>> A.rank(), nullity, X.rank()
            (2, 3, 3)
            >>> A * X
            nmod_mat(3, 5, [0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L], 23)
            >>> print(X)
            [ 1,  2,  3, 0, 0]
            [21, 20, 19, 0, 0]
            [ 1,  0,  0, 0, 0]
            [ 0,  1,  0, 0, 0]
            [ 0,  0,  1, 0, 0]

        """
        cdef nmod_mat res
        res = nmod_mat.__new__(nmod_mat)
        nmod_mat_init(res.val, nmod_mat_ncols(self.val), nmod_mat_ncols(self.val), self.val.mod.n)
        nullity = nmod_mat_nullspace(res.val, self.val)
        return res, nullity


#----------------------------------------------------------------------------#
#                                                                            #
#                           Special functions                                #
#                                                                            #
#----------------------------------------------------------------------------#

def number_of_partitions(n):
    """
    Returns p(n), the number of partitions of n, as an fmpz.

        >>> [number_of_partitions(n) for n in range(8)]
        [fmpz(1), fmpz(1), fmpz(2), fmpz(3), fmpz(5), fmpz(7), fmpz(11), fmpz(15)]
        >>> number_of_partitions(100)
        fmpz(190569292)
        >>> len(str(number_of_partitions(10**9)))
        35219
    """
    cdef fmpz v = fmpz()
    arith_number_of_partitions(v.val, n)
    return v

def moebius_mu(n):
    """
    Evaluates the Moebius function mu(n), returning a Python
    integer.

        >>> [moebius_mu(n) for n in range(10)]
        [0, 1, -1, -1, 0, -1, 1, -1, 0, 0]
    """
    n = any_as_fmpz(n)
    if n is NotImplemented:
        raise TypeError("the input must be an integer")
    return arith_moebius_mu((<fmpz>n).val)

def divisor_sigma(n, k):
    """
    Evaluates the divisor sum function sigma_k(n), returning an fmpz.

        >>> divisor_sigma(60, 0)
        fmpz(12)
        >>> divisor_sigma(60, 1)
        fmpz(168)
        >>> divisor_sigma(60, 10)
        fmpz(605263138639095300)
    """
    cdef fmpz v = fmpz()
    n = any_as_fmpz(n)
    if n is NotImplemented:
        raise TypeError("the input must be an integer")
    arith_divisor_sigma(v.val, (<fmpz>n).val, k)
    return v

def euler_phi(n):
    """
    Evaluates the Euler totient function phi(n), returning an fmpz.

        >>> euler_phi(60)
        fmpz(16)
        >>> euler_phi(3**10)
        fmpz(39366)
    """
    cdef fmpz v = fmpz()
    n = any_as_fmpz(n)
    if n is NotImplemented:
        raise TypeError("the input must be an integer")
    arith_euler_phi(v.val, (<fmpz>n).val)
    return v

def bell_number(n):
    """
    Returns the nth Bell number B_n as an fmpz.

        >>> [bell_number(n) for n in range(5)]
        [fmpz(1), fmpz(1), fmpz(2), fmpz(5), fmpz(15)]
        >>> bell_number(50)
        fmpz(185724268771078270438257767181908917499221852770)
    """
    cdef fmpz v = fmpz()
    arith_bell_number(v.val, n)
    return v

def euler_number(n):
    """
    Returns the nth Euler number E_n as an fmpz.

        >>> [euler_number(n) for n in range(6)]
        [fmpz(1), fmpz(0), fmpz(-1), fmpz(0), fmpz(5), fmpz(0)]
        >>> euler_number(50)
        fmpz(-6053285248188621896314383785111649088103498225146815121)
    """
    cdef fmpz v = fmpz()
    arith_euler_number(v.val, n)
    return v

def bernoulli_number(n):
    """
    Returns the nth Bernoulli number B_n as an fmpq.

        >>> [bernoulli_number(n) for n in range(8)]
        [fmpq(1,1), fmpq(-1,2), fmpq(1,6), fmpq(0,1), fmpq(-1,30), fmpq(0,1), fmpq(1,42), fmpq(0,1)]
        >>> bernoulli_number(50)
        fmpq(495057205241079648212477525,66)
    """
    cdef fmpq v = fmpq()
    arith_bernoulli_number(v.val, n)
    return v

def stirling_number_1(n, k):
    """
    Returns the Stirling number of the first kind s(n, k) as an fmpz.

        >>> [stirling_number_1(5, k) for k in range(6)]
        [fmpz(0), fmpz(24), fmpz(-50), fmpz(35), fmpz(-10), fmpz(1)]
    """
    cdef fmpz v = fmpz()
    arith_stirling_number_1(v.val, n, k)
    return v

def stirling_number_2(n, k):
    """
    Returns the Stirling number of the second kind S(n, k) as an fmpz.

        >>> [stirling_number_2(5, k) for k in range(6)]
        [fmpz(0), fmpz(1), fmpz(15), fmpz(25), fmpz(10), fmpz(1)]
    """
    cdef fmpz v = fmpz()
    arith_stirling_number_2(v.val, n, k)
    return v

def harmonic_number(n):
    """
    Returns the harmonic number H_n as an fmpq.

        >>> [harmonic_number(n) for n in range(6)]
        [fmpq(0,1), fmpq(1,1), fmpq(3,2), fmpq(11,6), fmpq(25,12), fmpq(137,60)]
        >>> harmonic_number(50)
        fmpq(13943237577224054960759,3099044504245996706400)
    """
    cdef fmpq v = fmpq()
    arith_harmonic_number(v.val, n)
    return v

def bernoulli_polynomial(n):
    """
    Returns the nth Bernoulli polynomial B_n(x) as an fmpq_poly.

        >>> bernoulli_polynomial(2)
        fmpq_poly([1, -6, 6], 6)
    """
    cdef fmpq_poly v = fmpq_poly()
    arith_bernoulli_polynomial(v.val, n)
    return v

def euler_polynomial(n):
    """
    Returns the nth Euler polynomial E_n(x) as an fmpq_poly.

        >>> euler_polynomial(3)
        fmpq_poly([1, 0, -6, 4], 4)
    """
    cdef fmpq_poly v = fmpq_poly()
    arith_euler_polynomial(v.val, n)
    return v

def legendre_polynomial(n):
    """
    Returns the nth Legendre polynomial P_n(x) as an fmpq_poly.

        >>> legendre_polynomial(3)
        fmpq_poly([0, -3, 0, 5], 2)

    """
    cdef fmpq_poly v = fmpq_poly()
    arith_legendre_polynomial(v.val, n)
    return v

def chebyshev_t_polynomial(n):
    """
    Returns the nth Chebyshev polynomial of the first kind T_n(x)
    as an fmpz_poly.

        >>> chebyshev_t_polynomial(3)
        fmpz_poly([0, -3, 0, 4])

    """
    cdef fmpz_poly v = fmpz_poly()
    arith_chebyshev_t_polynomial(v.val, n)
    return v

def chebyshev_u_polynomial(n):
    """
    Returns the nth Chebyshev polynomial of the second kind U_n(x)
    as an fmpz_poly.

        >>> chebyshev_u_polynomial(3)
        fmpz_poly([0, -4, 0, 8])

    """
    cdef fmpz_poly v = fmpz_poly()
    arith_chebyshev_u_polynomial(v.val, n)
    return v

def cyclotomic_polynomial(n):
    """
    Returns the nth cyclotomic polynomial Phi_n(x) as an fmpz_poly.

        >>> cyclotomic_polynomial(12)
        fmpz_poly([1, 0, -1, 0, 1])

    """
    cdef fmpz_poly v = fmpz_poly()
    arith_cyclotomic_polynomial(v.val, n)
    return v
