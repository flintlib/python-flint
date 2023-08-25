from flint.flint_base.flint_context cimport thectx
cimport libc.stdlib

cdef class flint_elem:
    def __repr__(self):
        if thectx.pretty:
            return self.str()
        else:
            return self.repr()

    def __str__(self):
        return self.str()


cdef class flint_scalar(flint_elem):
    pass


cdef class flint_poly(flint_elem):
    """
    Base class for polynomials.
    """

    def __iter__(self):
        cdef long i, n
        n = self.length()
        for i in range(n):
            yield self[i]

    def coeffs(self):
        return list(self)

    def str(self, bint ascending=False):
        """
        Convert to a human-readable string (generic implementation for
        all polynomial types).

        If *ascending* is *True*, the monomials are output from low degree to
        high, otherwise from high to low.
        """
        coeffs = [str(c) for c in self]
        if not coeffs:
            return "0"
        s = []
        coeffs = enumerate(coeffs)
        if not ascending:
            coeffs = reversed(list(coeffs))
        for i, c in coeffs:
            if c == "0":
                continue
            else:
                if c.startswith("-") or (" " in c):
                    c = "(" + c + ")"
                if i == 0:
                    s.append("%s" % c)
                elif i == 1:
                    if c == "1":
                        s.append("x")
                    else:
                        s.append("%s*x" % c)
                else:
                    if c == "1":
                        s.append("x^%s" % i)
                    else:
                        s.append("%s*x^%s" % (c, i))
        return " + ".join(s)

    def roots(self):
        """
        Deprecated function.

        To recover roots of a polynomial, first convert to acb:

        acb_poly(input_poly).roots()
        """
        raise NotImplementedError('This method is no longer supported. To recover the complex roots first convert to acb_poly')

cdef class flint_mpoly_context(flint_elem):
    """
    Base class for multivariate ring contexts
    """

    def __cinit__(self):
       self._init = False

    def __init__(self, long nvars, names):
        assert nvars >= 1
        assert len(names) == nvars
        self.py_names = tuple(bytes(name, 'utf-8') for name in names)
        self.c_names = <char**>libc.stdlib.malloc(nvars * sizeof(char *))
        self._init = True
        for i in range(nvars):
            self.c_names[i] = self.py_names[i]

    def __dealloc__(self):
        if self._init:
            libc.stdlib.free(self.c_names)
        self._init = False

    def name(self, long i):
        assert i >= 0 and i < len(self.py_names)
        return self.py_names[i].decode('utf-8')

    def gens(self):
        return tuple(self.gen(i) for i in range(self.nvars()))



cdef class flint_mpoly(flint_elem):
    """
    Base class for multivariate polynomials.
    """

    def leading_coefficient(self):
        return self.coefficient(0)

    def __hash__(self):
        s = repr(self)
        return hash(s)


cdef class flint_series(flint_elem):
    """
    Base class for power series.
    """
    def __iter__(self):
        cdef long i, n
        n = self.length()
        for i in range(n):
            yield self[i]

    def coeffs(self):
        return list(self)


cdef class flint_mat(flint_elem):
    """
    Base class for matrices.
    """

    def repr(self):
        if thectx.pretty:
            return str(self)
        # XXX
        return "%s(%i, %i, [%s])" % (type(self).__name__,
            self.nrows(), self.ncols(), (", ".join(map(str, self.entries()))))

    def str(self, *args, **kwargs):
        tab = self.table()
        if len(tab) == 0 or len(tab[0]) == 0:
            return "[]"
        tab = [[r.str(*args, **kwargs) for r in row] for row in tab]
        widths = []
        for i in xrange(len(tab[0])):
            w = max([len(row[i]) for row in tab])
            widths.append(w)
        for i in xrange(len(tab)):
            tab[i] = [s.rjust(widths[j]) for j, s in enumerate(tab[i])]
            tab[i] = "[" + (", ".join(tab[i])) + "]"
        return "\n".join(tab)

    def entries(self):
        cdef long i, j, m, n
        m = self.nrows()
        n = self.ncols()
        L = [None] * (m * n)
        for i from 0 <= i < m:
            for j from 0 <= j < n:
                L[i*n + j] = self[i, j]
        return L

    def __iter__(self):
        cdef long i, j, m, n
        m = self.nrows()
        n = self.ncols()
        for i from 0 <= i < m:
            for j from 0 <= j < n:
                yield self[i, j]

    def table(self):
        cdef long i, m, n
        m = self.nrows()
        n = self.ncols()
        L = self.entries()
        return [L[i*n : (i+1)*n] for i in range(m)]

    # supports mpmath conversions
    tolist = table
