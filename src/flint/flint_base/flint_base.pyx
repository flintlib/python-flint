from flint.flintlib.flint cimport (
    FLINT_BITS as _FLINT_BITS,
    FLINT_VERSION as _FLINT_VERSION,
    __FLINT_RELEASE as _FLINT_RELEASE,
    slong
)
from flint.flintlib.mpoly cimport ordering_t
from flint.flint_base.flint_context cimport thectx
from flint.flint_base.flint_base cimport Ordering
from flint.utils.typecheck cimport typecheck
cimport libc.stdlib

from typing import Optional


FLINT_BITS = _FLINT_BITS
FLINT_VERSION = _FLINT_VERSION.decode("ascii")
FLINT_RELEASE = _FLINT_RELEASE


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
        """
        Returns the coefficients of ``self`` as a list

            >>> from flint import fmpz_poly
            >>> f = fmpz_poly([1,2,3,4,5])
            >>> f.coeffs()
            [1, 2, 3, 4, 5]
        """
        return list(self)

    def str(self, bint ascending=False, *args, **kwargs):
        """
        Convert to a human-readable string (generic implementation for
        all polynomial types).

        If *ascending* is *True*, the monomials are output from low degree to
        high, otherwise from high to low.
        """
        coeffs = [c.str(*args, **kwargs) for c in self]
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
        Computes all the roots in the base ring of the polynomial.
        Returns a list of all pairs (*v*, *m*) where *v* is the 
        integer root and *m* is the multiplicity of the root.

        To compute complex roots of a polynomial, instead use
        the `.complex_roots()` method, which is available on
        certain polynomial rings.

            >>> from flint import fmpz_poly
            >>> fmpz_poly([1, 2]).roots()
            []
            >>> fmpz_poly([2, 1]).roots()
            [(-2, 1)]
            >>> fmpz_poly([12, 7, 1]).roots()
            [(-3, 1), (-4, 1)]
            >>> (fmpz_poly([-5,1]) * fmpz_poly([-5,1]) * fmpz_poly([-3,1])).roots()
            [(3, 1), (5, 2)]
        """
        factor_fn = getattr(self, "factor", None)
        if not callable(factor_fn):
            raise NotImplementedError("Polynomial has no factor method, roots cannot be determined")

        roots = []
        factors = self.factor()
        for fac, m in factors[1]:
            if fac.degree() == fac[1] == 1:
                v = - fac[0]
                roots.append((v, m))
        return roots

    def complex_roots(self):
        raise AttributeError("Complex roots are not supported for this polynomial")


cdef class flint_mpoly_context(flint_elem):
    """
    Base class for multivariate ring contexts
    """

    _ctx_cache = None

    def __init__(self, int nvars, names):
        if nvars < 0:
            raise ValueError("Cannot have a negative amount of variables")
        elif len(names) != nvars:
            raise ValueError("Number of variables must match number of variable names")
        self.py_names = tuple(name.encode("ascii") if not isinstance(name, bytes) else name for name in names)
        self.c_names = <const char**> libc.stdlib.malloc(nvars * sizeof(const char *))
        for i in range(nvars):
            self.c_names[i] = self.py_names[i]

    def __dealloc__(self):
        libc.stdlib.free(self.c_names)
        self.c_names = NULL

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        return f"{self.__class__.__name__}({self.nvars()}, '{repr(self.ordering())}', {self.names()})"

    def name(self, long i):
        if not 0 <= i < len(self.py_names):
            raise IndexError("Variable name index out of range")
        return self.py_names[i].decode("ascii")

    def names(self):
        return tuple(name.decode("ascii") for name in self.py_names)

    def gens(self):
        return tuple(self.gen(i) for i in range(self.nvars()))

    def variable_to_index(self, var: Union[int, str]):
        """Convert a variable name string or possible index to its index in the context."""
        if isinstance(var, str):
            try:
                i = self.names().index(var)
            except ValueError:
                raise ValueError("Variable not in context")
        elif isinstance(var, int):
            if not 0 <= var < self.nvars():
                raise IndexError("Generator index out of range")
            i = var
        else:
            raise TypeError("Invalid variable type")

        return i

    @staticmethod
    def create_variable_names(slong nvars, names: str):
        """
        Create a tuple of variable names based on the comma separated `names` string.

        If `names` contains a single value, and `nvars` > 1, then the variables are numbered, e.g.

            >>> flint_mpoly_context.create_variable_names(3, "x")
            ('x0', 'x1', 'x2')

        """
        nametup = tuple(name.strip() for name in names.split(','))
        if len(nametup) != nvars:
            if len(nametup) == 1:
                nametup = tuple(nametup[0] + str(i) for i in range(nvars))
            else:
                raise ValueError("Number of variables does not equal number of names")
        return nametup

    @classmethod
    def get_context(cls, slong nvars=1, ordering=Ordering.lex, names: Optional[str] = "x", nametup: Optional[tuple] = None):
        """
        Retrieve a context via the number of variables, `nvars`, the ordering, `ordering`, and either a variable
        name string, `names`, or a tuple of variable names, `nametup`.
        """

        # A type hint of `ordering: Ordering` results in the error "TypeError: an integer is required" if a Ordering
        # object is not provided. This is pretty obtuse so we check it's type ourselves
        if not isinstance(ordering, Ordering):
            raise TypeError(f"`ordering` ('{ordering}') is not an instance of flint.Ordering")

        if nametup is not None:
            key = nvars, ordering, nametup
        elif nametup is None and names is not None:
            key = nvars, ordering, cls.create_variable_names(nvars, names)
        else:
            raise ValueError("Must provide either `names` or `nametup`")

        ctx = cls._ctx_cache.get(key)
        if ctx is None:
            ctx = cls._ctx_cache.setdefault(key, cls(*key))
        return ctx

    @classmethod
    def from_context(cls, ctx: flint_mpoly_context):
        return cls.get_context(
            nvars=ctx.nvars(),
            ordering=ctx.ordering(),
            names=None,
            nametup=ctx.names()
        )


cdef class flint_mpoly(flint_elem):
    """
    Base class for multivariate polynomials.
    """

    def leading_coefficient(self):
        return self.coefficient(0)

    def to_dict(self):
        return {self.monomial(i): self.coefficient(i) for i in range(len(self))}

    def __contains__(self, x):
        """
        Returns True if `self` contains a term with exponent vector `x` and a non-zero coefficient.

            >>> from flint import fmpq_mpoly_ctx, Ordering
            >>> ctx = fmpq_mpoly_ctx.get_context(2, Ordering.lex, 'x')
            >>> p = ctx.from_dict({(0, 1): 2, (1, 1): 3})
            >>> (1, 1) in p
            True
            >>> (5, 1) in p
            False

        """
        return bool(self[x])

    def __iter__(self):
        return iter(self.keys())

    def __pos__(self):
        return self

    def terms(self):
        """
        Return the exponent vectors and coefficient of each term.

            >>> from flint import fmpq_mpoly_ctx, Ordering
            >>> ctx = fmpq_mpoly_ctx.get_context(2, Ordering.lex, 'x')
            >>> f = ctx.from_dict({(0, 0): 1, (1, 0): 2, (0, 1): 3, (1, 1): 4})
            >>> list(f.terms())
            [((1, 1), 4), ((1, 0), 2), ((0, 1), 3), ((0, 0), 1)]

        """
        return zip(self.monoms(), self.coeffs())


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


cdef ordering_t ordering_py_to_c(ordering):  # Cython does not like an "Ordering" type hint here
    if not isinstance(ordering, Ordering):
        raise TypeError(f"`ordering` ('{ordering}') is not an instance of flint.Ordering")

    if ordering == Ordering.lex:
        return ordering_t.ORD_LEX
    elif ordering == Ordering.deglex:
        return ordering_t.ORD_DEGLEX
    elif ordering == Ordering.degrevlex:
        return ordering_t.ORD_DEGREVLEX


cdef ordering_c_to_py(ordering_t ordering):
    if ordering == ordering_t.ORD_LEX:
        return Ordering.lex
    elif ordering == ordering_t.ORD_DEGLEX:
        return Ordering.deglex
    elif ordering == ordering_t.ORD_DEGREVLEX:
        return Ordering.degrevlex
    else:
        raise ValueError("Unimplemented term order %d" % ordering)
