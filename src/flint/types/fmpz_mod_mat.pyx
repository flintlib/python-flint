from flint.utils.typecheck cimport (
    typecheck,
)
from flint.flintlib.functions.fmpz cimport (
    fmpz_struct,
    fmpz_t,
    fmpz_init_set,
    fmpz_set,
)
from flint.flintlib.types.fmpz_mod_mat_compat cimport (
    compat_fmpz_mod_mat_init,
    compat_fmpz_mod_mat_init_set,
    compat_fmpz_mod_mat_clear,
    compat_fmpz_mod_mat_set,
    compat_fmpz_mod_mat_nrows,
    compat_fmpz_mod_mat_ncols,
    compat_fmpz_mod_mat_entry,
    compat_fmpz_mod_mat_set_entry,
    compat_fmpz_mod_mat_one,
    compat_fmpz_mod_mat_equal,
    compat_fmpz_mod_mat_is_zero,
    compat_fmpz_mod_mat_neg,
    compat_fmpz_mod_mat_add,
    compat_fmpz_mod_mat_sub,
    compat_fmpz_mod_mat_mul,
    compat_fmpz_mod_mat_scalar_mul_fmpz,
    compat_fmpz_mod_mat_inv,
    compat_fmpz_mod_mat_transpose,
    compat_fmpz_mod_mat_solve,
    compat_fmpz_mod_mat_rref,
    compat_fmpz_mod_mat_charpoly,
    compat_fmpz_mod_mat_minpoly,
)

from flint.flint_base.flint_base cimport (
    flint_mat,
)
from flint.types.fmpz cimport (
    fmpz,
)
from flint.types.fmpz_mod cimport (
    fmpz_mod_ctx,
    fmpz_mod,
)
from flint.types.fmpz_mod_poly cimport (
    fmpz_mod_poly,
    fmpz_mod_poly_ctx,
)
from flint.types.fmpz_mat cimport (
    fmpz_mat,
)
from flint.types.nmod_mat cimport (
    nmod_mat,
)


cdef any_as_fmpz_mod_mat(x):
    if typecheck(x, fmpz_mod_mat):
        return x
    elif typecheck(x, nmod_mat):
        return fmpz_mod_mat(x)
    else:
        return NotImplemented


cdef any_as_fmpz_mod_mat_ctx(x, fmpz_mod_ctx ctx):
    y = any_as_fmpz_mod_mat(x)
    if y is not NotImplemented:
        return y
    if typecheck(x, fmpz_mat):
        return fmpz_mod_mat(x, ctx)
    return NotImplemented


cdef class fmpz_mod_mat(flint_mat):
    """
    The ``fmpz_mod_mat`` type represents dense matrices over ``Z/nZ`` for
    arbitrary ``n`` (not necessarily word-sized unlike ``nmod_mat``). Some
    operations may assume that n is a prime.
    """
    def __dealloc__(self):
        if self._initialized:
            compat_fmpz_mod_mat_clear(self.val, self.ctx.val)

    def __init__(self, *args):
        """Construct an ``fmpz_mod_mat`` matrix.

        >>> from flint import fmpz_mod_mat, fmpz_mod_ctx, fmpz_mat
        >>> ctx = fmpz_mod_ctx(7)
        >>> fmpz_mod_mat(2, 2, [1, 2, 3, 4], ctx)
        [1, 2]
        [3, 4]
        >>> fmpz_mod_mat([[1, 2], [3, 4]], ctx)
        [1, 2]
        [3, 4]
        >>> fmpz_mod_mat(2, 2, ctx)
        [0, 0]
        [0, 0]
        >>> zmat = fmpz_mat(2, 2, [1, 2, 3, 4])
        >>> fmpz_mod_mat(zmat, ctx)
        [1, 2]
        [3, 4]
        """
        # XXX: This logic should be moved to flint_mat.
        if len(args) >= 2 and isinstance(args[0], int) and isinstance(args[1], int):
            if len(args) >= 3 and isinstance(args[2], (list, tuple)):
                m, n, entries, *arg = args
                if len(entries) != m * n:
                    raise ValueError("fmpz_mod_mat: invalid number of entries")
                self._init_from_list(m, n, entries, arg)
            else:
                m, n, *arg = args
                self._init_empty(m, n, arg)

        elif (len(args) >= 1 and isinstance(args[0], (list, tuple))
              and all(isinstance(arg, (list, tuple)) for arg in args[0])):
            lol, *arg = args
            m = len(lol)
            if m == 0:
                n = 0
                entries = []
            else:
                n = len(lol[0])
                if any(len(row) != n for row in lol):
                    raise ValueError("fmpz_mod_mat: inconsistent row lengths")
                entries = [x for row in lol for x in row]
            self._init_from_list(m, n, entries, arg)

        # XXX: nmod_mat is not a subclass of flint_mat (it should be)
        elif len(args) >= 1 and isinstance(args[0], (flint_mat, nmod_mat)):
            M, *arg = args
            self._init_from_matrix(M, arg)

        else:
            raise TypeError("fmpz_mod_mat: invalid arguments")

    cdef fmpz_mod_ctx _parse_args(self, list args):
        """Parse the context argument."""
        if len(args) != 1:
            if len(args) == 0:
                raise TypeError("fmpz_mod_mat: missing modulus argument")
            else:
                raise TypeError("fmpz_mod_mat: too many arguments")
        arg = args[0]
        if not typecheck(arg, fmpz_mod_ctx):
            raise TypeError("fmpz_mod_mat: invalid modulus argument")
        return arg

    cdef void _init_empty_ctx(self, slong m, slong n, fmpz_mod_ctx ctx):
        """Initialize an empty matrix with a given modulus context."""
        self.ctx = ctx
        compat_fmpz_mod_mat_init(self.val, m, n, ctx.val)
        self._initialized = True

    cdef void _init_empty(self, slong m, slong n, list args):
        """Initialize an empty matrix using the constructor arguments."""
        cdef fmpz_mod_ctx ctx
        ctx = self._parse_args(args)
        self._init_empty_ctx(m, n, ctx)

    cdef void _init_from_list(self, slong m, slong n, list entries, list args):
        """Initialize a matrix from a list of entries."""
        cdef fmpz_mod_ctx ctx
        cdef fmpz_mod e
        ctx = self._parse_args(args)
        entries = [ctx.any_as_fmpz_mod(x) for x in entries]

        self._init_empty_ctx(m, n, ctx)
        for i in range(m):
            for j in range(n):
                val = ctx.any_as_fmpz_mod(entries[i*n + j])
                if val is NotImplemented:
                    raise TypeError("fmpz_mod_mat: incompatible entries")
                e = <fmpz_mod> val
                compat_fmpz_mod_mat_set_entry(self.val, i, j, e.val, ctx.val)

    # XXX: Should be possible to type this as flint_mat but nmod_mat is not a subclass
    # cdef void _init_from_matrix(self, flint_mat M, list args):
    cdef void _init_from_matrix(self, M, list args):
        """Initialize from another matrix."""
        cdef fmpz_mod_ctx ctx
        cdef fmpz_mod_mat N1

        if typecheck(M, fmpz_mod_mat):
            N1 = M
            if args:
                ctx = self._parse_args(args)
            else:
                ctx = N1.ctx
            if N1.ctx != ctx:
                raise TypeError("fmpz_mod_mat: incompatible moduli")
            compat_fmpz_mod_mat_init_set(self.val, N1.val, ctx.val)
            self.ctx = ctx
            self._initialized = True
        elif typecheck(M, fmpz_mat):
            # XXX: This is inefficient.
            self._init_from_list(M.nrows(), M.ncols(), M.entries(), args)
        elif typecheck(M, nmod_mat):
            m = M.modulus()
            if args:
                ctx = self._parse_args(args)
                if m != ctx.modulus():
                    raise TypeError("fmpz_mod_mat: incompatible moduli")
            else:
                ctx = fmpz_mod_ctx(m)
            # XXX: This is inefficient.
            entries = [int(x) for x in M.entries()]
            self._init_from_list(M.nrows(), M.ncols(), entries, [ctx])
        else:
            raise TypeError("fmpz_mod_mat: unrocognized matrix type")

    cdef fmpz_mod_mat _new(self, slong m, slong n, fmpz_mod_ctx ctx):
        """Create an initialized matrix of given shape and context."""
        cdef fmpz_mod_mat res
        res = fmpz_mod_mat.__new__(fmpz_mod_mat)
        compat_fmpz_mod_mat_init(res.val, m, n, ctx.val)
        res.ctx = ctx
        res._initialized = True
        return res

    cdef fmpz_mod_mat _newlike(self):
        """Create an uninitialized matrix of the same shape and context."""
        return self._new(self.nrows(), self.ncols(), self.ctx)

    cdef fmpz_mod_mat _copy(self):
        """Create a copy of the matrix."""
        cdef fmpz_mod_mat res
        res = self._newlike()
        compat_fmpz_mod_mat_set(res.val, self.val, self.ctx.val)
        return res

    cpdef slong nrows(self):
        """Return the number of rows."""
        return compat_fmpz_mod_mat_nrows(self.val, self.ctx.val)

    cpdef slong ncols(self):
        """Return the number of columns."""
        return compat_fmpz_mod_mat_ncols(self.val, self.ctx.val)

    def modulus(self):
        """Return the modulus."""
        cdef fmpz mod
        mod = fmpz.__new__(fmpz)
        fmpz_init_set(mod.val, self.ctx.val.n)
        return mod

    cdef fmpz_mod _getitem(self, slong i, slong j):
        """Retrieve an element as an ``fmpz_mod``."""
        cdef fmpz_struct * p_e
        cdef fmpz_mod e
        p_e = compat_fmpz_mod_mat_entry(self.val, i, j)
        e = fmpz_mod.__new__(fmpz_mod)
        fmpz_set(e.val, p_e)
        e.ctx = self.ctx
        return e

    cdef void _setitem(self, slong i, slong j, fmpz_t e):
        """Set an element from a raw ``fmpz_t``."""
        compat_fmpz_mod_mat_set_entry(self.val, i, j, e, self.ctx.val)

    def repr(self):
        """Return a representation string."""
        # XXX: Generic repr does not include modulus argument.
        m = self.nrows()
        n = self.ncols()
        mod = self.modulus()
        entries = (", ".join([x.repr() for x in self.entries()]))
        return f"fmpz_mod_mat({m}, {n}, [{entries}], {mod})"

    def entries(self):
        """Return a list of entries."""
        # XXX: Ideally move this to flint_mat
        cdef slong i, j, m, n
        m = self.nrows()
        n = self.ncols()
        L = [None] * (m * n)
        for i from 0 <= i < m:
            for j from 0 <= j < n:
                L[i*n + j] = self._getitem(i, j)
        return L

    def __getitem__(self, index):
        """Retrieve an element."""
        # XXX: Ideally move this to flint_mat
        cdef slong i, j
        i, j = index
        if i < 0 or i >= self.nrows() or j < 0 or j >= self.ncols():
            raise IndexError("index %i,%i exceeds matrix dimensions" % (i, j))
        return self._getitem(i, j)

    def __setitem__(self, index, value):
        """Set an element."""
        cdef slong i, j
        cdef fmpz_mod e
        if not (isinstance(index, tuple) and len(index) == 2):
            raise TypeError("fmpz_mod_mat indices must be 2-tuples")
        i, j = index
        if i < 0 or i >= self.nrows() or j < 0 or j >= self.ncols():
            raise IndexError("index %i,%i exceeds matrix dimensions" % (i, j))
        e = self.ctx.any_as_fmpz_mod(value)
        self._setitem(i, j, e.val)

    def __bool__(self):
        """Return ``True`` if the matrix has any nonzero entries."""
        cdef bint zero
        zero = compat_fmpz_mod_mat_is_zero(self.val, self.ctx.val)
        return not zero

    def __richcmp__(self, other, int op):
        """Compare two matrices."""
        cdef bint res

        if op != 2 and op != 3:
            raise TypeError("matrices cannot be ordered")

        other = any_as_fmpz_mod_mat_ctx(other, self.ctx)
        if other is NotImplemented:
            return other

        if (<fmpz_mod_mat>self).ctx != (<fmpz_mod_mat>other).ctx:
            res = 0
        else:
            res = compat_fmpz_mod_mat_equal((<fmpz_mod_mat>self).val, (<fmpz_mod_mat>other).val, self.ctx.val)

        if op == 2:
            return res
        if op == 3:
            return not res

    def __pos__(self):
        """``+M``: Return self (not a copy)."""
        return self

    def __neg__(self):
        """``-M``"""
        res = self._newlike()
        compat_fmpz_mod_mat_neg((<fmpz_mod_mat> res).val, self.val, self.ctx.val)
        return res

    def _as_fmpz_mod_mat(self, other, same_shape=True):
        """Convert to ``fmpz_mod_mat`` but raise if shape or modulus do not match."""
        cdef fmpz_mod_mat o
        other = any_as_fmpz_mod_mat_ctx(other, self.ctx)
        if other is NotImplemented:
            return NotImplemented
        o = other
        if same_shape and not (self.nrows() == o.nrows() and self.ncols() == o.ncols()):
            raise ValueError("Shape mismatch: cannot add matrices")
        if self.ctx != o.ctx:
            raise ValueError("fmpz_mod_mat: incompatible moduli")
        return o

    def _add(self, fmpz_mod_mat other):
        """Add two ``fmpz_mod_mat`` matrices."""
        res = self._newlike()
        compat_fmpz_mod_mat_add(res.val, self.val, other.val, self.ctx.val)
        return res

    def _sub(self, fmpz_mod_mat other):
        """Subtract two ``fmpz_mod_mat`` matrices."""
        res = self._newlike()
        compat_fmpz_mod_mat_sub(res.val, self.val, other.val, self.ctx.val)
        return res

    def _matmul(self, fmpz_mod_mat other):
        """Multiply two ``fmpz_mod_mat`` matrices."""
        if self.ncols() != other.nrows():
            raise ValueError("Shape mismatch: cannot multiply matrices")
        res = self._new(self.nrows(), other.ncols(), self.ctx)
        compat_fmpz_mod_mat_mul(res.val, self.val, other.val, self.ctx.val)
        return res

    def _scalarmul(self, fmpz_mod other):
        """Multiply an ``fmpz_mod_mat`` matrix by an ``fmpz_mod`` scalar."""
        res = self._newlike()
        compat_fmpz_mod_mat_scalar_mul_fmpz(res.val, self.val, other.val, self.ctx.val)
        return res

    def _pow(self, slong other):
        """Raise an ``fmpz_mod_mat`` matrix to an integer power."""
        cdef fmpz_mod_mat res, tmp

        if self.nrows() != self.ncols():
            raise ValueError("fmpz_mod_mat pow: matrix must be square")
        if other < 0:
            self = self.inv()
            other = -other

        res = self._newlike()
        compat_fmpz_mod_mat_one(res.val, self.ctx.val)

        tmp = self._copy()

        while other > 0:
            if other % 2 == 1:
                compat_fmpz_mod_mat_mul(res.val, res.val, tmp.val, self.ctx.val)
            compat_fmpz_mod_mat_mul(tmp.val, tmp.val, tmp.val, self.ctx.val)
            other //= 2

        return res

    def _div(self, fmpz_mod other):
        """Divide an ``fmpz_mod_mat`` matrix by an ``fmpz_mod`` scalar."""
        return self._scalarmul(other.inverse())

    def __add__(self, other):
        """``M + N``: Add two matrices."""
        other = self._as_fmpz_mod_mat(other)
        if other is NotImplemented:
            return NotImplemented
        return self._add(other)

    def __radd__(self, other):
        """``N + M``: Add two matrices."""
        other = self._as_fmpz_mod_mat(other)
        if other is NotImplemented:
            return NotImplemented
        return other._add(self)

    def __sub__(self, other):
        """``M - N``: Subtract two matrices."""
        other = self._as_fmpz_mod_mat(other)
        if other is NotImplemented:
            return NotImplemented
        return self._sub(other)

    def __rsub__(self, other):
        """``N - M``: Subtract two matrices."""
        other = self._as_fmpz_mod_mat(other)
        if other is NotImplemented:
            return NotImplemented
        return other._sub(self)

    def __mul__(self, other):
        """``M * N``: Multiply two matrices."""
        cdef fmpz_mod e
        other_mat = self._as_fmpz_mod_mat(other, same_shape=False)
        if other_mat is not NotImplemented:
            return self._matmul(other_mat)
        other_scalar = self.ctx.any_as_fmpz_mod(other)
        if other_scalar is not NotImplemented:
            e = other_scalar
            return self._scalarmul(e)
        return NotImplemented

    def __rmul__(self, other):
        """``N * M``: Multiply two matrices."""
        cdef fmpz_mod e
        other_mat = self._as_fmpz_mod_mat(other)
        if other_mat is not NotImplemented:
            return other_mat._matmul(self)
        other_scalar = self.ctx.any_as_fmpz_mod(other)
        if other_scalar is not NotImplemented:
            e = other_scalar
            return self._scalarmul(e)
        return NotImplemented

    def __pow__(self, other):
        """``M ** n``: Raise a matrix to an integer power."""
        if not isinstance(other, int):
            return NotImplemented
        return self._pow(other)

    def __truediv__(self, other):
        """``M / n``: Divide a matrix by a scalar."""
        cdef fmpz_mod e
        other_scalar = self.ctx.any_as_fmpz_mod(other)
        if other_scalar is not NotImplemented:
            e = other_scalar
            return self._div(e)
        return NotImplemented

    def inv(self):
        """Return the inverse of a matrix.

        >>> from flint import fmpz_mod_mat, fmpz_mod_ctx
        >>> ctx = fmpz_mod_ctx(7)
        >>> M = fmpz_mod_mat([[1, 2], [3, 4]], ctx)
        >>> M.inv()
        [5, 1]
        [5, 3]

        Assumes that the modulus is prime.
        """
        cdef fmpz_mod_mat res
        if self.nrows() != self.ncols():
            raise ValueError("fmpz_mod_mat inv: matrix must be square")
        res = self._newlike()
        r = compat_fmpz_mod_mat_inv(res.val, self.val, self.ctx.val)
        if r == 0:
            raise ZeroDivisionError("fmpz_mod_mat inv: matrix is not invertible")
        return res

    def det(self):
        """Return the determinant of a matrix.

        >>> from flint import fmpz_mod_mat, fmpz_mod_ctx
        >>> ctx = fmpz_mod_ctx(7)
        >>> M = fmpz_mod_mat([[1, 2], [3, 4]], ctx)
        >>> M.det()
        fmpz_mod(5, 7)

        """
        # XXX: No fmpz_mod_mat_det function...
        p = self.charpoly()
        p0 = p[0]
        if self.nrows() % 2:
            p0 = -p0
        return p0

    def charpoly(self):
        """Return the characteristic polynomial of a matrix.

        >>> from flint import fmpz_mod_mat, fmpz_mod_ctx
        >>> ctx = fmpz_mod_ctx(7)
        >>> M = fmpz_mod_mat([[1, 2], [3, 4]], ctx)
        >>> M.charpoly()
        x^2 + 2*x + 5

        """
        cdef fmpz_mod_poly_ctx pctx
        cdef fmpz_mod_poly res

        if self.nrows() != self.ncols():
            raise ValueError("fmpz_mod_mat charpoly: matrix must be square")

        pctx = fmpz_mod_poly_ctx(self.ctx)
        res = fmpz_mod_poly(0, pctx)
        compat_fmpz_mod_mat_charpoly(res.val, self.val, self.ctx.val)
        return res

    def minpoly(self):
        """Return the minimal polynomial of a matrix.

        >>> from flint import fmpz_mod_mat, fmpz_mod_ctx
        >>> ctx = fmpz_mod_ctx(7)
        >>> M = fmpz_mod_mat([[2, 1, 0], [0, 2, 0], [0, 0, 2]], ctx)
        >>> M.charpoly()
        x^3 + x^2 + 5*x + 6
        >>> M.minpoly()
        x^2 + 3*x + 4

        """
        cdef fmpz_mod_poly_ctx pctx
        cdef fmpz_mod_poly res

        if self.nrows() != self.ncols():
            raise ValueError("fmpz_mod_mat minpoly: matrix must be square")

        pctx = fmpz_mod_poly_ctx(self.ctx)
        res = fmpz_mod_poly(0, pctx)
        compat_fmpz_mod_mat_minpoly(res.val, self.val, self.ctx.val)
        return res

    def transpose(self):
        """Return the transpose of a matrix.

        >>> from flint import fmpz_mod_mat, fmpz_mod_ctx
        >>> ctx = fmpz_mod_ctx(7)
        >>> M = fmpz_mod_mat([[1, 2], [3, 4]], ctx)
        >>> M
        [1, 2]
        [3, 4]
        >>> M.transpose()
        [1, 3]
        [2, 4]
        """
        cdef fmpz_mod_mat res
        res = self._new(self.ncols(), self.nrows(), self.ctx)
        compat_fmpz_mod_mat_transpose(res.val, self.val, self.ctx.val)
        return res

    def solve(self, rhs):
        """Solve a linear system.

        >>> from flint import fmpz_mod_mat, fmpz_mod_ctx
        >>> ctx = fmpz_mod_ctx(7)
        >>> M = fmpz_mod_mat([[1, 2], [3, 4]], ctx)
        >>> rhs = fmpz_mod_mat([[5], [6]], ctx)
        >>> M.solve(rhs)
        [3]
        [1]

        The matrix must be square and invertible.

        Assumes that the modulus is prime.
        """
        cdef bint success
        cdef fmpz_mod_mat res

        rhs = any_as_fmpz_mod_mat(rhs)

        if rhs is NotImplemented:
            raise TypeError("fmpz_mod_mat solve: invalid rhs")
        if self.nrows() != self.ncols():
            raise ValueError("fmpz_mod_mat solve: matrix must be square")
        if self.nrows() != rhs.nrows():
            raise ValueError("fmpz_mod_mat solve: shape mismatch")

        res = self._new(rhs.nrows(), rhs.ncols(), self.ctx)
        success = compat_fmpz_mod_mat_solve(res.val, self.val, (<fmpz_mod_mat> rhs).val, self.ctx.val)

        if not success:
            raise ZeroDivisionError("Singular matrix in solve")

        return res

    def rank(self):
        """Return the rank of a matrix.

        >>> from flint import fmpz_mod_mat, fmpz_mod_ctx
        >>> ctx = fmpz_mod_ctx(11)
        >>> M = fmpz_mod_mat([[1, 2, 3], [4, 5, 6], [7, 8, 9]], ctx)
        >>> M.rank()
        2

        Assumes that the modulus is prime.
        """
        return self.rref()[1]

    def rref(self, inplace=False):
        """Return the reduced row echelon form of a matrix and the rank.

        >>> from flint import fmpz_mod_mat, fmpz_mod_ctx
        >>> ctx = fmpz_mod_ctx(7)
        >>> M = fmpz_mod_mat([[1, 2, 3], [4, 5, 6]], ctx)
        >>> Mr, r = M.rref()
        >>> Mr
        [1, 0, 6]
        [0, 1, 2]
        >>> r
        2

        If ``inplace`` is ``True``, the matrix is modified in place.

        Assumes that the modulus is prime.
        """
        cdef fmpz_mod_mat res
        cdef slong r
        if inplace:
            res = self
        else:
            res = self._copy()
        r = compat_fmpz_mod_mat_rref(res.val, res.ctx.val)
        return (res, r)
