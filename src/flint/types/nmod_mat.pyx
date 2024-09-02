cimport cython

from flint.flintlib.flint cimport ulong, mp_limb_t
from flint.flintlib.nmod cimport nmod_t

from flint.flintlib.nmod_poly cimport (
    nmod_poly_init,
)

from flint.flintlib.fmpz_mat cimport fmpz_mat_nrows, fmpz_mat_ncols
from flint.flintlib.fmpz_mat cimport fmpz_mat_get_nmod_mat

from flint.flintlib.nmod_mat cimport (
    nmod_mat_struct,
    nmod_mat_init,
    nmod_mat_init_set,
    nmod_mat_clear,
    nmod_mat_nrows,
    nmod_mat_ncols,
    nmod_mat_is_square,
    nmod_mat_entry,
    nmod_mat_set_entry,
    nmod_mat_equal,
    nmod_mat_is_zero,
    nmod_mat_transpose,
    nmod_mat_scalar_mul,
    nmod_mat_neg,
    nmod_mat_add,
    nmod_mat_sub,
    nmod_mat_mul,
    nmod_mat_pow,
    nmod_mat_inv,
    nmod_mat_solve,
    nmod_mat_nullspace,
    nmod_mat_det,
    nmod_mat_rref,
    nmod_mat_rank,
    nmod_mat_charpoly,
    nmod_mat_minpoly,
    nmod_mat_randtest,
)

from flint.utils.typecheck cimport typecheck
from flint.types.fmpz_mat cimport any_as_fmpz_mat
from flint.types.fmpz_mat cimport fmpz_mat
from flint.types.nmod cimport nmod, any_as_nmod_ctx
from flint.types.nmod_poly cimport nmod_poly, nmod_poly_new_init, any_as_nmod_poly_ctx
from flint.pyflint cimport global_random_state
from flint.flint_base.flint_context cimport thectx

from flint.flint_base.flint_base cimport flint_mat

from flint.utils.flint_exceptions import DomainError


ctx = thectx


cdef nmod_mat new_nmod_mat_init(nmod_ctx ctx, ulong m, ulong n):
    """New initialized nmod_mat of size m x n with context ctx."""
    cdef nmod_mat r = nmod_mat.__new__(nmod_mat)
    nmod_mat_init(r.val, m, n, ctx.mod.n)
    r.ctx = ctx
    return r


cdef nmod_mat new_nmod_mat_copy(nmod_mat other):
    """New copy of nmod_mat other."""
    cdef nmod_mat r = nmod_mat.__new__(nmod_mat)
    nmod_mat_init_set(r.val, other.val)
    r.ctx = other.ctx
    return r


cdef any_as_nmod_mat(obj, nmod_ctx ctx):
    """Convert obj to nmod_mat or return NotImplemented."""
    cdef nmod_mat r
    if typecheck(obj, nmod_mat):
        return obj

    x = any_as_fmpz_mat(obj)
    if x is not NotImplemented:
        r = new_nmod_mat_init(ctx,
                              fmpz_mat_nrows((<fmpz_mat>x).val),
                              fmpz_mat_ncols((<fmpz_mat>x).val))
        fmpz_mat_get_nmod_mat(r.val, (<fmpz_mat>x).val)
        return r

    return NotImplemented


cdef class nmod_mat(flint_mat):
    """
    The nmod_mat type represents dense matrices over Z/nZ for word-size n (see
    fmpz_mod_mat for larger moduli).

    Some operations may assume that n is a prime.
    """

#    cdef nmod_mat_t val

    def __dealloc__(self):
        nmod_mat_clear(self.val)

    @cython.embedsignature(False)
    def __init__(self, *args):
        cdef long m, n, i, j
        cdef mp_limb_t mod
        cdef nmod_ctx ctx

        if len(args) == 1:
            val = args[0]
            if typecheck(val, nmod_mat):
                nmod_mat_init_set(self.val, (<nmod_mat>val).val)
                self.ctx = (<nmod_mat>val).ctx
                return

        mod = args[-1]
        args = args[:-1]
        if mod == 0:
            raise ValueError("modulus must be nonzero")

        ctx = any_as_nmod_ctx(mod)
        self.ctx = ctx

        if len(args) == 1:
            val = args[0]
            if typecheck(val, fmpz_mat):
                nmod_mat_init(self.val, fmpz_mat_nrows((<fmpz_mat>val).val),
                              fmpz_mat_ncols((<fmpz_mat>val).val), ctx.mod.n)
                fmpz_mat_get_nmod_mat(self.val, (<fmpz_mat>val).val)
            elif isinstance(val, (list, tuple)):
                m = len(val)
                n = 0
                if m != 0:
                    if not isinstance(val[0], (list, tuple)):
                        raise TypeError("single input to nmod_mat must be a list of lists")
                    n = len(val[0])
                    for i from 1 <= i < m:
                        if len(val[i]) != n:
                            raise ValueError("input rows have different lengths")
                nmod_mat_init(self.val, m, n, ctx.mod.n)
                for i from 0 <= i < m:
                    row = val[i]
                    for j from 0 <= j < n:
                        x = nmod(row[j], ctx)        # XXX: slow
                        self.val.rows[i][j] = (<nmod>x).val
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
                    x = nmod(entries[i*n + j], ctx)         # XXX: slow
                    self.val.rows[i][j] = (<nmod>x).val
        else:
            raise TypeError("nmod_mat: expected 1-3 arguments plus modulus")

    def __bool__(self):
        return not nmod_mat_is_zero(self.val)

    def __richcmp__(s, t, int op):
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

    def repr(self):
        m = self.nrows()
        n = self.ncols()
        entries = ', '.join(map(str, self.entries()))
        return f"nmod_mat({m}, {n}, [{entries}], {self.modulus()})"

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

    def __getitem__(self, index):
        cdef long i, j
        cdef nmod x
        i, j = index
        if i < 0 or i >= self.nrows() or j < 0 or j >= self.ncols():
            raise IndexError("index %i,%i exceeds matrix dimensions" % (i, j))
        x = nmod(nmod_mat_entry(self.val, i, j), self.modulus())  # XXX: slow
        return x

    def __setitem__(self, index, value):
        cdef long i, j
        cdef mp_limb_t v
        i, j = index
        if i < 0 or i >= self.nrows() or j < 0 or j >= self.ncols():
            raise IndexError("index %i,%i exceeds matrix dimensions" % (i, j))
        if self.ctx.any_as_nmod(&v, value):
            nmod_mat_set_entry(self.val, i, j, v)
        else:
            raise TypeError("cannot set item of type %s" % type(value))

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
        sv = &(<nmod_mat>s).val[0]
        t = any_as_nmod_mat(t, s.ctx)
        if t is NotImplemented:
            return t
        tv = &(<nmod_mat>t).val[0]
        if sv.mod.n != tv.mod.n:
            raise ValueError("cannot add nmod_mats with different moduli")
        if sv.r != tv.r or sv.c != tv.c:
            raise ValueError("incompatible shapes for matrix addition")
        r = new_nmod_mat_init(s.ctx, sv.r, sv.c)
        nmod_mat_add(r.val, sv, tv)
        return r

    def __radd__(s, t):
        cdef nmod_mat r
        cdef nmod_mat_struct *sv
        cdef nmod_mat_struct *tv
        sv = &(<nmod_mat>s).val[0]
        t = any_as_nmod_mat(t, s.ctx)
        if t is NotImplemented:
            return t
        tv = &(<nmod_mat>t).val[0]
        if sv.mod.n != tv.mod.n:
            raise ValueError("cannot add nmod_mats with different moduli")
        if sv.r != tv.r or sv.c != tv.c:
            raise ValueError("incompatible shapes for matrix addition")
        r = new_nmod_mat_init(s.ctx, sv.r, sv.c)
        nmod_mat_add(r.val, sv, tv)
        return r

    def __sub__(s, t):
        cdef nmod_mat r
        cdef nmod_mat_struct *sv
        cdef nmod_mat_struct *tv
        sv = &(<nmod_mat>s).val[0]
        t = any_as_nmod_mat(t, s.ctx)
        if t is NotImplemented:
            return t
        tv = &(<nmod_mat>t).val[0]
        if sv.mod.n != tv.mod.n:
            raise ValueError("cannot subtract nmod_mats with different moduli")
        if sv.r != tv.r or sv.c != tv.c:
            raise ValueError("incompatible shapes for matrix subtraction")
        r = new_nmod_mat_init(s.ctx, sv.r, sv.c)
        nmod_mat_sub(r.val, sv, tv)
        return r

    def __rsub__(s, t):
        cdef nmod_mat r
        cdef nmod_mat_struct *sv
        cdef nmod_mat_struct *tv
        sv = &(<nmod_mat>s).val[0]
        t = any_as_nmod_mat(t, s.ctx)
        if t is NotImplemented:
            return t
        tv = &(<nmod_mat>t).val[0]
        if sv.mod.n != tv.mod.n:
            raise ValueError("cannot subtract nmod_mats with different moduli")
        if sv.r != tv.r or sv.c != tv.c:
            raise ValueError("incompatible shapes for matrix subtraction")
        r = new_nmod_mat_init(s.ctx, sv.r, sv.c)
        nmod_mat_sub(r.val, tv, sv)
        return r

    cdef __mul_nmod(self, mp_limb_t c):
        cdef nmod_mat r
        r = new_nmod_mat_init(self.ctx, self.val.r, self.val.c)
        nmod_mat_scalar_mul(r.val, self.val, c)
        return r

    def __mul__(s, t):
        cdef nmod_mat r
        cdef nmod_mat_struct *sv
        cdef nmod_mat_struct *tv
        cdef mp_limb_t c
        sv = &(<nmod_mat>s).val[0]
        u = any_as_nmod_mat(t, s.ctx)
        if u is NotImplemented:
            if s.ctx.any_as_nmod(&c, t):
                return (<nmod_mat>s).__mul_nmod(c)
            return NotImplemented
        tv = &(<nmod_mat>u).val[0]
        if sv.mod.n != tv.mod.n:
            raise ValueError("cannot multiply nmod_mats with different moduli")
        if sv.c != tv.r:
            raise ValueError("incompatible shapes for matrix multiplication")
        r = new_nmod_mat_init(s.ctx, sv.r, tv.c)
        nmod_mat_mul(r.val, sv, tv)
        return r

    def __rmul__(s, t):
        cdef nmod_mat_struct *sv
        cdef mp_limb_t c
        sv = &(<nmod_mat>s).val[0]
        if s.ctx.any_as_nmod(&c, t):
            return (<nmod_mat>s).__mul_nmod(c)
        u = any_as_nmod_mat(t, s.ctx)
        if u is NotImplemented:
            return u
        return u * s

    def __pow__(self, e, m):
        cdef nmod_mat t
        cdef ulong ee
        if not self.nrows() == self.ncols():
            raise ValueError("matrix must be square")
        if m is not None:
            raise NotImplementedError("modular matrix exponentiation")
        if e < 0:
            if not self.ctx._is_prime:
                raise DomainError("negative matrix power needs prime modulus")
            self = self.inv()
            e = -e
        ee = e
        t = new_nmod_mat_copy(self)
        nmod_mat_pow(t.val, t.val, ee)
        return t

    @staticmethod
    def _div_(nmod_mat s, t):
        cdef mp_limb_t v
        if not s.ctx.any_as_nmod(&v, t):
            return NotImplemented
        t = nmod(v, s.val.mod.n)
        try:
            tinv = ~t
        except ZeroDivisionError:
            # XXX: Maybe nmod.__invert__ should raise DomainError instead?
            if t == 0:
                raise ZeroDivisionError("division by zero")
            else:
                raise DomainError("nmod_mat division: modulus must be prime")
        return s * tinv

    def __truediv__(s, t):
        return nmod_mat._div_(s, t)

    def __div__(s, t):
        return nmod_mat._div_(s, t)

    def det(self):
        """
        Returns the determinant of self as an nmod.

            >>> nmod_mat(2,2,[1,2,3,4],17).det()
            15

        """
        if not nmod_mat_is_square(self.val):
            raise ValueError("matrix must be square")
        return nmod(nmod_mat_det(self.val), self.modulus())

    def inv(self):
        """
        Returns the inverse of self.

            >>> from flint import nmod_mat
            >>> A = nmod_mat(2,2,[1,2,3,4],17)
            >>> A.inv()
            [15, 1]
            [10, 8]

        """
        cdef nmod_mat u

        if not nmod_mat_is_square(self.val):
            raise ValueError("matrix must be square")
        if not self.ctx._is_prime:
            raise DomainError("nmod_mat inv: modulus must be prime")

        u = new_nmod_mat_init(self.ctx, nmod_mat_nrows(self.val),
                              nmod_mat_ncols(self.val))
        if not nmod_mat_inv(u.val, self.val):
            raise ZeroDivisionError("matrix is singular")
        return u

    def transpose(self):
        """
        Returns the transpose of self.

            >>> nmod_mat(2,3,range(6),7).transpose()
            [0, 3]
            [1, 4]
            [2, 5]
        """
        cdef nmod_mat u
        u = new_nmod_mat_init(self.ctx, nmod_mat_ncols(self.val),
                              nmod_mat_nrows(self.val))
        nmod_mat_transpose(u.val, self.val)
        return u

    def solve(self, other):
        """
        Given self = A and other = B, returns a matrix X such
        that A*X = B, assuming that self is square and invertible.

            >>> A = nmod_mat(2, 2, [1,4,8,3], 11)
            >>> B = nmod_mat(2, 3, range(6), 11)
            >>> X = A.solve(B)
            >>> X
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
        t = any_as_nmod_mat(other, self.ctx)
        if t is NotImplemented:
            raise TypeError("cannot convert input to nmod_mat")
        if (nmod_mat_nrows(self.val) != nmod_mat_ncols(self.val) or
            nmod_mat_nrows(self.val) != nmod_mat_nrows((<nmod_mat>t).val)):
            raise ValueError("need a square system and compatible right hand side")
        if not self.ctx._is_prime:
            raise DomainError("nmod_mat solve: modulus must be prime")

        u = new_nmod_mat_init(self.ctx, nmod_mat_nrows((<nmod_mat>t).val),
                              nmod_mat_ncols((<nmod_mat>t).val))
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
            ([1, 0, 22]
            [0, 1,  2]
            [0, 0,  0], 2)
            >>> A.rref(inplace=True)
            ([1, 0, 22]
            [0, 1,  2]
            [0, 0,  0], 2)
            >>> A
            [1, 0, 22]
            [0, 1,  2]
            [0, 0,  0]

        """
        if not self.ctx._is_prime:
            raise DomainError("rref only works for prime moduli")

        if inplace:
            res = self
        else:
            res = new_nmod_mat_copy(self)
        rank = nmod_mat_rref((<nmod_mat>res).val)
        return res, rank

    def rank(self):
        """Return the rank of a matrix.

        >>> from flint import nmod_mat
        >>> M = nmod_mat([[1, 2], [3, 4]], 11)
        >>> M.rank()
        2
        """
        if not self.ctx._is_prime:
            raise DomainError("rank only works for prime moduli")
        return nmod_mat_rank(self.val)

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
            [0, 0, 0, 0, 0]
            [0, 0, 0, 0, 0]
            [0, 0, 0, 0, 0]
            >>> X
            [ 1,  2,  3, 0, 0]
            [21, 20, 19, 0, 0]
            [ 1,  0,  0, 0, 0]
            [ 0,  1,  0, 0, 0]
            [ 0,  0,  1, 0, 0]

        """
        cdef nmod_mat res
        res = new_nmod_mat_init(self.ctx, nmod_mat_ncols(self.val),
                                          nmod_mat_ncols(self.val))
        nullity = nmod_mat_nullspace(res.val, self.val)
        return res, nullity

    def charpoly(self):
        """Return the characteristic polynomial of a matrix.

        >>> from flint import nmod_mat
        >>> M = nmod_mat([[1, 2], [3, 4]], 11)
        >>> M.charpoly()
        x^2 + 6*x + 9

        """
        cdef nmod_poly res

        if self.nrows() != self.ncols():
            raise ValueError("fmpz_mod_mat charpoly: matrix must be square")

        # XXX: don't create a new context for the polynomial
        res = nmod_poly_new_init(any_as_nmod_poly_ctx(self.ctx.mod.n))
        nmod_mat_charpoly(res.val, self.val)
        return res

    def minpoly(self):
        """Return the minimal polynomial of a matrix.

        >>> from flint import nmod_mat
        >>> M = nmod_mat([[2, 1, 0], [0, 2, 0], [0, 0, 2]], 7)
        >>> M.charpoly()
        x^3 + x^2 + 5*x + 6
        >>> M.minpoly()
        x^2 + 3*x + 4

        """
        cdef nmod_poly res

        if self.nrows() != self.ncols():
            raise ValueError("fmpz_mod_mat minpoly: matrix must be square")
        if not self.ctx._is_prime:
            raise DomainError("minpoly only works for prime moduli")

        # XXX: don't create a new context for the polynomial
        res = nmod_poly_new_init(any_as_nmod_poly_ctx(self.ctx.mod.n))
        nmod_mat_minpoly(res.val, self.val)
        return res
