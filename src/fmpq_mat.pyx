cdef any_as_fmpq_mat(obj):
    if typecheck(obj, fmpq_mat):
        return obj
    if typecheck(obj, fmpz_mat):
        return fmpq_mat(obj)
    return NotImplemented

cdef class fmpq_mat(flint_mat):
    """
    Represents a dense matrix over the rational numbers.

        >>> A = fmpq_mat(3,3,[1,3,5,2,4,6,fmpq(2,3),2,4])
        >>> A.inv()
        [-3,  3/2, 3/2]
        [ 3, -1/2,  -3]
        [-1,    0, 3/2]
        >>> A.inv() * A
        [1, 0, 0]
        [0, 1, 0]
        [0, 0, 1]

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
                fmpq_mat_init(self.val, fmpq_mat_nrows((<fmpq_mat>val).val),
                                        fmpq_mat_ncols((<fmpq_mat>val).val))
                fmpq_mat_set(self.val, (<fmpq_mat>val).val)
            elif typecheck(val, fmpz_mat):
                fmpq_mat_init(self.val, fmpz_mat_nrows((<fmpz_mat>val).val),
                                        fmpz_mat_ncols((<fmpz_mat>val).val))
                fmpq_mat_set_fmpz_mat(self.val, (<fmpz_mat>val).val)
            elif isinstance(val, (list, tuple)):
                m = len(val)
                n = 0
                if m != 0:
                    if not isinstance(val[0], (list, tuple)):
                        raise TypeError("single input to fmpq_mat must be a list of lists")
                    n = len(val[0])
                    for i from 1 <= i < m:
                        if len(val[i]) != n:
                            raise ValueError("input rows have different lengths")
                fmpq_mat_init(self.val, m, n)
                for i from 0 <= i < m:
                    row = val[i]
                    for j from 0 <= j < n:
                        x = fmpq(row[j])
                        fmpq_set(fmpq_mat_entry(self.val, i, j), (<fmpq>x).val)
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

    def det(self):
        """
        Returns the determinant of *self* as an *fmpq*.

            >>> (fmpq_mat(2,2,[1,2,3,4]) / 5).det()
            -2/25
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

    @staticmethod
    def _div_(fmpq_mat s, t):
        t = any_as_fmpq(t)
        if t is NotImplemented:
            return t
        return s * (1 / t)

    def __truediv__(s, t):
        return fmpq_mat._div_(s, t)

    def __div__(s, t):
        return fmpq_mat._div_(s, t)

    def inv(self):
        """
        Returns the inverse matrix of *self*.

            >>> (fmpq_mat([[1,2],[3,4]]) / 5).inv()
            [ -10,    5]
            [15/2, -5/2]
            >>> (fmpq_mat([[1,2],[3,6]]) / 5).inv()
            Traceback (most recent call last):
              ...
            ZeroDivisionError: matrix is singular
        """
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
        Returns the transpose of *self*.
        
            >>> fmpq_mat(2,3,range(6)).transpose()
            [0, 3]
            [1, 4]
            [2, 5]
        """
        cdef fmpq_mat u
        u = fmpq_mat.__new__(fmpq_mat)
        fmpq_mat_init(u.val, fmpq_mat_ncols(self.val), fmpq_mat_nrows(self.val))
        fmpq_mat_transpose(u.val, self.val)
        return u

    def solve(self, other, algorithm=None):
        """
        Given matrices *A* and *B* represented by *self* and *other*,
        returns an *fmpq_mat* *X* such that `AX = B`, assuming that
        *A* is square and invertible.

        Algorithm can *None* for a default choice, or "fflu" or "dixon"
        (faster for large matrices).

            >>> A = fmpq_mat(2, 2, [1,4,8,3])
            >>> B = fmpq_mat(2, 3, range(6))
            >>> X = A.solve(B)
            >>> X
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
        if algorithm is None:
            if fmpq_mat_nrows(self.val) < 25:
                result = fmpq_mat_solve_fraction_free(u.val, self.val, (<fmpq_mat>t).val)
            else:
                result = fmpq_mat_solve_dixon(u.val, self.val, (<fmpq_mat>t).val)
        elif algorithm == "fflu":
            result = fmpq_mat_solve_fraction_free(u.val, self.val, (<fmpq_mat>t).val)
        elif algorithm == "dixon":
            result = fmpq_mat_solve_dixon(u.val, self.val, (<fmpq_mat>t).val)
        else:
            raise ValueError("unknown algorithm")
        if not result:
            raise ZeroDivisionError("singular matrix in solve()")
        return u

    def rref(self, inplace=False):
        """
        Computes the reduced row echelon form (rref) of *self*,
        either returning a new copy or modifying self in-place.
        Returns (*rref*, *rank*).

            >>> A = fmpq_mat(3,3,range(9))
            >>> A.rref()
            ([1, 0, -1]
            [0, 1,  2]
            [0, 0,  0], 2)
            >>> A.rref(inplace=True)
            ([1, 0, -1]
            [0, 1,  2]
            [0, 0,  0], 2)
            >>> A
            [1, 0, -1]
            [0, 1,  2]
            [0, 0,  0]

        """
        if inplace:
            res = self
        else:
            res = fmpq_mat.__new__(fmpq_mat)
            fmpq_mat_init((<fmpq_mat>res).val, fmpq_mat_nrows(self.val), fmpq_mat_ncols(self.val))
        rank = fmpq_mat_rref((<fmpq_mat>res).val, self.val)
        return res, rank

    @classmethod
    def hilbert(cls, long n, long m):
        """
        Returns the *n* by *m* truncated Hilbert matrix.

            >>> fmpq_mat.hilbert(2,3)
            [  1, 1/2, 1/3]
            [1/2, 1/3, 1/4]
        """
        cdef fmpq_mat u
        u = fmpq_mat(n, m)
        fmpq_mat_hilbert_matrix(u.val)
        return u

    def numer_denom(self):
        """
        Returns (*A*, *d*) where *A* is an *fmpz_mat* and *d* is an
        *fmpz* representing the minimal denominator such that
        *A* times *d* equals *self*.

            >>> A, d = fmpq_mat.hilbert(3,3).numer_denom()
            >>> A
            [60, 30, 20]
            [30, 20, 15]
            [20, 15, 12]
            >>> d
            60
        """
        cdef fmpz_mat num
        cdef fmpz_den
        num = fmpz_mat(self.nrows(), self.ncols())
        den = fmpz()
        fmpq_mat_get_fmpz_mat_matwise(num.val, den.val, self.val)
        return num, den

    def charpoly(self):
        cdef fmpq_poly u
        u = fmpq_poly.__new__(fmpq_poly)
        fmpq_poly_init(u.val)
        fmpq_mat_charpoly(u.val, self.val)
        return u

    def minpoly(self):
        cdef fmpq_poly u
        u = fmpq_poly.__new__(fmpq_poly)
        fmpq_poly_init(u.val)
        fmpq_mat_minpoly(u.val, self.val)
        return u

    def __pow__(self, n, z):
        cdef fmpq_mat v
        assert z is None
        n = int(n)
        if n == 0:
            r, c = self.nrows(), self.ncols()
            assert r == c
            v = fmpq_mat(r, c)
            fmpq_mat_one(v.val)
            return v
        if n == 1:
            return self
        if n == 2:
            return self * self
        if n < 0:
            return self.inv() ** (-n)
        v = self ** (n // 2)
        v = v * v
        if n % 2:
            v *= self
        return v


