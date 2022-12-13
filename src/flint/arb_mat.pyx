cdef arb_mat_coerce_operands(x, y):
    if typecheck(x, arb_mat):
        if isinstance(y, (fmpz_mat, fmpq_mat)):
            return x, arb_mat(y)
        if isinstance(y, (int, long, float, fmpz, fmpq, arb)):
            return x, arb_mat(x.nrows(), x.ncols(), y)
        if isinstance(y, (complex, acb)):
            return acb_mat(x), acb_mat(x.nrows(), x.ncols(), y)
    elif typecheck(y, arb_mat):
        if isinstance(x, (fmpz_mat, fmpq_mat)):
            return arb_mat(x), y
        if isinstance(y, (int, long, float, fmpz, fmpq, arb)):
            return arb_mat(y.nrows(), y.ncols(), x), y
        if isinstance(y, (complex, acb)):
            return acb_mat(y.nrows(), y.ncols(), x), acb_mat(y)
    return NotImplemented, NotImplemented

cdef arb_mat_coerce_scalar(x, y):
    if isinstance(y, (int, long, float, fmpz, fmpq, arb)):
        return x, any_as_arb(y)
    if isinstance(y, (complex, acb)):
        return acb_mat(x), any_as_acb(y)
    return NotImplemented, NotImplemented

cdef class arb_mat(flint_mat):
    """
    Represents a matrix over the real numbers.

        >>> A = arb_mat([[1,2],[3,4]]) ** 2 / 5
        >>> A
        [[1.40000000000000 +/- 3.12e-16],                2.00000000000000]
        [               3.00000000000000, [4.40000000000000 +/- 1.43e-15]]
        >>> print(A.str(5, radius=False))
        [1.4000, 2.0000]
        [3.0000, 4.4000]

    """

    cdef arb_mat_t val

    def __cinit__(self):
        arb_mat_init(self.val, 0, 0)

    def __dealloc__(self):
        arb_mat_clear(self.val)

    @classmethod
    def convert_operand(cls, x):
        """
        Attempts to convert *x* to an *arb_mat*, returning NotImplemented
        if unsuccessful.
        """
        if typecheck(x, cls):
            return x
        if typecheck(x, fmpz_mat) or typecheck(x, fmpq_mat):
            return cls(x)
        return NotImplemented

    @classmethod
    def convert(cls, x):
        """
        Attempts to convert *x* to an *arb_mat*, raising TypeError if
        unsuccessful.
        """
        x = cls.convert_operand(x)
        if x is NotImplemented:
            raise TypeError("unable to convert type %s to type %s" % (type(x), cls))
        return x

    @cython.embedsignature(False)
    def __init__(self, *args):
        cdef long m, n, i, j
        if len(args) == 1:
            val = args[0]
            if typecheck(val, arb_mat):
                m = arb_mat_nrows((<arb_mat>val).val)
                n = arb_mat_ncols((<arb_mat>val).val)
                arb_mat_init(self.val, m, n)
                arb_mat_set(self.val, (<arb_mat>val).val)
            elif typecheck(val, fmpz_mat):
                m = fmpz_mat_nrows((<fmpz_mat>val).val)
                n = fmpz_mat_ncols((<fmpz_mat>val).val)
                arb_mat_init(self.val, m, n)
                arb_mat_set_fmpz_mat(self.val, (<fmpz_mat>val).val)
            elif typecheck(val, fmpq_mat):
                m = fmpq_mat_nrows((<fmpq_mat>val).val)
                n = fmpq_mat_ncols((<fmpq_mat>val).val)
                arb_mat_init(self.val, m, n)
                arb_mat_set_fmpq_mat(self.val, (<fmpq_mat>val).val, getprec())
            elif isinstance(val, (list, tuple)):
                m = len(val)
                n = 0
                if m != 0:
                    if not isinstance(val[0], (list, tuple)):
                        raise TypeError("single input to arb_mat must be a list of lists")
                    n = len(val[0])
                    for i from 1 <= i < m:
                        if len(val[i]) != n:
                            raise ValueError("input rows have different lengths")
                arb_mat_init(self.val, m, n)
                for i from 0 <= i < m:
                    row = val[i]
                    for j from 0 <= j < n:
                        x = arb(row[j])
                        arb_set(arb_mat_entry(self.val, i, j), (<arb>x).val)
            elif hasattr(val, "rows"):   # allows conversion from mpmath matrices
                m = val.rows
                n = val.cols
                arb_mat_init(self.val, m, n)
                for i from 0 <= i < m:
                    for j from 0 <= j < n:
                        x = arb(val[i,j])
                        arb_set(arb_mat_entry(self.val, i, j), (<arb>x).val)
            else:
                raise TypeError("cannot create arb_mat from input of type %s" % type(val))
        elif len(args) == 2:
            m, n = args
            arb_mat_init(self.val, m, n)
        elif len(args) == 3:
            m, n, entries = args
            arb_mat_init(self.val, m, n)
            if isinstance(entries, (int, long, float, fmpz, fmpq, arb)):
                c = entries
                entries = [0] * (m * n)
                for i in range(min(m,n)):
                    entries[i*n + i] = c
            else:
                entries = list(entries)
            if len(entries) != m*n:
                raise ValueError("list of entries has the wrong length")
            for i from 0 <= i < m:
                for j from 0 <= j < n:
                    x = arb(entries[i*n + j])
                    arb_set(arb_mat_entry(self.val, i, j), (<arb>x).val)
        else:
            raise ValueError("arb_mat: expected 1-3 arguments")

    def __nonzero__(self):
        raise NotImplementedError

    cpdef long nrows(s):
        """
        Returns the number of rows of *s*.
        """
        return arb_mat_nrows(s.val)

    cpdef long ncols(s):
        """
        Returns the number of columns of *s*.
        """
        return arb_mat_ncols(s.val)

    def __getitem__(self, index):
        cdef long i, j
        cdef arb x
        i, j = index
        if i < 0 or i >= self.nrows() or j < 0 or j >= self.ncols():
            raise ValueError("index %i,%i exceeds matrix dimensions" % (i, j))
        x = arb.__new__(arb)
        arb_set(x.val, arb_mat_entry(self.val, i, j))
        return x

    def __setitem__(self, index, value):
        cdef long i, j
        i, j = index
        if i < 0 or i >= self.nrows() or j < 0 or j >= self.ncols():
            raise ValueError("index %i,%i exceeds matrix dimensions" % (i, j))
        c = any_as_arb(value)
        arb_set(arb_mat_entry(self.val, i, j), (<arb>c).val)

    def transpose(s):
        """
        Returns the transpose of *s*.
        """
        cdef arb_mat u
        u = arb_mat.__new__(arb_mat)
        arb_mat_init(u.val, arb_mat_ncols(s.val), arb_mat_nrows(s.val))
        arb_mat_transpose(u.val, s.val)
        return u

    def det(s):
        """
        Returns the determinant of the square matrix *s* as an *arb*.

            >>> A = arb_mat(3, 3, range(9))
            >>> showgood(lambda: A.det(), dps=25)    # exact singular
            0
            >>> A[2,2] = 10
            >>> showgood(lambda: A.det(), dps=25)
            -6.000000000000000000000000
            >>> showgood(lambda: (A * A).det())
            36.0000000000000
            >>> print(arb_mat(0, 0).det())
            1.00000000000000
        """
        cdef arb d
        if arb_mat_nrows(s.val) != arb_mat_ncols(s.val):
            raise ValueError("matrix must be square")
        d = arb.__new__(arb)
        arb_mat_det(d.val, s.val, getprec())
        return d

    def __pos__(s):
        cdef arb_mat u
        u = arb_mat.__new__(arb_mat)
        arb_mat_init(u.val, arb_mat_nrows(s.val), arb_mat_ncols(s.val))
        arb_mat_set(u.val, s.val)   # round?
        return u

    def __neg__(s):
        cdef arb_mat u
        u = arb_mat.__new__(arb_mat)
        arb_mat_init(u.val, arb_mat_nrows(s.val), arb_mat_ncols(s.val))
        arb_mat_neg(u.val, s.val)   # round?
        return u

    def __add__(s, t):
        cdef long m, n
        if type(s) is type(t):
            m = (<arb_mat>s).nrows()
            n = (<arb_mat>s).ncols()
            if m != (<arb_mat>t).nrows() or n != (<arb_mat>t).ncols():
                raise ValueError("incompatible shapes for matrix addition")
            u = arb_mat.__new__(arb_mat)
            arb_mat_init((<arb_mat>u).val, m, n)
            arb_mat_add((<arb_mat>u).val, (<arb_mat>s).val, (<arb_mat>t).val, getprec())
            return u
        s, t = arb_mat_coerce_operands(s, t)
        if s is NotImplemented:
            return s
        return s + t

    def __sub__(s, t):
        cdef long m, n
        if type(s) is type(t):
            m = (<arb_mat>s).nrows()
            n = (<arb_mat>s).ncols()
            if m != (<arb_mat>t).nrows() or n != (<arb_mat>t).ncols():
                raise ValueError("incompatible shapes for matrix addition")
            u = arb_mat.__new__(arb_mat)
            arb_mat_init((<arb_mat>u).val, m, n)
            arb_mat_sub((<arb_mat>u).val, (<arb_mat>s).val, (<arb_mat>t).val, getprec())
            return u
        s, t = arb_mat_coerce_operands(s, t)
        if s is NotImplemented:
            return s
        return s - t

    def _scalar_mul_(s, arb t):
        cdef arb_mat u
        u = arb_mat.__new__(arb_mat)
        arb_mat_init(u.val, arb_mat_nrows(s.val), arb_mat_ncols(s.val))
        arb_mat_scalar_mul_arb(u.val, s.val, t.val, getprec())
        return u

    def __mul__(s, t):
        cdef arb_mat u
        if type(s) is type(t):
            if arb_mat_ncols((<arb_mat>s).val) != arb_mat_nrows((<arb_mat>t).val):
                raise ValueError("incompatible shapes for matrix multiplication")
            u = arb_mat.__new__(arb_mat)
            arb_mat_init(u.val, arb_mat_nrows((<arb_mat>s).val), arb_mat_ncols((<arb_mat>t).val))
            arb_mat_mul(u.val, (<arb_mat>s).val, (<arb_mat>t).val, getprec())
            return u
        if typecheck(s, arb_mat):
            c, d = arb_mat_coerce_scalar(s, t)
            if c is not NotImplemented:
                return c._scalar_mul_(d)
        else:
            d, c = arb_mat_coerce_scalar(t, s)
            if d is not NotImplemented:
                return d._scalar_mul_(c)
        s, t = arb_mat_coerce_operands(s, t)
        if s is NotImplemented:
            return s
        return s * t

    def _scalar_div_(s, arb t):
        cdef arb_mat u
        u = arb_mat.__new__(arb_mat)
        arb_mat_init(u.val, arb_mat_nrows(s.val), arb_mat_ncols(s.val))
        arb_mat_scalar_div_arb(u.val, s.val, t.val, getprec())
        return u

    @staticmethod
    def _div_(s, t):
        cdef arb_mat u
        if typecheck(s, arb_mat):
            s, t = arb_mat_coerce_scalar(s, t)
            if s is NotImplemented:
                return s
            return s._scalar_div_(t)
        return NotImplemented

    def __truediv__(s, t):
        return arb_mat._div_(s, t)

    def __div__(s, t):
        return arb_mat._div_(s, t)

    def __pow__(s, e, m):
        cdef arb_mat u
        cdef ulong exp
        cdef long n
        if not typecheck(s, arb_mat):
            return NotImplemented
        exp = e
        n = arb_mat_nrows((<arb_mat>s).val)
        if n != arb_mat_ncols((<arb_mat>s).val):
            raise ValueError("matrix must be square")
        if m is not None:
            raise NotImplementedError("modular matrix exponentiation")
        u = arb_mat.__new__(arb_mat)
        arb_mat_init((<arb_mat>u).val, n, n)
        arb_mat_pow_ui((<arb_mat>u).val, (<arb_mat>s).val, exp, getprec())
        return u

    def inv(s, bint nonstop=False):
        """
        Returns the inverse matrix of the square matrix *s*.

        If *s* is numerically singular, raises :exc:`ZeroDivisionError`
        unless *nonstop* is set in which case a matrix with NaN entries
        is returned.

            >>> A = arb_mat(2, 2, [1, 5, 2, 4])
            >>> print(A * A.inv())
            [[1.00000000000000 +/- 6.11e-16],                  [+/- 3.34e-16]]
            [                 [+/- 4.45e-16], [1.00000000000000 +/- 5.56e-16]]
            >>> A = arb_mat(2, 2, [1, 5, 2, 10])
            >>> A.inv()
            Traceback (most recent call last):
              ...
            ZeroDivisionError: matrix is singular
            >>> A.inv(nonstop=True)
            [nan, nan]
            [nan, nan]

        """
        cdef arb_mat u
        if arb_mat_nrows(s.val) != arb_mat_ncols(s.val):
            raise ValueError("matrix must be square")
        u = arb_mat.__new__(arb_mat)
        arb_mat_init(u.val, arb_mat_nrows(s.val), arb_mat_ncols(s.val))
        if not arb_mat_inv(u.val, s.val, getprec()):
            if nonstop:
                for i from 0 <= i < arb_mat_nrows(u.val):
                    for j from 0 <= j < arb_mat_nrows(u.val):
                        arb_indeterminate(arb_mat_entry(u.val, i, j))
            else:
                raise ZeroDivisionError("matrix is singular")
        return u

    def solve(s, t, bint nonstop=False, algorithm=None):
        """
        Solves `AX = B` where *A* is a square matrix given by *s* and
        `B` is a matrix given by *t*.

        If *A* is numerically singular, raises :exc:`ZeroDivisionError`
        unless *nonstop* is set in which case a matrix with NaN entries
        is returned.

            >>> A = arb_mat(2, 2, [1, 2, 3, 4])
            >>> X = arb_mat(2, 3, range(6))
            >>> B = A * X
            >>> print(A.solve(B))
            [                 [+/- 4.74e-15], [1.00000000000000 +/- 4.78e-15], [2.00000000000000 +/- 8.52e-15]]
            [[3.00000000000000 +/- 3.56e-15], [4.00000000000000 +/- 3.59e-15], [5.00000000000000 +/- 6.28e-15]]
            >>> arb_mat([[1,1],[0,0]]).solve(arb_mat(2,3))
            Traceback (most recent call last):
              ...
            ZeroDivisionError: singular matrix in solve()
            >>> arb_mat([[1,1],[0,0]]).solve(arb_mat(2,3), nonstop=True)
            [nan, nan, nan]
            [nan, nan, nan]

        The optional *algorithm* can be None (default), "lu", or "precond".
        It can also be set to "approx" in which case an approximate
        floating-point solution (warning: without error bounds!) is returned.
        """
        cdef arb_mat u
        cdef bint res
        cdef long i, j
        t = arb_mat.convert(t)
        if (arb_mat_nrows(s.val) != arb_mat_ncols(s.val) or
            arb_mat_nrows(s.val) != arb_mat_nrows((<arb_mat>t).val)):
            raise ValueError("need a square system and compatible right hand side")
        u = arb_mat.__new__(arb_mat)
        arb_mat_init(u.val, arb_mat_nrows((<arb_mat>t).val), arb_mat_ncols((<arb_mat>t).val))
        if algorithm is None:
            res = arb_mat_solve(u.val, s.val, (<arb_mat>t).val, getprec())
        elif algorithm == 'lu':
            res = arb_mat_solve_lu(u.val, s.val, (<arb_mat>t).val, getprec())
        elif algorithm == 'precond':
            res = arb_mat_solve_precond(u.val, s.val, (<arb_mat>t).val, getprec())
        elif algorithm == "approx":
            res = arb_mat_approx_solve(u.val, s.val, (<arb_mat>t).val, getprec())
        else:
            raise ValueError("unknown algorithm")
        if not res:
            if nonstop:
                for i from 0 <= i < arb_mat_nrows(u.val):
                    for j from 0 <= j < arb_mat_ncols(u.val):
                        arb_indeterminate(arb_mat_entry(u.val, i, j))
            else:
                raise ZeroDivisionError("singular matrix in solve()")
        return u

    def exp(s):
        """
        Returns the matrix exponential of *s*.

            >>> print(arb_mat(2, 2, [1, 4, -2, 1]).exp())
            [ [-2.58607310345045 +/- 5.06e-15],  [1.18429895089106 +/- 1.15e-15]]
            [[-0.592149475445530 +/- 5.73e-16], [-2.58607310345045 +/- 5.06e-15]]
        """
        cdef arb_mat u
        if arb_mat_nrows(s.val) != arb_mat_ncols(s.val):
            raise ValueError("matrix must be square")
        u = arb_mat.__new__(arb_mat)
        arb_mat_init(u.val, arb_mat_nrows(s.val), arb_mat_ncols(s.val))
        arb_mat_exp(u.val, s.val, getprec())
        return u

    def charpoly(s):
        """
        Returns the characteristic polynomial of *s* as an *arb_poly*.

            >>> print(arb_mat(2, 2, [1, 1, 1, 0]).charpoly())
            1.00000000000000*x^2 + (-1.00000000000000)*x + (-1.00000000000000)
        """
        cdef arb_poly u
        if arb_mat_nrows(s.val) != arb_mat_ncols(s.val):
            raise ValueError("matrix must be square")
        u = arb_poly.__new__(arb_poly)
        arb_mat_charpoly(u.val, s.val, getprec())
        return u

    def mid(s):
        """
        Returns the matrix consisting of the midpoints of the entries of *s*.

            >>> arb_mat([["1.5 +/- 0.1", 3]]).mid()
            [1.50000000000000, 3.00000000000000]
        """
        cdef arb_mat u
        u = arb_mat.__new__(arb_mat)
        arb_mat_init(u.val, arb_mat_nrows(s.val), arb_mat_ncols(s.val))
        arb_mat_get_mid(u.val, s.val)
        return u

    def trace(s):
        """
        Returns the trace of the square matrix *s* as an *arb*.

            >>> arb_mat([[3,4],[5,7]]).trace()
            10.0000000000000
        """
        cdef arb d
        if arb_mat_nrows(s.val) != arb_mat_ncols(s.val):
            raise ValueError("matrix must be square")
        d = arb.__new__(arb)
        arb_mat_trace(d.val, s.val, getprec())
        return d

    @classmethod
    def hilbert(cls, long n, long m):
        """
        Returns the *n* by *m* truncated Hilbert matrix.

            >>> arb_mat.hilbert(6,2)
            [                1.00000000000000,                0.500000000000000]
            [               0.500000000000000, [0.333333333333333 +/- 3.71e-16]]
            [[0.333333333333333 +/- 3.71e-16],                0.250000000000000]
            [               0.250000000000000, [0.200000000000000 +/- 4.45e-17]]
            [[0.200000000000000 +/- 4.45e-17], [0.166666666666667 +/- 3.71e-16]]
            [[0.166666666666667 +/- 3.71e-16], [0.142857142857143 +/- 1.79e-16]]
        """
        cdef arb_mat u
        assert n >= 0
        assert m >= 0
        u = arb_mat.__new__(arb_mat)
        arb_mat_init(u.val, n, m)
        arb_mat_hilbert(u.val, getprec())
        return u

    @classmethod
    def pascal(cls, long n, long m, int triangular=0):
        """
        Returns the *n* by *m* truncated Pascal matrix. If *triangular*
        is 0, the symmetric version of this matrix is returned; if
        *triangular* is -1 or 1, the lower or upper triangular version
        of the Pascal matrix is returned.

            >>> arb_mat.pascal(3, 4)
            [1.00000000000000, 1.00000000000000, 1.00000000000000, 1.00000000000000]
            [1.00000000000000, 2.00000000000000, 3.00000000000000, 4.00000000000000]
            [1.00000000000000, 3.00000000000000, 6.00000000000000, 10.0000000000000]
            >>> arb_mat.pascal(3, 4, 1)
            [1.00000000000000, 1.00000000000000, 1.00000000000000, 1.00000000000000]
            [               0, 1.00000000000000, 2.00000000000000, 3.00000000000000]
            [               0,                0, 1.00000000000000, 3.00000000000000]
            >>> arb_mat.pascal(3, 4, -1)
            [1.00000000000000,                0,                0, 0]
            [1.00000000000000, 1.00000000000000,                0, 0]
            [1.00000000000000, 2.00000000000000, 1.00000000000000, 0]

        """
        cdef arb_mat u
        assert n >= 0
        assert m >= 0
        u = arb_mat.__new__(arb_mat)
        arb_mat_init(u.val, n, m)
        arb_mat_pascal(u.val, triangular, getprec())
        return u

    @classmethod
    def stirling(cls, long n, long m, int kind=0):
        """
        Returns the *n* by *m* truncated Stirling matrix. The
        parameter *kind* can be 0 for unsigned Stirling numbers of the
        first kind, 1 for signed Stirling numbers of the first kind,
        and 2 for Stirling numbers of the second kind.

            >>> arb_mat.stirling(5, 4)
            [1.00000000000000,                0,                0,                0]
            [               0, 1.00000000000000,                0,                0]
            [               0, 1.00000000000000, 1.00000000000000,                0]
            [               0, 2.00000000000000, 3.00000000000000, 1.00000000000000]
            [               0, 6.00000000000000, 11.0000000000000, 6.00000000000000]
            >>> arb_mat.stirling(5, 4, 1)
            [1.00000000000000,                 0,                 0,                 0]
            [               0,  1.00000000000000,                 0,                 0]
            [               0, -1.00000000000000,  1.00000000000000,                 0]
            [               0,  2.00000000000000, -3.00000000000000,  1.00000000000000]
            [               0, -6.00000000000000,  11.0000000000000, -6.00000000000000]
            >>> arb_mat.stirling(5, 4, 2)
            [1.00000000000000,                0,                0,                0]
            [               0, 1.00000000000000,                0,                0]
            [               0, 1.00000000000000, 1.00000000000000,                0]
            [               0, 1.00000000000000, 3.00000000000000, 1.00000000000000]
            [               0, 1.00000000000000, 7.00000000000000, 6.00000000000000]
        """
        cdef arb_mat u
        assert n >= 0
        assert m >= 0
        if not 0 <= kind <= 2:
            raise ValueError("expected kind = 0, 1 or 2")
        u = arb_mat.__new__(arb_mat)
        arb_mat_init(u.val, n, m)
        arb_mat_stirling(u.val, kind, getprec())
        return u

    @classmethod
    def dct(cls, long n, long m=-1):
        """
        Returns the size *n* by *n* DCT matrix (optionally a separate
        number of columns *m* can be given in which case the periodic
        extension of the smaller dimension is used).

            >>> print(arb_mat.dct(4).str(4))
            [              0.5000,                0.5000,                0.5000,                0.5000]
            [[0.6533 +/- 1.86e-5],  [0.2706 +/- 1.96e-6], [-0.2706 +/- 1.96e-6], [-0.6533 +/- 1.86e-5]]
            [   [0.5000 +/- 3e-9],    [-0.5000 +/- 3e-9],    [-0.5000 +/- 3e-9],     [0.5000 +/- 3e-9]]
            [[0.2706 +/- 1.96e-6], [-0.6533 +/- 1.86e-5],  [0.6533 +/- 1.86e-5], [-0.2706 +/- 1.96e-6]]

        """
        cdef arb_mat u
        if m < 0:
            m = n
        assert n >= 0
        u = arb_mat.__new__(arb_mat)
        arb_mat_init(u.val, n, m)
        arb_mat_dct(u.val, 0, getprec())
        return u

    def overlaps(s, arb_mat t):
        """
        Returns whether *s* and *t* overlap (in the sense of balls).

            >>> A = arb_mat([[1,2],[3,4]])
            >>> ((A / 3) * 3).overlaps(A)
            True
            >>> ((A / 3) * 3 + 0.0001).overlaps(A)
            False
        """
        return bool(arb_mat_overlaps(s.val, t.val))

    def contains(s, t):
        """
        Returns whether *t* is contained in *s* (in the sense of balls).

            >>> A = arb_mat([[1,2],[3,4]])
            >>> ((A / 3) * 3).contains(A)
            True
            >>> A.contains((A / 3) * 3)
            False
            >>> ((A / 3) * 3).contains(fmpz_mat([[1,2],[3,4]]))
            True
            >>> ((A / 3) * 3).contains(fmpz_mat([[1,2],[3,5]]))
            False
            >>> (A / 3).contains(fmpq_mat([[1,2],[3,4]]) / 3)
            True
        """
        if isinstance(t, arb_mat):
            return bool(arb_mat_contains(s.val, (<arb_mat>t).val))
        if isinstance(t, fmpz_mat):
            return bool(arb_mat_contains_fmpz_mat(s.val, (<fmpz_mat>t).val))
        if isinstance(t, fmpq_mat):
            return bool(arb_mat_contains_fmpq_mat(s.val, (<fmpq_mat>t).val))
        raise TypeError("expected a matrix of compatible type")

    def chop(s, tol):
        """
        Returns a copy of *s* where entries that are bounded by *tol* in
        magnitude have been replaced by exact zeros.

            >>> print(arb_mat.stirling(4, 4).inv().str(5, radius=False))
            [1.0000,       0,       0,      0]
            [     0,  1.0000,   0e-14,  0e-15]
            [     0, -1.0000,  1.0000,  0e-15]
            [     0,  1.0000, -3.0000, 1.0000]
            >>> print(arb_mat.stirling(4, 4).inv().chop(1e-6).str(5, radius=False))
            [1.0000,       0,       0,      0]
            [     0,  1.0000,       0,      0]
            [     0, -1.0000,  1.0000,      0]
            [     0,  1.0000, -3.0000, 1.0000]

        """
        cdef arb_mat u
        cdef arb b
        cdef long i, j, n, m
        u = arb_mat(s)
        n = s.nrows()
        m = s.ncols()
        b = arb(tol)
        arb_get_mag_lower(arb_radref(b.val), b.val)
        arf_zero(arb_midref(b.val))
        for i from 0 <= i < n:
            for j from 0 <= j < m:
                # and arb_contains_zero(...)?
                if arb_contains(b.val, arb_mat_entry(u.val, i, j)):
                    arb_zero(arb_mat_entry(u.val, i, j))
        return u

    def __richcmp__(s, t, int op):
        cdef int stype, ttype
        cdef bint res
        if not (op == 2 or op == 3):
            raise ValueError("comparing matrices")
        if type(s) is not type(t):
            s, t = arb_mat_coerce_operands(s, t)
            if s is NotImplemented:
                return s
        if op == 2:
            res = arb_mat_eq((<arb_mat>s).val, (<arb_mat>t).val)
        else:
            res = arb_mat_ne((<arb_mat>s).val, (<arb_mat>t).val)
        return res

    def eig(s, *args, **kwargs):
        r"""
        Computes eigenvalues and/or eigenvectors of this matrix.
        This is just a wrapper for :meth:`.acb_mat.eig`; see the
        documentation for that method for details.
        """
        return acb_mat(s).eig(*args, **kwargs)
