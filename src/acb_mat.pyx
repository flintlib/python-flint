cdef acb_mat_coerce_operands(x, y):
    if typecheck(x, acb_mat):
        if isinstance(y, (fmpz_mat, fmpq_mat, arb_mat)):
            return x, acb_mat(y)
        if isinstance(y, (int, long, float, complex, fmpz, fmpq, arb, acb)):
            return x, acb_mat(x.nrows(), x.ncols(), y)
    elif typecheck(y, acb_mat):
        if isinstance(x, (fmpz_mat, fmpq_mat, arb_mat)):
            return acb_mat(x), y
        if isinstance(y, (int, long, float, complex, fmpz, fmpq, arb, acb)):
            return acb_mat(y.nrows(), y.ncols(), x), y
    return NotImplemented, NotImplemented

cdef acb_mat_coerce_scalar(x, y):
    if isinstance(y, (int, long, float, complex, fmpz, fmpq, arb, acb)):
        return x, any_as_acb(y)
    return NotImplemented, NotImplemented

cdef class acb_mat(flint_mat):
    """
    Represents a matrix over the complex numbers.

        >>> A = acb_mat([[1,2+1j],[3,4]]) ** 2 / 5
        >>> print(A.str(5, radius=False))
        [1.4000 + 0.60000j,  2.0000 + 1.0000j]
        [           3.0000, 4.4000 + 0.60000j]

    """
    cdef acb_mat_t val

    def __cinit__(self):
        acb_mat_init(self.val, 0, 0)

    def __dealloc__(self):
        acb_mat_clear(self.val)

    @classmethod
    def convert_operand(cls, x):
        """
        Attempts to convert *x* to an *acb_mat*, returning NotImplemented
        if unsuccessful.
        """
        if typecheck(x, cls):
            return x
        if typecheck(x, fmpz_mat) or typecheck(x, fmpq_mat) or typecheck(x, arb_mat):
            return cls(x)
        return NotImplemented

    @classmethod
    def convert(cls, x):
        """
        Attempts to convert *x* to an *acb_mat*, raising TypeError if
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
            if typecheck(val, acb_mat):
                m = acb_mat_nrows((<acb_mat>val).val)
                n = acb_mat_ncols((<acb_mat>val).val)
                acb_mat_init(self.val, m, n)
                acb_mat_set(self.val, (<acb_mat>val).val)
            elif typecheck(val, arb_mat):
                m = arb_mat_nrows((<arb_mat>val).val)
                n = arb_mat_ncols((<arb_mat>val).val)
                acb_mat_init(self.val, m, n)
                for i in range(m):
                    for j in range(n):
                        acb_set_arb(acb_mat_entry(self.val, i, j), arb_mat_entry((<arb_mat>val).val, i, j))
            elif typecheck(val, fmpz_mat):
                m = fmpz_mat_nrows((<fmpz_mat>val).val)
                n = fmpz_mat_ncols((<fmpz_mat>val).val)
                acb_mat_init(self.val, m, n)
                acb_mat_set_fmpz_mat(self.val, (<fmpz_mat>val).val)
            elif typecheck(val, fmpq_mat):
                m = fmpq_mat_nrows((<fmpq_mat>val).val)
                n = fmpq_mat_ncols((<fmpq_mat>val).val)
                acb_mat_init(self.val, m, n)
                acb_mat_set_fmpq_mat(self.val, (<fmpq_mat>val).val, getprec())
            elif isinstance(val, (list, tuple)):
                m = len(val)
                n = 0
                if m != 0:
                    if not isinstance(val[0], (list, tuple)):
                        raise TypeError("single input to acb_mat must be a list of lists")
                    n = len(val[0])
                    for i from 1 <= i < m:
                        if len(val[i]) != n:
                            raise ValueError("input rows have different lengths")
                acb_mat_init(self.val, m, n)
                for i from 0 <= i < m:
                    row = val[i]
                    for j from 0 <= j < n:
                        x = acb(row[j])
                        acb_set(acb_mat_entry(self.val, i, j), (<acb>x).val)
            elif hasattr(val, "rows"):   # allows conversion from mpmath matrices
                m = val.rows
                n = val.cols
                acb_mat_init(self.val, m, n)
                for i from 0 <= i < m:
                    for j from 0 <= j < n:
                        x = acb(val[i,j])
                        acb_set(acb_mat_entry(self.val, i, j), (<acb>x).val)
            else:
                raise TypeError("cannot create acb_mat from input of type %s" % type(val))
        elif len(args) == 2:
            m, n = args
            acb_mat_init(self.val, m, n)
        elif len(args) == 3:
            m, n, entries = args
            acb_mat_init(self.val, m, n)
            if isinstance(entries, (int, long, float, complex, fmpz, fmpq, arb, acb)):
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
                    x = acb(entries[i*n + j])
                    acb_set(acb_mat_entry(self.val, i, j), (<acb>x).val)
        else:
            raise ValueError("acb_mat: expected 1-3 arguments")

    def __nonzero__(self):
        raise NotImplementedError

    cpdef long nrows(s):
        """
        Returns the number of rows of *s*.
        """
        return acb_mat_nrows(s.val)

    cpdef long ncols(s):
        """
        Returns the number of columns of *s*.
        """
        return acb_mat_ncols(s.val)

    def __getitem__(self, index):
        cdef long i, j
        cdef acb x
        i, j = index
        if i < 0 or i >= self.nrows() or j < 0 or j >= self.ncols():
            raise ValueError("index %i,%i exceeds matrix dimensions" % (i, j))
        x = acb.__new__(acb)
        acb_set(x.val, acb_mat_entry(self.val, i, j))
        return x

    def __setitem__(self, index, value):
        cdef long i, j
        i, j = index
        if i < 0 or i >= self.nrows() or j < 0 or j >= self.ncols():
            raise ValueError("index %i,%i exceeds matrix dimensions" % (i, j))
        c = any_as_acb(value)
        acb_set(acb_mat_entry(self.val, i, j), (<acb>c).val)

    def transpose(s):
        """
        Returns the transpose of *s*.
        """
        cdef acb_mat u
        u = acb_mat.__new__(acb_mat)
        acb_mat_init(u.val, acb_mat_ncols(s.val), acb_mat_nrows(s.val))
        acb_mat_transpose(u.val, s.val)
        return u

    def conjugate(s):
        """
        Returns the entrywise conjugate (not the conjugate transpose) of *s*.

            >>> acb_mat([[1-2j, 1+3j]]).conjugate()
            [1.00000000000000 + 2.00000000000000j, 1.00000000000000 - 3.00000000000000j]
        """
        cdef acb_mat u
        u = acb_mat.__new__(acb_mat)
        acb_mat_init(u.val, acb_mat_nrows(s.val), acb_mat_ncols(s.val))
        acb_mat_conjugate(u.val, s.val)
        return u

    def det(s):
        """
        Returns the determinant of the square matrix *s* as an *acb*.

            >>> A = acb_mat(3, 3, range(9))
            >>> showgood(lambda: A.det(), dps=25)    # singular
            0
            >>> A[2,2] = 10
            >>> showgood(lambda: A.det(), dps=25)
            -6.000000000000000000000000
            >>> showgood(lambda: (A * A).det())
            36.0000000000000
            >>> print(acb_mat(0, 0).det())
            1.00000000000000
        """
        cdef acb d
        if acb_mat_nrows(s.val) != acb_mat_ncols(s.val):
            raise ValueError("matrix must be square")
        d = acb.__new__(acb)
        acb_mat_det(d.val, s.val, getprec())
        return d

    def __pos__(s):
        cdef acb_mat u
        u = acb_mat.__new__(acb_mat)
        acb_mat_init(u.val, acb_mat_nrows(s.val), acb_mat_ncols(s.val))
        acb_mat_set(u.val, s.val)   # round?
        return u

    def __neg__(s):
        cdef acb_mat u
        u = acb_mat.__new__(acb_mat)
        acb_mat_init(u.val, acb_mat_nrows(s.val), acb_mat_ncols(s.val))
        acb_mat_neg(u.val, s.val)   # round?
        return u

    def __add__(s, t):
        cdef long m, n
        if type(s) is type(t):
            m = (<acb_mat>s).nrows()
            n = (<acb_mat>s).ncols()
            if m != (<acb_mat>t).nrows() or n != (<acb_mat>t).ncols():
                raise ValueError("incompatible shapes for matrix addition")
            u = acb_mat.__new__(acb_mat)
            acb_mat_init((<acb_mat>u).val, m, n)
            acb_mat_add((<acb_mat>u).val, (<acb_mat>s).val, (<acb_mat>t).val, getprec())
            return u
        s, t = acb_mat_coerce_operands(s, t)
        if s is NotImplemented:
            return s
        return s + t

    def __sub__(s, t):
        cdef long m, n
        if type(s) is type(t):
            m = (<acb_mat>s).nrows()
            n = (<acb_mat>s).ncols()
            if m != (<acb_mat>t).nrows() or n != (<acb_mat>t).ncols():
                raise ValueError("incompatible shapes for matrix addition")
            u = acb_mat.__new__(acb_mat)
            acb_mat_init((<acb_mat>u).val, m, n)
            acb_mat_sub((<acb_mat>u).val, (<acb_mat>s).val, (<acb_mat>t).val, getprec())
            return u
        s, t = acb_mat_coerce_operands(s, t)
        if s is NotImplemented:
            return s
        return s - t

    def _scalar_mul_(s, acb t):
        cdef acb_mat u
        u = acb_mat.__new__(acb_mat)
        acb_mat_init(u.val, acb_mat_nrows(s.val), acb_mat_ncols(s.val))
        acb_mat_scalar_mul_acb(u.val, s.val, t.val, getprec())
        return u

    def __mul__(s, t):
        cdef acb_mat u
        if type(s) is type(t):
            if acb_mat_ncols((<acb_mat>s).val) != acb_mat_nrows((<acb_mat>t).val):
                raise ValueError("incompatible shapes for matrix multiplication")
            u = acb_mat.__new__(acb_mat)
            acb_mat_init(u.val, acb_mat_nrows((<acb_mat>s).val), acb_mat_ncols((<acb_mat>t).val))
            acb_mat_mul(u.val, (<acb_mat>s).val, (<acb_mat>t).val, getprec())
            return u
        if typecheck(s, acb_mat):
            c, d = acb_mat_coerce_scalar(s, t)
            if c is not NotImplemented:
                return c._scalar_mul_(d)
        else:
            d, c = acb_mat_coerce_scalar(t, s)
            if d is not NotImplemented:
                return d._scalar_mul_(c)
        s, t = acb_mat_coerce_operands(s, t)
        if s is NotImplemented:
            return s
        return s * t

    def _scalar_div_(s, acb t):
        cdef acb_mat u
        u = acb_mat.__new__(acb_mat)
        acb_mat_init(u.val, acb_mat_nrows(s.val), acb_mat_ncols(s.val))
        acb_mat_scalar_div_acb(u.val, s.val, t.val, getprec())
        return u

    @staticmethod
    def _div_(s, t):
        cdef acb_mat u
        if typecheck(s, acb_mat):
            s, t = acb_mat_coerce_scalar(s, t)
            if s is NotImplemented:
                return s
            return s._scalar_div_(t)
        return NotImplemented

    def __truediv__(s, t):
        return acb_mat._div_(s, t)

    def __div__(s, t):
        return acb_mat._div_(s, t)

    def __pow__(s, e, m):
        cdef acb_mat u
        cdef ulong exp
        cdef long n
        if not typecheck(s, acb_mat):
            return NotImplemented
        exp = e
        n = acb_mat_nrows((<acb_mat>s).val)
        if n != acb_mat_ncols((<acb_mat>s).val):
            raise ValueError("matrix must be square")
        if m is not None:
            raise NotImplementedError("modular matrix exponentiation")
        u = acb_mat.__new__(acb_mat)
        acb_mat_init((<acb_mat>u).val, n, n)
        acb_mat_pow_ui((<acb_mat>u).val, (<acb_mat>s).val, exp, getprec())
        return u

    def inv(s, bint nonstop=False):
        """
        Returns the inverse matrix of the square matrix *s*.

        If *s* is numerically singular, raises :exc:`ZeroDivisionError`
        unless *nonstop* is set in which case a matrix with NaN entries
        is returned.

            >>> A = acb_mat(2, 2, [1, 5, 2, 4])
            >>> print(A * A.inv())
            [[1.00000000000000 +/- 6.11e-16],                  [+/- 3.34e-16]]
            [                 [+/- 4.45e-16], [1.00000000000000 +/- 5.56e-16]]
            >>> A = acb_mat(2, 2, [1, 5, 2, 10])
            >>> A.inv()
            Traceback (most recent call last):
              ...
            ZeroDivisionError: matrix is singular
            >>> A.inv(nonstop=True)
            [nan + nanj, nan + nanj]
            [nan + nanj, nan + nanj]

        """
        cdef acb_mat u
        cdef long i, j
        if acb_mat_nrows(s.val) != acb_mat_ncols(s.val):
            raise ValueError("matrix must be square")
        u = acb_mat.__new__(acb_mat)
        acb_mat_init(u.val, acb_mat_nrows(s.val), acb_mat_ncols(s.val))
        if not acb_mat_inv(u.val, s.val, getprec()):
            if nonstop:
                for i from 0 <= i < acb_mat_nrows(u.val):
                    for j from 0 <= j < acb_mat_nrows(u.val):
                        acb_indeterminate(acb_mat_entry(u.val, i, j))
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

            >>> A = acb_mat(2, 2, [1, 2, 3, 4])
            >>> X = acb_mat(2, 3, range(6))
            >>> B = A * X
            >>> print(A.solve(B))
            [                 [+/- 4.74e-15], [1.00000000000000 +/- 4.78e-15], [2.00000000000000 +/- 8.52e-15]]
            [[3.00000000000000 +/- 3.56e-15], [4.00000000000000 +/- 3.59e-15], [5.00000000000000 +/- 6.28e-15]]
            >>> acb_mat([[1,1],[0,0]]).solve(acb_mat(2,3), nonstop=True)
            [nan + nanj, nan + nanj, nan + nanj]
            [nan + nanj, nan + nanj, nan + nanj]

        The optional *algorithm* can be None (default), "lu", or "precond".
        It can also be set to "approx" in which case an approximate
        floating-point solution (warning: without error bounds!) is returned.
        """
        cdef acb_mat u
        cdef bint res
        cdef long i, j
        t = acb_mat.convert(t)
        if (acb_mat_nrows(s.val) != acb_mat_ncols(s.val) or
            acb_mat_nrows(s.val) != acb_mat_nrows((<acb_mat>t).val)):
            raise ValueError("need a square system and compatible right hand side")
        u = acb_mat.__new__(acb_mat)
        acb_mat_init(u.val, acb_mat_nrows((<acb_mat>t).val), acb_mat_ncols((<acb_mat>t).val))
        if algorithm is None:
            res = acb_mat_solve(u.val, s.val, (<acb_mat>t).val, getprec())
        elif algorithm == 'lu':
            res = acb_mat_solve_lu(u.val, s.val, (<acb_mat>t).val, getprec())
        elif algorithm == 'precond':
            res = acb_mat_solve_precond(u.val, s.val, (<acb_mat>t).val, getprec())
        elif algorithm == 'approx':
            res = acb_mat_approx_solve(u.val, s.val, (<acb_mat>t).val, getprec())
        else:
            raise ValueError("unknown algorithm")
        if not res:
            if nonstop:
                for i from 0 <= i < acb_mat_nrows(u.val):
                    for j from 0 <= j < acb_mat_ncols(u.val):
                        acb_indeterminate(acb_mat_entry(u.val, i, j))
            else:
                raise ZeroDivisionError("singular matrix in solve()")
        return u

    def exp(s):
        """
        Returns the matrix exponential of *s*.

            >>> print(acb_mat(2, 2, [1, 4, -2, 1]).exp())
            [ [-2.58607310345045 +/- 5.06e-15],  [1.18429895089106 +/- 1.15e-15]]
            [[-0.592149475445530 +/- 5.73e-16], [-2.58607310345045 +/- 5.06e-15]]
        """
        cdef acb_mat u
        if acb_mat_nrows(s.val) != acb_mat_ncols(s.val):
            raise ValueError("matrix must be square")
        u = acb_mat.__new__(acb_mat)
        acb_mat_init(u.val, acb_mat_nrows(s.val), acb_mat_ncols(s.val))
        acb_mat_exp(u.val, s.val, getprec())
        return u

    def charpoly(s):
        """
        Returns the characteristic polynomial of *s* as an *acb_poly*.

            >>> print(acb_mat(2, 2, [1, 1, 1, 0]).charpoly())
            1.00000000000000*x^2 + (-1.00000000000000)*x + (-1.00000000000000)
        """
        cdef acb_poly u
        if acb_mat_nrows(s.val) != acb_mat_ncols(s.val):
            raise ValueError("matrix must be square")
        u = acb_poly.__new__(acb_poly)
        acb_mat_charpoly(u.val, s.val, getprec())
        return u

    def mid(s):
        """
        Returns the matrix consisting of the midpoints of the entries of *s*.

            >>> acb_mat([["1.5 +/- 0.1", 3]]).mid()
            [1.50000000000000, 3.00000000000000]
        """
        cdef acb_mat u
        u = acb_mat.__new__(acb_mat)
        acb_mat_init(u.val, acb_mat_nrows(s.val), acb_mat_ncols(s.val))
        acb_mat_get_mid(u.val, s.val)
        return u

    def trace(s):
        """
        Returns the trace of the square matrix *s* as an *acb*.

            >>> acb_mat([[3,4],[5,7]]).trace()
            10.0000000000000
        """
        cdef acb d
        if acb_mat_nrows(s.val) != acb_mat_ncols(s.val):
            raise ValueError("matrix must be square")
        d = acb.__new__(acb)
        acb_mat_trace(d.val, s.val, getprec())
        return d

    @classmethod
    def dft(cls, long n, long m=-1):
        """
        Returns the size *n* by *n* DFT matrix (optionally a separate
        number of columns *m* can be given in which case the periodic
        extension of the smaller dimension is used).

            >>> print(acb_mat.dft(3).str(5, radius=False))
            [0.57735,             0.57735,             0.57735]
            [0.57735, -0.28868 - 0.50000j, -0.28868 + 0.50000j]
            [0.57735, -0.28868 + 0.50000j, -0.28868 - 0.50000j]

        """
        cdef acb_mat u
        if m < 0:
            m = n
        assert n >= 0
        u = acb_mat.__new__(acb_mat)
        acb_mat_init(u.val, n, m)
        acb_mat_dft(u.val, 0, getprec())
        return u

    def overlaps(s, acb_mat t):
        """
        Returns whether *s* and *t* overlap (in the sense of balls).

            >>> A = acb_mat([[1,2],[3,4]])
            >>> ((A / 3) * 3).overlaps(A)
            True
            >>> ((A / 3) * 3 + 0.0001).overlaps(A)
            False
        """
        return bool(acb_mat_overlaps(s.val, t.val))

    def contains(s, t):
        """
        Returns whether *t* is contained in *s* (in the sense of balls).

            >>> A = acb_mat([[1,2],[3,4]])
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
            >>> ((A / 3) * 3).contains(arb_mat([[1,2],[3,4]]))
            True
        """
        if isinstance(t, acb_mat):
            return bool(acb_mat_contains(s.val, (<acb_mat>t).val))
        if isinstance(t, arb_mat):
            return bool(acb_mat_contains(s.val, (<acb_mat>(acb_mat(t))).val))
        if isinstance(t, fmpz_mat):
            return bool(acb_mat_contains_fmpz_mat(s.val, (<fmpz_mat>t).val))
        if isinstance(t, fmpq_mat):
            return bool(acb_mat_contains_fmpq_mat(s.val, (<fmpq_mat>t).val))
        raise TypeError("expected a matrix of compatible type")

    def chop(s, tol):
        """
        Returns a copy of *s* where real and imaginary parts of entries
        that are bounded by *tol* in magnitude have been replaced by
        exact zeros.

            >>> A = acb_mat([[1, 1+1e-20j], [1e-20+1j, 1e-20+1e-20j]])
            >>> A.chop(1e-6)
            [ 1.00000000000000, 1.00000000000000]
            [1.00000000000000j,                0]

        """
        cdef acb_mat u
        cdef arb b
        cdef long i, j, n, m
        u = acb_mat(s)
        n = s.nrows()
        m = s.ncols()
        b = arb(tol)
        arb_get_mag_lower(arb_radref(b.val), b.val)
        arf_zero(arb_midref(b.val))
        for i from 0 <= i < n:
            for j from 0 <= j < m:
                # and arb_contains_zero(...)?
                if arb_contains(b.val, acb_realref(acb_mat_entry(u.val, i, j))):
                    arb_zero(acb_realref(acb_mat_entry(u.val, i, j)))
                if arb_contains(b.val, acb_imagref(acb_mat_entry(u.val, i, j))):
                    arb_zero(acb_imagref(acb_mat_entry(u.val, i, j)))
        return u

    def __richcmp__(s, t, int op):
        cdef int stype, ttype
        cdef bint res
        if not (op == 2 or op == 3):
            raise ValueError("comparing matrices")
        if type(s) is not type(t):
            s, t = acb_mat_coerce_operands(s, t)
            if s is NotImplemented:
                return s
        if op == 2:
            res = acb_mat_eq((<acb_mat>s).val, (<acb_mat>t).val)
        else:
            res = acb_mat_ne((<acb_mat>s).val, (<acb_mat>t).val)
        return res

    @property
    def real(s):
        """
        Entrywise real part of this matrix as an *arb_mat*.

            >>> print(acb_mat.dft(3).real.str(5, radius=False))
            [0.57735,  0.57735,  0.57735]
            [0.57735, -0.28868, -0.28868]
            [0.57735, -0.28868, -0.28868]
        """
        cdef arb_mat u
        cdef long i, j, n, m
        n = s.nrows()
        m = s.ncols()
        u = arb_mat(n, m)
        for i from 0 <= i < n:
            for j from 0 <= j < m:
                acb_get_real(arb_mat_entry(u.val, i, j), acb_mat_entry(s.val, i, j))
        return u

    @property
    def imag(s):
        """
        Entrywise imaginary part of this matrix as an *arb_mat*.

            >>> print(acb_mat.dft(3).imag.str(5, radius=False))
            [0,        0,        0]
            [0, -0.50000,  0.50000]
            [0,  0.50000, -0.50000]
        """
        cdef arb_mat u
        cdef long i, j, n, m
        n = s.nrows()
        m = s.ncols()
        u = arb_mat(n, m)
        for i from 0 <= i < n:
            for j from 0 <= j < m:
                acb_get_imag(arb_mat_entry(u.val, i, j), acb_mat_entry(s.val, i, j))
        return u
