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

    def eig(s, bint left=False, bint right=False, multiple=False, algorithm=None, tol=None, long maxiter=0, bint nonstop=False):
        r"""
        Computes eigenvalues and optionally eigenvectors of this matrix.
        Returns either *E*, (*E*, *L*), (*E*, *R*) or (*E*, *L*, *R*)
        depending on whether the flags *left* and *right* are set,
        where *E* is a list of the eigenvalues, *L* is a matrix
        with the left eigenvectors as rows, and *R* is a matrix
        with the right eigenvectors as columns.

        The *algorithm* can be "rump", "vdhoeven_mourrain", or *None*
        to use a default algorithm. Typically "rump" is slower and more
        accurate while "vdhoeven_mourrain" (the current default)
        is faster and less accurate.

            >>> A = acb_mat([[2,3,5],[7,11,13],[17,19,23]])
            >>> for c in A.eig(): print(c)
            ...
            [1.105299634957 +/- 6.34e-13] + [+/- 1.83e-13]j
            [-1.917027627441 +/- 2.64e-13] + [+/- 1.83e-13]j
            [36.811727992483 +/- 6.97e-13] + [+/- 1.83e-13]j
            >>> for c in A.eig(algorithm="rump"): print(c)
            ...
            [1.10529963495745 +/- 4.71e-15] + [+/- 2.92e-15]j
            [-1.91702762744092 +/- 8.45e-15] + [+/- 3.86e-15]j
            [36.8117279924835 +/- 4.72e-14] + [+/- 9.07e-15]j

        With the left and right eigenvector matrices, a complete
        diagonalization of the matrix is produced:

            >>> A = acb_mat([[2,3,5],[7,11,13],[17,19,23]])
            >>> E, L, R = A.eig(left=True, right=True)
            >>> D = acb_mat(3,3)
            >>> for i in range(3): D[i,i] = E[i]
            ...
            >>> (L*A*R - D).contains(acb_mat(3,3))
            True
            >>> (R*D*L - A).contains(acb_mat(3,3))
            True

        Ill-conditioned or large matrices may require high precision
        to isolate the eigenvalues::

            >>> sum(acb_mat(arb_mat.hilbert(20,20)).eig())
            Traceback (most recent call last):
              ...
            ValueError: failed to isolate eigenvalues (try higher prec, multiple=True for multiple eigenvalues, or nonstop=True to avoid the exception)
            >>> sum(acb_mat(arb_mat.hilbert(20,20)).eig(nonstop=True))
            nan + nanj
            >>> showgood(lambda: sum(acb_mat(arb_mat.hilbert(20,20)).eig(nonstop=True)), parts=False)
            2.47967321036454 + 0e-55j

        With default options, the method only succeeds if all eigenvalues can be
        isolated. Multiple (overlapping) eigenvalues can be handled by
        setting *multiple* = *True*.

            >>> acb_mat.dft(4).eig()
            Traceback (most recent call last):
              ...
            ValueError: failed to isolate eigenvalues (try higher prec, multiple=True for multiple eigenvalues, or nonstop=True to avoid the exception)
            >>> acb_mat.dft(4).eig(nonstop=True)
            [nan + nanj, nan + nanj, nan + nanj, nan + nanj]
            >>> acb_mat.dft(4).eig(multiple=True)
            [[-1.0000000000000 +/- 2.26e-15] + [+/- 1.23e-15]j, [+/- 4.96e-16] + [-1.00000000000000 +/- 3.72e-16]j, [1.00000000000000 +/- 4.98e-16] + [+/- 3.42e-16]j, [1.00000000000000 +/- 4.98e-16] + [+/- 3.42e-16]j]

        At this time, computing the eigenvectors is not supported
        with multiple eigenvalues:

            >>> acb_mat.dft(4).eig(multiple=True, right=True)
            Traceback (most recent call last):
              ...
            NotImplementedError: eigenvectors not supported with multiple=True

        The *algorithm* can also be set to "approx" to compute
        approximate eigenvalues and/or eigenvectors without error bounds.

            >>> for c in acb_mat.dft(4).eig(algorithm="approx"): print(c.str(radius=False))
            ...
            -0.999999999999999 - 7.85046229341892e-17j
            -2.35513868802566e-16 - 1.00000000000000j
            1.00000000000000 - 6.64346650360854e-17j
            0.999999999999999 - 5.14675360671472e-17j

        If *algorithm* is set to "approx", then *multiple* has
        no effect, and both eigenvalues and eigenvectors can be computed
        regardless of overlap.

            >>> E = acb_mat.dft(4).eig(algorithm="approx")
            >>> E, R = acb_mat.dft(4).eig(right=True, algorithm="approx")
            >>> E, L = acb_mat.dft(4).eig(left=True, algorithm="approx")
            >>> E, L, R = acb_mat.dft(4).eig(left=True, right=True, algorithm="approx")

        """
        cdef acb_mat E, L, R
        cdef acb_mat_struct *LP
        cdef acb_mat_struct *RP
        cdef mag_struct * magp
        cdef long i, n, prec
        cdef int success
        cdef mag_t tolm
        n = s.nrows()
        if n != s.ncols():
            raise ValueError("matrix must be square")
        prec = getprec()
        E = acb_mat(1, n)
        if left:
            L = acb_mat(n, n)
            LP = &(L.val[0])
        else:
            LP = NULL
        if right or algorithm != "approx":
            R = acb_mat(n, n)
            RP = &(R.val[0])
        else:
            RP = NULL
        if tol is not None:
            mag_init(tolm)
            tol = arb(tol)
            arb_get_mag(tolm, (<arb>tol).val)
            magp = tolm
        else:
            magp = NULL
        if n != 0:
            if algorithm == "approx":
                acb_mat_approx_eig_qr(acb_mat_entry(E.val, 0, 0),
                    LP, RP,  s.val, magp, maxiter, getprec())
            else:
                acb_mat_approx_eig_qr(acb_mat_entry(E.val, 0, 0),
                    NULL, RP, s.val, magp, maxiter, getprec())
                if multiple:
                    if left or right:
                        raise NotImplementedError("eigenvectors not supported with multiple=True")
                    if algorithm == "rump":
                        success = acb_mat_eig_multiple_rump(acb_mat_entry(E.val, 0, 0),
                            s.val, acb_mat_entry(E.val, 0, 0), RP, prec)
                    else:
                        success = acb_mat_eig_multiple(acb_mat_entry(E.val, 0, 0),
                            s.val, acb_mat_entry(E.val, 0, 0), RP, prec)
                else:
                    if algorithm == "rump":
                        success = acb_mat_eig_simple_rump(acb_mat_entry(E.val, 0, 0),
                            LP, RP, s.val, acb_mat_entry(E.val, 0, 0), RP, prec)
                    elif algorithm == "vdhoeven_mourrain":
                        success = acb_mat_eig_simple_vdhoeven_mourrain(acb_mat_entry(E.val, 0, 0),
                            LP, RP, s.val, acb_mat_entry(E.val, 0, 0), RP, prec)
                    else:
                        success = acb_mat_eig_simple(acb_mat_entry(E.val, 0, 0),
                            LP, RP, s.val, acb_mat_entry(E.val, 0, 0), RP, prec)
                if not (nonstop or success):
                    raise ValueError("failed to isolate eigenvalues (try higher prec, multiple=True for multiple eigenvalues, or nonstop=True to avoid the exception)")
        if tol is not None:
            mag_clear(tolm)
        Elist = [acb() for i in range(n)]
        for i in range(n):
            acb_swap((<acb>(Elist [i])).val, acb_mat_entry(E.val, 0, i))
        if not left and not right:
            return Elist
        if left and right:
            return Elist, L, R
        if left:
            return Elist, L
        return Elist, R
