cdef any_as_fmpz_mat(obj):
    if typecheck(obj, fmpz_mat):
        return obj
    return NotImplemented

cdef class fmpz_mat(flint_mat):
    """
    The *fmpz_mat* type represents dense matrices over the integers.

    An *fmpz_mat* can be constructed empty, from a list of entries in
    row-major order, from a list of lists, or from an existing matrix::

        >>> fmpz_mat(2, 4)
        [0, 0, 0, 0]
        [0, 0, 0, 0]
        >>> A = fmpz_mat(2, 4, range(8))
        >>> A
        [0, 1, 2, 3]
        [4, 5, 6, 7]
        >>> fmpz_mat(A) == A
        True
        >>> fmpz_mat([[1,2,3],[4,5,6]])
        [1, 2, 3]
        [4, 5, 6]

    Entries can be accessed and set::

        >>> A[0,1]
        1
        >>> A[1,3] = 8
        >>> A
        [0, 1, 2, 3]
        [4, 5, 6, 8]

    Arithmetic operations are supported::

        >>> A * 3
        [ 0,  3,  6,  9]
        [12, 15, 18, 24]
        >>> A + (-A)
        [0, 0, 0, 0]
        [0, 0, 0, 0]
        >>> A * fmpz_mat(4, 3, range(12))
        [ 42,  48,  54]
        [123, 146, 169]
        >>> A * fmpz_mat(3, 3, range(9))
        Traceback (most recent call last):
          ...
        ValueError: incompatible shapes for matrix multiplication

    Disabling pretty-printing::

        >>> ctx.pretty = False
        >>> fmpz_mat([[1,2,3],[4,5,6]])
        fmpz_mat(2, 3, [1, 2, 3, 4, 5, 6])
        >>> ctx.pretty = True

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
            elif isinstance(val, (list, tuple)):
                m = len(val)
                n = 0
                if m != 0:
                    if not isinstance(val[0], (list, tuple)):
                        raise TypeError("single input to fmpz_mat must be a list of lists")
                    n = len(val[0])
                    for i from 1 <= i < m:
                        if len(val[i]) != n:
                            raise ValueError("input rows have different lengths")
                fmpz_mat_init(self.val, m, n)
                for i from 0 <= i < m:
                    row = val[i]
                    for j from 0 <= j < n:
                        x = fmpz(row[j])
                        fmpz_set(fmpz_mat_entry(self.val, i, j), (<fmpz>x).val)
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
        Returns the number of columns of *self*.
        """
        return fmpz_mat_ncols(self.val)

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

    def det(self):
        """
        Returns the determinant of *self* as an *fmpz*.

            >>> A = fmpz_mat(3, 3, range(9))
            >>> A.det()
            0
            >>> A[2,2] = 10
            >>> A.det()
            -6
            >>> (A * A).det()
            36
            >>> fmpz_mat(0, 0).det()
            1

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

    @staticmethod
    def _div_(fmpz_mat s, t):
        return s * (1 / fmpq(t))

    def __truediv__(s, t):
        return fmpz_mat._div_(s, t)

    def __div__(s, t):
        return fmpz_mat._div_(s, t)

    def __pow__(self, e, m):
        cdef fmpz_mat t
        cdef ulong ee
        if not typecheck(self, fmpz_mat):
            return NotImplemented
        if not fmpz_mat_is_square((<fmpz_mat>self).val):
            raise ValueError("matrix must be square")
        if m is not None:
            raise NotImplementedError("modular matrix exponentiation")
        ee = e
        t = fmpz_mat(self)   # XXX
        fmpz_mat_pow(t.val, t.val, ee)
        return t

    @classmethod
    def hadamard(cls, ulong n):
        """
        Attempts to construct a Hadamard matrix of size *n*.
        Raises :exc:`ValueError` if no such Hadamard matrix exists.
        The method can also raise :exc:`ValueError` if it fails to
        construct a Hadamard matrix of the specified size, even if one exists.
        It always succeeds if *n* is a power of two.

            >>> fmpz_mat.hadamard(1)
            [1]
            >>> fmpz_mat.hadamard(2)
            [1,  1]
            [1, -1]
            >>> fmpz_mat.hadamard(4)
            [1,  1,  1,  1]
            [1, -1,  1, -1]
            [1,  1, -1, -1]
            [1, -1, -1,  1]
            >>> fmpz_mat.hadamard(12)
            [ 1, -1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1]
            [-1, -1,  1, -1,  1, -1,  1, -1,  1, -1,  1, -1]
            [ 1,  1,  1, -1,  1,  1, -1, -1, -1, -1,  1,  1]
            [ 1, -1, -1, -1,  1, -1, -1,  1, -1,  1,  1, -1]
            [ 1,  1,  1,  1,  1, -1,  1,  1, -1, -1, -1, -1]
            [ 1, -1,  1, -1, -1, -1,  1, -1, -1,  1, -1,  1]
            [ 1,  1, -1, -1,  1,  1,  1, -1,  1,  1, -1, -1]
            [ 1, -1, -1,  1,  1, -1, -1, -1,  1, -1, -1,  1]
            [ 1,  1, -1, -1, -1, -1,  1,  1,  1, -1,  1,  1]
            [ 1, -1, -1,  1, -1,  1,  1, -1, -1, -1,  1, -1]
            [ 1,  1,  1,  1, -1, -1, -1, -1,  1,  1,  1, -1]
            [ 1, -1,  1, -1, -1,  1, -1,  1,  1, -1, -1, -1]
            >>> fmpz_mat.hadamard(10)
            Traceback (most recent call last):
              ...
            ValueError: unable to construct Hadamard matrix of size 10

        """
        cdef fmpz_mat res = fmpz_mat(n, n)
        if fmpz_mat_hadamard(res.val):
            return res
        else:
            raise ValueError("unable to construct Hadamard matrix of size %i" % n)

    def is_hadamard(self):
        """
        Determines whether *self* is a Hadamard matrix.

            >>> (fmpz_mat.hadamard(20) * -1).is_hadamard()
            True
            >>> (fmpz_mat.hadamard(20) * 2).is_hadamard()
            False
        """
        return bool(fmpz_mat_is_hadamard(self.val))

    @classmethod
    def randtest(cls, ulong m, ulong n, ulong bits):
        """
        Returns a random (*m*, *n*) matrix with non-uniformly chosen
        entries up to the specified number of bits in size. Small and
        zero entries are generated with increased probability.

            >>> fmpz_mat.randtest(3, 2, 100)   # doctest: +SKIP
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
        Returns a random (*m*, *n*) matrix with uniformly chosen entries up
        to the specified number of bits in size.

            >>> fmpz_mat.randbits(3, 2, 100)   # doctest: +SKIP
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
        Returns a random sparse (*m*, *n*) matrix of the specified rank
        with entries up to the specified number of bits in size.

            >>> fmpz_mat.randrank(3,6,2,20)   # doctest: +SKIP
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
        Returns the rank of *self*.

            >>> A = fmpz_mat(3, 3, range(9))
            >>> A.rank()
            2
            >>> A[2,2] = 10
            >>> A.rank()
            3
        """
        return fmpz_mat_rank(self.val)

    def inv(self, bint integer=False):
        """
        Returns the inverse matrix of *self*, which by default will
        be of type *fmpq_mat*::

            >>> fmpz_mat(3,3,[1,2,4,0,1,1,2,-1,0]).inv()
            [-1/3,  4/3,  2/3]
            [-2/3,  8/3,  1/3]
            [ 2/3, -5/3, -1/3]

        If *integer* is set, returns the inverse as an *fmpz_mat*,
        or raises an exception if the matrix is not invertible
        over the integers.

            >>> fmpz_mat([[5,1],[19,4]]).inv(integer=True)
            [  4, -1]
            [-19,  5]
            >>> fmpz_mat([[5,1],[19,5]]).inv(integer=True)
            Traceback (most recent call last):
              ...
            ValueError: matrix is not invertible over the integers

        """
        cdef fmpz_mat_t tmp
        cdef fmpz_mat v
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
            if integer:
                if not fmpz_is_pm1(den):
                    raise ValueError("matrix is not invertible over the integers")
                v = fmpz_mat.__new__(fmpz_mat)
                fmpz_mat_init_set(v.val, tmp)
                return v
            else:
                u = fmpq_mat.__new__(fmpq_mat)
                fmpq_mat_init(u.val, fmpz_mat_nrows(self.val), fmpz_mat_ncols(self.val))
                fmpq_mat_set_fmpz_mat_div_fmpz(u.val, tmp, den)
                return u
        finally:
            fmpz_clear(den)
            fmpz_mat_clear(tmp)

    def transpose(self):
        """
        Returns the transpose of *self*.

            >>> fmpz_mat(2,3,range(6)).transpose()
            [0, 3]
            [1, 4]
            [2, 5]
        """
        cdef fmpz_mat u
        u = fmpz_mat.__new__(fmpz_mat)
        fmpz_mat_init(u.val, fmpz_mat_ncols(self.val), fmpz_mat_nrows(self.val))
        fmpz_mat_transpose(u.val, self.val)
        return u

    def solve(self, other, bint integer=False):
        """
        Given matrices *A* and *B* represented by *self* and *other*,
        returns an *fmpq_mat* *X* such that `AX = B`, assuming that
        *A* is square and invertible.

        If *integer* is *True*, returns an *fmpz_mat*, solving the
        system only if the system matrix is invertible over the integers.
        (Warning: solving with *integer* set to *True* is
        currently slow for large matrices.)

            >>> A = fmpz_mat(2, 2, [1,4,8,3])
            >>> B = fmpz_mat(2, 3, range(6))
            >>> A.solve(B)
            [12/29, 13/29, 14/29]
            [-3/29,  4/29, 11/29]
            >>> A.solve(B, integer=True)
            Traceback (most recent call last):
              ...
            ValueError: matrix is not invertible over the integers
            >>> fmpz_mat([[1,2], [3,5]]).solve(B, integer=True)
            [ 6,  3, 0]
            [-3, -1, 1]
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
        cdef fmpq_mat v
        cdef fmpz d
        cdef int result
        t = any_as_fmpz_mat(other)
        if t is NotImplemented:
            raise TypeError("cannot convert input to fmpz_mat")
        if (fmpz_mat_nrows(self.val) != fmpz_mat_ncols(self.val) or
            fmpz_mat_nrows(self.val) != fmpz_mat_nrows((<fmpz_mat>t).val)):
            raise ValueError("need a square system and compatible right hand side")
        if not integer:
            v = fmpq_mat(fmpz_mat_nrows((<fmpz_mat>t).val), fmpz_mat_ncols((<fmpz_mat>t).val))
            result = fmpq_mat_solve_fmpz_mat(v.val, self.val, (<fmpz_mat>t).val)
            if not result:
                raise ZeroDivisionError("singular matrix in solve()")
            return v
        else:
            u = fmpz_mat.__new__(fmpz_mat)
            fmpz_mat_init(u.val, fmpz_mat_nrows((<fmpz_mat>t).val),
                fmpz_mat_ncols((<fmpz_mat>t).val))
            d = fmpz.__new__(fmpz)
            result = fmpz_mat_solve(u.val, d.val, self.val, (<fmpz_mat>t).val)
            if not fmpz_is_pm1(d.val):
                raise ValueError("matrix is not invertible over the integers")
            u *= d
            if not result:
                raise ZeroDivisionError("singular matrix in solve()")
            return u

    def rref(self, inplace=False):
        """
        Computes the reduced row echelon form (rref) of *self*,
        either returning a new copy or modifying self in-place.
        Returns (*rref*, *denominator*, *rank*).

            >>> ctx.pretty = False
            >>> A = fmpz_mat(3,3,range(9))
            >>> A.rref()
            (fmpz_mat(3, 3, [3, 0, -3, 0, 3, 6, 0, 0, 0]), fmpz(3), 2)
            >>> A.rref(inplace=True)
            (fmpz_mat(3, 3, [3, 0, -3, 0, 3, 6, 0, 0, 0]), fmpz(3), 2)
            >>> A
            fmpz_mat(3, 3, [3, 0, -3, 0, 3, 6, 0, 0, 0])
            >>> ctx.pretty = True

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
        Computes a basis for the nullspace of the matrix *A* represented
        by *self*. Returns (*X*, *nullity*) where nullity is the rank of
        the nullspace of *A* and *X* is a matrix whose first *nullity*
        columns are linearly independent, and such that `AX = 0`.

            >>> A = fmpz_mat(3,5,range(1,16))
            >>> X, nullity = A.nullspace()
            >>> A.rank(), nullity, X.rank()
            (2, 3, 3)
            >>> A * X
            [0, 0, 0, 0, 0]
            [0, 0, 0, 0, 0]
            [0, 0, 0, 0, 0]
            >>> X
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

    def lll(self, bint transform=False, double delta=0.99, double eta=0.51, rep="zbasis", gram="approx"):
        r"""
        Returns the LLL reduction of *self*, optionally along with
        a transformation matrix.

            >>> M = fmpz_mat([[11,17],[13,19]])
            >>> M.lll()
            [ 2, 2]
            [-3, 3]
            >>> L, T = M.lll(transform=True)
            >>> T * M == L
            True

        """
        cdef fmpz_mat u, v
        cdef fmpz_lll_t ctx
        cdef long i
        cdef int rt
        cdef int gt
        if rep == "zbasis":
            rt = 1
        elif rep == "gram":
            rt = 0
        else:
            raise ValueError("rep must be 'zbasis' or 'gram'")
        if gram == "approx":
            gt = 0
        elif gram == "exact":
            gt = 1
        else:
            raise ValueError("gram must be 'approx' or 'exact'")
        fmpz_lll_context_init(ctx, delta, eta, rt, gt)
        u = fmpz_mat(self)
        if transform:
            v = fmpz_mat(self.nrows(), self.nrows())
            for 0 <= i < self.nrows():
                v[i,i] = 1
            fmpz_lll(u.val, v.val, ctx)
            return u, v
        else:
            fmpz_lll(u.val, NULL, ctx)
            return u

    def hnf(self, bint transform=False):
        """
        Returns the Hermite normal form of *self*, optionally
        with a transformation matrix.

            >>> A = fmpz_mat(3,4,range(12))
            >>> A.hnf()
            [4, 0, -4, -8]
            [0, 1,  2,  3]
            [0, 0,  0,  0]
            >>> H, T = A.hnf(transform=True)
            >>> H == T * A
            True
        """
        cdef fmpz_mat H
        cdef fmpz_mat U
        H = fmpz_mat(self.nrows(), self.ncols())
        if transform:
            U = fmpz_mat(self.nrows(), self.nrows())
            fmpz_mat_hnf_transform(H.val, U.val, self.val)
            return H, U
        else:
            fmpz_mat_hnf(H.val, self.val)
            return H

    def is_hnf(self):
        """
        Determines whether *self* is in Hermite normal form.

            >>> fmpz_mat(3,4,range(12)).is_hnf()
            False
            >>> fmpz_mat(3,4,range(12)).hnf().is_hnf()
            True
        """
        return bool(fmpz_mat_is_in_hnf(self.val))

    def snf(self):
        """
        Returns the Smith normal form of *self*.

            >>> A = fmpz_mat(3,4,range(12))
            >>> A.snf()
            [1, 0, 0, 0]
            [0, 4, 0, 0]
            [0, 0, 0, 0]
        """
        cdef fmpz_mat H
        H = fmpz_mat(self.nrows(), self.ncols())
        fmpz_mat_snf(H.val, self.val)
        return H

    def is_snf(self):
        """
        Determines whether *self* is in Smith normal form.

            >>> fmpz_mat(3,4,range(12)).is_snf()
            False
            >>> fmpz_mat(3,4,range(12)).snf().is_snf()
            True
        """
        return bool(fmpz_mat_is_in_snf(self.val))

    def charpoly(self):
        cdef fmpz_poly u
        u = fmpz_poly.__new__(fmpz_poly)
        fmpz_poly_init(u.val)
        fmpz_mat_charpoly(u.val, self.val)
        return u

    def minpoly(self):
        cdef fmpz_poly u
        u = fmpz_poly.__new__(fmpz_poly)
        fmpz_poly_init(u.val)
        fmpz_mat_minpoly(u.val, self.val)
        return u

