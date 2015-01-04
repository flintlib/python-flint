cdef class arb_mat(flint_mat):
    """
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
            else:
                raise TypeError("cannot create arb_mat from input of type %s" % type(val))
        elif len(args) == 2:
            m, n = args
            arb_mat_init(self.val, m, n)
        elif len(args) == 3:
            m, n, entries = args
            arb_mat_init(self.val, m, n)
            entries = list(entries)
            if len(entries) != m*n:
                raise ValueError("list of entries has the wrong length")
            for i from 0 <= i < m:
                for j from 0 <= j < n:
                    x = arb(entries[i*n + j])   # slow
                    arb_set(arb_mat_entry(self.val, i, j), (<arb>x).val)
        else:
            raise ValueError("arb_mat: expected 1-3 arguments")

    def __nonzero__(self):
        raise NotImplementedError

    def __richcmp__(s, t, int op):
        return NotImplementedError

    cpdef long nrows(self):
        """
        Returns the number of rows of self.
        """
        return arb_mat_nrows(self.val)

    cpdef long ncols(self):
        """
        Returns the number of columns of self.
        """
        return arb_mat_ncols(self.val)

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
        c = arb(value)  # XXX
        arb_set(arb_mat_entry(self.val, i, j), (<arb>c).val)

    def det(s):
        """
        Returns the determinant of s as an arb.

        If the matrix is singular, the result will be an interval
        containing zero. As currently implemented, the width of this interval
        will not generally converge to zero as the precision increases.

            >>> A = arb_mat(3, 3, range(9))
            >>> showgood(lambda: A.det(), dps=25)    # singular
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
        s = arb_mat.convert_operand(s)
        if s is NotImplemented:
            return s
        t = arb_mat.convert_operand(t)
        if t is NotImplemented:
            return t
        m = (<arb_mat>s).nrows()
        n = (<arb_mat>s).ncols()
        if m != (<arb_mat>t).nrows() or n != (<arb_mat>t).ncols():
            raise ValueError("incompatible shapes for matrix addition")
        u = arb_mat.__new__(arb_mat)
        arb_mat_init((<arb_mat>u).val, m, n)
        arb_mat_add((<arb_mat>u).val, (<arb_mat>s).val, (<arb_mat>t).val, getprec())
        return u

    def __sub__(s, t):
        cdef long m, n
        s = arb_mat.convert_operand(s)
        if s is NotImplemented:
            return s
        t = arb_mat.convert_operand(t)
        if t is NotImplemented:
            return t
        m = (<arb_mat>s).nrows()
        n = (<arb_mat>s).ncols()
        if m != (<arb_mat>t).nrows() or n != (<arb_mat>t).ncols():
            raise ValueError("incompatible shapes for matrix subtraction")
        u = arb_mat.__new__(arb_mat)
        arb_mat_init((<arb_mat>u).val, m, n)
        arb_mat_sub((<arb_mat>u).val, (<arb_mat>s).val, (<arb_mat>t).val, getprec())
        return u

    cdef _mul_fmpz_(s, fmpz t):
        cdef arb_mat u
        u = arb_mat.__new__(arb_mat)
        arb_mat_init(u.val, arb_mat_nrows(s.val), arb_mat_ncols(s.val))
        arb_mat_scalar_mul_fmpz(u.val, s.val, t.val, getprec())
        return u

    cdef _mul_arb_(s, arb t):
        cdef arb_mat u
        u = arb_mat.__new__(arb_mat)
        arb_mat_init(u.val, arb_mat_nrows(s.val), arb_mat_ncols(s.val))
        arb_mat_scalar_mul_arb(u.val, s.val, t.val, getprec())
        return u

    cdef _mul_arb_mat_(s, arb_mat t):
        cdef arb_mat u
        if arb_mat_ncols(s.val) != arb_mat_nrows(t.val):
            raise ValueError("incompatible shapes for matrix multiplication")
        u = arb_mat.__new__(arb_mat)
        arb_mat_init(u.val, arb_mat_nrows(s.val), arb_mat_ncols(t.val))
        arb_mat_mul(u.val, s.val, t.val, getprec())
        return u

    def __mul__(s, t):
        if typecheck(s, arb_mat):
            if typecheck(t, arb_mat):
                return (<arb_mat>s)._mul_arb_mat_(t)
            if typecheck(t, fmpz):
                return (<arb_mat>s)._mul_fmpz_(t)
            if typecheck(t, arb):
                return (<arb_mat>s)._mul_arb_(t)
            if typecheck(t, acb):
                return acb_mat(s) * t
            if typecheck(t, (int, long)):
                return (<arb_mat>s)._mul_fmpz_(fmpz(t))
        else:
            if typecheck(s, fmpz):
                return (<arb_mat>t)._mul_fmpz_(s)
            if typecheck(s, arb):
                return (<arb_mat>t)._mul_arb_(s)
            if typecheck(s, acb):
                return acb_mat(t) * s
            if typecheck(s, (int, long)):
                return (<arb_mat>t)._mul_fmpz_(fmpz(s))
        s = arb_mat.convert_operand(s)
        if s is NotImplemented:
            return s
        t = arb_mat.convert_operand(t)
        if t is NotImplemented:
            return t
        return (<arb_mat>s)._mul_arb_mat(<arb_mat>t)

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

    def __invert__(s):
        cdef arb_mat u
        if arb_mat_nrows(s.val) != arb_mat_ncols(s.val):
            raise ValueError("matrix must be square")
        u = arb_mat.__new__(arb_mat)
        arb_mat_init(u.val, arb_mat_nrows(s.val), arb_mat_ncols(s.val))
        if not arb_mat_inv(u.val, s.val, getprec()):
            raise ZeroDivisionError("matrix is singular")
        return u

    # ???
    inv = __invert__

    def solve(s, t):
        """
        Solves `AX = B` where *A* is a square matrix given by *s* and
        `B` is a matrix given by *t*.
        Raises :exc:`ZeroDivisionError` if *A* is numerically singular.
        """
        cdef arb_mat u
        t = arb_mat.convert(t)
        if (arb_mat_nrows(s.val) != arb_mat_ncols(s.val) or
            arb_mat_nrows(s.val) != arb_mat_nrows((<arb_mat>t).val)):
            raise ValueError("need a square system and compatible right hand side")
        u = arb_mat.__new__(arb_mat)
        arb_mat_init(u.val, arb_mat_nrows((<arb_mat>t).val), arb_mat_ncols((<arb_mat>t).val))
        if not arb_mat_solve(u.val, s.val, (<arb_mat>t).val, getprec()):
            raise ZeroDivisionError("singular matrix in solve()")
        return u

    def exp(s):
        """
        Computes the matrix exponential of *s*.
        """
        cdef arb_mat u
        if arb_mat_nrows(s.val) != arb_mat_ncols(s.val):
            raise ValueError("matrix must be square")
        u = arb_mat.__new__(arb_mat)
        arb_mat_init(u.val, arb_mat_nrows(s.val), arb_mat_ncols(s.val))
        arb_mat_exp(u.val, s.val, getprec())
        return u

