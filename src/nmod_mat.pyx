cdef any_as_nmod_mat(obj, nmod_t mod):
    cdef nmod_mat r
    cdef mp_limb_t v
    if typecheck(obj, nmod_mat):
        return obj
    x = any_as_fmpz_mat(obj)
    if x is not NotImplemented:
        r = nmod_mat.__new__(nmod_mat)
        nmod_mat_init(r.val, fmpz_mat_nrows((<fmpz_mat>x).val),
                             fmpz_mat_ncols((<fmpz_mat>x).val), mod.n)
        fmpz_mat_get_nmod_mat(r.val, (<fmpz_mat>x).val)
        return r
    return NotImplemented

cdef class nmod_mat:
    """
    The nmod_poly type represents dense matrices over Z/nZ for
    word-size n. Some operations may assume that n is a prime.
    """

    cdef nmod_mat_t val

    def __dealloc__(self):
        nmod_mat_clear(self.val)

    @cython.embedsignature(False)
    def __init__(self, *args):
        cdef long m, n, i, j
        cdef mp_limb_t mod
        if len(args) == 1:
            val = args[0]
            if typecheck(val, nmod_mat):
                nmod_mat_init_set(self.val, (<nmod_mat>val).val)
                return
        mod = args[-1]
        args = args[:-1]
        if mod == 0:
            raise ValueError("modulus must be nonzero")
        if len(args) == 1:
            val = args[0]
            if typecheck(val, fmpz_mat):
                nmod_mat_init(self.val, fmpz_mat_nrows((<fmpz_mat>val).val),
                    fmpz_mat_ncols((<fmpz_mat>val).val), mod)
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
                nmod_mat_init(self.val, m, n, mod)
                for i from 0 <= i < m:
                    row = val[i]
                    for j from 0 <= j < n:
                        x = nmod(row[j], mod)
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
                    x = nmod(entries[i*n + j], mod)         # XXX: slow
                    self.val.rows[i][j] = (<nmod>x).val
        else:
            raise ValueError("nmod_mat: expected 1-3 arguments plus modulus")

    def __nonzero__(self):
        return not nmod_mat_is_zero(self.val)

    def __richcmp__(s, t, int op):
        cdef mp_limb_t v
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

    def __repr__(self):
        if ctx.pretty:
            return str(self)
        return "nmod_mat(%i, %i, %s, %i)" % (self.nrows(), self.ncols(),
            [int(c) for c in self.entries()], self.modulus())

    def __str__(self):
        return matrix_to_str(self.table())

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

    def table(self):
        cdef long i, m, n
        m = self.nrows()
        n = self.ncols()
        L = self.entries()
        return [L[i*n:(i+1)*n] for i in range(m)]

    def __getitem__(self, index):
        cdef long i, j
        cdef nmod x
        i, j = index
        if i < 0 or i >= self.nrows() or j < 0 or j >= self.ncols():
            raise ValueError("index %i,%i exceeds matrix dimensions" % (i, j))
        x = nmod(nmod_mat_entry(self.val, i, j), self.modulus()) # XXX: slow
        return x

    def __setitem__(self, index, value):
        cdef long i, j
        cdef mp_limb_t v
        i, j = index
        if i < 0 or i >= self.nrows() or j < 0 or j >= self.ncols():
            raise ValueError("index %i,%i exceeds matrix dimensions" % (i, j))
        if any_as_nmod(&v, value, self.val.mod):
            self.val.rows[i][j] = v
        else:
            raise ValueError("cannot set item of type %s" % type(value))

    def det(self):
        """
        Returns the determinant of self as an nmod.

            >>> nmod_mat(2,2,[1,2,3,4],17).det()
            15

        """
        if not nmod_mat_is_square(self.val):
            raise ValueError("matrix must be square")
        return nmod(nmod_mat_det(self.val), self.modulus())

    def rank(self):
        return nmod_mat_rank(self.val)

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
        if typecheck(s, nmod_mat):
            sv = &(<nmod_mat>s).val[0]
            t = any_as_nmod_mat(t, sv.mod)
            if t is NotImplemented:
                return t
            tv = &(<nmod_mat>t).val[0]
        else:
            tv = &(<nmod_mat>t).val[0]
            s = any_as_nmod_mat(s, tv.mod)
            if s is NotImplemented:
                return s
            sv = &(<nmod_mat>s).val[0]
        if sv.mod.n != tv.mod.n:
            raise ValueError("cannot add nmod_mats with different moduli")
        if sv.r != tv.r or sv.c != tv.c:
            raise ValueError("incompatible shapes for matrix addition")
        r = nmod_mat.__new__(nmod_mat)
        nmod_mat_init(r.val, sv.r, sv.c, sv.mod.n)
        nmod_mat_add(r.val, sv, tv)
        return r

    def __sub__(s, t):
        cdef nmod_mat r
        cdef nmod_mat_struct *sv
        cdef nmod_mat_struct *tv
        if typecheck(s, nmod_mat):
            sv = &(<nmod_mat>s).val[0]
            t = any_as_nmod_mat(t, sv.mod)
            if t is NotImplemented:
                return t
            tv = &(<nmod_mat>t).val[0]
        else:
            tv = &(<nmod_mat>t).val[0]
            s = any_as_nmod_mat(s, tv.mod)
            if s is NotImplemented:
                return s
            sv = &(<nmod_mat>s).val[0]
        if sv.mod.n != tv.mod.n:
            raise ValueError("cannot subtract nmod_mats with different moduli")
        if sv.r != tv.r or sv.c != tv.c:
            raise ValueError("incompatible shapes for matrix subtraction")
        r = nmod_mat.__new__(nmod_mat)
        nmod_mat_init(r.val, sv.r, sv.c, sv.mod.n)
        nmod_mat_sub(r.val, sv, tv)
        return r

    cdef __mul_nmod(self, mp_limb_t c):
        cdef nmod_mat r = nmod_mat.__new__(nmod_mat)
        nmod_mat_init(r.val, self.val.r, self.val.c, self.val.mod.n)
        nmod_mat_scalar_mul(r.val, self.val, c)
        return r

    def __mul__(s, t):
        cdef nmod_mat r
        cdef nmod_mat_struct *sv
        cdef nmod_mat_struct *tv
        cdef mp_limb_t c
        if typecheck(s, nmod_mat):
            sv = &(<nmod_mat>s).val[0]
            u = any_as_nmod_mat(t, sv.mod)
            if u is NotImplemented:
                if any_as_nmod(&c, t, sv.mod):
                    return (<nmod_mat>s).__mul_nmod(c)
                return NotImplemented
            tv = &(<nmod_mat>u).val[0]
        else:
            tv = &(<nmod_mat>t).val[0]
            u = any_as_nmod_mat(s, tv.mod)
            if u is NotImplemented:
                if any_as_nmod(&c, s, tv.mod):
                    return (<nmod_mat>t).__mul_nmod(c)
                return NotImplemented
            sv = &(<nmod_mat>u).val[0]
        if sv.mod.n != tv.mod.n:
            raise ValueError("cannot multiply nmod_mats with different moduli")
        if sv.c != tv.r:
            raise ValueError("incompatible shapes for matrix multiplication")
        r = nmod_mat.__new__(nmod_mat)
        nmod_mat_init(r.val, sv.r, tv.c, sv.mod.n)
        nmod_mat_mul(r.val, sv, tv)
        return r

    @staticmethod
    def _div_(nmod_mat s, t):
        cdef mp_limb_t v
        if not any_as_nmod(&v, t, s.val.mod):
            return NotImplemented
        t = nmod(v, s.val.mod.n)
        return s * (~t)

    def __truediv__(s, t):
        return nmod_mat._div_(s, t)

    def __div__(s, t):
        return nmod_mat._div_(s, t)

    def inv(self):
        cdef nmod_mat u
        if not nmod_mat_is_square(self.val):
            raise ValueError("matrix must be square")
        u = nmod_mat.__new__(nmod_mat)
        nmod_mat_init(u.val, nmod_mat_nrows(self.val),
            nmod_mat_ncols(self.val), self.val.mod.n)
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
        u = nmod_mat.__new__(nmod_mat)
        nmod_mat_init(u.val, nmod_mat_ncols(self.val),
            nmod_mat_nrows(self.val), self.val.mod.n)
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
        t = any_as_nmod_mat(other, self.val.mod)
        if t is NotImplemented:
            raise TypeError("cannot convert input to nmod_mat")
        if (nmod_mat_nrows(self.val) != nmod_mat_ncols(self.val) or
            nmod_mat_nrows(self.val) != nmod_mat_nrows((<nmod_mat>t).val)):
            raise ValueError("need a square system and compatible right hand side")
        u = nmod_mat.__new__(nmod_mat)
        nmod_mat_init(u.val, nmod_mat_nrows((<nmod_mat>t).val),
            nmod_mat_ncols((<nmod_mat>t).val), self.val.mod.n)
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
        if inplace:
            res = self
        else:
            res = nmod_mat.__new__(nmod_mat)
            nmod_mat_init_set((<nmod_mat>res).val, self.val)
        rank = nmod_mat_rref((<nmod_mat>res).val)
        return res, rank

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
        res = nmod_mat.__new__(nmod_mat)
        nmod_mat_init(res.val, nmod_mat_ncols(self.val), nmod_mat_ncols(self.val), self.val.mod.n)
        nullity = nmod_mat_nullspace(res.val, self.val)
        return res, nullity

