cdef any_as_nmod_poly(obj, nmod_t mod):
    cdef nmod_poly r
    cdef mp_limb_t v
    # XXX: should check that modulus is the same here, and not all over the place
    if typecheck(obj, nmod_poly):
        return obj
    if any_as_nmod(&v, obj, mod):
        r = nmod_poly.__new__(nmod_poly)
        nmod_poly_init(r.val, mod.n)
        nmod_poly_set_coeff_ui(r.val, 0, v)
        return r
    x = any_as_fmpz_poly(obj)
    if x is not NotImplemented:
        r = nmod_poly.__new__(nmod_poly)
        nmod_poly_init(r.val, mod.n)   # XXX: create flint _nmod_poly_set_modulus for this?
        fmpz_poly_get_nmod_poly(r.val, (<fmpz_poly>x).val)
        return r
    return NotImplemented

cdef nmod_poly_set_list(nmod_poly_t poly, list val):
    cdef long i, n
    cdef nmod_t mod
    cdef mp_limb_t v
    nmod_init(&mod, nmod_poly_modulus(poly)) # XXX
    n = PyList_GET_SIZE(<PyObject*>val)
    nmod_poly_fit_length(poly, n)
    for i from 0 <= i < n:
        c = val[i]
        if any_as_nmod(&v, val[i], mod):
            nmod_poly_set_coeff_ui(poly, i, v)
        else:
            raise TypeError("unsupported coefficient in list")

cdef class nmod_poly(flint_poly):
    """
    The nmod_poly type represents dense univariate polynomials
    over Z/nZ for word-size n.

        >>> a = nmod_poly([5,1,10,14,8], 7)
        >>> a
        x^4 + 3*x^2 + x + 5
        >>> -a
        6*x^4 + 4*x^2 + 6*x + 2
        >>> ctx.pretty = False
        >>> list(nmod_poly(list(range(3)), 2))
        [nmod(0, 2), nmod(1, 2)]
        >>> nmod_poly([1, 2, 3], 23) ** 3
        nmod_poly([1, 6, 21, 21, 17, 8, 4], 23)
        >>> divmod(nmod_poly([2,0,1,1,6],23), nmod_poly([3,5,7],23))
        (nmod_poly([4, 0, 14], 23), nmod_poly([13, 3], 23))
        >>> ctx.pretty = True

    """

    cdef nmod_poly_t val

    #def __cinit__(self):

    def __dealloc__(self):
        nmod_poly_clear(self.val)

    def __init__(self, val=None, ulong mod=0):
        cdef ulong m2
        if typecheck(val, nmod_poly):
            m2 = nmod_poly_modulus((<nmod_poly>val).val)
            if m2 != mod:
                raise ValueError("different moduli!")
            nmod_poly_init(self.val, m2)
            nmod_poly_set(self.val, (<nmod_poly>val).val)
        else:
            if mod == 0:
                raise ValueError("a nonzero modulus is required")
            nmod_poly_init(self.val, mod)
            if typecheck(val, fmpz_poly):
                fmpz_poly_get_nmod_poly(self.val, (<fmpz_poly>val).val)
            elif typecheck(val, list):
                nmod_poly_set_list(self.val, val)
            else:
                raise TypeError("cannot create nmod_poly from input of type %s", type(val))

    def __len__(self):
        return nmod_poly_length(self.val)

    cpdef long length(self):
        return nmod_poly_length(self.val)

    cpdef long degree(self):
        return nmod_poly_degree(self.val)

    cpdef mp_limb_t modulus(self):
        return nmod_poly_modulus(self.val)

    def __richcmp__(s, t, int op):
        cdef mp_limb_t v
        cdef bint res
        if op != 2 and op != 3:
            raise TypeError("nmod_polyss cannot be ordered")
        if typecheck(s, nmod_poly) and typecheck(t, nmod_poly):
            if (<nmod_poly>s).val.mod.n != (<nmod_poly>t).val.mod.n:
                res = False
            else:
                res = nmod_poly_equal((<nmod_poly>s).val, (<nmod_poly>t).val)
            if op == 2:
                return res
            if op == 3:
                return not res
        return NotImplemented

    def __iter__(self):
        cdef long i, n
        n = self.length()
        for i from 0 <= i < n:
            yield self[i]

    def coeffs(self):
        cdef long i, n
        cdef list L
        cdef mp_limb_t m
        n = self.length()
        m = self.modulus()
        L = [nmod(0,m) for i in range(n)]   # XXX: speed up
        for i from 0 <= i < n:
            (<nmod>(L[i])).val = nmod_poly_get_coeff_ui(self.val, i)
        return L

    def repr(self):
        return "nmod_poly(%s, %s)" % ([int(c) for c in self.coeffs()], self.modulus())

    def __getitem__(self, long i):
        cdef nmod x
        x = nmod(0, self.modulus())
        if i < 0:
            return x
        x.val = nmod_poly_get_coeff_ui(self.val, i)
        return x

    def __setitem__(self, long i, x):
        cdef mp_limb_t v
        if i < 0:
            raise ValueError("cannot assign to index < 0 of polynomial")
        if any_as_nmod(&v, x, self.val.mod):
            nmod_poly_set_coeff_ui(self.val, i, v)
        else:
            raise TypeError("cannot set element of type %s" % type(x))

    def __nonzero__(self):
        return not nmod_poly_is_zero(self.val)

    def __call__(self, other):
        cdef mp_limb_t c
        if any_as_nmod(&c, other, self.val.mod):
            v = nmod(0, self.modulus())
            (<nmod>v).val = nmod_poly_evaluate_nmod(self.val, c)
            return v
        t = any_as_nmod_poly(other, self.val.mod)
        if t is not NotImplemented:
            r = nmod_poly.__new__(nmod_poly)
            nmod_poly_init_preinv((<nmod_poly>r).val, self.val.mod.n, self.val.mod.ninv)
            nmod_poly_compose((<nmod_poly>r).val, self.val, (<nmod_poly>t).val)
            return r
        raise TypeError("cannot call nmod_poly with input of type %s", type(other))

    def __pos__(self):
        return self

    def __neg__(self):
        cdef nmod_poly r = nmod_poly.__new__(nmod_poly)
        nmod_poly_init_preinv(r.val, self.val.mod.n, self.val.mod.ninv)
        nmod_poly_neg(r.val, self.val)
        return r

    def __add__(s, t):
        cdef nmod_poly r
        if typecheck(s, nmod_poly):
            t = any_as_nmod_poly(t, (<nmod_poly>s).val.mod)
            if t is NotImplemented:
                return t
        else:
            s = any_as_nmod_poly(s, (<nmod_poly>t).val.mod)
            if s is NotImplemented:
                return s
        if (<nmod_poly>s).val.mod.n != (<nmod_poly>t).val.mod.n:
            raise ValueError("cannot add nmod_polys with different moduli")
        r = nmod_poly.__new__(nmod_poly)
        nmod_poly_init_preinv(r.val, (<nmod_poly>t).val.mod.n, (<nmod_poly>t).val.mod.ninv)
        nmod_poly_add(r.val, (<nmod_poly>s).val, (<nmod_poly>t).val)
        return r

    def __sub__(s, t):
        cdef nmod_poly r
        if typecheck(s, nmod_poly):
            t = any_as_nmod_poly(t, (<nmod_poly>s).val.mod)
            if t is NotImplemented:
                return t
        else:
            s = any_as_nmod_poly(s, (<nmod_poly>t).val.mod)
            if s is NotImplemented:
                return s
        if (<nmod_poly>s).val.mod.n != (<nmod_poly>t).val.mod.n:
            raise ValueError("cannot subtract nmod_polys with different moduli")
        r = nmod_poly.__new__(nmod_poly)
        nmod_poly_init_preinv(r.val, (<nmod_poly>t).val.mod.n, (<nmod_poly>t).val.mod.ninv)
        nmod_poly_sub(r.val, (<nmod_poly>s).val, (<nmod_poly>t).val)
        return r

    def __mul__(s, t):
        cdef nmod_poly r
        if typecheck(s, nmod_poly):
            t = any_as_nmod_poly(t, (<nmod_poly>s).val.mod)
            if t is NotImplemented:
                return t
        else:
            s = any_as_nmod_poly(s, (<nmod_poly>t).val.mod)
            if s is NotImplemented:
                return s
        if (<nmod_poly>s).val.mod.n != (<nmod_poly>t).val.mod.n:
            raise ValueError("cannot multiply nmod_polys with different moduli")
        r = nmod_poly.__new__(nmod_poly)
        nmod_poly_init_preinv(r.val, (<nmod_poly>t).val.mod.n, (<nmod_poly>t).val.mod.ninv)
        nmod_poly_mul(r.val, (<nmod_poly>s).val, (<nmod_poly>t).val)
        return r

    # TODO: __div__, __truediv__

    def __floordiv__(s, t):
        cdef nmod_poly r
        if typecheck(s, nmod_poly):
            t = any_as_nmod_poly(t, (<nmod_poly>s).val.mod)
            if t is NotImplemented:
                return t
        else:
            s = any_as_nmod_poly(s, (<nmod_poly>t).val.mod)
            if s is NotImplemented:
                return s
        if (<nmod_poly>s).val.mod.n != (<nmod_poly>t).val.mod.n:
            raise ValueError("cannot divide nmod_polys with different moduli")
        if nmod_poly_is_zero((<nmod_poly>t).val):
            raise ZeroDivisionError("polynomial division by zero")
        r = nmod_poly.__new__(nmod_poly)
        nmod_poly_init_preinv(r.val, (<nmod_poly>t).val.mod.n, (<nmod_poly>t).val.mod.ninv)
        nmod_poly_div(r.val, (<nmod_poly>s).val, (<nmod_poly>t).val)
        return r

    def __divmod__(s, t):
        cdef nmod_poly P, Q
        if typecheck(s, nmod_poly):
            t = any_as_nmod_poly(t, (<nmod_poly>s).val.mod)
            if t is NotImplemented:
                return t
        else:
            s = any_as_nmod_poly(s, (<nmod_poly>t).val.mod)
            if s is NotImplemented:
                return s
        if (<nmod_poly>s).val.mod.n != (<nmod_poly>t).val.mod.n:
            raise ValueError("cannot divide nmod_polys with different moduli")
        if nmod_poly_is_zero((<nmod_poly>t).val):
            raise ZeroDivisionError("polynomial division by zero")
        P = nmod_poly.__new__(nmod_poly)
        Q = nmod_poly.__new__(nmod_poly)
        nmod_poly_init_preinv(P.val, (<nmod_poly>t).val.mod.n, (<nmod_poly>t).val.mod.ninv)
        nmod_poly_init_preinv(Q.val, (<nmod_poly>t).val.mod.n, (<nmod_poly>t).val.mod.ninv)
        nmod_poly_divrem(P.val, Q.val, (<nmod_poly>s).val, (<nmod_poly>t).val)
        return P, Q

    def __mod__(s, t):
        return divmod(s, t)[1]      # XXX

    def __pow__(nmod_poly self, ulong exp, mod):
        cdef nmod_poly res
        if mod is not None:
            raise NotImplementedError("nmod_poly modular exponentiation")
        res = nmod_poly.__new__(nmod_poly)
        nmod_poly_init_preinv(res.val, (<nmod_poly>self).val.mod.n, (<nmod_poly>self).val.mod.ninv)
        nmod_poly_pow(res.val, self.val, exp)
        return res

    def gcd(self, other):
        """
        Returns the monic greatest common divisor of self and other.

            >>> A = nmod_poly([1,2,3,4], 7); B = nmod_poly([4,1,5], 7)
            >>> (A * B).gcd(B) * 5
            5*x^2 + x + 4

        """
        cdef nmod_poly res
        other = any_as_nmod_poly(other, (<nmod_poly>self).val.mod)
        if other is NotImplemented:
            raise TypeError("cannot convert input to nmod_poly")
        if self.val.mod.n != (<nmod_poly>other).val.mod.n:
            raise ValueError("moduli must be the same")
        res = nmod_poly.__new__(nmod_poly)
        nmod_poly_init_preinv(res.val, self.val.mod.n, self.val.mod.ninv)
        nmod_poly_gcd(res.val, self.val, (<nmod_poly>other).val)
        return res

    def factor(self, algorithm=None):
        """
        Factors self into irreducible factors, returning a tuple
        (c, factors) where c is the leading coefficient and
        factors is a list of (poly, exp) pairs with all polynomials
        monic.

            >>> nmod_poly(list(range(10)), 3).factor()
            (2, [(x, 1), (x + 2, 7)])
            >>> nmod_poly(list(range(10)), 19).factor()
            (9, [(x, 1), (x^4 + 15*x^3 + 2*x^2 + 7*x + 3, 1), (x^4 + 7*x^3 + 12*x^2 + 15*x + 12, 1)])
            >>> nmod_poly(list(range(10)), 53).factor()
            (9, [(x, 1), (x^8 + 48*x^7 + 42*x^6 + 36*x^5 + 30*x^4 + 24*x^3 + 18*x^2 + 12*x + 6, 1)])

        Algorithm can be None (default), 'berlekamp', or
        'cantor-zassenhaus'.

            >>> nmod_poly([3,2,1,2,3], 7).factor(algorithm='berlekamp')
            (3, [(x + 2, 1), (x + 4, 1), (x^2 + 4*x + 1, 1)])
            >>> nmod_poly([3,2,1,2,3], 7).factor(algorithm='cantor-zassenhaus')
            (3, [(x + 4, 1), (x + 2, 1), (x^2 + 4*x + 1, 1)])

        """
        cdef nmod_poly_factor_t fac
        cdef mp_limb_t lead
        cdef int i
        nmod_poly_factor_init(fac)
        if algorithm == 'berlekamp':
            lead = nmod_poly_factor_with_berlekamp(fac, self.val)
        elif algorithm == 'cantor-zassenhaus':
            lead = nmod_poly_factor_with_cantor_zassenhaus(fac, self.val)
        else:
            lead = nmod_poly_factor(fac, self.val)
        res = [None] * fac.num
        for 0 <= i < fac.num:
            u = nmod_poly.__new__(nmod_poly)
            nmod_poly_init_preinv((<nmod_poly>u).val,
                (<nmod_poly>self).val.mod.n, (<nmod_poly>self).val.mod.ninv)
            nmod_poly_set((<nmod_poly>u).val, &fac.p[i])
            exp = fac.exp[i]
            res[i] = (u, exp)
        c = nmod.__new__(nmod)
        (<nmod>c).mod = self.val.mod
        (<nmod>c).val = lead
        nmod_poly_factor_clear(fac)
        return c, res

