cdef any_as_fmpz_poly(x):
    cdef fmpz_poly res
    if typecheck(x, fmpz_poly):
        return x
    elif typecheck(x, fmpz):
        res = fmpz_poly.__new__(fmpz_poly)
        fmpz_poly_set_fmpz(res.val, (<fmpz>x).val)
        return res
    elif PyInt_Check(<PyObject*>x):
        res = fmpz_poly.__new__(fmpz_poly)
        fmpz_poly_set_si(res.val, PyInt_AS_LONG(<PyObject*>x))
        return res
    elif PyLong_Check(<PyObject*>x):
        res = fmpz_poly.__new__(fmpz_poly)
        t = fmpz(x)
        fmpz_poly_set_fmpz(res.val, (<fmpz>t).val)
        return res
    return NotImplemented

cdef fmpz_poly_set_list(fmpz_poly_t poly, list val):
    cdef long i, n
    cdef fmpz_t x
    n = PyList_GET_SIZE(<PyObject*>val)
    fmpz_poly_fit_length(poly, n)
    fmpz_init(x)
    for i from 0 <= i < n:
        if typecheck(val[i], fmpz):
            fmpz_poly_set_coeff_fmpz(poly, i, (<fmpz>(val[i])).val)
        elif fmpz_set_python(x, val[i]):
            fmpz_poly_set_coeff_fmpz(poly, i, x)
        else:
            fmpz_clear(x)
            raise TypeError("unsupported coefficient in list")
    fmpz_clear(x)

cdef class fmpz_poly(flint_poly):
    """
    The fmpz_poly type represents dense univariate polynomials over
    the integers.

        >>> fmpz_poly([1,2,3]) ** 3
        fmpz_poly([1, 6, 21, 44, 63, 54, 27])
        >>> divmod(fmpz_poly([2,0,1,1,6]), fmpz_poly([3,5,7]))
        (fmpz_poly([]), fmpz_poly([2, 0, 1, 1, 6]))
    """

    cdef fmpz_poly_t val

    def __cinit__(self):
        fmpz_poly_init(self.val)

    def __dealloc__(self):
        fmpz_poly_clear(self.val)

    def __init__(self, val=None):
        if val is not None:
            if typecheck(val, fmpz_poly):
                fmpz_poly_set(self.val, (<fmpz_poly>val).val)
            elif isinstance(val, list):
                fmpz_poly_set_list(self.val, val)
            else:
                raise TypeError("cannot create fmpz_poly from input of type %s", type(val))

    def __len__(self):
        return fmpz_poly_length(self.val)

    cpdef long length(self):
        return fmpz_poly_length(self.val)

    cpdef long degree(self):
        return fmpz_poly_degree(self.val)

    def __richcmp__(self, other, int op):
        cdef bint r
        if op != 2 and op != 3:
            raise TypeError("polynomials cannot be ordered")
        self = any_as_fmpz_poly(self)
        if self is NotImplemented:
            return self
        other = any_as_fmpz_poly(other)
        if other is NotImplemented:
            return other
        r = fmpz_poly_equal((<fmpz_poly>self).val, (<fmpz_poly>other).val)
        if op == 3:
            r = not r
        return r

    def __getitem__(self, long i):
        cdef fmpz x
        x = fmpz()
        if i < 0:
            return x
        fmpz_poly_get_coeff_fmpz(x.val, self.val, i)
        return x

    def __setitem__(self, long i, x):
        if i < 0:
            raise ValueError("cannot assign to index < 0 of polynomial")
        v = fmpz(x)  # XXX
        fmpz_poly_set_coeff_fmpz(self.val, i, (<fmpz>v).val)

    def repr(self):
        return "fmpz_poly([%s])" % (", ".join(map(str, self.coeffs())))

    def __nonzero__(self):
        return not fmpz_poly_is_zero(self.val)

    def __call__(self, other):
        t = any_as_fmpz(other)
        if t is not NotImplemented:
            v = fmpz.__new__(fmpz)
            fmpz_poly_evaluate_fmpz((<fmpz>v).val, self.val, (<fmpz>t).val)
            return v
        t = any_as_fmpz_poly(other)
        if t is not NotImplemented:
            v = fmpz_poly.__new__(fmpz_poly)
            fmpz_poly_compose((<fmpz_poly>v).val, self.val, (<fmpz_poly>t).val)
            return v
        t = any_as_fmpq(other)
        if t is not NotImplemented:
            return fmpq_poly(self)(t)
        t = any_as_fmpq_poly(other)
        if t is not NotImplemented:
            return fmpq_poly(self)(t)
        raise TypeError("cannot call fmpz_poly with input of type %s", type(other))

    def __pos__(self):
        return self

    def __neg__(self):
        cdef fmpz_poly res = fmpz_poly.__new__(fmpz_poly)
        fmpz_poly_neg(res.val, self.val)
        return res

    def __add__(self, other):
        cdef fmpz_poly res
        self = any_as_fmpz_poly(self)
        if self is NotImplemented:
            return self
        other = any_as_fmpz_poly(other)
        if other is NotImplemented:
            return other
        res = fmpz_poly.__new__(fmpz_poly)
        fmpz_poly_add(res.val, (<fmpz_poly>self).val, (<fmpz_poly>other).val)
        return res

    def __sub__(self, other):
        cdef fmpz_poly res
        self = any_as_fmpz_poly(self)
        if self is NotImplemented:
            return self
        other = any_as_fmpz_poly(other)
        if other is NotImplemented:
            return other
        res = fmpz_poly.__new__(fmpz_poly)
        fmpz_poly_sub(res.val, (<fmpz_poly>self).val, (<fmpz_poly>other).val)
        return res

    def __mul__(self, other):
        cdef fmpz_poly res
        self = any_as_fmpz_poly(self)
        if self is NotImplemented:
            return self
        other = any_as_fmpz_poly(other)
        if other is NotImplemented:
            return other
        res = fmpz_poly.__new__(fmpz_poly)
        fmpz_poly_mul(res.val, (<fmpz_poly>self).val, (<fmpz_poly>other).val)
        return res

    def __floordiv__(self, other):
        cdef fmpz_poly res
        self = any_as_fmpz_poly(self)
        if self is NotImplemented:
            return self
        other = any_as_fmpz_poly(other)
        if other is NotImplemented:
            return other
        if fmpz_poly_is_zero((<fmpz_poly>other).val):
            raise ZeroDivisionError("fmpz_poly division by 0")
        res = fmpz_poly.__new__(fmpz_poly)
        fmpz_poly_div(res.val, (<fmpz_poly>self).val, (<fmpz_poly>other).val)
        return res

    def __mod__(self, other):
        cdef fmpz_poly res
        self = any_as_fmpz_poly(self)
        if self is NotImplemented:
            return self
        other = any_as_fmpz_poly(other)
        if other is NotImplemented:
            return other
        if fmpz_poly_is_zero((<fmpz_poly>other).val):
            raise ZeroDivisionError("fmpz_poly division by 0")
        res = fmpz_poly.__new__(fmpz_poly)
        fmpz_poly_rem(res.val, (<fmpz_poly>self).val, (<fmpz_poly>other).val)
        return res

    def __divmod__(self, other):
        cdef fmpz_poly P, Q
        self = any_as_fmpz_poly(self)
        if self is NotImplemented:
            return self
        other = any_as_fmpz_poly(other)
        if other is NotImplemented:
            return other
        if fmpz_poly_is_zero((<fmpz_poly>other).val):
            raise ZeroDivisionError("fmpz_poly divmod by 0")
        P = fmpz_poly.__new__(fmpz_poly)
        Q = fmpz_poly.__new__(fmpz_poly)
        fmpz_poly_divrem(P.val, Q.val, (<fmpz_poly>self).val, (<fmpz_poly>other).val)
        return P, Q

    def __pow__(fmpz_poly self, ulong exp, mod):
        cdef fmpz_poly res
        if mod is not None:
            raise NotImplementedError("fmpz_poly modular exponentiation")
        res = fmpz_poly.__new__(fmpz_poly)
        fmpz_poly_pow(res.val, self.val, exp)
        return res

    def gcd(self, other):
        """
        Returns the greatest common divisor of self and other.

            >>> A = fmpz_poly([2,0,1,0,5]); B = fmpz_poly([2,3,4])
            >>> (A*B).gcd(B)
            fmpz_poly([2, 3, 4])

        """
        cdef fmpz_poly res
        other = any_as_fmpz_poly(other)
        if other is NotImplemented:
            raise TypeError("cannot convert input to fmpz_poly")
        res = fmpz_poly.__new__(fmpz_poly)
        fmpz_poly_gcd(res.val, self.val, (<fmpz_poly>other).val)
        return res

    def factor(self):
        """
        Factors self into irreducible factors, returning a tuple
        (c, factors) where c is the content of the coefficients and
        factors is a list of (poly, exp) pairs.

            >>> (-73 * fmpz_poly([1,2,3]) ** 3 * fmpz_poly([5,6,7,8,9]) ** 8).factor()
            (fmpz(-73), [(fmpz_poly([1, 2, 3]), 3), (fmpz_poly([5, 6, 7, 8, 9]), 8)])
            >>> chebyshev_t_polynomial(6).factor()
            (fmpz(1), [(fmpz_poly([-1, 0, 2]), 1), (fmpz_poly([1, 0, -16, 0, 16]), 1)])
            >>> (fmpz_poly([-1,1])**100).factor()
            (fmpz(1), [(fmpz_poly([-1, 1]), 100)])
            >>> fmpz_poly([1,2,3,4,5,6]).factor()
            (fmpz(1), [(fmpz_poly([1, 2, 3, 4, 5, 6]), 1)])

        """
        cdef fmpz_poly_factor_t fac
        cdef int i
        fmpz_poly_factor_init(fac)
        fmpz_poly_factor_zassenhaus(fac, self.val)
        res = [0] * fac.num
        for 0 <= i < fac.num:
            u = fmpz_poly.__new__(fmpz_poly)
            fmpz_poly_set((<fmpz_poly>u).val, &fac.p[i])
            exp = fac.exp[i]
            res[i] = (u, exp)
        c = fmpz.__new__(fmpz)
        fmpz_set((<fmpz>c).val, &fac.c)
        fmpz_poly_factor_clear(fac)
        return c, res

    @staticmethod
    def cyclotomic_ui(ulong n):
        """
        Returns the nth cyclotomic polynomial Phi_n(x) as an fmpz_poly.

            >>> fmpz_poly.cyclotomic_ui(12)
            fmpz_poly([1, 0, -1, 0, 1])

        """
        u = fmpz_poly.__new__(fmpz_poly)
        fmpz_poly_cyclotomic((<fmpz_poly>u).val, n)
        return u

    @staticmethod
    def cos_minpoly_ui(ulong n):
        u = fmpz_poly.__new__(fmpz_poly)
        fmpz_poly_cos_minpoly((<fmpz_poly>u).val, n)
        return u

    @staticmethod
    def swinnerton_dyer_ui(ulong n):
        if n > 20:
            raise OverflowError("that's way too large...")
        u = fmpz_poly.__new__(fmpz_poly)
        fmpz_poly_swinnerton_dyer((<fmpz_poly>u).val, n)
        return u

