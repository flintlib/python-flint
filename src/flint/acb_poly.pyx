cdef acb_poly_coerce_operands(x, y):
    if typecheck(x, acb_poly):
        if isinstance(y, (int, long, float, complex, fmpz, fmpq, arb, acb, fmpz_poly, fmpq_poly, arb_poly)):
            return x, acb_poly(y)
    else:
        if isinstance(x, (int, long, float, complex, fmpz, fmpq, arb, acb, fmpz_poly, fmpq_poly, arb_poly)):
            return acb_poly(x), y
    return NotImplemented, NotImplemented

cdef acb_poly_set_list(acb_poly_t poly, list val, long prec):
    cdef long i, n
    cdef acb_t x
    n = PyList_GET_SIZE(<PyObject*>val)
    acb_poly_fit_length(poly, n)
    acb_init(x)
    for i in range(n):
        if typecheck(val[i], acb):
            acb_poly_set_coeff_acb(poly, i, (<acb>(val[i])).val)
        elif acb_set_python(x, val[i], 1):
            acb_poly_set_coeff_acb(poly, i, x)
        else:
            acb_clear(x)
            raise TypeError("unsupported coefficient type for acb_poly")
    acb_clear(x)

cdef int acb_vec_check_tol(acb_srcptr vec, long n, tol) except -1:
    cdef acb t
    tol = arb(tol)
    for i in range(n):
        t = acb.__new__(acb)
        acb_set(t.val, &vec[i])
        if not (abs(t.rad()) <= tol):
            return 0
    return 1

cdef class acb_poly(flint_poly):

    cdef acb_poly_t val

    def __cinit__(self):
        acb_poly_init(self.val)

    def __dealloc__(self):
        acb_poly_clear(self.val)

    def __init__(self, val=None):
        if val is not None:
            if typecheck(val, acb_poly):
                acb_poly_set(self.val, (<acb_poly>val).val)
            elif typecheck(val, arb_poly):
                acb_poly_set_arb_poly(self.val, (<arb_poly>val).val)
            elif typecheck(val, fmpz_poly):
                acb_poly_set_fmpz_poly(self.val, (<fmpz_poly>val).val, getprec())
            elif typecheck(val, fmpq_poly):
                acb_poly_set_fmpq_poly(self.val, (<fmpq_poly>val).val, getprec())
            elif typecheck(val, list):
                acb_poly_set_list(self.val, val, getprec())
            else:
                val = any_as_acb(val)
                acb_poly_set_acb(self.val, (<acb>val).val)

    def __len__(self):
        return acb_poly_length(self.val)

    cpdef long length(self):
        return acb_poly_length(self.val)

    cpdef long degree(self):
        return acb_poly_degree(self.val)

    def __getitem__(self, long i):
        cdef acb x
        x = acb()
        if i >= 0:
            acb_poly_get_coeff_acb(x.val, self.val, i)
        return x

    def __setitem__(self, long i, x):
        if i < 0:
            raise ValueError("cannot assign to index < 0 of polynomial")
        if not typecheck(x, acb):
            x = acb(x)
        acb_poly_set_coeff_acb(self.val, i, (<acb>x).val)

    def repr(self):
        return "acb_poly([%s])" % (", ".join(map(str, self)))

    @classmethod
    def from_roots(cls, roots):
        """
        Constructs the monic polynomial whose roots are the given complex numbers.

            >>> acb_poly.from_roots(range(4))
            1.00000000000000*x^4 + (-6.00000000000000)*x^3 + 11.0000000000000*x^2 + (-6.00000000000000)*x
        """
        cdef acb_ptr xs
        cdef long i, n
        roots = [any_as_acb(x) for x in roots]
        n = len(roots)
        xs = <acb_ptr> libc.stdlib.malloc(sizeof(acb_struct) * n)
        for i in range(n):
            xs[i] = (<acb>(roots[i])).val[0]
        u = acb_poly.__new__(acb_poly)
        acb_poly_product_roots((<acb_poly>u).val, xs, n, getprec())
        libc.stdlib.free(xs)
        return u

    def evaluate(self, xs, algorithm='fast'):
        """
        Multipoint evaluation: evaluates *self* at the list of
        points *xs*. The algorithm can be 'iter' or 'fast'. The 'fast'
        algorithm is asymptotically fast, but has worse numerical stability.

        Note: for ordinary single-point evaluation, just call the polynomial
        with the point as the argument.
        """
        cdef acb_ptr xsv, ysv
        cdef long i, n
        assert algorithm in ('fast', 'iter')
        xs = [any_as_acb(x) for x in xs]
        n = len(xs)
        if n == 0:
            return []
        ys = [acb.__new__(acb) for i in range(n)]
        xsv = <acb_ptr> libc.stdlib.malloc(sizeof(acb_struct) * n)
        ysv = <acb_ptr> libc.stdlib.malloc(sizeof(acb_struct) * n)
        for i in range(n):
            xsv[i] = (<acb>(xs[i])).val[0]
            ysv[i] = (<acb>(ys[i])).val[0]
        if algorithm == 'fast':
            acb_poly_evaluate_vec_fast(ysv, (<acb_poly>self).val, xsv, n, getprec())
        else:
            acb_poly_evaluate_vec_iter(ysv, (<acb_poly>self).val, xsv, n, getprec())
        for i in range(n):  # move values
            (<acb>(ys[i])).val[0] = ysv[i]
        libc.stdlib.free(xsv)
        libc.stdlib.free(ysv)
        return ys

    @classmethod
    def interpolate(cls, xs, ys, algorithm='fast'):
        """
        Constructs the unique interpolating polynomial of length at most *n*
        taking the values *ys* when evaluated at the *n* distinct points *xs*.
        Algorithm can be 'newton', 'barycentric' or 'fast'. The 'fast'
        algorithm is asymptotically fast, but has worse numerical stability.
        """
        cdef acb_ptr xsv, ysv
        cdef long i, n
        assert algorithm in ('fast', 'newton', 'barycentric')
        xs = [any_as_acb(x) for x in xs]
        ys = [any_as_acb(y) for y in ys]
        n = len(xs)
        if len(ys) != n:
            raise ValueError("xs and ys must have the same length")
        xsv = <acb_ptr> libc.stdlib.malloc(sizeof(acb_struct) * n)
        ysv = <acb_ptr> libc.stdlib.malloc(sizeof(acb_struct) * n)
        for i in range(n):
            xsv[i] = (<acb>(xs[i])).val[0]
            ysv[i] = (<acb>(ys[i])).val[0]
        u = acb_poly.__new__(acb_poly)
        if algorithm == 'fast':
            acb_poly_interpolate_fast((<acb_poly>u).val, xsv, ysv, n, getprec())
        elif algorithm == 'newton':
            acb_poly_interpolate_newton((<acb_poly>u).val, xsv, ysv, n, getprec())
        else:
            acb_poly_interpolate_barycentric((<acb_poly>u).val, xsv, ysv, n, getprec())
        libc.stdlib.free(xsv)
        libc.stdlib.free(ysv)
        return u

    def derivative(self):
        u = acb_poly.__new__(acb_poly)
        acb_poly_derivative((<acb_poly>u).val, (<acb_poly>self).val, getprec())
        return u

    def integral(self):
        u = acb_poly.__new__(acb_poly)
        acb_poly_integral((<acb_poly>u).val, (<acb_poly>self).val, getprec())
        return u

    def __pos__(self):
        return self # ?

    def __neg__(s):
        u = acb_poly.__new__(acb_poly)
        acb_poly_neg((<acb_poly>u).val, (<acb_poly>s).val)
        return u

    def __add__(s, t):
        if type(s) is type(t):
            u = acb_poly.__new__(acb_poly)
            acb_poly_add((<acb_poly>u).val, (<acb_poly>s).val, (<acb_poly>t).val, getprec())
            return u
        s, t = acb_poly_coerce_operands(s, t)
        if s is NotImplemented:
            return s
        return s + t

    def __sub__(s, t):
        if type(s) is type(t):
            u = acb_poly.__new__(acb_poly)
            acb_poly_sub((<acb_poly>u).val, (<acb_poly>s).val, (<acb_poly>t).val, getprec())
            return u
        s, t = acb_poly_coerce_operands(s, t)
        if s is NotImplemented:
            return s
        return s - t

    def __mul__(s, t):
        if type(s) is type(t):
            u = acb_poly.__new__(acb_poly)
            acb_poly_mul((<acb_poly>u).val, (<acb_poly>s).val, (<acb_poly>t).val, getprec())
            return u
        s, t = acb_poly_coerce_operands(s, t)
        if s is NotImplemented:
            return s
        return s * t

    def __floordiv__(s, t):
        if type(s) is type(t):
            q = acb_poly.__new__(acb_poly)
            r = acb_poly.__new__(acb_poly)
            if acb_poly_divrem((<acb_poly>q).val, (<acb_poly>r).val,
                    (<acb_poly>s).val, (<acb_poly>t).val, getprec()):
                return q
            else:
                raise ZeroDivisionError("acb_poly leading coefficient must be nonzero")
        s, t = acb_poly_coerce_operands(s, t)
        if s is NotImplemented:
            return s
        return s // t

    def __mod__(s, t):
        if type(s) is type(t):
            q = acb_poly.__new__(acb_poly)
            r = acb_poly.__new__(acb_poly)
            if acb_poly_divrem((<acb_poly>q).val, (<acb_poly>r).val,
                    (<acb_poly>s).val, (<acb_poly>t).val, getprec()):
                return r
            else:
                raise ZeroDivisionError("acb_poly leading coefficient must be nonzero")
        s, t = acb_poly_coerce_operands(s, t)
        if s is NotImplemented:
            return s
        return s % t

    def __divmod__(s, t):
        if type(s) is type(t):
            q = acb_poly.__new__(acb_poly)
            r = acb_poly.__new__(acb_poly)
            if acb_poly_divrem((<acb_poly>q).val, (<acb_poly>r).val,
                    (<acb_poly>s).val, (<acb_poly>t).val, getprec()):
                return q, r
            else:
                raise ZeroDivisionError("acb_poly leading coefficient must be nonzero")
        s, t = acb_poly_coerce_operands(s, t)
        if s is NotImplemented:
            return s
        return divmod(s, t)

    def __pow__(acb_poly s, ulong exp, mod):
        if mod is not None:
            raise NotImplementedError("acb_poly modular exponentiation")
        u = acb_poly.__new__(acb_poly)
        acb_poly_pow_ui((<acb_poly>u).val, (<acb_poly>s).val, exp, getprec())
        return u

    def __call__(s, t):
        if typecheck(t, acb_poly):
            u = acb_poly.__new__(acb_poly)
            acb_poly_compose((<acb_poly>u).val, (<acb_poly>s).val, (<acb_poly>t).val, getprec())
            return u
        if typecheck(t, acb):
            u = acb.__new__(acb)
            acb_poly_evaluate((<acb>u).val, (<acb_poly>s).val, (<acb>t).val, getprec())
            return u
        if isinstance(t, (int, long, float, complex, fmpz, fmpq, arb)):
            return s(acb(t))
        if isinstance(t, (fmpz_poly, fmpq_poly, arb_poly)):
            return s(acb_poly(t))
        raise TypeError("cannot call acb_poly with input of type %s", type(t))

    def unique_fmpz_poly(self):
        u = fmpz_poly.__new__(fmpz_poly)
        if acb_poly_get_unique_fmpz_poly((<fmpz_poly>u).val, self.val):
            return u
        else:
            return None

    def roots(s, tol=None, maxprec=None):
        """
        Attempts to isolate all the complex roots of *s*.
        If *tol* is specified, the roots are further refined
        to at least the requested tolerance.
        The input polynomial must be squarefree and sufficiently
        accurate. Raises an exception if unsuccessful.

            >>> for c in acb_poly.from_roots([1,2,3,4,5]).roots(1e-10): print(c)
            ...
            [1.00000000000000 +/- 7.71e-18] + [+/- 7.59e-18]j
            [2.00000000000000 +/- 1.18e-16] + [+/- 1.16e-16]j
            [3.00000000000000 +/- 4.24e-16] + [+/- 4.22e-16]j
            [4.00000000000000 +/- 7.21e-16] + [+/- 7.00e-16]j
            [5.00000000000000 +/- 2.41e-16] + [+/- 2.24e-16]j
            >>> for c in acb_poly.from_roots([1,2,2,4,5]).roots(1e-10): print(c)
            ...
            Traceback (most recent call last):
              ...
            ValueError: roots() failed to converge: insufficient precision, or squareful input

        """
        cdef long prec, initial_prec, target_prec, isolated, maxiter, deg, i
        cdef acb_ptr roots
        cdef acb_poly_t tmp
        deg = s.degree()
        if deg < 0:
            return []

        roots = _acb_vec_init(deg)
        acb_poly_init(tmp)
        pyroots = []

        try:
            initial_prec = 32
            prec = initial_prec

            if maxprec is None:
                maxprec = 2 * getprec()

            while 1:
                maxiter = min(max(deg, 32), prec)
                acb_poly_set_round(tmp, s.val, prec)

                if prec == initial_prec:
                    isolated = acb_poly_find_roots(roots, tmp, NULL, maxiter, prec)
                else:
                    isolated = acb_poly_find_roots(roots, tmp, roots, maxiter, prec)

                if isolated == deg:
                    if tol is None:
                        break
                    elif acb_vec_check_tol(roots, deg, tol) == 1:
                        break

                prec *= 2

                if prec > maxprec:
                    raise ValueError("roots() failed to converge: insufficient precision, or squareful input")

            _acb_vec_sort_pretty(roots, deg)
            pyroots = [acb.__new__(acb) for i in range(deg)]
            for i in range(deg):
                acb_set((<acb>(pyroots[i])).val, &roots[i])

        finally:
            acb_poly_clear(tmp)
            _acb_vec_clear(roots, deg)

        return pyroots

    def root_bound(self):
        """Returns an upper bound for the absolute value of
        the roots of self."""
        cdef mag_t m
        mag_init(m)
        r = arb()
        acb_poly_root_bound_fujiwara(m, self.val)
        arf_set_mag(arb_midref((<arb>r).val), m)
        mag_clear(m)
        return r

