cdef arb_poly_coerce_operands(x, y):
    if typecheck(x, arb_poly):
        if isinstance(y, (int, long, float, fmpz, fmpq, arb, fmpz_poly, fmpq_poly)):
            return x, arb_poly(y)
        if isinstance(y, (complex, acb)):
            return acb_poly(x), acb_poly(y)
    else:
        if isinstance(x, (int, long, float, fmpz, fmpq, arb, fmpz_poly, fmpq_poly)):
            return arb_poly(x), y
        if isinstance(x, (complex, acb)):
            return acb_poly(x), acb_poly(y)
    return NotImplemented, NotImplemented

cdef arb_poly_set_list(arb_poly_t poly, list val, long prec):
    cdef long i, n
    cdef arb_t x
    n = PyList_GET_SIZE(<PyObject*>val)
    arb_poly_fit_length(poly, n)
    arb_init(x)
    for i in range(n):
        if typecheck(val[i], arb):
            arb_poly_set_coeff_arb(poly, i, (<arb>(val[i])).val)
        elif arb_set_python(x, val[i], 1):
            arb_poly_set_coeff_arb(poly, i, x)
        else:
            arb_clear(x)
            raise TypeError("unsupported coefficient type for arb_poly")
    arb_clear(x)

cdef class arb_poly(flint_poly):

    cdef arb_poly_t val

    def __cinit__(self):
        arb_poly_init(self.val)

    def __dealloc__(self):
        arb_poly_clear(self.val)

    def __init__(self, val=None):
        if val is not None:
            if typecheck(val, arb_poly):
                arb_poly_set(self.val, (<arb_poly>val).val)
            elif typecheck(val, fmpz_poly):
                arb_poly_set_fmpz_poly(self.val, (<fmpz_poly>val).val, getprec())
            elif typecheck(val, fmpq_poly):
                arb_poly_set_fmpq_poly(self.val, (<fmpq_poly>val).val, getprec())
            elif typecheck(val, list):
                arb_poly_set_list(self.val, val, getprec())
            else:
                val = any_as_arb(val)
                arb_poly_set_arb(self.val, (<arb>val).val)

    def __len__(self):
        return arb_poly_length(self.val)

    cpdef long length(self):
        return arb_poly_length(self.val)

    cpdef long degree(self):
        return arb_poly_degree(self.val)

    def __getitem__(self, long i):
        cdef arb x
        x = arb()
        if i >= 0:
            arb_poly_get_coeff_arb(x.val, self.val, i)
        return x

    def __setitem__(self, long i, x):
        if i < 0:
            raise ValueError("cannot assign to index < 0 of polynomial")
        if not typecheck(x, arb):
            x = arb(x)
        arb_poly_set_coeff_arb(self.val, i, (<arb>x).val)

    def repr(self):
        return "arb_poly([%s])" % (", ".join(map(str, self)))

    @classmethod
    def from_roots(cls, roots):
        """
        Constructs the monic polynomial whose roots are the given real numbers.

            >>> arb_poly.from_roots(range(4))
            1.00000000000000*x^4 + (-6.00000000000000)*x^3 + 11.0000000000000*x^2 + (-6.00000000000000)*x

        There is currently no dedicated method to construct a real polynomial
        from complex conjugate roots (use :meth:`.acb_poly.from_roots`).
        """
        cdef arb_ptr xs
        cdef long i, n
        roots = [any_as_arb(x) for x in roots]
        n = len(roots)
        xs = <arb_ptr> libc.stdlib.malloc(sizeof(arb_struct) * n)
        for i in range(n):
            xs[i] = (<arb>(roots[i])).val[0]
        u = arb_poly.__new__(arb_poly)
        arb_poly_product_roots((<arb_poly>u).val, xs, n, getprec())
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
        cdef arb_ptr xsv, ysv
        cdef long i, n
        assert algorithm in ('fast', 'iter')
        xs = [any_as_arb(x) for x in xs]
        n = len(xs)
        if n == 0:
            return []
        ys = [arb.__new__(arb) for i in range(n)]
        xsv = <arb_ptr> libc.stdlib.malloc(sizeof(arb_struct) * n)
        ysv = <arb_ptr> libc.stdlib.malloc(sizeof(arb_struct) * n)
        for i in range(n):
            xsv[i] = (<arb>(xs[i])).val[0]
            ysv[i] = (<arb>(ys[i])).val[0]
        if algorithm == 'fast':
            arb_poly_evaluate_vec_fast(ysv, (<arb_poly>self).val, xsv, n, getprec())
        else:
            arb_poly_evaluate_vec_iter(ysv, (<arb_poly>self).val, xsv, n, getprec())
        for i in range(n):  # move values
            (<arb>(ys[i])).val[0] = ysv[i]
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
        cdef arb_ptr xsv, ysv
        cdef long i, n
        assert algorithm in ('fast', 'newton', 'barycentric')
        xs = [any_as_arb(x) for x in xs]
        ys = [any_as_arb(y) for y in ys]
        n = len(xs)
        if len(ys) != n:
            raise ValueError("xs and ys must have the same length")
        xsv = <arb_ptr> libc.stdlib.malloc(sizeof(arb_struct) * n)
        ysv = <arb_ptr> libc.stdlib.malloc(sizeof(arb_struct) * n)
        for i in range(n):
            xsv[i] = (<arb>(xs[i])).val[0]
            ysv[i] = (<arb>(ys[i])).val[0]
        u = arb_poly.__new__(arb_poly)
        if algorithm == 'fast':
            arb_poly_interpolate_fast((<arb_poly>u).val, xsv, ysv, n, getprec())
        elif algorithm == 'newton':
            arb_poly_interpolate_newton((<arb_poly>u).val, xsv, ysv, n, getprec())
        else:
            arb_poly_interpolate_barycentric((<arb_poly>u).val, xsv, ysv, n, getprec())
        libc.stdlib.free(xsv)
        libc.stdlib.free(ysv)
        return u

    def derivative(self):
        u = arb_poly.__new__(arb_poly)
        arb_poly_derivative((<arb_poly>u).val, (<arb_poly>self).val, getprec())
        return u

    def integral(self):
        u = arb_poly.__new__(arb_poly)
        arb_poly_integral((<arb_poly>u).val, (<arb_poly>self).val, getprec())
        return u

    def __pos__(self):
        return self # ?

    def __neg__(s):
        u = arb_poly.__new__(arb_poly)
        arb_poly_neg((<arb_poly>u).val, (<arb_poly>s).val)
        return u

    def __add__(s, t):
        if type(s) is type(t):
            u = arb_poly.__new__(arb_poly)
            arb_poly_add((<arb_poly>u).val, (<arb_poly>s).val, (<arb_poly>t).val, getprec())
            return u
        s, t = arb_poly_coerce_operands(s, t)
        if s is NotImplemented:
            return s
        return s + t

    def __sub__(s, t):
        if type(s) is type(t):
            u = arb_poly.__new__(arb_poly)
            arb_poly_sub((<arb_poly>u).val, (<arb_poly>s).val, (<arb_poly>t).val, getprec())
            return u
        s, t = arb_poly_coerce_operands(s, t)
        if s is NotImplemented:
            return s
        return s - t

    def __mul__(s, t):
        if type(s) is type(t):
            u = arb_poly.__new__(arb_poly)
            arb_poly_mul((<arb_poly>u).val, (<arb_poly>s).val, (<arb_poly>t).val, getprec())
            return u
        s, t = arb_poly_coerce_operands(s, t)
        if s is NotImplemented:
            return s
        return s * t

    def __floordiv__(s, t):
        if type(s) is type(t):
            q = arb_poly.__new__(arb_poly)
            r = arb_poly.__new__(arb_poly)
            if arb_poly_divrem((<arb_poly>q).val, (<arb_poly>r).val,
                    (<arb_poly>s).val, (<arb_poly>t).val, getprec()):
                return q
            else:
                raise ZeroDivisionError("arb_poly leading coefficient must be nonzero")
        s, t = arb_poly_coerce_operands(s, t)
        if s is NotImplemented:
            return s
        return s // t

    def __mod__(s, t):
        if type(s) is type(t):
            q = arb_poly.__new__(arb_poly)
            r = arb_poly.__new__(arb_poly)
            if arb_poly_divrem((<arb_poly>q).val, (<arb_poly>r).val,
                    (<arb_poly>s).val, (<arb_poly>t).val, getprec()):
                return r
            else:
                raise ZeroDivisionError("arb_poly leading coefficient must be nonzero")
        s, t = arb_poly_coerce_operands(s, t)
        if s is NotImplemented:
            return s
        return s % t

    def __divmod__(s, t):
        if type(s) is type(t):
            q = arb_poly.__new__(arb_poly)
            r = arb_poly.__new__(arb_poly)
            if arb_poly_divrem((<arb_poly>q).val, (<arb_poly>r).val,
                    (<arb_poly>s).val, (<arb_poly>t).val, getprec()):
                return q, r
            else:
                raise ZeroDivisionError("arb_poly leading coefficient must be nonzero")
        s, t = arb_poly_coerce_operands(s, t)
        if s is NotImplemented:
            return s
        return divmod(s, t)

    def __pow__(arb_poly s, ulong exp, mod):
        if mod is not None:
            raise NotImplementedError("arb_poly modular exponentiation")
        u = arb_poly.__new__(arb_poly)
        arb_poly_pow_ui((<arb_poly>u).val, (<arb_poly>s).val, exp, getprec())
        return u

    def __call__(s, t):
        if typecheck(t, arb_poly):
            u = arb_poly.__new__(arb_poly)
            arb_poly_compose((<arb_poly>u).val, (<arb_poly>s).val, (<arb_poly>t).val, getprec())
            return u
        if typecheck(t, arb):
            u = arb.__new__(arb)
            arb_poly_evaluate((<arb>u).val, (<arb_poly>s).val, (<arb>t).val, getprec())
            return u
        if typecheck(t, acb):
            u = acb.__new__(acb)
            arb_poly_evaluate_acb((<acb>u).val, (<arb_poly>s).val, (<acb>t).val, getprec())
            return u
        if isinstance(t, (int, long, float, fmpz, fmpq)):
            return s(arb(t))
        if isinstance(t, (fmpz_poly, fmpq_poly)):
            return s(arb_poly(t))
        if isinstance(t, (complex)):
            return s(acb(t))
        raise TypeError("cannot call arb_poly with input of type %s", type(t))

    def unique_fmpz_poly(self):
        u = fmpz_poly.__new__(fmpz_poly)
        if arb_poly_get_unique_fmpz_poly((<fmpz_poly>u).val, self.val):
            return u
        else:
            return None

