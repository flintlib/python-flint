def goodness(x, bint bothcomplex=True, metric=None):
    if metric is not None:
        x = metric(x)
    if isinstance(x, arb):
        return arb_rel_accuracy_bits((<arb>x).val)
    if isinstance(x, acb):
        if bothcomplex:
            return min(goodness(x.real), goodness(x.imag))
        else:
            return acb_rel_accuracy_bits((<acb>x).val)
    if isinstance(x, (tuple, list, arb_mat, acb_mat, arb_poly, acb_poly, arb_series, acb_series)):
        return min(goodness(y, bothcomplex) for y in x)
    raise TypeError("must have arb or acb")

def goodstr(x):
    if isinstance(x, tuple):
        return "(" + ", ".join(goodstr(y) for y in x) + ")"
    if isinstance(x, list):
        return "[" + ", ".join(goodstr(y) for y in x) + "]"
    if isinstance(x, arb):
        return x.str(radius=False)
    if isinstance(x, acb):
        return x.str(radius=False)
    raise TypeError("must have arb or acb")

def good(func, long prec=0, long maxprec=0, long dps=0,
        long maxdps=0, long padding=10, bint verbose=False, bint show=False, bint bothcomplex=True, metric=None):
    cdef long orig, morebits, acc

    if dps > 0:
        prec = dps_to_prec(dps)
    if maxdps > 0:
        maxprec = dps_to_prec(maxdps)

    if prec == 0:
        prec = ctx.prec
    if maxprec == 0:
        maxprec = 10 * prec + 100

    if metric == "abssum":
        metric = lambda L: sum(abs(c) for c in L)

    # for printing
    if dps == 0:
        dps = ctx.dps

    orig = ctx.prec
    morebits = max(20, prec / 5)

    if verbose:
        print "prec = %i, maxprec = %i" % (prec, maxprec)
    try:
        ctx.prec = prec * 1.01 + 2 * padding
        while ctx.prec < maxprec:
            if verbose:
                print "eval prec = %i" % ctx.prec
            v = func()
            acc = goodness(v, bothcomplex, metric)
            if verbose:
                print "good bits = %i" % acc
            if acc > prec + padding:
                if show:
                    ctx.dps = dps
                    print(goodstr(v))
                    return
                else:
                    return v
            ctx.prec += morebits
            morebits *= 2
    finally:
        ctx.prec = orig
    raise ValueError("no convergence (maxprec=%i, try higher maxprec)" % maxprec)

def showgood(func, **kwargs):
    kwargs["show"] = True
    good(func, **kwargs)


def divisor_sigma(n, k):
    """
    Evaluates the divisor sum function sigma_k(n), returning an fmpz.

        >>> divisor_sigma(60, 0)
        12
        >>> divisor_sigma(60, 1)
        168
        >>> divisor_sigma(60, 10)
        605263138639095300
    """
    cdef fmpz v = fmpz()
    n = any_as_fmpz(n)
    if n is NotImplemented:
        raise TypeError("the input must be an integer")
    arith_divisor_sigma(v.val, (<fmpz>n).val, k)
    return v

def euler_phi(n):
    """
    Evaluates the Euler totient function phi(n), returning an fmpz.

        >>> euler_phi(60)
        16
        >>> euler_phi(3**10)
        39366
    """
    cdef fmpz v = fmpz()
    n = any_as_fmpz(n)
    if n is NotImplemented:
        raise TypeError("the input must be an integer")
    arith_euler_phi(v.val, (<fmpz>n).val)
    return v

def bell_number(n):
    """
    Returns the nth Bell number B_n as an fmpz.

        >>> [bell_number(n) for n in range(5)]
        [1, 1, 2, 5, 15]
        >>> bell_number(50)
        185724268771078270438257767181908917499221852770
    """
    cdef fmpz v = fmpz()
    arith_bell_number(v.val, n)
    return v

def euler_number(n):
    """
    Returns the nth Euler number E_n as an fmpz.

        >>> [euler_number(n) for n in range(6)]
        [1, 0, -1, 0, 5, 0]
        >>> euler_number(50)
        -6053285248188621896314383785111649088103498225146815121
    """
    cdef fmpz v = fmpz()
    arith_euler_number(v.val, n)
    return v

def stirling_number_1(n, k):
    """
    Returns the Stirling number of the first kind s(n, k) as an fmpz.

        >>> [stirling_number_1(5, k) for k in range(6)]
        [0, 24, -50, 35, -10, 1]
    """
    cdef fmpz v = fmpz()
    arith_stirling_number_1(v.val, n, k)
    return v

def stirling_number_2(n, k):
    """
    Returns the Stirling number of the second kind S(n, k) as an fmpz.

        >>> [stirling_number_2(5, k) for k in range(6)]
        [0, 1, 15, 25, 10, 1]
    """
    cdef fmpz v = fmpz()
    arith_stirling_number_2(v.val, n, k)
    return v

def bernoulli_polynomial(n):
    """
    Returns the nth Bernoulli polynomial B_n(x) as an fmpq_poly.

        >>> bernoulli_polynomial(2)
        x^2 + (-1)*x + 1/6
    """
    cdef fmpq_poly v = fmpq_poly()
    arith_bernoulli_polynomial(v.val, n)
    return v

def euler_polynomial(n):
    """
    Returns the nth Euler polynomial E_n(x) as an fmpq_poly.

        >>> euler_polynomial(3)
        x^3 + (-3/2)*x^2 + 1/4
    """
    cdef fmpq_poly v = fmpq_poly()
    arith_euler_polynomial(v.val, n)
    return v

def legendre_polynomial(n):
    """
    Returns the nth Legendre polynomial P_n(x) as an fmpq_poly.

        >>> legendre_polynomial(3)
        5/2*x^3 + (-3/2)*x

    """
    cdef fmpq_poly v = fmpq_poly()
    arith_legendre_polynomial(v.val, n)
    return v

def chebyshev_t_polynomial(n):
    """
    Returns the nth Chebyshev polynomial of the first kind T_n(x)
    as an fmpz_poly.

        >>> chebyshev_t_polynomial(3)
        4*x^3 + (-3)*x

    """
    cdef fmpz_poly v = fmpz_poly()
    arith_chebyshev_t_polynomial(v.val, n)
    return v

def chebyshev_u_polynomial(n):
    """
    Returns the nth Chebyshev polynomial of the second kind U_n(x)
    as an fmpz_poly.

        >>> chebyshev_u_polynomial(3)
        8*x^3 + (-4)*x

    """
    cdef fmpz_poly v = fmpz_poly()
    arith_chebyshev_u_polynomial(v.val, n)
    return v

