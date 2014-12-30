def goodness(x):
    if isinstance(x, tuple) or isinstance(x, list):
        return min(goodness(y) for y in x)
    assert isinstance(x, arb)
    return arb_rel_accuracy_bits((<arb>x).val)

def goodstr(x):
    if isinstance(x, tuple):
        return "(" + ", ".join(goodstr(y) for y in x) + ")"
    if isinstance(x, list):
        return "[" + ", ".join(goodstr(y) for y in x) + "]"
    assert isinstance(x, arb)
    return str(x.mid())

def good(func, long prec=0, long maxprec=0, long dps=0,
        long maxdps=0, long padding=10, bint verbose=False, bint show=False):
    cdef long orig, morebits, acc

    if dps > 0:
        prec = dps_to_prec(dps)
    if maxdps > 0:
        maxprec = dps_to_prec(maxdps)

    if prec == 0:
        prec = ctx.prec
    if maxprec == 0:
        maxprec = 10 * prec + 100

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
            acc = goodness(v)
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
    raise ValueError("no convergence")

def showgood(func, **kwargs):
    kwargs["show"] = True
    good(func, **kwargs)


def divisor_sigma(n, k):
    """
    Evaluates the divisor sum function sigma_k(n), returning an fmpz.

        >>> divisor_sigma(60, 0)
        fmpz(12)
        >>> divisor_sigma(60, 1)
        fmpz(168)
        >>> divisor_sigma(60, 10)
        fmpz(605263138639095300)
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
        fmpz(16)
        >>> euler_phi(3**10)
        fmpz(39366)
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
        [fmpz(1), fmpz(1), fmpz(2), fmpz(5), fmpz(15)]
        >>> bell_number(50)
        fmpz(185724268771078270438257767181908917499221852770)
    """
    cdef fmpz v = fmpz()
    arith_bell_number(v.val, n)
    return v

def euler_number(n):
    """
    Returns the nth Euler number E_n as an fmpz.

        >>> [euler_number(n) for n in range(6)]
        [fmpz(1), fmpz(0), fmpz(-1), fmpz(0), fmpz(5), fmpz(0)]
        >>> euler_number(50)
        fmpz(-6053285248188621896314383785111649088103498225146815121)
    """
    cdef fmpz v = fmpz()
    arith_euler_number(v.val, n)
    return v

def stirling_number_1(n, k):
    """
    Returns the Stirling number of the first kind s(n, k) as an fmpz.

        >>> [stirling_number_1(5, k) for k in range(6)]
        [fmpz(0), fmpz(24), fmpz(-50), fmpz(35), fmpz(-10), fmpz(1)]
    """
    cdef fmpz v = fmpz()
    arith_stirling_number_1(v.val, n, k)
    return v

def stirling_number_2(n, k):
    """
    Returns the Stirling number of the second kind S(n, k) as an fmpz.

        >>> [stirling_number_2(5, k) for k in range(6)]
        [fmpz(0), fmpz(1), fmpz(15), fmpz(25), fmpz(10), fmpz(1)]
    """
    cdef fmpz v = fmpz()
    arith_stirling_number_2(v.val, n, k)
    return v

def bernoulli_polynomial(n):
    """
    Returns the nth Bernoulli polynomial B_n(x) as an fmpq_poly.

        >>> bernoulli_polynomial(2)
        fmpq_poly([1, -6, 6], 6)
    """
    cdef fmpq_poly v = fmpq_poly()
    arith_bernoulli_polynomial(v.val, n)
    return v

def euler_polynomial(n):
    """
    Returns the nth Euler polynomial E_n(x) as an fmpq_poly.

        >>> euler_polynomial(3)
        fmpq_poly([1, 0, -6, 4], 4)
    """
    cdef fmpq_poly v = fmpq_poly()
    arith_euler_polynomial(v.val, n)
    return v

def legendre_polynomial(n):
    """
    Returns the nth Legendre polynomial P_n(x) as an fmpq_poly.

        >>> legendre_polynomial(3)
        fmpq_poly([0, -3, 0, 5], 2)

    """
    cdef fmpq_poly v = fmpq_poly()
    arith_legendre_polynomial(v.val, n)
    return v

def chebyshev_t_polynomial(n):
    """
    Returns the nth Chebyshev polynomial of the first kind T_n(x)
    as an fmpz_poly.

        >>> chebyshev_t_polynomial(3)
        fmpz_poly([0, -3, 0, 4])

    """
    cdef fmpz_poly v = fmpz_poly()
    arith_chebyshev_t_polynomial(v.val, n)
    return v

def chebyshev_u_polynomial(n):
    """
    Returns the nth Chebyshev polynomial of the second kind U_n(x)
    as an fmpz_poly.

        >>> chebyshev_u_polynomial(3)
        fmpz_poly([0, -4, 0, 8])

    """
    cdef fmpz_poly v = fmpz_poly()
    arith_chebyshev_u_polynomial(v.val, n)
    return v

