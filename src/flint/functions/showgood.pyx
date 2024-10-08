from flint.flint_base.flint_context cimport thectx
from flint.utils.conversion cimport dps_to_prec
from flint.types.arb_mat cimport arb_mat
from flint.types.acb_mat cimport acb_mat
from flint.types.arb_poly cimport arb_poly
from flint.types.acb_poly cimport acb_poly
from flint.types.arb cimport arb
from flint.types.acb cimport acb
from flint.types.arb_series cimport arb_series
from flint.types.acb_series cimport acb_series
from flint.flintlib.functions.arb cimport *
from flint.flintlib.functions.acb cimport *

ctx = thectx

cdef slong __goodness(x, bint parts, metric):
    if metric is not None:
        x = metric(x)
    if isinstance(x, arb):
        return arb_rel_accuracy_bits((<arb>x).val)
    if isinstance(x, acb):
        if parts:
            return min(__goodness(x.real, parts, metric), __goodness(x.imag, parts, metric))
        else:
            return acb_rel_accuracy_bits((<acb>x).val)
    if isinstance(x, (tuple, list, arb_mat, acb_mat, arb_poly, acb_poly, arb_series, acb_series)):
        return min(__goodness(y, parts, metric) for y in x)
    raise TypeError("must have arb or acb")

cdef goodstr(x):
    if isinstance(x, tuple):
        return "(" + ", ".join(goodstr(y) for y in x) + ")"
    if isinstance(x, list):
        return "[" + ", ".join(goodstr(y) for y in x) + "]"
    if isinstance(x, (arb, acb, arb_mat, acb_mat, arb_poly, acb_poly, arb_series, acb_series)):
        return x.str(radius=False)
    raise TypeError("must be a element/tuple/list of arb, acb, arb_mat, acb_mat, arb_poly, acb_poly, arb_series, or acb_series")


def good(func, slong prec=0, slong maxprec=0, slong dps=0,
         slong maxdps=0, slong padding=10, bint verbose=False,
         bint show=False, bint parts=True, metric=None):
    """
    Evaluates *func*, automatically increasing the precision to get
    a result accurate to the current working precision (or the
    precision specified by *prec* or *dps*).

        >>> from flint import arb
        >>> good(lambda: (arb.pi() + arb("1e-100")).sin())
        Traceback (most recent call last):
          ...
        ValueError: no convergence (maxprec=630, try higher maxprec)
        >>> good(lambda: (arb.pi() + arb("1e-100")).sin(), maxprec=1000)
        [-1.00000000000000e-100 +/- 3e-119]

    The function *func* can return an *arb*, an *acb*, or a composite
    object such as a tuple or a matrix. By default all real and
    imaginary parts of all components must be accurate. This means
    that convergence is not possible in case of inexact zeros.
    This behavior can be overridden by setting *parts* to *False*.

        >>> from flint import acb
        >>> good(lambda: (acb(0,-1) ** 0.5) ** 2)
        Traceback (most recent call last):
          ...
        ValueError: no convergence (maxprec=630, try higher maxprec)
        >>> good(lambda: (acb(0,-1) ** 0.5) ** 2, parts=False)
        [+/- 4.50e-22] + [-1.00000000000000 +/- 3e-20]j


    """
    cdef slong orig, morebits, acc

    if dps > 0:
        prec = dps_to_prec(dps)
    if maxdps > 0:
        maxprec = dps_to_prec(maxdps)

    if prec == 0:
        prec = ctx.prec
    if maxprec == 0:
        maxprec = 10 * prec + 100

    if metric == "abssum":
        def metric(L): return sum(abs(c) for c in L)

    # for printing
    if dps == 0:
        dps = ctx.dps

    orig = ctx.prec
    morebits = max(20, prec // 5)

    if verbose:
        print("prec = %i, maxprec = %i" % (prec, maxprec))
    try:
        ctx.prec = prec * 1.01 + 2 * padding
        while ctx.prec < maxprec:
            if verbose:
                print("eval prec = %i" % ctx.prec)
            v = func()
            acc = __goodness(v, parts, metric)
            if verbose:
                print("good bits = %i" % acc)
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
    """
    Evaluates *func* accurately with :func:`good`, printing the decimal
    value of the result (without an explicit radius) instead
    of returning it.

        >>> from flint import arb
        >>> showgood(lambda: arb.pi())
        3.14159265358979
        >>> showgood(lambda: arb.pi(), dps=50)
        3.1415926535897932384626433832795028841971693993751
    """
    kwargs["show"] = True
    good(func, **kwargs)
