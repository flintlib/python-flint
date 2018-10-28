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
