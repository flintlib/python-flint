cdef int any_as_nmod(mp_limb_t * val, obj, nmod_t mod) except -1:
    cdef int success
    cdef fmpz_t t
    if typecheck(obj, nmod):
        if (<nmod>obj).mod.n != mod.n:
            raise ValueError("cannot coerce integers mod n with different n")
        val[0] = (<nmod>obj).val
        return 1
    z = any_as_fmpz(obj)
    if z is not NotImplemented:
        val[0] = fmpz_fdiv_ui((<fmpz>z).val, mod.n)
        return 1
    q = any_as_fmpq(obj)
    if q is not NotImplemented:
        fmpz_init(t)
        fmpz_set_ui(t, mod.n)
        success = fmpq_mod_fmpz(t, (<fmpq>q).val, t)
        val[0] = fmpz_get_ui(t)
        fmpz_clear(t)
        if not success:
            raise ZeroDivisionError("%s does not exist mod %i!" % (q, mod.n))
        return 1
    return 0

cdef class nmod(flint_scalar):
    """
    The nmod type represents elements of Z/nZ for word-size n.

        >>> nmod(10,17) * 2
        3

    """

    cdef mp_limb_t val
    cdef nmod_t mod

    def __init__(self, val, mod):
        cdef mp_limb_t m
        m = mod
        nmod_init(&self.mod, m)
        if not any_as_nmod(&self.val, val, self.mod):
            raise TypeError("cannot create nmod from object of type %s" % type(val))

    def repr(self):
        return "nmod(%s, %s)" % (self.val, self.mod.n)

    def str(self):
        return str(int(self.val))

    def __int__(self):
        return int(self.val)

    def __long__(self):
        return self.val

    def modulus(self):
        return self.mod.n

    def __richcmp__(s, t, int op):
        cdef mp_limb_t v
        cdef bint res
        if op != 2 and op != 3:
            raise TypeError("nmods cannot be ordered")
        if typecheck(s, nmod) and typecheck(t, nmod):
            res = ((<nmod>s).val == (<nmod>t).val) and \
                  ((<nmod>s).mod.n == (<nmod>t).mod.n)
            if op == 2:
                return res
            else:
                return not res
        return NotImplemented

    def __nonzero__(self):
        return self.val != 0

    def __pos__(self):
        return self

    def __neg__(self):
        cdef nmod r = nmod.__new__(nmod)
        r.mod = self.mod
        r.val = nmod_neg(self.val, self.mod)
        return r

    def __add__(s, t):
        cdef nmod r
        cdef mp_limb_t val
        if not typecheck(s, nmod):
            s, t = t, s
        if any_as_nmod(&val, t, (<nmod>s).mod):
            r = nmod.__new__(nmod)
            r.mod = (<nmod>s).mod
            r.val = nmod_add(val, (<nmod>s).val, r.mod)
            return r
        return NotImplemented

    def __sub__(s, t):
        cdef nmod r
        cdef mp_limb_t val
        if typecheck(s, nmod):
            if any_as_nmod(&val, t, (<nmod>s).mod):
                r = nmod.__new__(nmod)
                r.mod = (<nmod>s).mod
                r.val = nmod_sub((<nmod>s).val, val, r.mod)
                return r
        else:
            if any_as_nmod(&val, s, (<nmod>t).mod):
                r = nmod.__new__(nmod)
                r.mod = (<nmod>t).mod
                r.val = nmod_sub(val, (<nmod>t).val, r.mod)
                return r
        return NotImplemented

    def __mul__(s, t):
        cdef nmod r
        cdef mp_limb_t val
        if not typecheck(s, nmod):
            s, t = t, s
        if any_as_nmod(&val, t, (<nmod>s).mod):
            r = nmod.__new__(nmod)
            r.mod = (<nmod>s).mod
            r.val = nmod_mul(val, (<nmod>s).val, r.mod)
            return r
        return NotImplemented

    @staticmethod
    def _div_(s, t):
        cdef nmod r
        cdef mp_limb_t sval, tval, x
        cdef nmod_t mod
        if typecheck(s, nmod):
            mod = (<nmod>s).mod
            sval = (<nmod>s).val
            if not any_as_nmod(&tval, t, mod):
                return NotImplemented
        else:
            mod = (<nmod>t).mod
            tval = (<nmod>t).val
            if not any_as_nmod(&sval, s, mod):
                return NotImplemented
        if tval == 0:
            raise ZeroDivisionError("%s is not invertible mod %s" % (tval, mod.n))
        # XXX: check invertibility?
        x = nmod_div(sval, tval, mod)
        if x == 0:
            raise ZeroDivisionError("%s is not invertible mod %s" % (tval, mod.n))
        r = nmod.__new__(nmod)
        r.mod = mod
        r.val = x
        return r

    def __truediv__(s, t):
        return nmod._div_(s, t)

    def __div__(s, t):
        return nmod._div_(s, t)

    def __invert__(self):
        return (1 / self)   # XXX: speed up

