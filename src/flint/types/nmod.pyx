from flint.flint_base.flint_base cimport flint_scalar
from flint.utils.typecheck cimport typecheck
from flint.types.fmpq cimport any_as_fmpq
from flint.types.fmpz cimport any_as_fmpz
from flint.types.fmpz cimport fmpz
from flint.types.fmpq cimport fmpq

from flint.flintlib.types.flint cimport ulong, nmod_t, mp_limb_t
from flint.flintlib.functions.fmpz cimport fmpz_t
from flint.flintlib.functions.nmod cimport (
    nmod_init,
    nmod_pow_fmpz,
    nmod_neg,
    nmod_add,
    nmod_sub,
    nmod_mul,
)
from flint.flintlib.functions.fmpz cimport fmpz_fdiv_ui, fmpz_init, fmpz_clear
from flint.flintlib.functions.fmpz cimport fmpz_set_ui, fmpz_get_ui
from flint.flintlib.functions.fmpq cimport fmpq_mod_fmpz
from flint.flintlib.functions.ulong_extras cimport n_gcdinv, n_sqrtmod

from flint.utils.flint_exceptions import DomainError

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

    def modulus(self):
        return self.mod.n

    def __richcmp__(s, t, int op):
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
        elif typecheck(s, nmod) and typecheck(t, int):
            res = s.val == (t % s.mod.n)
            if op == 2:
                return res
            else:
                return not res
        return NotImplemented

    def __hash__(self):
        return hash((int(self.val), self.modulus))

    def __bool__(self):
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
        if any_as_nmod(&val, t, (<nmod>s).mod):
            r = nmod.__new__(nmod)
            r.mod = (<nmod>s).mod
            r.val = nmod_add(val, (<nmod>s).val, r.mod)
            return r
        return NotImplemented

    def __radd__(s, t):
        cdef nmod r
        cdef mp_limb_t val
        if any_as_nmod(&val, t, (<nmod>s).mod):
            r = nmod.__new__(nmod)
            r.mod = (<nmod>s).mod
            r.val = nmod_add((<nmod>s).val, val, r.mod)
            return r
        return NotImplemented

    def __sub__(s, t):
        cdef nmod r
        cdef mp_limb_t val
        if any_as_nmod(&val, t, (<nmod>s).mod):
            r = nmod.__new__(nmod)
            r.mod = (<nmod>s).mod
            r.val = nmod_sub((<nmod>s).val, val, r.mod)
            return r
        return NotImplemented

    def __rsub__(s, t):
        cdef nmod r
        cdef mp_limb_t val
        if any_as_nmod(&val, t, (<nmod>s).mod):
            r = nmod.__new__(nmod)
            r.mod = (<nmod>s).mod
            r.val = nmod_sub(val, (<nmod>s).val, r.mod)
            return r
        return NotImplemented

    def __mul__(s, t):
        cdef nmod r
        cdef mp_limb_t val
        if any_as_nmod(&val, t, (<nmod>s).mod):
            r = nmod.__new__(nmod)
            r.mod = (<nmod>s).mod
            r.val = nmod_mul(val, (<nmod>s).val, r.mod)
            return r
        return NotImplemented

    def __rmul__(s, t):
        cdef nmod r
        cdef mp_limb_t val
        if any_as_nmod(&val, t, (<nmod>s).mod):
            r = nmod.__new__(nmod)
            r.mod = (<nmod>s).mod
            r.val = nmod_mul((<nmod>s).val, val, r.mod)
            return r
        return NotImplemented

    @staticmethod
    def _div_(s, t):
        cdef nmod r
        cdef mp_limb_t sval, tval
        cdef nmod_t mod
        cdef ulong tinvval

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
        if not s:
            return s

        g = n_gcdinv(&tinvval, <ulong>tval, <ulong>mod.n)
        if g != 1:
            raise ZeroDivisionError("%s is not invertible mod %s" % (tval, mod.n))

        r = nmod.__new__(nmod)
        r.mod = mod
        r.val = nmod_mul(sval, <mp_limb_t>tinvval, mod)
        return r

    def __truediv__(s, t):
        return nmod._div_(s, t)

    def __rtruediv__(s, t):
        return nmod._div_(t, s)

    def __invert__(self):
        cdef nmod r
        cdef ulong g, inv, sval
        sval = <ulong>(<nmod>self).val
        g = n_gcdinv(&inv, sval, self.mod.n)
        if g != 1:
            raise ZeroDivisionError("%s is not invertible mod %s" % (sval, self.mod.n))
        r = nmod.__new__(nmod)
        r.mod = self.mod
        r.val = <mp_limb_t>inv
        return r

    def __pow__(self, exp, modulus=None):
        cdef nmod r
        cdef mp_limb_t rval, mod
        cdef ulong g, rinv

        if modulus is not None:
            raise TypeError("three-argument pow() not supported by nmod")

        e = any_as_fmpz(exp)
        if e is NotImplemented:
            return NotImplemented

        rval = (<nmod>self).val
        mod = (<nmod>self).mod.n

        # XXX: It is not clear that it is necessary to special case negative
        # exponents here. The nmod_pow_fmpz function seems to handle this fine
        # but the Flint docs say that the exponent must be nonnegative.
        if e < 0:
            g = n_gcdinv(&rinv, <ulong>rval, <ulong>mod)
            if g != 1:
                raise ZeroDivisionError("%s is not invertible mod %s" % (rval, mod))
            rval = <mp_limb_t>rinv
            e = -e

        r = nmod.__new__(nmod)
        r.mod = self.mod
        r.val = nmod_pow_fmpz(rval, (<fmpz>e).val, self.mod)
        return r

    def sqrt(self):
        """
        Return the square root of this nmod or raise an exception.

            >>> s = nmod(10, 13).sqrt()
            >>> s
            6
            >>> s * s
            10
            >>> nmod(11, 13).sqrt()
            Traceback (most recent call last):
                ...
            flint.utils.flint_exceptions.DomainError: no square root exists for 11 mod 13

        The modulus must be prime.

        """
        cdef nmod r
        cdef mp_limb_t val
        r = nmod.__new__(nmod)
        r.mod = self.mod

        if self.val == 0:
            return r

        val = n_sqrtmod(self.val, self.mod.n)
        if val == 0:
            raise DomainError("no square root exists for %s mod %s" % (self.val, self.mod.n))

        r.val = val
        return r
