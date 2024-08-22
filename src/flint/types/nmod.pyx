from flint.flint_base.flint_base cimport flint_scalar
from flint.utils.typecheck cimport typecheck
from flint.types.fmpq cimport any_as_fmpq
from flint.types.fmpz cimport any_as_fmpz
from flint.types.fmpz cimport fmpz
from flint.types.fmpq cimport fmpq

from flint.flintlib.flint cimport ulong
from flint.flintlib.fmpz cimport fmpz_t
from flint.flintlib.nmod cimport nmod_pow_fmpz
from flint.flintlib.nmod_vec cimport *
from flint.flintlib.fmpz cimport fmpz_fdiv_ui, fmpz_init, fmpz_clear
from flint.flintlib.fmpz cimport fmpz_set_ui, fmpz_get_ui
from flint.flintlib.fmpq cimport fmpq_mod_fmpz
from flint.flintlib.ulong_extras cimport n_gcdinv, n_is_prime, n_sqrtmod

from flint.utils.flint_exceptions import DomainError


#cdef int any_as_nmod(mp_limb_t * val, obj, nmod_t mod) except -1:
#    return mod.ctx.any_as_nmod(val, obj)


_nmod_ctx_cache = {}


cdef nmod_ctx any_as_nmod_ctx(obj):
    """Convert an int to an nmod_ctx."""
    if typecheck(obj, nmod_ctx):
        return obj
    if typecheck(obj, int):
        ctx = _nmod_ctx_cache.get(obj)
        if ctx is None:
            ctx = nmod_ctx(obj)
            _nmod_ctx_cache[obj] = ctx
        return ctx
    return NotImplemented


cdef class nmod_ctx:
    """
    Context object for creating :class:`~.nmod` initalised 
    with modulus :math:`N`.

        >>> nmod_ctx(17)
        nmod_ctx(17)

    """
    def __init__(self, mod):
        cdef mp_limb_t m
        m = mod
        nmod_init(&self.mod, m)
        self._is_prime = n_is_prime(m)

    def __eq__(self, other):
        # XXX: If we could ensure uniqueness of nmod_ctx for given modulus then
        # we would need to implement __eq__ and __hash__ at all...
        #
        # It isn't possible to ensure uniqueness in __new__ like it is in
        # Python because we can't return an existing object from __new__. What
        # we could do though is make it so that __init__ raises an error and
        # use a static method .new() to create new objects.
        if self is other:
            return True
        if not typecheck(other, nmod_ctx):
            return NotImplemented
        return self.mod.n == (<nmod_ctx>other).mod.n

    def __repr__(self):
        return f"nmod_ctx({self.modulus()})"

    cdef int any_as_nmod(self, mp_limb_t * val, obj) except -1:
        """Convert an object to an nmod element."""
        cdef int success
        cdef fmpz_t t
        if typecheck(obj, nmod):
            if (<nmod>obj).ctx != self:
                raise ValueError("cannot coerce integers mod n with different n")
            val[0] = (<nmod>obj).val
            return 1
        z = any_as_fmpz(obj)
        if z is not NotImplemented:
            val[0] = fmpz_fdiv_ui((<fmpz>z).val, self.mod.n)
            return 1
        q = any_as_fmpq(obj)
        if q is not NotImplemented:
            fmpz_init(t)
            fmpz_set_ui(t, self.mod.n)
            success = fmpq_mod_fmpz(t, (<fmpq>q).val, t)
            val[0] = fmpz_get_ui(t)
            fmpz_clear(t)
            if not success:
                raise ZeroDivisionError("%s does not exist mod %i!" % (q, self.mod.n))
            return 1
        return 0

    def modulus(self):
        """Get the modulus of the context.

        >>> ctx = nmod_ctx(17)
        >>> ctx.modulus()
        17

        """
        return fmpz(self.mod.n)

    def is_prime(self):
        """Check if the modulus is prime.

        >>> ctx = nmod_ctx(17)
        >>> ctx.is_prime()
        True

        """
        return self._is_prime

    def zero(self):
        """Return the zero element of the context.

        >>> ctx = nmod_ctx(17)
        >>> ctx.zero()
        0

        """
        return self(0)

    def one(self):
        """Return the one element of the context.

        >>> ctx = nmod_ctx(17)
        >>> ctx.one()
        1

        """
        return self(1)

    def __hash__(self):
        return hash(self.mod)

    def __eq__(self, other):
        if typecheck(other, nmod_ctx):
            return self.mod.n == (<nmod_ctx>other).mod.n
        else:
            return NotImplemented

    def __str__(self):
        return f"Context for nmod with modulus: {self.modulus()}"

    def __repr__(self):
        return f"nmod_ctx({self.modulus()})"

    cdef nmod _new(self, mp_limb_t * val):
        cdef nmod r = nmod.__new__(nmod)
        r.val = val[0]
        r.ctx = self
        return r

    def __call__(self, val):
        """Create an nmod element from an integer.

        >>> ctx = nmod_ctx(17)
        >>> ctx(10)
        10

        """
        cdef mp_limb_t v
        v = val
        return self._new(&v)


cdef class nmod(flint_scalar):
    """
    The nmod type represents elements of Z/nZ for word-size n.

        >>> nmod(10,17) * 2
        3

    """
    def __init__(self, val, mod):
        ctx = any_as_nmod_ctx(mod)
        if ctx is NotImplemented:
            raise TypeError("Invalid context/modulus for nmod: %s" % mod)
        if not ctx.any_as_nmod(&self.val, val):
            raise TypeError("cannot create nmod from object of type %s" % type(val))
        self.ctx = ctx

    def repr(self):
        return "nmod(%s, %s)" % (self.val, self.ctx.mod.n)

    def str(self):
        return str(int(self.val))

    def __int__(self):
        return int(self.val)

    def modulus(self):
        return self.ctx.mod.n

    def __richcmp__(s, t, int op):
        cdef bint res
        if op != 2 and op != 3:
            raise TypeError("nmods cannot be ordered")
        if typecheck(s, nmod) and typecheck(t, nmod):
            res = ((<nmod>s).val == (<nmod>t).val) and \
                  ((<nmod>s).ctx.mod.n == (<nmod>t).ctx.mod.n)
            if op == 2:
                return res
            else:
                return not res
        elif typecheck(s, nmod) and typecheck(t, int):
            res = s.val == (t % s.ctx.mod.n)
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
        r.ctx = self.ctx
        r.val = nmod_neg(self.val, self.ctx.mod)
        return r

    def __add__(s, t):
        cdef nmod r, s2
        cdef mp_limb_t val
        s2 = s
        if s2.ctx.any_as_nmod(&val, t):
            r = nmod.__new__(nmod)
            r.ctx = s2.ctx
            r.val = nmod_add(val, s2.val, s2.ctx.mod)
            return r
        return NotImplemented

    def __radd__(s, t):
        cdef nmod r, s2
        cdef mp_limb_t val
        s2 = s
        if s2.ctx.any_as_nmod(&val, t):
            r = nmod.__new__(nmod)
            r.ctx = s2.ctx
            r.val = nmod_add(s2.val, val, s2.ctx.mod)
            return r
        return NotImplemented

    def __sub__(s, t):
        cdef nmod r, s2
        cdef mp_limb_t val
        s2 = s
        if s2.ctx.any_as_nmod(&val, t):
            r = nmod.__new__(nmod)
            r.ctx = s2.ctx
            r.val = nmod_sub(s2.val, val, s2.ctx.mod)
            return r
        return NotImplemented

    def __rsub__(s, t):
        cdef nmod r
        cdef mp_limb_t val
        s2 = s
        if s2.ctx.any_as_nmod(&val, t):
            r = nmod.__new__(nmod)
            r.ctx = s2.ctx
            r.val = nmod_sub(val, s2.val, s2.ctx.mod)
            return r
        return NotImplemented

    def __mul__(s, t):
        cdef nmod r, s2
        cdef mp_limb_t val
        s2 = s
        if s2.ctx.any_as_nmod(&val, t):
            r = nmod.__new__(nmod)
            r.ctx = s2.ctx
            r.val = nmod_mul(val, s2.val, s2.ctx.mod)
            return r
        return NotImplemented

    def __rmul__(s, t):
        cdef nmod r, s2
        cdef mp_limb_t val
        s2 = s
        if s2.ctx.any_as_nmod(&val, t):
            r = nmod.__new__(nmod)
            r.ctx = s2.ctx
            r.val = nmod_mul(s2.val, val, s2.ctx.mod)
            return r
        return NotImplemented

    @staticmethod
    def _div_(s, t):
        cdef nmod r, s2, t2
        cdef mp_limb_t sval, tval
        cdef nmod_ctx ctx
        cdef ulong tinvval

        if typecheck(s, nmod):
            s2 = s
            ctx = s2.ctx
            sval = s2.val
            if not ctx.any_as_nmod(&tval, t):
                return NotImplemented
        else:
            t2 = t
            ctx = t2.ctx
            tval = t2.val
            if not ctx.any_as_nmod(&sval, s):
                return NotImplemented

        if tval == 0:
            raise ZeroDivisionError("%s is not invertible mod %s" % (tval, ctx.mod.n))
        if not s:
            return s

        g = n_gcdinv(&tinvval, <ulong>tval, <ulong>ctx.mod.n)
        if g != 1:
            raise ZeroDivisionError("%s is not invertible mod %s" % (tval, ctx.mod.n))

        r = nmod.__new__(nmod)
        r.ctx = ctx
        r.val = nmod_mul(sval, <mp_limb_t>tinvval, ctx.mod)
        return r

    def __truediv__(s, t):
        return nmod._div_(s, t)

    def __rtruediv__(s, t):
        return nmod._div_(t, s)

    def __invert__(self):
        cdef nmod r, s
        cdef nmod_ctx ctx
        cdef ulong g, inv, sval
        s = self
        ctx = s.ctx
        sval = <ulong>s.val
        g = n_gcdinv(&inv, sval, ctx.mod.n)
        if g != 1:
            raise ZeroDivisionError("%s is not invertible mod %s" % (sval, ctx.mod.n))
        r = nmod.__new__(nmod)
        r.ctx = ctx
        r.val = <mp_limb_t>inv
        return r

    def __pow__(self, exp, modulus=None):
        cdef nmod r, s
        cdef mp_limb_t rval, mod
        cdef ulong g, rinv

        if modulus is not None:
            raise TypeError("three-argument pow() not supported by nmod")

        e = any_as_fmpz(exp)
        if e is NotImplemented:
            return NotImplemented

        s = self
        ctx = s.ctx
        rval = s.val
        mod = ctx.mod.n

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
        r.ctx = ctx
        r.val = nmod_pow_fmpz(rval, (<fmpz>e).val, ctx.mod)
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
        r.ctx = self.ctx

        if self.val == 0:
            return r

        val = n_sqrtmod(self.val, self.ctx.mod.n)
        if val == 0:
            raise DomainError("no square root exists for %s mod %s" % (self.val, self.ctx.mod.n))

        r.val = val
        return r
