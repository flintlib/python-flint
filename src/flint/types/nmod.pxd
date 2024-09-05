cimport cython

from flint.flintlib.types.flint cimport mp_limb_t, ulong
from flint.flintlib.types.fmpz cimport fmpz_t
from flint.flintlib.types.nmod cimport nmod_t

from flint.flintlib.functions.nmod cimport nmod_init
from flint.flintlib.functions.ulong_extras cimport n_is_prime
from flint.flintlib.functions.fmpq cimport fmpq_mod_fmpz
from flint.flintlib.functions.fmpz cimport (
    fmpz_t,
    fmpz_fdiv_ui,
    fmpz_init,
    fmpz_clear,
    fmpz_set_ui,
    fmpz_get_ui,
)

from flint.flint_base.flint_base cimport flint_scalar
from flint.utils.typecheck cimport typecheck

from flint.types.fmpz cimport fmpz, any_as_fmpz
from flint.types.fmpq cimport fmpq, any_as_fmpq


cdef dict _nmod_ctx_cache


@cython.no_gc
cdef class nmod_ctx:
    cdef nmod_t mod
    cdef bint _is_prime

    @staticmethod
    cdef inline nmod_ctx any_as_nmod_ctx(obj):
        """Convert an ``nmod_ctx`` or ``int`` to an ``nmod_ctx``."""
        if typecheck(obj, nmod_ctx):
            return obj
        if typecheck(obj, int):
            return nmod_ctx._get_ctx(obj)
        raise TypeError("Invalid context/modulus for nmod: %s" % obj)

    @staticmethod
    cdef inline nmod_ctx _get_ctx(int mod):
        """Retrieve an nmod context from the cache or create a new one."""
        ctx = _nmod_ctx_cache.get(mod)
        if ctx is None:
            ctx = _nmod_ctx_cache.setdefault(mod, nmod_ctx._new_ctx(mod))
        return ctx

    @staticmethod
    cdef inline nmod_ctx _new_ctx(ulong mod):
        """Create a new nmod context."""
        cdef nmod_ctx ctx = nmod_ctx.__new__(nmod_ctx)
        nmod_init(&ctx.mod, mod)
        ctx._is_prime = n_is_prime(mod)
        return ctx

    @cython.final
    cdef inline int any_as_nmod(nmod_ctx ctx, mp_limb_t * val, obj) except -1:
        """Convert an object to an nmod element."""
        cdef int success
        cdef fmpz_t t
        if typecheck(obj, nmod):
            if (<nmod>obj).ctx.mod.n != ctx.mod.n:
                raise ValueError("cannot coerce integers mod n with different n")
            val[0] = (<nmod>obj).val
            return 1
        z = any_as_fmpz(obj)
        if z is not NotImplemented:
            val[0] = fmpz_fdiv_ui((<fmpz>z).val, ctx.mod.n)
            return 1
        q = any_as_fmpq(obj)
        if q is not NotImplemented:
            fmpz_init(t)
            fmpz_set_ui(t, ctx.mod.n)
            success = fmpq_mod_fmpz(t, (<fmpq>q).val, t)
            val[0] = fmpz_get_ui(t)
            fmpz_clear(t)
            if not success:
                raise ZeroDivisionError("%s does not exist mod %i!" % (q, ctx.mod.n))
            return 1
        return 0

    @cython.final
    cdef inline nmod new_nmod(self):
        cdef nmod r = nmod.__new__(nmod)
        r.ctx = self
        return r


@cython.no_gc
cdef class nmod(flint_scalar):
    cdef mp_limb_t val
    cdef nmod_ctx ctx
