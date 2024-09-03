cimport cython

from cpython.list cimport PyList_GET_SIZE

from flint.utils.typecheck cimport typecheck
from flint.flint_base.flint_base cimport flint_poly

from flint.flintlib.flint cimport mp_limb_t, ulong
from flint.flintlib.fmpz_poly cimport fmpz_poly_get_nmod_poly

from flint.flintlib.nmod cimport nmod_t
from flint.flintlib.nmod_poly cimport (
    nmod_poly_t,
    nmod_poly_init_preinv,
    nmod_poly_fit_length,
    nmod_poly_set_coeff_ui,
)

from flint.types.fmpz_poly cimport fmpz_poly, any_as_fmpz_poly
from flint.types.nmod cimport nmod_ctx, nmod


cdef dict _nmod_poly_ctx_cache = {}


@cython.no_gc
cdef class nmod_poly_ctx:
    cdef nmod_t mod
    cdef bint _is_prime
    cdef nmod_ctx scalar_ctx

    @staticmethod
    cdef inline nmod_poly_ctx any_as_nmod_poly_ctx(obj):
        """Convert an ``nmod_poly_ctx`` or ``int`` to an ``nmod_poly_ctx``."""
        if typecheck(obj, nmod_poly_ctx):
            return obj
        if typecheck(obj, int):
            return nmod_poly_ctx._get_ctx(obj)
        raise TypeError("Invalid context/modulus for nmod_poly: %s" % obj)

    @staticmethod
    cdef inline nmod_poly_ctx _get_ctx(int mod):
        """Retrieve an nmod_poly context from the cache or create a new one."""
        ctx = _nmod_poly_ctx_cache.get(mod)
        if ctx is None:
            ctx = _nmod_poly_ctx_cache.setdefault(mod, nmod_poly_ctx._new_ctx(mod))
        return ctx

    @staticmethod
    cdef inline nmod_poly_ctx _new_ctx(ulong mod):
        """Create a new nmod_poly context."""
        cdef nmod_ctx scalar_ctx
        cdef nmod_poly_ctx ctx
        scalar_ctx = nmod_ctx.new(mod)

        ctx = nmod_poly_ctx.__new__(nmod_poly_ctx)
        ctx.mod = scalar_ctx.mod
        ctx._is_prime = scalar_ctx._is_prime
        ctx.scalar_ctx = scalar_ctx

        return ctx

    @cython.final
    cdef inline int any_as_nmod(self, mp_limb_t * val, obj) except -1:
        return self.scalar_ctx.any_as_nmod(val, obj)

    @cython.final
    cdef inline any_as_nmod_poly(self, obj):
        cdef nmod_poly r
        cdef mp_limb_t v
        # XXX: should check that modulus is the same here, and not all over the place
        if typecheck(obj, nmod_poly):
            return obj
        if self.any_as_nmod(&v, obj):
            r = self.new_nmod_poly()
            nmod_poly_set_coeff_ui(r.val, 0, v)
            return r
        x = any_as_fmpz_poly(obj)
        if x is not NotImplemented:
            r = self.new_nmod_poly()
            fmpz_poly_get_nmod_poly(r.val, (<fmpz_poly>x).val)
            return r
        return NotImplemented

    @cython.final
    cdef inline nmod new_nmod(self):
        return self.scalar_ctx.new_nmod()

    @cython.final
    cdef inline nmod_poly new_nmod_poly(self):
        cdef nmod_poly p = nmod_poly.__new__(nmod_poly)
        nmod_poly_init_preinv(p.val, self.mod.n, self.mod.ninv)
        p.ctx = self
        return p

    @cython.final
    cdef inline nmod_poly_set_list(self, nmod_poly_t poly, list val):
        cdef long i, n
        cdef mp_limb_t v
        n = PyList_GET_SIZE(val)
        nmod_poly_fit_length(poly, n)
        for i from 0 <= i < n:
            if self.any_as_nmod(&v, val[i]):
                nmod_poly_set_coeff_ui(poly, i, v)
            else:
                raise TypeError("unsupported coefficient in list")


@cython.no_gc
cdef class nmod_poly(flint_poly):
    cdef nmod_poly_t val
    cdef nmod_poly_ctx ctx

    cpdef long length(self)
    cpdef long degree(self)
    cpdef mp_limb_t modulus(self)
