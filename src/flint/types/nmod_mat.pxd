cimport cython

from flint.utils.typecheck cimport typecheck
from flint.flint_base.flint_base cimport flint_mat

from flint.flintlib.flint cimport mp_limb_t, ulong
from flint.flintlib.fmpz_mat cimport (
    fmpz_mat_nrows,
    fmpz_mat_ncols,
    fmpz_mat_get_nmod_mat,
)
from flint.flintlib.nmod cimport nmod_t
from flint.flintlib.nmod_mat cimport (
    nmod_mat_t,
    nmod_mat_init,
    nmod_mat_init_set,
)

from flint.types.fmpz cimport fmpz
from flint.types.fmpz_mat cimport fmpz_mat, any_as_fmpz_mat
from flint.types.nmod cimport nmod_ctx, nmod
from flint.types.nmod_poly cimport nmod_poly_ctx, nmod_poly


cdef dict _nmod_mat_ctx_cache


@cython.no_gc
cdef class nmod_mat_ctx:
    cdef nmod_t mod
    cdef bint _is_prime
    cdef nmod_ctx scalar_ctx
    cdef nmod_poly_ctx poly_ctx

    @staticmethod
    cdef inline nmod_mat_ctx any_as_nmod_mat_ctx(obj):
        """Convert an ``nmod_mat_ctx`` or ``int`` to an ``nmod_mat_ctx``."""
        if typecheck(obj, nmod_mat_ctx):
            return obj
        if typecheck(obj, int):
            return nmod_mat_ctx._get_ctx(obj)
        elif typecheck(obj, fmpz):
            return nmod_mat_ctx._get_ctx(int(obj))
        raise TypeError("nmod_mat: expected last argument to be an nmod_mat_ctx or an integer")

    @staticmethod
    cdef inline nmod_mat_ctx _get_ctx(int mod):
        """Retrieve an nmod_mat context from the cache or create a new one."""
        ctx = _nmod_mat_ctx_cache.get(mod)
        if ctx is None:
            ctx = _nmod_mat_ctx_cache.setdefault(mod, nmod_mat_ctx._new_ctx(mod))
        return ctx

    @staticmethod
    cdef inline nmod_mat_ctx _new_ctx(ulong mod):
        """Create a new nmod_mat context."""
        cdef nmod_ctx scalar_ctx
        cdef nmod_poly_ctx poly_ctx
        cdef nmod_mat_ctx ctx

        poly_ctx = nmod_poly_ctx.new(mod)
        scalar_ctx = poly_ctx.scalar_ctx

        ctx = nmod_mat_ctx.__new__(nmod_mat_ctx)
        ctx.mod = scalar_ctx.mod
        ctx._is_prime = scalar_ctx._is_prime
        ctx.scalar_ctx = scalar_ctx
        ctx.poly_ctx = poly_ctx

        return ctx

    @cython.final
    cdef inline int any_as_nmod(self, mp_limb_t * val, obj) except -1:
        return self.scalar_ctx.any_as_nmod(val, obj)

    @cython.final
    cdef inline any_as_nmod_mat(self, obj):
        """Convert obj to nmod_mat or return NotImplemented."""
        cdef nmod_mat r

        if typecheck(obj, nmod_mat):
            return obj

        x = any_as_fmpz_mat(obj)
        if x is not NotImplemented:
            r = self.new_nmod_mat(fmpz_mat_nrows((<fmpz_mat>x).val),
                                  fmpz_mat_ncols((<fmpz_mat>x).val))
            fmpz_mat_get_nmod_mat(r.val, (<fmpz_mat>x).val)
            return r

        return NotImplemented

    @cython.final
    cdef inline nmod new_nmod(self):
        return self.scalar_ctx.new_nmod()

    @cython.final
    cdef inline nmod_poly new_nmod_poly(self):
        return self.poly_ctx.new_nmod_poly()

    @cython.final
    cdef inline nmod_mat new_nmod_mat(self, ulong m, ulong n):
        """New initialized nmod_mat of size m x n with context ctx."""
        cdef nmod_mat r = nmod_mat.__new__(nmod_mat)
        nmod_mat_init(r.val, m, n, self.mod.n)
        r.ctx = self
        return r

    @cython.final
    cdef inline nmod_mat new_nmod_mat_copy(self, nmod_mat other):
        """New copy of nmod_mat other."""
        cdef nmod_mat r = nmod_mat.__new__(nmod_mat)
        nmod_mat_init_set(r.val, other.val)
        r.ctx = other.ctx
        return r


@cython.no_gc
cdef class nmod_mat(flint_mat):
    cdef nmod_mat_t val
    cdef nmod_mat_ctx ctx

    cpdef long nrows(self)
    cpdef long ncols(self)
    cpdef mp_limb_t modulus(self)
    cdef __mul_nmod(self, mp_limb_t c)
