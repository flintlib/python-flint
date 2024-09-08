from flint.flintlib.functions.fq_default cimport *
from flint.types.fmpz cimport fmpz
from flint.flint_base.flint_base cimport flint_scalar

cpdef enum fq_default_type:
    """
    Enum for the fq_default context types:

    - 1. `fq_default_ctx.FQ_ZECH`: Use `fq_zech_t`,
    - 2. `fq_default_ctx.FQ_NMOD`: Use `fq_nmod_t`,
    - 3. `fq_default_ctx.FQ`: Use `fq_t`.
    - 4. `fq_default_ctx.NMOD`: Use `nmod` for degree = 1,
    - 5. `fq_default_ctx.FMPZ_MOD`: Use `fmpz_mod` for degree = 1.

    These can be manually selected, or type: `fq_default_ctx.DEFAULT` can be used
    for the implementation to be automatically decided by Flint (default),
    """
    DEFAULT = 0
    FQ_ZECH = 1
    FQ_NMOD = 2
    FQ = 3
    NMOD = 4
    FMPZ_MOD = 5

cdef class fq_default_ctx:
    cdef fq_default_ctx_t val
    cdef readonly char *var
    cdef bint _initialized

    cdef new_ctype_fq_default(self)
    cdef set_list_as_fq_default(self, fq_default_t val, obj)
    cdef set_any_scalar_as_fq_default(self, fq_default_t fq_ele, obj)
    cdef set_any_as_fq_default(self, fq_default_t val, obj)
    cdef any_as_fq_default(self, obj)

    cdef _set_from_order(self, p, d, var, fq_type=*, check_prime=*)
    cdef _set_from_modulus(self, modulus, var, fq_type=*, check_prime=*, check_modulus=*)

cdef class fq_default(flint_scalar):
    cdef fq_default_ctx ctx
    cdef fq_default_t val

    # Arithmetic for flint_scalar base class
    cpdef fq_default _neg_(fq_default self)
    cpdef fq_default _add_(fq_default self, fq_default other)
    cpdef fq_default _mul_(fq_default self, fq_default other)
    cpdef fq_default _sub_(fq_default self, fq_default other)
    cpdef fq_default _rsub_(fq_default self, fq_default other)
    cpdef fq_default _div_(fq_default self, fq_default other)
    cpdef fq_default _rdiv_(fq_default self, fq_default other)
    cpdef fq_default _invert_(fq_default self)
