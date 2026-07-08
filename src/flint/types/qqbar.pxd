# src/flint/types/qqbar.pxd
# src/flint/types/qqbar.pxd
from flint.flint_base.flint_base cimport flint_ctx, flint_scalar

# Note: In python-flint's repo layout, flint_base is flat in the flint/ folder, 
# so we cimport directly from 'flint.flint_base' rather than subdirectories.

cdef extern from "flint/qqbar.h":
    ctypedef struct qqbar_struct:
        # Leaving this block completely blank or using a mock field
        int _dummy
    ctypedef qqbar_struct qqbar_t[1]

    # Core lifetime and initialization C signatures
    void qqbar_init(qqbar_t res)
    void qqbar_clear(qqbar_t res)
    void qqbar_set(qqbar_t res, const qqbar_t x)

cdef class qqbar(flint_scalar):
    cdef qqbar_t val
    cpdef int is_rational(self)
