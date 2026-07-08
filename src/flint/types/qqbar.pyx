# src/flint/types/qqbar.pyx
from flint.flint_base.flint_base cimport flint_scalar
cdef class qqbar(flint_scalar):
    """
    Python wrapper class for FLINT's qqbar_t (Algebraic Numbers Fields)
    """
    def __cinit__(self):
        # Always safe-allocate raw C pointers inside __cinit__
        qqbar_init(self.val)

    def __dealloc__(self):
        # Prevent memory leaks when Python cleans up the object instance
        qqbar_clear(self.val)

    def __init__(self, val=None):
        if val is not None:
            if isinstance(val, qqbar):
                qqbar_set(self.val, (<qqbar>val).val)
            else:
                raise TypeError("Cannot initialize qqbar type from input")

    cpdef int is_rational(self):
        # Example method shell pointing to underlying properties
        return 0

    def __repr__(self):
        return "Algebraic Number (qqbar)"