# src/flint/types/qqbar.pxd
from flint.flint_base.flint_base cimport flint_ctx, flint_scalar
from flint.flintlib.types.flint cimport slong
from flint.flintlib.functions.fmpz cimport fmpz_t
from flint.flintlib.functions.fmpq cimport fmpq_t


cdef extern from "flint/qqbar.h" nogil:
    ctypedef struct qqbar_struct:
        pass
    ctypedef qqbar_struct qqbar_t[1]

    # Lifecycle methods
    void qqbar_init(qqbar_t res)
    void qqbar_clear(qqbar_t res)
    void qqbar_set(qqbar_t res, const qqbar_t x)
    void qqbar_set_si(qqbar_t res, slong x)
    
    # Printing methods
    void qqbar_print(const qqbar_t x)
    void qqbar_printn(const qqbar_t x, slong n)

    # Interoperability setters
    void qqbar_set_fmpz(qqbar_t res, const fmpz_t x)
    void qqbar_set_fmpq(qqbar_t res, const fmpq_t x)

    # Arithmetic operations
    void qqbar_add(qqbar_t res, const qqbar_t x, const qqbar_t y)
    void qqbar_sub(qqbar_t res, const qqbar_t x, const qqbar_t y)
    void qqbar_mul(qqbar_t res, const qqbar_t x, const qqbar_t y)
    void qqbar_div(qqbar_t res, const qqbar_t x, const qqbar_t y)

cdef class qqbar(flint_scalar):
    cdef qqbar_t val