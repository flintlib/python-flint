from flint.flintlib.types.flint cimport fmpz_struct, slong, ulong, nn_ptr

cdef extern from "flint/arf.h":
    ctypedef enum arf_rnd_t:
        ARF_RND_DOWN
        ARF_RND_NEAR
        ARF_RND_FLOOR
        ARF_RND_CEIL
        ARF_RND_UP

cdef extern from "flint/arf_types.h":

    cdef const int ARF_NOPTR_LIMBS

    ctypedef struct mantissa_noptr_struct:
        ulong d[ARF_NOPTR_LIMBS]

    ctypedef struct mantissa_ptr_struct:
        slong alloc
        nn_ptr d

    ctypedef union mantissa_struct:
        mantissa_noptr_struct noptr
        mantissa_ptr_struct ptr

    ctypedef struct arf_struct:
        fmpz_struct exp
        slong size
        mantissa_struct d

    ctypedef arf_struct arf_t[1]
    ctypedef arf_struct * arf_ptr
    ctypedef const arf_struct * arf_srcptr
