# _flint.pxd
#
# Define fundamental types and constants

cdef extern from "Python.h":
    ctypedef void PyObject

cdef enum:
    FMPZ_UNKNOWN = 0
    FMPZ_REF = 1
    FMPZ_TMP = 2


#
# Note: ulong and slong are used throughout Flint/Arb. They are expected to be
# 32 bit unsigned and signed integer types on a 32 bit system and 64 bit on a
# 64 bit system. We denote them as unsigned long and long here which would be
# incorrect on 64 bit Windows but the definition here does not matter because
# their actual sizes will be determined by the values from gmp.h and
# flint/flint.h. Their size in bits (32 or 64) is recorded in the FLINT_BITS
# macro which is defined in flint/flint.h.
#

cdef extern from "gmp.h":
    ctypedef unsigned long ulong
    ctypedef unsigned long mp_limb_t
    ctypedef long mp_size_t
    ctypedef long mp_exp_t
    ctypedef mp_limb_t* mp_ptr
    ctypedef mp_limb_t* mp_srcptr
    ctypedef unsigned long mp_bitcnt_t

cdef extern from "flint/fmpz.h":
    ctypedef long slong
    ctypedef ulong flint_bitcnt_t
    ctypedef ulong * nn_ptr
    ctypedef const ulong * nn_srcptr

ctypedef slong fmpz_struct
ctypedef fmpz_struct fmpz_t[1]

cdef extern from "flint/fmpz.h":
    # Macros
    int COEFF_IS_MPZ(fmpz_struct x)

cdef extern from *:
    """
    #if __FLINT_RELEASE < 30200 /* Flint < 3.2.0 */

    /* Functions renamed in Flint 3.2.0 */
    #define flint_rand_init flint_randinit
    #define flint_rand_clear flint_randclear

    #endif
    """

cdef extern from "flint/flint.h":
    # These defines are needed to work around a Cython bug.
    # Otherwise sizeof(ulong) will give the wrong size on 64 bit Windows.
    # https://github.com/cython/cython/issues/6339
    """
    #define SIZEOF_ULONG sizeof(ulong)
    #define SIZEOF_SLONG sizeof(slong)
    """
    int SIZEOF_ULONG
    int SIZEOF_SLONG

    ctypedef struct __FLINT_FILE:
        pass
    ctypedef __FLINT_FILE FLINT_FILE

    const char * FLINT_VERSION
    const int __FLINT_RELEASE

    const int FLINT_BITS

    ctypedef void * flint_rand_t
    void flint_rand_init(flint_rand_t state)
    void flint_rand_clear(flint_rand_t state)

    void flint_set_num_threads(long)
    long flint_get_num_threads()

    void flint_cleanup()

    ctypedef struct nmod_t:
        mp_limb_t n
        mp_limb_t ninv
        flint_bitcnt_t norm

    ctypedef struct fmpq:
        fmpz_struct num
        fmpz_struct den

cdef extern from *:
    """
    /* FLINT_BITS is not known until C compile time. We need to check if long
     * or long long matches FLINT_BITS to know which CPython function to call.
     */
    #if FLINT_BITS == 32 && LONG_MAX == 2147483647
    #define pylong_as_slong PyLong_AsLongAndOverflow
    #elif FLINT_BITS == 64 && LLONG_MAX == 9223372036854775807
    #define pylong_as_slong PyLong_AsLongLongAndOverflow
    #else
    #error FLINT_BITS does not match width of long or long long.
    #endif
    """
    slong pylong_as_slong(PyObject *pylong, int *overflow)
