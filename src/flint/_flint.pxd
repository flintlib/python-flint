# _flint.pxd
#
# Define the contents of the Python, GMP, Flint and Arb headers.

cdef extern from "Python.h":
    ctypedef void PyObject
#    ctypedef void PyTypeObject
#     ctypedef long Py_ssize_t
#     int PyObject_TypeCheck(object, PyTypeObject*)
#     int PyInt_Check(PyObject *o)
#     int PyLong_Check(PyObject *o)
#     long PyInt_AS_LONG(PyObject *io)
#     double PyFloat_AS_DOUBLE(PyObject *io)
#     Py_ssize_t PyList_GET_SIZE(PyObject *list)
#     long PyLong_AsLongAndOverflow(PyObject *pylong, int *overflow)
#     long long PyLong_AsLongLongAndOverflow(PyObject *pylong, int *overflow)

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

ctypedef slong fmpz_struct

cdef extern from "flint/flint.h":
    const int FLINT_BITS
    ctypedef void * flint_rand_t
    void flint_randinit(flint_rand_t state)
    void flint_randclear(flint_rand_t state)
    void flint_set_num_threads(long)
    long flint_get_num_threads()
    void flint_cleanup()

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

from flint.flintlib.nmod_poly cimport nmod_poly_t
from flint.flintlib.nmod_mat cimport nmod_mat_t
from flint.flintlib.fmpz cimport fmpz_t
from flint.flintlib.fmpz_poly cimport fmpz_poly_t, fmpz_poly_struct
from flint.flintlib.fmpz_mat cimport fmpz_mat_struct, fmpz_mat_t
from flint.flintlib.fmpq cimport fmpq_t, fmpq_struct
from flint.flintlib.fmpq_poly cimport fmpq_poly_struct, fmpq_poly_t
from flint.flintlib.fmpq_mat cimport fmpq_mat_t
from flint.flintlib.mag cimport mag_struct, mag_t, mag_ptr, mag_srcptr

from flint.flintlib.arf cimport arf_struct, arf_t, arf_ptr, arf_srcptr, arf_rnd_t
from flint.flintlib.arb cimport arb_struct, arb_ptr, arb_srcptr, arb_t
from flint.flintlib.arb cimport arb_midref, arb_radref
from flint.flintlib.acb cimport acb_struct, acb_ptr, acb_srcptr, acb_t
from flint.flintlib.acb cimport acb_realref, acb_imagref
from flint.flintlib.arb_poly cimport arb_poly_struct, arb_poly_t
from flint.flintlib.arb_mat cimport arb_mat_struct, arb_mat_t
# from flint.flintlib.acb_poly cimport acb_poly_struct, acb_poly_t
# from flint.flintlib.acb_mat cimport acb_mat_struct, acb_mat_t
# from flint.flintlib.dirichlet cimport dirichlet_group_struct, dirichlet_group_t
# from flint.flintlib.dirichlet cimport dirichlet_char_struct, dirichlet_char_t
# from flint.flintlib.mpoly cimport ordering_t, mpoly_ctx_struct, mpoly_ctx_t



"""
cdef extern from "flint/fmpz_mpoly_factor.h":

    ctypedef struct fmpz_mpoly_factor_struct:
        fmpz_t content
        fmpz_mpoly_struct * poly
        fmpz_struct * exp
        slong length
        slong alloc

    ctypedef fmpz_mpoly_factor_struct fmpz_mpoly_factor_t[1]


    void fmpz_mpoly_factor_init(fmpz_mpoly_factor_t fac, const fmpz_mpoly_ctx_t ctx)

    void fmpz_mpoly_factor_clear(fmpz_mpoly_factor_t fac, const fmpz_mpoly_ctx_t ctx)

    int fmpz_mpoly_factor(fmpz_mpoly_factor_t fac, const fmpz_mpoly_t A, int full, const fmpz_mpoly_ctx_t ctx)
"""

