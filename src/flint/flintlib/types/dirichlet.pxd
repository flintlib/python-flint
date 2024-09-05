from flint.flintlib.types.flint cimport ulong
from flint.flintlib.functions.nmod cimport nmod_t


cdef extern from "flint/dirichlet.h":

    ctypedef struct dirichlet_group_struct:
        ulong q
        ulong q_even
        nmod_t mod
        ulong rad_q
        ulong phi_q
        long neven
        long num
        ulong expo
        void * P
        ulong * generators
        ulong * PHI
    ctypedef dirichlet_group_struct dirichlet_group_t[1]

    ctypedef struct dirichlet_char_struct:
        ulong n
        ulong * log
    ctypedef dirichlet_char_struct dirichlet_char_t[1]

    cdef ulong DIRICHLET_CHI_NULL
