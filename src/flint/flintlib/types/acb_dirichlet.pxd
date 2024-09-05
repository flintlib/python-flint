from flint.flintlib.types.flint cimport ulong, slong
from flint.flintlib.types.arb cimport mag_struct
from flint.flintlib.types.acb cimport acb_t, acb_ptr, acb_struct

cdef extern from "flint/acb_dirichlet.h":

    ctypedef struct acb_dirichlet_roots_struct:
        ulong order
        ulong reduced_order
        acb_t z
        slong size
        slong depth
        acb_ptr * Z
        int use_pow

    ctypedef acb_dirichlet_roots_struct acb_dirichlet_roots_t[1]

    ctypedef struct acb_dirichlet_hurwitz_precomp_struct:
        acb_struct s
        mag_struct err
        acb_ptr coeffs
        int deflate
        slong A
        slong N
        slong K

    ctypedef acb_dirichlet_hurwitz_precomp_struct acb_dirichlet_hurwitz_precomp_t[1]
