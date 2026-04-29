from flint.flintlib.types.flint cimport slong, ulong
from flint.flintlib.types.fq_nmod cimport fq_nmod_ctx_struct


cdef extern from "flint/fq_zech.h":

    ctypedef struct fq_zech_struct:
        ulong value

    ctypedef fq_zech_struct fq_zech_t[1]

    ctypedef struct fq_zech_ctx_struct:
        ulong qm1              # q - 1
        ulong qm1o2            # (q - 1) / 2 or 1 when p == 2
        ulong qm1opm1          # (q - 1) / (p - 1)
        ulong p
        double ppre
        ulong prime_root       # primitive root for prime subfield
        ulong * zech_log_table
        ulong * prime_field_table
        ulong * eval_table

        fq_nmod_ctx_struct * fq_nmod_ctx
        int owns_fq_nmod_ctx
        int is_conway  # whether field was generated using Flint Conway tables (assures primitivity)

    ctypedef fq_zech_ctx_struct fq_zech_ctx_t[1]

    ctypedef struct fq_zech_mat_struct:
        fq_zech_struct * entries
        slong r
        slong s
        fq_zech_struct ** rows

    ctypedef fq_zech_mat_struct fq_zech_mat_t[1]

    ctypedef struct fq_zech_poly_struct:
        fq_zech_struct * coeffs
        slong alloc
        slong length

    ctypedef fq_zech_poly_struct fq_zech_poly_t[1]

    ctypedef struct fq_zech_poly_factor_struct:
        fq_zech_poly_struct * poly
        slong * exp
        slong num
        slong alloc

    ctypedef fq_zech_poly_factor_struct fq_zech_poly_factor_t[1]
