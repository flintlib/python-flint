cdef extern from "flint/gr.h":

    cdef enum truth_t_enum:
        T_TRUE
        T_FALSE
        T_UNKOWN

    ctypedef truth_t_enum truth_t

    ctypedef struct gr_stream_struct:
        FLINT_FILE * fp
        char * s
        slong len
        slong alloc

    ctypedef gr_stream_struct gr_stream_t[1]

    ctypedef int (*gr_funcptr)(void)

    cdef const int GR_CTX_STRUCT_DATA_BYTES

    ctypedef int (*gr_funcptr)(void)

    cdef struct gr_ctx_struct:
        char data[GR_CTX_STRUCT_DATA_BYTES]
        ulong which_ring
        slong sizeof_elem
        gr_funcptr * methods
        ulong size_limit

    ctypedef gr_ctx_struct gr_ctx_t[1]

    ctypedef void * gr_ptr
    ctypedef const void * gr_srcptr
    ctypedef void * gr_ctx_ptr

    ctypedef struct gr_vec_struct:
        gr_ptr entries
        slong alloc
        slong length

    ctypedef gr_vec_struct gr_vec_t[1]

    ctypedef struct gr_mat_struct:
        gr_ptr entries
        slong r
        slong c
        gr_ptr * rows

    ctypedef gr_mat_struct gr_mat_t[1]

    ctypedef struct gr_poly_struct:
        gr_ptr coeffs
        slong alloc
        slong length

    ctypedef gr_poly_struct gr_poly_t[1]

    ctypedef struct gr_mpoly_struct:
        gr_ptr coeffs
        ulong * exps
        slong length
        flint_bitcnt_t bits  # number of bits per exponent
        slong coeffs_alloc   # abs size in ulong units
        slong exps_alloc     # abs size in ulong units

    ctypedef gr_mpoly_struct gr_mpoly_t[1]
