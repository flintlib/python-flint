from flint.flintlib.flint cimport ulong
from flint.flintlib.nmod_vec cimport nmod_t

cdef extern from "dirichlet.h":
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

    ulong dirichlet_group_size(const dirichlet_group_t G)
    void dirichlet_group_init(dirichlet_group_t G, ulong q)
    void dirichlet_group_clear(dirichlet_group_t G)
    ulong dirichlet_number_primitive(const dirichlet_group_t G)

    void dirichlet_char_init(dirichlet_char_t x, const dirichlet_group_t G)
    void dirichlet_char_clear(dirichlet_char_t x)
    void dirichlet_char_print(const dirichlet_group_t G, const dirichlet_char_t x)

    void dirichlet_char_set(dirichlet_char_t x, const dirichlet_group_t G, const dirichlet_char_t y)
    int dirichlet_char_eq(const dirichlet_char_t x, const dirichlet_char_t y)
    int dirichlet_parity_char(const dirichlet_group_t G, const dirichlet_char_t x)
    ulong dirichlet_conductor_char(const dirichlet_group_t G, const dirichlet_char_t x)
    ulong dirichlet_order_char(const dirichlet_group_t G, const dirichlet_char_t x)

    void dirichlet_char_log(dirichlet_char_t x, const dirichlet_group_t G, ulong m)
    ulong dirichlet_char_exp(const dirichlet_group_t G, const dirichlet_char_t x)
    ulong _dirichlet_char_exp(dirichlet_char_t x, const dirichlet_group_t G)
    void dirichlet_char_index(dirichlet_char_t x, const dirichlet_group_t G, ulong j)
    ulong dirichlet_index_char(const dirichlet_group_t G, const dirichlet_char_t x)
    void dirichlet_char_one(dirichlet_char_t x, const dirichlet_group_t G)
    void dirichlet_char_first_primitive(dirichlet_char_t x, const dirichlet_group_t G)
    int dirichlet_char_next(dirichlet_char_t x, const dirichlet_group_t G)
    int dirichlet_char_next_primitive(dirichlet_char_t x, const dirichlet_group_t G)
    void dirichlet_char_mul(dirichlet_char_t c, const dirichlet_group_t G, const dirichlet_char_t a, const dirichlet_char_t b)
    void dirichlet_char_pow(dirichlet_char_t c, const dirichlet_group_t G, const dirichlet_char_t a, ulong n)
    void dirichlet_char_lower(dirichlet_char_t y, const dirichlet_group_t H, const dirichlet_char_t x, const dirichlet_group_t G)
    void dirichlet_char_lift(dirichlet_char_t y, const dirichlet_group_t G, const dirichlet_char_t x, const dirichlet_group_t H)

    cdef ulong DIRICHLET_CHI_NULL

    ulong dirichlet_pairing(const dirichlet_group_t G, ulong m, ulong n)
    ulong dirichlet_pairing_char(const dirichlet_group_t G, const dirichlet_char_t a, const dirichlet_char_t b)

    int dirichlet_char_is_principal(const dirichlet_group_t G, const dirichlet_char_t chi)
    int dirichlet_char_is_real(const dirichlet_group_t G, const dirichlet_char_t chi)
    int dirichlet_char_is_primitive(const dirichlet_group_t G, const dirichlet_char_t chi)
    ulong dirichlet_chi(const dirichlet_group_t G, const dirichlet_char_t chi, ulong n)
