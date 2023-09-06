from flint.flintlib.acb cimport acb_t, acb_ptr
from flint.flintlib.dirichlet cimport dirichlet_group_t, dirichlet_char_t
from flint.flintlib.flint cimport ulong
from flint.flintlib.acb_poly cimport acb_poly_t
from flint.flintlib.fmpz cimport fmpz_t
from flint.flintlib.arb cimport arb_t


cdef extern from "acb_dirichlet.h":
    void acb_dirichlet_eta(acb_t res, const acb_t s, long prec)
    void acb_dirichlet_chi(acb_t res, const dirichlet_group_t G, const dirichlet_char_t chi, ulong n, long prec)

    void acb_dirichlet_l(acb_t res, const acb_t s, const dirichlet_group_t G, const dirichlet_char_t chi, long prec)
    void acb_dirichlet_hardy_z(acb_ptr res, const acb_t t, const dirichlet_group_t G, const dirichlet_char_t chi, long len, long prec)
    void acb_dirichlet_l_series(acb_poly_t res, const acb_poly_t s, const dirichlet_group_t G, const dirichlet_char_t chi, int deflate, long len, long prec)

    void acb_dirichlet_stieltjes(acb_t res, const fmpz_t n, const acb_t a, long prec)

    void acb_dirichlet_gram_point(arb_t res, const fmpz_t n, const dirichlet_group_t G, const dirichlet_char_t chi, long prec)
    void acb_dirichlet_zeta_zeros(acb_ptr res, const fmpz_t n, long len, long prec)
    void acb_dirichlet_zeta_nzeros(arb_t res, const arb_t t, long prec)
    void acb_dirichlet_backlund_s(arb_t res, const arb_t t, long prec)
    void acb_dirichlet_zeta_zero(acb_t res, const fmpz_t n, long prec)
    void acb_dirichlet_zeta_zeros(acb_ptr res, const fmpz_t n, long len, long prec)
