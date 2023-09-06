from flint._flint cimport mp_limb_t
from flint.flintlib.nmod_poly cimport nmod_poly_struct, nmod_poly_t

cdef extern from "flint/nmod_poly_factor.h":
    ctypedef struct nmod_poly_factor_struct:
        nmod_poly_struct *p
        long *exp
        long num
        long alloc
    ctypedef nmod_poly_factor_struct nmod_poly_factor_t[1]

    int nmod_poly_is_irreducible(nmod_poly_t f)
    mp_limb_t nmod_poly_factor_with_berlekamp(nmod_poly_factor_t result, nmod_poly_t poly)
    mp_limb_t nmod_poly_factor_with_cantor_zassenhaus(nmod_poly_factor_t result, nmod_poly_t poly)
    mp_limb_t nmod_poly_factor(nmod_poly_factor_t result, nmod_poly_t input)
    void nmod_poly_factor_init(nmod_poly_factor_t fac)
    void nmod_poly_factor_clear(nmod_poly_factor_t fac)
