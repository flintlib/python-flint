from flint.flintlib.flint cimport mp_limb_t, mp_bitcnt_t
from flint.flintlib.nmod cimport nmod_t

cdef extern from "flint/nmod_vec.h":
    void nmod_init(nmod_t * mod, mp_limb_t n)
    mp_limb_t nmod_add(mp_limb_t a, mp_limb_t b, nmod_t mod)
    mp_limb_t nmod_sub(mp_limb_t a, mp_limb_t b, nmod_t mod)
    mp_limb_t nmod_neg(mp_limb_t a, nmod_t mod)
    mp_limb_t nmod_mul(mp_limb_t a, mp_limb_t b, nmod_t mod)
    mp_limb_t nmod_div(mp_limb_t a, mp_limb_t b, nmod_t mod)
