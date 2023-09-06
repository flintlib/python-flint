from flint._flint cimport mp_limb_t, mp_bitcnt_t

cdef extern from "flint/nmod.h":
    ctypedef struct nmod_t:
       mp_limb_t n
       mp_limb_t ninv
       mp_bitcnt_t norm
    void nmod_init(nmod_t * mod, mp_limb_t n)
    mp_limb_t nmod_add(mp_limb_t a, mp_limb_t b, nmod_t mod)
    mp_limb_t nmod_sub(mp_limb_t a, mp_limb_t b, nmod_t mod)
    mp_limb_t nmod_neg(mp_limb_t a, nmod_t mod)
    mp_limb_t nmod_mul(mp_limb_t a, mp_limb_t b, nmod_t mod)
    mp_limb_t nmod_div(mp_limb_t a, mp_limb_t b, nmod_t mod)
