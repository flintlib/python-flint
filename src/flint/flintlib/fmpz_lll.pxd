from flint.flintlib.fmpz_mat cimport fmpz_mat_t

cdef extern from "flint/fmpz_lll.h":
    ctypedef struct fmpz_lll_struct:
        double delta
        double eta
        int rt
        int gt

    ctypedef fmpz_lll_struct fmpz_lll_t[1]

    void fmpz_lll_context_init(fmpz_lll_t fl, double delta, double eta, int rt, int gt)
    void fmpz_lll(fmpz_mat_t B, fmpz_mat_t U, const fmpz_lll_t fl)
