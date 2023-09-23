from flint.flintlib.flint cimport ulong, slong
from flint.flintlib.fmpz cimport fmpz_t
from flint.flintlib.arf cimport arf_t
from flint.flintlib.arb cimport arb_t

cdef extern from "flint/partitions.h":
# from here on is parsed
    void partitions_rademacher_bound(arf_t b, const fmpz_t n, ulong N)
    void partitions_hrr_sum_arb(arb_t x, const fmpz_t n, slong N0, slong N, int use_doubles)
    void partitions_fmpz_fmpz(fmpz_t p, const fmpz_t n, int use_doubles)
    void partitions_fmpz_ui(fmpz_t p, ulong n)
    void partitions_fmpz_ui_using_doubles(fmpz_t p, ulong n)
    void partitions_leading_fmpz(arb_t res, const fmpz_t n, slong prec)
