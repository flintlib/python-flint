from flint.flintlib.types.arb cimport arb_t
from flint.flintlib.types.arf cimport arf_t
from flint.flintlib.types.flint cimport fmpz_t, slong, ulong



cdef extern from "flint/partitions.h":
    void partitions_rademacher_bound(arf_t b, const fmpz_t n, ulong N)
    void partitions_hrr_sum_arb(arb_t x, const fmpz_t n, slong N0, slong N, int use_doubles)
    void partitions_fmpz_fmpz(fmpz_t p, const fmpz_t n, int use_doubles)
    void partitions_fmpz_ui(fmpz_t p, ulong n)
    void partitions_fmpz_ui_using_doubles(fmpz_t p, ulong n)
    void partitions_leading_fmpz(arb_t res, const fmpz_t n, slong prec)
