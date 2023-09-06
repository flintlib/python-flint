from flint.flintlib.flint cimport ulong
from flint.flintlib.fmpq cimport fmpq_t

cdef extern from "bernoulli.h":
    void bernoulli_fmpq_ui(fmpq_t, ulong)
    void bernoulli_cache_compute(long n)
