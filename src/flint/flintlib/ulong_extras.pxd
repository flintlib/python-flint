from flint._flint cimport ulong

cdef extern from "flint/ulong_extras.h":
    ulong n_gcd(ulong n, ulong k)
    int n_is_prime(ulong n)
