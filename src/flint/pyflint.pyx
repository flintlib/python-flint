"""
Python wrapper for FLINT and Arb.
"""

# cimport flint
cimport libc.stdlib

from flint.flint_base.flint_context cimport thectx

cdef flint_rand_t global_random_state
flint_rand_init(global_random_state)

ctx = thectx
