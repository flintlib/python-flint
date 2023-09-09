"""
Python wrapper for FLINT and Arb.
"""

cimport flint
cimport libc.stdlib
cimport cython

from flint.flint_base.flint_context cimport thectx

cdef flint_rand_t global_random_state
flint_randinit(global_random_state)

ctx = thectx
