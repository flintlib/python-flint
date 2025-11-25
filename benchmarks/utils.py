"""Utility functions for benchmarking."""

import inspect
import sys


def get_all_benchmarks(modulename):
    """Returns a list of all functions in the current module named 'benchmark_*'."""
    return [
        obj
        for name, obj in inspect.getmembers(sys.modules[modulename])
        if inspect.isfunction(obj) and name.startswith("benchmark_")
    ]
