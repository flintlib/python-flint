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


def run_all_benchmarks(modulename):
    """Runs all benchmarks in the current file and returns their times."""
    return {
        benchmark.__name__[len("benchmark_"):]: benchmark()
        for benchmark in get_all_benchmarks(modulename)
    }

def print_benchmark_results(results):
    """Pretty-prints a set of benchmark results."""
    print(results)
    longest_name = max(len(name) for name in results.keys())
    for name, time in results.items():
        print(f"{name:<{longest_name+2}} | {time:>10.5f}")
