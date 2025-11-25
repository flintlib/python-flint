"""Benchmarks for basic operations of a number of flint types.

These benchmarks are written using pyperf, and can be run by 

python benchmarks/simple_benchmarks.py

Since each benchmark is very short ("microbenchmarks"), they may
require an increased number of samples to produce statistically
significant results (via the --values flag).
"""

import pyperf
from utils import get_all_benchmarks

NUM_EXECUTIONS = 10000000


def benchmark_arb_addition(runner: pyperf.Runner):
    """Simple benchmark for adding two arbs."""
    runner.timeit(
        name="arb addition",
        setup=[
            "from flint import arb",
            "a = arb(1, 2)",
            "b = arb.pi()",
        ],
        stmt="a + b"
    )


def benchmark_arb_multiplication(runner: pyperf.Runner):
    """Simple benchmark for multiplying two arbs."""
    runner.timeit(
        name="arb multiplication",
        setup=[
            "from flint import arb",
            "a = arb(1, 2)",
            "b = arb.pi()",
        ],
        stmt="a * b"
    )


def benchmark_arb_contains(runner: pyperf.Runner):
    """Simple benchmark for comparing two arbs."""
    runner.timeit(
        name="arb contains",
        setup=[
            "from flint import arb",
            "a = arb(1, 2)",
            "b = arb.pi()",
        ],
        stmt=[
            "a in b",
            "b in a",
        ]
    )


def benchmark_fmpz_addition(runner: pyperf.Runner):
    """Simple benchmark for adding two fmpzs."""
    runner.timeit(
        name="fmpz addition",
        setup=[
            "from flint import fmpz",
            "a = fmpz(1)",
            "b = fmpz(6)",
        ],
        stmt="a + b"
    )

def benchmark_fmpz_multiplication(runner: pyperf.Runner):
    """Simple benchmark for multiplying two fmpzs."""
    runner.timeit(
        name="fmpz multiplication",
        setup=[
            "from flint import fmpz",
            "a = fmpz(1)",
            "b = fmpz(6)",
        ],
        stmt="a * b"
    )


def benchmark_fmpz_eq(runner: pyperf.Runner):
    """Simple benchmark for equality on two fmpzs."""
    runner.timeit(
        name="fmpz equality",
        setup=[
            "from flint import fmpz",
            "a = fmpz(1)",
            "b = fmpz(6)",
        ],
        stmt=["a == b", "b == a"]
    )


def benchmark_fmpq_addition(runner: pyperf.Runner):
    """Simple benchmark for adding two fmpqs."""
    runner.timeit(
        name="fmpq addition",
        setup=[
            "from flint import fmpq",
            "a = fmpq(1, 2)",
            "b = fmpq(15, 7)",
        ],
        stmt="a + b"
    )


def benchmark_fmpq_multiplication(runner: pyperf.Runner):
    """Simple benchmark for multiplying two fmpqs."""
    runner.timeit(
        name="fmpq multiplication",
        setup=[
            "from flint import fmpq",
            "a = fmpq(1, 2)",
            "b = fmpq(15, 7)",
        ],
        stmt="a * b"
    )


def benchmark_fmpq_eq(runner: pyperf.Runner):
    """Simple benchmark for equality on two fmpqs."""
    runner.timeit(
        name="fmpq equality",
        setup=[
            "from flint import fmpq",
            "a = fmpq(1, 2)",
            "b = fmpq(15, 7)",
        ],
        stmt=["a == b", "b == a"]
    )


def benchmark_acb_addition(runner: pyperf.Runner):
    """Simple benchmark for adding two acbs."""
    runner.timeit(
        name="acb addition",
        setup=[
            "from flint import acb",
            "a = acb(1 + 3j)",
            "b = acb.pi()",
        ],
        stmt="a + b"
    )


def benchmark_acb_multiplication(runner: pyperf.Runner):
    """Simple benchmark for multiplying two acbs."""
    runner.timeit(
        name="acb multiplication",
        setup=[
            "from flint import acb",
            "a = acb(1 + 3j)",
            "b = acb.pi()",
        ],
        stmt="a * b"
    )


def benchmark_acb_eq(runner: pyperf.Runner):
    """Simple benchmark for containment on two acbs."""
    runner.timeit(
        name="acb contains",
        setup=[
            "from flint import acb",
            "a = acb(1 + 3j)",
            "b = acb.pi()",
        ],
        stmt=["a in b", "b in a"]
    )

def main():
    """Run all the benchmarks."""
    runner = pyperf.Runner()
    benchmarks = get_all_benchmarks(__name__)
    for benchmark in benchmarks:
        benchmark(runner)

if __name__ == "__main__":
    main()
