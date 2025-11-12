"""Benchmarks for basic operations of a number of flint types."""

import timeit
from flint import arb, fmpz, fmpq, acb
from utils import run_all_benchmarks, print_benchmark_results

NUM_EXECUTIONS = 10000000


def benchmark_arb_addition(num_executions = NUM_EXECUTIONS):
    """Simple benchmark for adding two arbs."""
    a = arb(1, 2)
    b = arb.pi()

    return timeit.timeit(lambda: a + b, number=num_executions)


def benchmark_arb_multiplication(num_executions = NUM_EXECUTIONS):
    """Simple benchmark for multiplying two arbs."""
    a = arb(1, 2)
    b = arb.pi()

    return timeit.timeit(lambda: a * b, number=num_executions)


def benchmark_arb_contains(num_executions = NUM_EXECUTIONS):
    """Simple benchmark for comparing two arbs."""
    a = arb(1, 2)
    b = arb.pi()

    def to_time():
        first = a in b
        second = a in a
        return first, second

    return timeit.timeit(to_time, number=num_executions)


def benchmark_fmpz_addition(num_executions = NUM_EXECUTIONS):
    """Simple benchmark for adding two fmpzs."""
    a = fmpz(1)
    b = fmpz(6)

    return timeit.timeit(lambda: a + b, number=num_executions)


def benchmark_fmpz_multiplication(num_executions = NUM_EXECUTIONS):
    """Simple benchmark for multiplying two fmpzs."""
    a = fmpz(1)
    b = fmpz(6)

    return timeit.timeit(lambda: a * b, number=num_executions)


def benchmark_fmpz_eq(num_executions = NUM_EXECUTIONS):
    """Simple benchmark for equality on two fmpzs."""
    a = fmpz(1)
    b = fmpz(6)

    def to_time():
        first = a == b
        second = a == a
        return first, second

    return timeit.timeit(to_time, number=num_executions)


def benchmark_fmpq_addition(num_executions = NUM_EXECUTIONS):
    """Simple benchmark for adding two fmpqs."""
    a = fmpq(1, 2)
    b = fmpq(15, 7)

    return timeit.timeit(lambda: a + b, number=num_executions)


def benchmark_fmpq_multiplication(num_executions = NUM_EXECUTIONS):
    """Simple benchmark for multiplying two fmpqs."""
    a = fmpq(1, 2)
    b = fmpq(15, 7)

    return timeit.timeit(lambda: a * b, number=num_executions)


def benchmark_fmpq_eq(num_executions = NUM_EXECUTIONS):
    """Simple benchmark for equality on two fmpqs."""
    a = fmpq(1, 2)
    b = fmpq(15, 7)

    def to_time():
        first = a == b
        second = a == a
        return first, second

    return timeit.timeit(to_time, number=num_executions)


def benchmark_acb_addition(num_executions = NUM_EXECUTIONS):
    """Simple benchmark for adding two acbs."""
    a = acb(1 + 3j)
    b = acb.pi()

    return timeit.timeit(lambda: a + b, number=num_executions)


def benchmark_acb_multiplication(num_executions = NUM_EXECUTIONS):
    """Simple benchmark for multiplying two acbs."""
    a = acb(1 + 3j)
    b = acb.pi()

    return timeit.timeit(lambda: a * b, number=num_executions)


def benchmark_acb_eq(num_executions = NUM_EXECUTIONS):
    """Simple benchmark for equality on two acbs."""
    a = acb(1 + 3j)
    b = acb.pi()

    def to_time():
        first = a in b
        second = a in a
        return first, second

    return timeit.timeit(to_time, number=num_executions)

if __name__ == "__main__":
    print_benchmark_results(run_all_benchmarks(__name__))