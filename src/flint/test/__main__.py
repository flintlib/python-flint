#
# Run the python-flint test-suite as
#
#    python -m flint.test
#

import sys
import doctest
import traceback
import argparse

import flint
from flint.test.test_all import all_tests
from flint.test.test_docstrings import find_doctests


def run_tests(verbose=None):
    """Run all tests

    This is usually called by

        $ python -m flint.flint
    """
    # Show the limited output by default
    if verbose is None:
        verbose = True

    total = 0
    failed = 0

    for test in all_tests:

        if verbose:
            print(f'{test.__name__}...', end='', flush=True)

        try:
            test()
        except Exception as e:
            print(f'Error in {test.__name__}')
            if verbose:
                traceback.print_exc()
            failed += 1
        else:
            if verbose:
                print('OK')

        total += 1

    return failed, total


def run_doctests(tests, verbose=False):
    runner = doctest.DocTestRunner()
    for module, test_set in tests:
        if test_set:
            print(f"{module}...", end="" if not verbose else "\n", flush=True)
            for test in test_set:
                if verbose:
                    print("\tTesting:", test.name)
                runner.run(test)
            print("OK")

    return runner.summarize()


def run_all_tests(tests=True, doctests=True, verbose=None):

    success = True

    if tests:
        print("Running tests...")
        t_failed, t_total = run_tests(verbose=verbose)

    if doctests:
        print("Running doctests...")
        d_failed, d_total = run_doctests(find_doctests(flint), verbose=verbose)

    if tests:
        if t_failed:
            print(f'flint.test: {t_failed} of {t_total} tests failed')
            success = False
        else:
            print(f'flint.test: all {t_total} tests passed!')

    if doctests:
        if d_failed:
            print(f'flint.test: {d_failed} of {d_total} doctests failed')
            success = False
        else:
            print(f'flint.test: all {d_total} doctests passed!')

    return success


def main(*args):
    """Run the python-flint test-suite"""

    parser = argparse.ArgumentParser(description="Run the python-flint test-suite")
    parser.add_argument("--quiet", "-q", action="store_true", help="less verbose output")
    parser.add_argument("--verbose", "-v", action="store_true", help="more verbose output")
    parser.add_argument("--tests", "-t", action="store_true", help="run tests")
    parser.add_argument("--doctests", "-d", action="store_true", help="run doctests")
    args = parser.parse_args(args)

    if not args.tests and not args.doctests:
        # Default is run all tests:
        tests = True
        doctests = True
    else:
        # Either --tests or --doctests was specified
        tests = args.tests
        doctests = args.doctests

    # Default is show output from tests but keep the doctests quiet.
    if args.verbose:
        verbose = True
    elif args.quiet:
        verbose = False
    else:
        verbose = None

    success = run_all_tests(tests=tests, doctests=doctests, verbose=verbose)

    if not success:
        print("----------------------------------------")
        print("!!!FAILED!!!: Something is wrong with your installation of python-flint!")
        print("----------------------------------------")
        return 1
    else:
        print("----------------------------------------")
        print("OK: Your installation of python-flint seems to be working just fine!")
        print("----------------------------------------")
        return 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:]))
