import doctest
import importlib
import pkgutil
import pytest
import re

import flint

dunder_test_regex = re.compile(r'^(.*?)__test__\..*?\.(.*) \(line (\d+)\)$')


def find_doctests(module):
    finder = doctest.DocTestFinder()
    tests = []
    for module_info in pkgutil.walk_packages(module.__path__, flint.__name__ + "."):
        try:
            module = importlib.import_module(module_info.name)

            res = []
            for test in filter(lambda x: bool(x.examples), finder.find(module)):
                m = dunder_test_regex.match(test.name)
                if m is not None:
                    groups = m.groups()
                    test.name = groups[0] + groups[1]
                    test.lineno = int(groups[2])
                    res.append(test)

            tests.append((module_info.name, res))

        except Exception as e:
            print(f"Error importing {module_info.name}: {e}")
    return tests


class PyTestDocTestRunner(doctest.DocTestRunner):
    def report_failure(self, out, test, example, got):
        pytest.fail(
            "\n".join([
                f"Failed: {test.name}, line: {test.lineno}",
                "Failed example:",
                f"\t{example.source.strip()}",
                "Expected:",
                f"\t{example.want.strip()}",
                "Got:",
                f"\t{got.strip()}"
            ]),
            pytrace=False,
        )


runner = PyTestDocTestRunner()

@pytest.mark.parametrize("module,test", [(module, test) for module, test_set in find_doctests(flint) for test in test_set])
def test_docstrings(module, test):
    runner.run(test)
