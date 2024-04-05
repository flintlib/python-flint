#!/bin/bash
#
# Note: cython's Cython/Coverage.py fails for pyx files that are included in
# other pyx files. This gives the following error:
#
#   $ coverage report -m
#   Plugin 'Cython.Coverage.Plugin' did not provide a file reporter for
#   '.../python-flint/src/flint/fmpz.pyx'.
#
# A patch to the file is needed:
#
#  --- Coverage.py.backup   2022-12-09 17:36:35.387690467 +0000
#  +++ Coverage.py 2022-12-09 17:08:06.282516837 +0000
#  @@ -172,7 +172,9 @@ class Plugin(CoveragePlugin):
#          else:
#              c_file, _ = self._find_source_files(filename)
#              if not c_file:
# -                return None
# +                c_file = os.path.join(os.path.dirname(filename), 'pyflint.c')
# +                if not os.path.exists(c_file):
# +                    return None
#              rel_file_path, code = self._read_source_lines(c_file, filename)
#              if code is None:
#                  return None  # no source found
#
#

set -o errexit

source bin/activate

export PYTHON_FLINT_COVERAGE=true

# Force a rebuild of everything with coverage tracing enabled:
# touch src/flint/flintlib/*

python setup.py build_ext --inplace

coverage run -m flint.test $@

#coverage report -m
coverage html
