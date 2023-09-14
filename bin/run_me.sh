#!/usr/bin/env bash

echo "================================================================================"
echo "Ref: https://github.com/flintlib/python-flint/issues/60#issuecomment-1682728914"
echo
echo "This folder contains useful scripts for development. To begin, set up"
echo "environment variables by running \`source bin/activate.fish\` (not executing the"
echo "script!), then build dependencies such as gmp and flint by"
echo "\`bin/build_dependencies_unix.fish\`."
echo
echo "Now to build \`python-flint\`, simply run \`bin/build_inplace.fish\`, which"
echo "places the extension modules right at \`src/flint\`. You can import by e.g."
echo
echo "    >>> import flint"
echo "    >>> print(flint.types.nmod.nmod(17, 5))"
echo
echo "Or run the tests from the command line:"
echo
echo "    $ python -m flint.test"
echo "================================================================================"
