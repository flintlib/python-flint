Setup
===============================================================================

First install both FLINT and Arb (currently, the git versions are required).
See:

* https://github.com/fredrik-johansson/flint2/
* https://github.com/fredrik-johansson/arb/

Install build dependencies::

    sudo apt-get install cython python-dev

Build and install python-flint from source::

    python setup.py build_ext
    sudo python setup.py install

Run the test suite::

    python test/test.py

Import using::

    >>> from flint import *

Additional paths
----------------

The FLINT and Arb header files and library files (libflint.so and libarb.so)
must be available at compile time. If they are in a nonstandard location
(for example, if they have been built but not installed),
use a command such as the following to build::

    python ./setup.py build_ext --include-dirs=/home/fredrik/src/flint2:/home/fredrik/src/arb --library-dirs=/home/fredrik/src/flint2:/home/fredrik/src/arb

Likewise, before starting the Python interpreter, tell the linker
where to find the library files using something like::

    export LD_LIBRARY_PATH=/home/fredrik/src/flint2:/home/fredrik/src/arb:$LD_LIBRARY_PATH

You may also have to install the CPimport file::

    sudo cp /home/fredrik/src/flint2/qadic/CPimport.txt /usr/local/share/flint/CPimport.txt

