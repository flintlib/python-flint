Setup
===============================================================================

First install both FLINT (version 2.5 or later) and Arb (version 2.16 or later).
See:

* http://flintlib.org/
* http://arblib.org/

Python-FLINT is available on PyPI, the Python Package Index
(https://pypi.org/project/python-flint/).
The latest release can be installed using::

    pip install python-flint

Python-FLINT is also available on conda-forge
(https://anaconda.org/conda-forge/python-flint).
It can be installed using::

    conda install -c conda-forge python-flint

Python-FLINT can also be installed from a local git checkout or a source archive
as follows::

    pip install .

To build Python-FLINT manually, you may first have to install
some build dependencies::

    sudo apt-get install cython python-dev

Then run::

    python setup.py build_ext
    sudo python setup.py install

Run the test suite::

    python test/test.py

Build the documentation::

    cd doc
    make html
    cd ..

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

