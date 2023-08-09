Setup
===============================================================================

First install both FLINT (version 2.9 or later) and Arb (version 2.23 or later).
See:

* http://flintlib.org/
* http://arblib.org/

Python-FLINT is available on PyPI, the Python Package Index
(https://pypi.org/project/python-flint/).
The latest release can be installed using::

    pip install python-flint

Binary wheels are provided for Windows amd64, Linux (manylinux 2_17) x86_64,
macOS x86_64 and macOS arm64. For other platforms, pip will attempt to build
Python-FLINT from source which requires a C compiler and the FLINT and Arb
header files and library files (libflint.so and libarb.so) to be available as
well as the Python development headers and Cython and numpy.

Python-FLINT is also available on conda-forge for Linux and macOS.
(https://anaconda.org/conda-forge/python-flint).
It can be installed using::

    conda install -c conda-forge python-flint

Python-FLINT can also be installed from a local git checkout or a source archive
as follows::

    pip install .

To build Python-FLINT manually, you first need to install some build
dependencies::

    pip install Cython numpy

Then run::

    python setup.py build_ext
    python setup.py install

Run the test suite::

    python -m flint.test

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

Build all dependencies from source
----------------------------------

From a VCS checkout, to build python-flint and all dependencies from source,
using the exact versions that are tested in CI and used for the binary PyPI
wheels, run the following in a unix shell::

    source bin/activate
    bin/build_dependencies_unix.sh

The script will download and build GMP, MPFR, FLINT and Arb and build them all
in a ``.local`` directory. The ``bin/activate`` script sets the appropriate
path environment variables for C headers and libraries which is needed for
the ``build_dependencies_unix.sh`` script to work. After running the script,
you can then build Python-FLINT in place with::

    python setup.py build_ext --in-place

and run the test suite with::

    python -m flint.test

This way of building Python-FLINT depends on the ``bin/activate`` script to
locate the shared libraries at runtime. The script will also set ``PYTHONPATH``
so that the in-place build of Python-FLINT can be imported.

These steps will also work under MinGW with the mingw64 toolchain, but you
should first run::

    echo '[build]' > setup.cfg
    echo 'compiler = mingw32' >> setup.cfg

    # Install the mingw-w64 toolchain
    pacman -S --noconfirm mingw-w64-x86_64-gcc m4 make mingw-w64-x86_64-tools-git

To change the versions of the dependencies that are built, edit the
``bin/build_variables.sh`` script.
