Install/Build from source
=========================

.. note::
   The instructions here are for building ``python-flint`` from source. For
   most users it is recommended to install prebuilt binaries from ``PyPI`` or
   ``conda-forge`` instead. The instructions here are only needed if a binary
   is not available for the platform. See :ref:`install_pip_conda`.


.. _simple_build_instructions:

Simple build instructions
-------------------------

The simple explanation of how to build ``python-flint`` from source is that
there are two steps:

- Install ``FLINT >= 3.0`` (see :ref:`installing_the_dependencies` below).
- Run ``pip install --no-binary python-flint python-flint``.

For example on Ubuntu 24.04 (but not older versions of Ubuntu) and when installing
``python-flint >= 0.7.0`` these two steps are::

    sudo apt-get install libflint-dev
    pip install --no-binary python-flint python-flint

The first command installs ``FLINT 3.0.1`` system-wide. With the second command
``pip`` will download the source code for the latest release of
``python-flint`` from PyPI, build it and install it into the active Python
environment. When building, ``pip`` will create a temporary isolated build
environment and will install the Python build dependencies (``Cython``,
``meson``, ...) into this environment so it is not necessary to install them
manually before running ``pip install``.

If you have the source code locally then you can build and install with::

    pip install path/to/python-flint-directory-or-archive

After installing from source it is recommended to run the tests to check that
everything is working correctly as described in :ref:`test_installation`.

The remainder of this page provides more detailed instructions for:

- :ref:`supported_versions`.
- :ref:`building_from_source`.
- :ref:`building_older_versions`.
- :ref:`installing_the_dependencies`.
- :ref:`building_on_windows`.
- :ref:`non_standard_location`.
- :ref:`editable_install`.

.. note::
   If you have more than one Python environment in your system then you need to
   ensure that you are installing ``python-flint`` into the correct one. This
   may require using the full path to ``pip`` or something like ``python3 -m
   pip`` or you may need to activate the environment first.


.. _supported_versions:

Compatibility and supported versions
------------------------------------

.. note::
   The compatibility information here is mostly only relevant when building
   ``python-flint`` from source. For most users it is recommended to install
   pre-built binaries from ``PyPI`` or ``conda-forge``. See
   :ref:`install_pip_conda`.

Generally each release of python-flint will be compatible with a range of
Python versions as described in `SPEC 0
<https://scientific-python.org/specs/spec-0000/>`_. Since python-flint 0.5.0
the minimum supported FLINT version is ``3.0`` and each release of python-flint
supports all versions of ``FLINT >= 3.0`` available at the time of release.

Compatible versions (note that 0.7.0 is not yet released):

.. list-table:: python-flint compatibility
   :header-rows: 1

   * - python-flint
     - Release date
     - CPython
     - FLINT
     - Cython
   * - 0.7.0
     - Not yet
     - 3.10-3.13
     - 3.0-3.2?
     - 3.0-3.1?
   * - 0.6.0
     - 1st Feb 2024
     - 3.9-3.12
     - 3.0
     - 3.0
   * - 0.5.0
     - 22nd Oct 2023
     - 3.9-3.12
     - 3.0
     - 3.0
   * - 0.4.0
     - 8th Aug 2023
     - 3.9-3.11
     - ``2.9.0`` (``Arb 2.23.0``)
     - 3.0
   * - 0.3.0
     - 7th Dec 2018
     - older Python versions
     - ``< 3.0``
     - ``< 3.0``

The minimum and maximum versions of Python represent the versions that are
tested in CI and for which binaries are provided on PyPI. It is likely that
``python-flint`` will work with other versions of Python (particularly older
Python versions) but this is not tested. It is possible that ``conda-forge``
may provide binaries for other versions of Python.

The minimum versions of Cython and FLINT are needed because it is known that
python-flint will not even build with older versions of these libraries. The
maximum versions of all dependencies are speculative and are based on the
versions that are known to work at the time of release. It is possible that
newer versions of Cython and FLINT will work but unlikely. During the year
following the release of ``python-flint 0.4.0`` every non-patch release of
Cython, FLINT, or CPython has required changes to the ``python-flint`` source
code to be able to build at all. In particular the following releases of
Cython, FLINT and CPython had changes that would prevent building all versions
of ``python-flint`` existing at the time:

- Flint 3.0 (Arb and Flint merged, lots of changes)
- Flint 3.1 (Function signature for ``fmpz_mod_mat`` changed)
- Flint 3.2 (``flint_randinit`` function name changed)
- Cython 3.0 (Handling of dunders changed)
- Cython 3.1 (Removal of ``PyInt_*`` functions)
- CPython 3.12 (Removal of distutils)

It is expected then that any future untested ``3.x`` version of Cython, FLINT,
or CPython will not be compatible with past versions of ``python-flint`` which
is why the table above lists the versions that were known to work at the time
of release.

As of python-flint 0.7.0, CPython ``3.13t`` free-threaded builds are tested in
CI but wheels are not provided on PyPI. There are no known issues related to
using python-flint in a `PEP 703 <https://peps.python.org/pep-0703/>`_
free-threaded build but it is likely that mutating objects shared by multiple
threads is not safe.

It is also possible to build and use python-flint for PyPy. Other Python
implementations may work but are not tested.


.. _building_from_source:

Installing python-flint from source
-----------------------------------

.. note::
   The instructions here are for building ``python-flint`` from source. For
   most users it is recommended to install prebuilt binaries from ``PyPI`` or
   ``conda-forge`` instead. The instructions here are only needed if a binary
   is not available for the platform. See :ref:`install_pip_conda`.

   Also if you are working on ``python-flint`` itself then it is not
   recommended to install the package as described here. Instead see the
   :ref:`development_workflow` page for how to work on ``python-flint``.

The source code for ``python-flint`` is available on `GitHub
<https://github.com/flintlib/python-flint/tags>`_ and source distributions can
be downloaded from PyPI.

To build from source you must first install the dependencies (see
:ref:`installing_the_dependencies` below for instructions). Once the
dependencies are installed the following command will download the
``python-flint`` source code from PyPI, then build and install it into the
active Python environment::

    pip install python-flint

This will try to install a binary first but will otherwise download, build and
install the latest release of ``python-flint`` from PyPI. If you definitely
want to build from source then you can use the ``--no-binary`` option::

    pip install --no-binary python-flint python-flint

To install a specific version of ``python-flint`` from PyPI use e.g.::

    pip install python-flint==0.7.0a4

To build and install the latest ``python-flint`` from git master you can
use::

    pip install git+https://github.com/flintlib/python-flint.git@master

If you already have the source code downloaded or checked out from git, you can
``cd`` in and build and install with::

    pip install .

Alternatively if you would like to build a wheel you can use
``pypa/build`` (first ``pip install build``)::

    python -m build

Note that wheels built in this way will not include the dependencies (unlike
those distributed on PyPI) and cannot generally be installed on other systems.

Since ``python-flint 0.7.0`` the build system is ``meson`` and the build
requirements and version constraints are listed in ``pyproject.toml``. When
using build isolation the build requirements are installed in a temporary
virtual environment and so it should not be necessary to install them in the
active Python environment before running ``pip install``.

To build without build isolation with ``python-flint >= 0.7.0`` the 
dependencies should first be installed in the active Python environment::

    pip install Cython==3.0 meson meson-python ninja
    pip install --no-build-isolation .

The ``meson`` build system will detect the versions of ``FLINT`` and Cython
installed in the system and will fail if they are not versions that were known
to be compatible at the time of the release of ``python-flint``. To build
against new, untested versions of ``FLINT`` or Cython you can pass the
``-Dflint_version_check=false`` option to the ``meson`` build system::

    pip install --config-settings=setup-args="-Dflint_version_check=false" .

This is useful for testing new versions of ``FLINT`` with ``python-flint`` for
example if you want to build ``python-flint`` against the latest git version of
``FLINT``. See :ref:`supported_versions` above for the versions of ``FLINT``
and Cython that are supported by each version of ``python-flint``.


.. _building_older_versions:

Installing older versions from source
-------------------------------------

For ``python-flint < 0.6.0`` the source distribution did not include
``pyproject.toml`` and did not list the build requirements. Also for
``python-flint < 0.7.0`` the build requirements were different and there were
no version constraints listed on the dependencies. An list of the build
requirements for older versions of ``python-flint`` is given above in
:ref:`supported_versions`.

For ``python-flint < 0.7.0`` you will need to install the build requirements
manually, pin the version of Cython, and disable build isolation::

    pip install Cython==3.0 setuptools numpy
    pip install --no-build-isolation .

For ``python-flint < 0.4.0`` older versions of Cython are needed (``<= 0.29``).
If the build fails during the Cython step then it is likely that a different
version of Cython is needed.


.. _installing_the_dependencies:

Installing the dependencies
---------------------------

.. note::
    It is not necessary to install the dependencies manually if you install
    from PyPI or conda-forge as is recommended. When installing with ``conda``
    the packages for the dependencies will also be installed from conda-forge
    automatically. The binaries on PyPI are built with the dependencies bundled
    in the wheel so that they do not need to be installed separately.

    The following instructions are only for when building ``python-flint`` from
    source if needed because a binary is not available for your platform. See
    :ref:`install_pip_conda`.

The dependencies for building ``python-flint`` have changed over time. See
:ref:`supported_versions` above for the versions of the dependencies that are
supported by each version of ``python-flint``.

As of ``python-flint 0.7.0`` the runtime dependencies are Python and FLINT (at
least version 3.0) and the build-time dependencies are a C compiler,
``Cython``, ``meson``, ``meson-python`` and ``ninja``. Commands shown above
such as ``pip install .`` will install dependencies like ``Cython``, ``meson``
etc automatically. If you already have Python and a C compiler then what needs
to be installed before building ``python-flint`` is ``FLINT``.

At the time of writing, few Linux distributions provide ``FLINT >= 3.0`` in
their package repositories but for example on ``Ubuntu 24.04`` (but not any
earlier Ubuntu versions) you can install ``FLINT 3.0.1`` with::

    sudo apt-get install libflint-dev

On MacOS you can install FLINT from homebrew with::

    brew install flint

Other package managers may also provide ``FLINT`` but make sure that it is at
least version ``3.0``.

Once ``FLINT`` is installed it should be possible to build ``python-flint``
with any of the commands shown above e.g.::

    pip install .

If it is not possible to install FLINT from a package manager then you need to
install GMP and MPFR and then build FLINT. You may still be able to install GMP
and MPFR from a package manager for example on Ubuntu::

    sudo apt-get install libgmp-dev libmpfr-dev

The python-flint git repo has a script `bin/install_flint_ubuntu.sh
<https://github.com/flintlib/python-flint/blob/master/bin/install_flint_ubuntu.sh>`_
that uses ``apt-get`` to install all dependencies needed to build ``FLINT``,
then builds ``FLINT`` from git using a specified git ref, and then installs
``FLINT`` system-wide::

    bin/install_flint_ubuntu.sh v3.0.1  # version 3.0.1
    bin/install_flint_ubuntu.sh main    # latest git

The script can be adapted for other Linux distributions or MacOS to use
something other than ``apt-get`` to install dependencies.

If the whole stack needs to be built from source then download the source for
all three (`GMP <https://gmplib.org/#DOWNLOAD>`_, `MPFR
<https://www.mpfr.org/mpfr-current/>`_, `FLINT
<https://flintlib.org/downloads.html>`_) and build each with the standard::

    ./configure
    make
    make install

Adapt the ``configure`` commands as needed. Once these are installed you should
again be able to install ``python-flint`` with::

    pip install .


.. _building_on_windows:

Installing from source on Windows
---------------------------------

.. note::
   Building from source is not the recommended way for most users to install
   ``python-flint``, especially on Windows. For most users it is recommended to
   use the binaries from ``PyPI`` or ``conda-forge`` except in cases where a
   binary is not available for the platform. See :ref:`install_pip_conda`.

The instructions in :ref:`installing_the_dependencies` above are for Unix-like systems
(e.g. Linux or MacOS). On Windows the dependencies can be built in a similar
way using MSYS2 or under WSL. It is also possible to build ``python-flint`` and
its dependencies using MSVC but we do not currently provide instructions for
this. The `conda-forge recipe
<https://github.com/conda-forge/python-flint-feedstock>`_ for ``python-flint``
builds on Windows using MSVC.

The `MSYS2 <https://www.msys2.org/>`_ project provides a Unix-like environment
for Windows and a package manager that can be used to install the dependencies.
The git repo for ``python-flint`` has a script `bin/cibw_before_all_windows.sh
<https://github.com/flintlib/python-flint/blob/master/bin/cibw_before_all_windows.sh>`_
that installs the dependencies under MSYS2 and builds ``GMP``, ``MPFR``,
``FLINT``. This script is used for building the Windows binaries for PyPI. We
use the ``MinGW64`` (``mingw-w64-x86_64``) toolchain for building on Windows
rather than MSVC because it makes it possible to have a fat build of ``GMP``
(``--enable-fat``) which bundles micro-architecture specific optimisations for
``x86_64`` in a redistributable shared library. This is important for
performance on modern ``x86_64`` CPUs and is not possible if building ``GMP``
with MSVC. Since we need to use ``MinGW64`` for building ``GMP`` it is simplest
to use it for building ``MPFR``, ``FLINT`` and ``python-flint`` as well and
means that the same Unix-style build scripts can be used for all platforms.

The ``python-flint`` project does not have much experience using MSVC. Possibly
it would be better to build ``GMP`` using ``MinGW64`` and then build ``MPFR``,
``FLINT`` and ``python-flint`` using MSVC. It is also possible that it would be
better to build ``GMP``, ``MPFR``, ``FLINT`` using MinGW64 and then build
``python-flint`` using MSVC. Someone with more experience with MSVC would need
to help with this. We would welcome contributions that explain how to build
``python-flint`` and its dependencies using MSVC and/or that improve the build
process for distributed binaries on Windows.


.. _non_standard_location:

Using ``FLINT`` from a non-standard location
--------------------------------------------

.. note::
    This section is only relevant when building ``python-flint`` from source.
    For most users it is recommended to use the binaries from ``PyPI`` or
    ``conda-forge``. See :ref:`install_pip_conda`. The instructions here are
    also not needed if you have installed ``FLINT`` and its dependencies
    system-wide (e.g. using a package manager like ``apt-get`` or ``brew``).

If you have installed ``FLINT`` in a non-standard location then you will need
to instruct the ``python-flint`` build system where to find it and ensure that
the ``FLINT`` shared library can be found at runtime.

Since ``python-flint 0.7.0`` the build system is `meson
<https://mesonbuild.com/>`_ and uses `pkg-config
<https://www.freedesktop.org/wiki/Software/pkg-config/>`_ to find the
dependencies ``FLINT``, ``GMP`` and ``MPFR``. If these are installed in a
non-standard location then you can set the ``PKG_CONFIG_PATH`` environment
variable to point to the directory containing the ``.pc`` files for these
libraries. For example if you have installed ``FLINT`` in ``~/.local`` then you
can set the environment variable like this::

    export PKG_CONFIG_PATH=$(pwd)/.local/lib/pkgconfig

Note that in some systems the ``lib/pkgconfig`` directory may be in a different
location e.g. ``lib64/pkgconfig``. It is also possible to pass the path to the
``pkg-config`` files to the ``meson-python`` build backend. For example if
building with ``pip``::

    pip install \
        --config-settings=setup-args="--pkg-config-path=$(pwd)/.local/lib/pkgconfig" \
        python-flint

Setting the path to the ``pkg-config`` files in this way will allow the
``python-flint`` build system to find the ``FLINT`` library at build time. At
runtime the ``GMP``, ``MPFR`` and ``FLINT`` shared libraries must be in a
location where the dynamic linker can find them. On Linux the environment
variable ``LD_LIBRARY_PATH`` can be used to add the directory containing the
shared libraries to the search path. On MacOS the environment variable is
``DYLD_LIBRARY_PATH`` and on Windows it is ``PATH``. For example on Linux if
``FLINT`` is installed in ``~/.local/lib`` then you can set the environment
variable::

    export LD_LIBRARY_PATH=$(pwd)/.local/lib

Using the environment variable like this means that it needs to be set every
time you run Python and use ``python-flint`` (the git repo provides ``source
bin/activate`` for doing this). A better option on Unix-like systems is to
install ``RPATH`` entries into the ``python-flint`` extension modules. On some
platforms this is done automatically by the ``meson`` build system but on
others it needs to be enabled explicitly. This can be done by passing the
``-Dadd_flint_rpath=true`` option to the ``meson`` build system::

    pip install \
        --config-settings=setup-args="--pkg-config-path=$(pwd)/.local/lib/pkgconfig" \
        --config-settings=setup-args="-Dadd_flint_rpath=true" \
        python-flint

For versions of ``python-flint`` before ``0.7.0`` the build system is
``setuptools`` (or ``numpy.distutils`` for ``Python < 3.12``). In this case
``pkg-config`` is not used. The following environment variables can be used to
set the location of the ``FLINT`` and other shared libraries at build time or
runtime::

    C_INCLUDE_PATH=$(pwd)/.local/include  # build-time
    LIBRARY_PATH=$(pwd)/.local/lib        # build-time
    LDFLAGS=-Wl,-rpath=$(pwd)/.local/lib  # build-time Linux or MacOS
    LD_LIBRARY_PATH=$(pwd)/.local/lib     # run-time Linux
    DYLD_LIBRARY_PATH=$(pwd)/.local/lib   # run-time MacOS
    PATH=$(pwd)/.local/bin:$PATH          # run-time Windows

A future improvement for ``python-flint`` could be if the meson build system
could build all dependencies (``GMP``, ``MPFR``, ``FLINT``) as shared libraries
and bundle them into ``python-flint`` although `this is not currently possible
with meson-python
<https://github.com/mesonbuild/meson-python/discussions/410>`_. Otherwise
perhaps it could be possible to link ``FLINT`` and the other libraries
statically into ``python-flint``.


.. _editable_install:

Installing in editable mode
---------------------------

.. note::
   For working on ``python-flint`` itself it is not recommended to install the
   package into the active Python environment. Instead the development workflow
   uses ``spin`` and ``meson`` to manage a local build of ``python-flint``. See
   the :ref:`development_workflow` page for more information on how to develop
   ``python-flint``.

If you are building and testing ``python-flint`` while working on another
project then it may be useful to install ``python-flint`` in editable mode.
This allows making changes to the code of ``python-flint`` and seeing the
changes reflected in the other environment without needing to reinstall
``python-flint`` each time. This might be useful for example if you are using
``git bisect`` to find a change in ``python-flint`` (although it will not work
if you go back to versions before ``0.7.0``).

Since ``0.7.0`` it is possible to install ``python-flint`` as a
`meson-python editable install
<https://meson-python.readthedocs.io/en/latest/how-to-guides/editable-installs.html>`_.
To install ``python-flint`` in editable mode, first install ``FLINT`` and
then::

    git clone https://github.com/flintlib/python-flint.git
    cd python-flint
    pip install meson meson-python cython ninja
    pip install --no-build-isolation --editable .
    python -m flint.test  # recommended if you have made changes

This requires ``--no-build-isolation`` so that the build directory is not
deleted after install. Once installed in editable mode, each time Python is
restarted and ``python-flint`` is imported (``import flint``) an import hook
will check if the source code has changed and if so will rebuild the extension
modules and update the Python files. The rebuild uses ``meson`` for fast,
parallel, incremental rebuilds. Note that for the rebuild to happen and for the
changes to take effect it is necessary to start a new Python process e.g. by
running ``python`` again or by restarting the Jupyter kernel.

If you have installed ``FLINT`` in a non-standard location then you should set
the ``pkg-config`` path as described in :ref:`non_standard_location` above::

    pip install
        --no-build-isolation \
        --config-settings=setup-args="--pkg-config-path=$(pwd)/.local/lib/pkgconfig" \
        --editable .

To fully remove the editable install you can run::

    pip uninstall python-flint

and then delete the ``build`` directory that was created in the root of the
``python-flint`` git repo.
