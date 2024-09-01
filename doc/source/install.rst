Install
=======

.. _install_pip_conda:

Install with pip or conda
-------------------------

The recommended way to install ``python-flint`` for general use is to install a
prebuilt binary package. Python-FLINT is available on PyPI, the `Python Package
Index <https://pypi.org/project/python-flint/>`_. The latest release can be
installed using::

    pip install python-flint

.. note::
    If you have more than one Python environment on your system then you should
    ensure that you are installing ``python-flint`` into the correct Python
    environment. This may require using the full path to ``pip`` or using
    something like ``python3 -m pip`` to ensure that ``pip`` is run from the
    correct Python environment.

Python-FLINT is also available from `conda-forge
<https://anaconda.org/conda-forge/python-flint>`_ and can be installed with
``conda``::

    conda install -c conda-forge python-flint

Both PyPI and conda-forge provide prebuilt binary packages for ``python-flint``
for common platforms such as Windows, MacOS and Linux for a range of different
architectures and Python versions. If binaries are available for your platform
then installing from PyPI or conda-forge is simplest and does not require any
other dependencies or configuration. It is also possible that other software
distribution channels may provide binaries for ``python-flint``.

A specific version of ``python-flint`` can be installed with ``pip`` by
specifying the version number, for example::

    pip install python-flint==0.6.0

After installing ``python-flint`` you can test your installation as described
in :ref:`test_installation` below.


Fully supported platforms
-------------------------

Generally each release of python-flint will be compatible with a range of
Python versions as described in `SPEC 0
<https://scientific-python.org/specs/spec-0000/>`_. At the time of writing, the
current release of ``python-flint`` is ``0.6.0`` and binaries are provided for
Python 3.9, 3.10 3.11 and 3.12 for the following platforms:

- Windows 64-bit (``x86_64``)
- MacOS 64-bit Intel and 64-bit ARM (i.e. Apple Silicon)
- Linux 64-bit (``x86_64``, ``manylinux``)

A WASM build of ``python-flint`` is also available for use in the browser
via `Pyodide <https://pyodide.org/en/stable/usage/packages-in-pyodide.html>`_.

You can see which binaries are provided for the latest release of
``python-flint`` `from PyPI <https://pypi.org/project/python-flint/#files>`_
and which are provided `from conda-forge
<https://anaconda.org/conda-forge/python-flint>`_. A list of supported versions
of Python and other ``python-flint`` dependencies can be found at
:ref:`supported_versions`.


Platforms without binaries
--------------------------

There are many other platforms on which ``python-flint`` works fine but for
which binaries are not provided. If a binary is not available for your platform
then you may be able to build from source as described in
:ref:`simple_build_instructions`.

Notably, at the time of writing the following platforms do *not* have binaries
available but ``python-flint`` should work if built from source:

- Linux aarch64 (``conda-forge`` has binaries but ``PyPI`` does not, see
  `gh-105 <https://github.com/flintlib/python-flint/issues/105>`_).
- non-glibc Linux distros (e.g. ``musllinux`` rather than ``manylinux``)
- PyPy

Binaries for Linux aarch64 will likely be added in future when the platform is
available for testing in CI.


Unsupported platforms
---------------------

It is *not* known or expected that ``python-flint`` will currently work on:

- Windows on ARM (has never been tested)
- Any 32-bit platform (previously worked but has not been tested for some
  time)
- Other Python implementations besides CPython and PyPy (e.g. GraalPython,
  Jython, IronPython)

Support for Windows on ARM will likely be added in future when the platform is
available for testing in CI.


.. _test_installation:

Testing your installation
-------------------------

However you install ``python-flint``, you can test the installation by running
the test suite (for ``python-flint >= 0.5.0``)::

    python -m flint.test

This does not take long and will run all the tests and doctests and should
hopefully show something like this::

    $ python -m flint.test --quiet
    Running tests...
    Running doctests...
    flint.test: all 54 tests passed!
    flint.test: all 4283 doctests passed!
    ----------------------------------------
    OK: Your installation of python-flint seems to be working just fine!

From Python you can instead run the same tests with::

    from flint.test.__main__ import main
    main()

If you have installed ``python-flint`` from PyPI or conda-forge as described
above and your installation passes the tests then you are ready to use
``python-flint``. If the tests fail (or do not complete) then please report the
problem on the `python-flint issue tracker
<https://github.com/flintlib/python-flint/issues>`_.

If the tests have all passed then ``python-flint`` is installed and ready to
use!
