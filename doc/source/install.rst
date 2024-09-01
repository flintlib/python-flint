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

At the time of writing the next release of ``python-flint`` will be ``0.7.0``
and binaries will be provided for Python 3.10, 3.11, 3.12 and 3.13 for the
following platforms:

- Windows 64-bit (``x86_64``)
- MacOS 64-bit Intel and 64-bit ARM (i.e. Apple Silicon)
- Linux 64-bit (``x86_64``, ``manylinux``)

Notably, at the time of writing the following platforms do *not* have binaries
available:

- Windows on ARM.
- Linux aarch64 (``conda-forge`` has binaries but ``PyPI`` does not, see
  `gh-105 <https://github.com/flintlib/python-flint/issues/105>`_).
- non-glibc Linux distros (``musllinux`` rather than ``manylinux``)

You can see which binaries are provided for the latest release of
``python-flint`` `from PyPI <https://pypi.org/project/python-flint/#files>`_
and which are provided `from conda-forge
<https://anaconda.org/conda-forge/python-flint>`_. A list of supported versions
of Python and other ``python-flint`` dependencies can be found at
:ref:`supported_versions`.

If a binary is not available for your platform then you can build from source
as described in :ref:`build_from_source`.


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
