Experimental C API
==================

python-flint provides an experimental version 1 C API for ``fmpz``. The API
may change between releases and is not yet a long-term ABI stability promise.

Cython consumers
----------------

Cython 3.1 or newer can use the installed declarations directly::

    from flint cimport fmpz, fmpz_add, fmpz_get_value

``fmpz_add(a, b)`` returns a new ``fmpz`` and operates entirely through the
python-flint capsule. Extensions using only this operation do not need FLINT
headers or a direct FLINT link. See ``examples/cython/capi``.

``fmpz_get_value(value)`` returns a read-only, borrowed pointer to FLINT
storage. It remains valid only while ``value`` is alive. Never clear or mutate
it, retain it independently of its owner, or use it after the owner is
released. Direct consumers must compile and link against the same compatible
FLINT ABI/build as python-flint. In particular, this is not safe with wheels
containing a private, renamed FLINT. See ``examples/cython/flint``.

C consumers
-----------

Include the installed header and import the capsule before using it::

    #include <python_flint/fmpz.h>

    if (PyFlint_Import() < 0)
        return NULL;
    result = PyFlint_FMPZ_Add(left, right);  /* new reference */

``PyFlint_FMPZ_GetValue`` has the same borrowed, read-only ownership rules as
the Cython wrapper. Both operations return ``NULL`` with a Python exception set
on failure, and reject objects that are not ``flint.fmpz`` instances.

Header location
---------------

Build systems can locate ``python_flint/fmpz.h`` using::

    import flint
    print(flint.get_include())

The public header and Cython declarations contain no details of python-flint's
private extension-object layout.
