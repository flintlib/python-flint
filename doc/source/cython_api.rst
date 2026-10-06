Experimental Cython API
=======================

python-flint provides an experimental Cython API for interoperating with
extension modules written in Cython. The API may grow or change between minor
python-flint releases. Within a minor release series, such as ``0.9.x``, it is
expected to remain ABI compatible.

The API currently exposes the ``fmpz`` extension type, FLINT's C ``fmpz``
type as ``flint_fmpz``, and functions for converting in both directions::

    from flint._capi.cython cimport (
        flint_fmpz,
        fmpz,
        fmpz_from_value,
        fmpz_get_value,
    )

``fmpz_get_value(value)`` returns a borrowed, read-only pointer to the
underlying FLINT ``fmpz`` value. Its return type is ``const flint_fmpz *``,
which Cython emits as ``const fmpz *`` and which can be passed directly to
read-only FLINT functions accepting an ``fmpz_t`` argument. The pointer is
valid only while ``value`` is alive. Consumers must not clear or mutate the
value through this pointer, keep the pointer after its owner is released, or
use it without holding the GIL.

``fmpz_from_value(value)`` returns a new Python ``fmpz`` containing a copy of
the supplied FLINT value. The caller retains ownership of the input, which only
needs to remain valid for the duration of the call.

These functions are intended to allow extensions linked against a separate
FLINT library to exchange values with python-flint, including when
python-flint's wheel contains bundled, renamed FLINT and GMP libraries. The
libraries do not need to be the same shared-library instances, but their FLINT
and GMP data ABIs must be compatible. Each library remains responsible for
allocating, mutating and clearing the values that it owns: borrowed values must
only be inspected, and values crossing into python-flint are copied by
``fmpz_from_value``.

See ``examples/cython`` for a complete extension that calls ``fmpz_bits`` on
the borrowed value.
