General concepts
===============================================================================

Importing
-----------------

The ``flint`` module exposes a set of distinctly-named types together
with a small number of top-level functions and objects.
Most functionality is provided as methods on the types. This means
that there should be no namespace conflicts with most user code,
with Python's builtin ``math`` and ``cmath`` modules, or with
packages such as ``gmpy``, ``numpy``, ``sympy`` and ``mpmath``.
For typical interactive use, it should therefore
generally be safe to ``import *``:

    >>> from flint import *
    >>> fmpq(3) / 2
    3/2

For non-interactive use, it is still good manners to use explicit
imports or preserve the ``flint`` namespace prefix::

    >>> import flint
    >>> flint.fmpq(3) / 2
    3/2

Global context
-----------------

Various settings are controlled by a global context object,
``flint.ctx``. Printing this object in the REPL shows the current
settings, with a brief explanation of each parameter::

    >>> from flint import ctx
    >>> ctx
    pretty = True      # pretty-print repr() output
    unicode = False    # use unicode characters in output
    prec = 53          # real/complex precision (in bits)
    dps = 15           # real/complex precision (in digits)
    cap = 10           # power series precision
    threads = 1        # max number of threads used internally

The user can mutate the properties directly, for example::

    >>> ctx.pretty = False
    >>> fmpq(3,2)
    fmpq(3,2)
    >>> ctx.pretty = True
    3/2

Calling ``ctx.default()`` restores the default settings.

The special method ``ctx.cleanup()`` frees up internal caches
used by MPFR, FLINT and Arb. The user does normally not have to
worry about this.

Types and methods
-----------------

As a general rule, C functions associated with a type in FLINT or Arb
are exposed as methods of the corresponding Python type.

For example, there is both an :meth:`.fmpq.bernoulli` (which computes
a Bernoulli number as an exact fraction) and :meth:`.arb.bernoulli`
(which computes a Bernoulli number as an approximate real number).

A function that transforms a single value to the same type
is usually an ordinary method of that type, for instance :meth:`.arb.exp`.
A function with a different signature can either provided as a
static method that takes all inputs as function arguments, or as a
method of the "primary" input, taking the other inputs
as arguments to the method (for example :meth:`.arb.bessel_j`).

When a method involves different types for inputs and outputs (or
just among the inputs), it will
typically be a method of the more "complex" type. For example, a matrix
type is more "complex" than the underlying scalar type, so
:meth:`.fmpz_mat.det` is a method of the matrix type, returning a scalar,
and not vice versa.

The method-based interface is intended to keep the code simple,
not to be aesthetically pleasing to mathematicians. A functional
top-level interface might be added in the future, allowing more idiomatic
mathematical notation (for example, :func:`exp` and
:func:`det` as regular functions).

Mutability
----------

Objects have immutable semantics. For example, the second line in::

    b = a
    a += c

leaves *b* unchanged.

However, mutation via direct element access is supported for matrices
and polynomials. Some methods also allow explicitly performing the
operation in-place. Civilized users will restrict their use of such
methods to the point in the code where the object is first constructed::

    def create_thing():   # ok
        a = thing()
        a.mutate()
        return a

Crashing and burning
---------------------------------------

Very little overflow checking is done ahead-of-time. Trying to compute an
object far too large to hold in memory (for example, the exact factorial
of `2^{64}-1`) will likely abort the process,
instead of raising an :exc:`OverflowError` or :exc:`MemoryError` that
can be caught at the Python level.

Input that is obviously *invalid* (for example a negative number passed
as a length) can also cause crashes or worse things to happen.
Ideally, bad input should be caught at the Python level and result in
appropriate exceptions being raised, but this is not yet done
systematically. At this time, users should assume that invalid
input leads to undefined behavior!

Inexact numbers and numerical evaluation
-----------------------------------------------------------------------

Real and complex numbers are represented by midpoint-radius intervals
(balls). All operations on real and complex numbers output intervals
representing rigorous error bounds. This also extends to polynomials
and matrices of real and complex numbers.

The working precision for real and complex arithmetic is controlled by the
global context object attributes :func:`ctx.prec` (in bits)
:func:`ctx.dps` (in decimal digits). Changing either attribute changes
the other to match.

Be careful about using Python float and complex literals as input.
Doing ``arb(0.1)`` actually gives an interval containing
the rational number

.. math ::

    3602879701896397 \times 2^{-55} = 0.1000000000000000055511151231257827021181583404541015625

which might not be what you want. Do ``arb("0.1")``, ``arb("1/10")``
or ``arb(fmpq(1,10))`` if
you want the correct decimal fraction. Small integers and
power-of-two denominators are still safe, for example ``arb(100.25)``.

Pointwise boolean predicates (such as the usual comparison operators)
involving inexact numbers return
*True* only if the predicate certainly is true (i.e. it holds for all
combinations of points that can be chosen from the set-valued inputs),
and return *False* if the
predicate either definitely is false or the truth cannot be determined.
To determine that a predicate is definitely false,
test both the predicate and the inverse predicate,
e.g. if either ``x < y`` or ``y <= x`` returns *True*, then the other
is definitely false; if both return *False*, then neither can be
determined from the available data.

The following convenience functions are provided for numerical evaluation
with adaptive working precision.

.. autofunction :: flint.good

.. autofunction :: flint.showgood

Power series
-----------------------------------------------------------------------

Power series objects track the precision (the number of known terms)
automatically.  The upper precision for power series is controlled by
``flint.ctx.cap``, with the default value 10.

    >>> fmpq_series([0,1]).exp()
    1 + x + 1/2*x^2 + 1/6*x^3 + 1/24*x^4 + 1/120*x^5 + 1/720*x^6 + 1/5040*x^7 + 1/40320*x^8 + 1/362880*x^9 + O(x^10)
    >>> ctx.cap = 4
    >>> fmpq_series([0,1]).exp()
    1 + x + 1/2*x^2 + 1/6*x^3 + O(x^4)
    >>> ctx.cap = 10
    >>> fmpq_series([0,1], prec=5).exp()
    1 + x + 1/2*x^2 + 1/6*x^3 + 1/24*x^4 + O(x^5)

    >>> ctx.cap = 3
    >>> ctx.dps = 10
    >>> arb_series([1,3,4]).exp()
    ([2.718281828 +/- 4.79e-10]) + ([8.154845485 +/- 4.36e-10])*x + ([23.10539554 +/- 2.25e-9])*x^2 + O(x^3)
    >>> ctx.default()
