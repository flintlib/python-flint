General concepts
===============================================================================

Types and methods
-----------------

(Preliminary content.)

As a general rule, C functions associated with a type in FLINT or Arb
are exposed as methods of the corresponding Python type.

For example, there is both an :meth:`.fmpq.bernoulli_ui` (which computes
a Bernoulli number as an exact fraction) and :meth:`.arb.bernoulli_ui`
(which computes a Bernoulli number as an approximate real number).

A function that transforms a single value to the same type
is usually an ordinary method of that type, for instance :meth:`.arb.exp`.
A function with a different signature can either provided as a classmethod
that takes all inputs as function arguments, or as a
method of the "primary" input, taking the other inputs
as arguments to the method (for example :meth:`.arb.rising_ui`).

When a method involves different types for inputs and outputs (or
just among the inputs), it will
typically be a method of the more "complex" type. For example, a matrix
type is more "complex" than the underlying scalar type, so
:meth:`.fmpz_mat.det` is a method of the matrix type, returning a scalar,
and not vice versa.
If the input type is not the same as the type that the method belongs
to, it will typically be a classmethod (for example
:meth:`.fmpq.bernoulli_ui`).

The method-based interface is intended to be precise and to simplify the implementation
somewhat, not primarily to be aesthetically pleasing. A functional
top-level interface might be added in the future, allowing more idiomatic
mathematical notation (for example, :func:`exp` and
:func:`det` as regular functions).
Such an interface will presumably just involve type-dispatch glue, leaving
all the interesting computational logic to the methods on the
various mathematical types.

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

    def mutated_thing():   # ok
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
Eventually, all such bad input should be caught at the Python level
and result in appropriate exceptions being raised. But until the code
matures, assume that such checks have not yet been implemented and that
invalid input leads to undefined behavior!

Inexact numbers
-----------------

Real and complex numbers are represented by intervals. All operations on
real and complex numbers output intervals representing rigorous error bounds.
This also extends to polynomials and matrices of real and complex numbers.

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
you want the correct decimal fraction. But small integers and
power-of-two denominators are safe, for example ``arb(100.25)``.

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
