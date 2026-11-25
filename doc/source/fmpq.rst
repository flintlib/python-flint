**fmpq** -- rational numbers
===============================================================================

An :class:`~flint.fmpq` represents an exact rational number. Its numerator
and denominator are arbitrary-precision integers in reduced form:
:math:`p/q` with :math:`q>0` and :math:`\gcd(p,q)=1`. Zero is represented
by :math:`0/1`.

Construction and arithmetic
---------------------------

For integers ``p`` and nonzero ``q``, ``fmpq(p, q)`` constructs
:math:`p/q` and reduces it automatically. A single integer or a string
such as ``"2/3"`` can also be used:

.. code-block:: python

    from flint import fmpq

    x = fmpq(6, -8)
    y = fmpq("2/3")
    z = fmpq(5)

These values are :math:`x=-3/4`, :math:`y=2/3` and :math:`z=5`.
The methods :meth:`~flint.fmpq.numer` and :meth:`~flint.fmpq.denom`
return the reduced numerator and positive denominator as
:class:`~flint.fmpz` objects. Thus ``x.numer()`` is :math:`-3` and
``x.denom()`` is :math:`4`.

The operators ``+``, ``-``, ``*`` and ``/`` perform exact rational
arithmetic. Integers can be mixed with :class:`~flint.fmpq` values;
division requires a nonzero divisor. For example, with

.. code-block:: python

    a = fmpq(2, 3)
    b = fmpq(5, 6)

the sum, product and quotient are

.. math::

    a+b = \frac{3}{2}, \qquad
    ab = \frac{5}{9}, \qquad
    \frac{a}{b} = \frac{4}{5}.

Harmonic numbers
----------------

For a nonnegative integer :math:`n` and a positive integer :math:`m`,
the generalized harmonic number is the exact rational sum

.. math::

    H_{n,m} = \sum_{k=1}^{n} \frac{1}{k^m}, \qquad H_{0,m} = 0.

The ordinary harmonic number is :math:`H_n = H_{n,1}`.
:meth:`~flint.fmpq.harmonic` computes :math:`H_n` as an :class:`~flint.fmpq`.

Consecutive values
~~~~~~~~~~~~~~~~~~

To compute all values from :math:`H_{0,m}` through :math:`H_{n,m}`, reuse
the preceding partial sum:

.. math::

    H_{k,m} = H_{k-1,m} + \frac{1}{k^m} \qquad (k \geq 1).

For nonnegative integer ``n`` and positive integer ``m``, the following
function returns these :math:`n+1` values in order:

.. code-block:: python

    from flint import fmpq

    def harmonic_prefix(n, m=1):
        values = [fmpq(0)]
        for k in range(1, n + 1):
            values.append(values[-1] + fmpq(1, k**m))
        return values

The constructor ``fmpq(1, k**m)`` forms each summand exactly. After
iteration :math:`k`, the new last entry is :math:`H_{k,m}`, since the
update adds precisely the next term of the defining sum.

For example, when :math:`m=2`, the first five entries are

.. math::

    (H_{0,2}, H_{1,2}, H_{2,2}, H_{3,2}, H_{4,2})
    = \left(0, 1, \frac{5}{4}, \frac{49}{36}, \frac{205}{144}\right).

Hurwitz zeta
~~~~~~~~~~~~

For an integer :math:`m>1`, splitting the convergent defining series of
the zeta function after its first :math:`n` terms gives

.. math::

    H_{n,m} = \zeta(m) - \zeta(m,n+1).

This is the Hurwitz zeta shift relation in
`DLMF 25.11.4 <https://dlmf.nist.gov/25.11.E4>`_ with :math:`a=1`.
With ``arb`` imported from ``flint``, the expression
``arb(m).zeta() - arb(m).zeta(n + 1)`` evaluates this difference using
:meth:`~flint.arb.zeta` in real-ball arithmetic.
The recurrence above computes an exact rational value.

At :math:`m=1`, each zeta term has a pole. Use
:meth:`~flint.fmpq.harmonic` or the recurrence for ordinary harmonic
numbers.

API reference
-------------

.. autoclass :: flint.fmpq
  :members:
  :inherited-members:
  :undoc-members:

