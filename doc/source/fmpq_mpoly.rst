**fmpq_mpoly** -- multivariate polynomials over the rational numbers
===============================================================================

A context fixes the variables and monomial ordering of a polynomial ring
over :math:`\mathbb{Q}`. For nonnegative integers :math:`i` and :math:`j`,
a dictionary entry ``(i, j): c`` represents the term :math:`c x^i y^j`
when the context's variables are ``x`` and ``y``, in that order.

For example, the polynomial

.. math::

    p(x, y) = \frac{1}{3} x y + \frac{1}{2} x + \frac{3}{4} y

can be constructed from the context's generators::

    >>> from flint import fmpq, fmpq_mpoly_ctx
    >>> ctx = fmpq_mpoly_ctx.get(("x", "y"), "lex")
    >>> x, y = ctx.gens()
    >>> p = x*y/3 + x/2 + 3*y/4
    >>> p
    1/3*x*y + 1/2*x + 3/4*y

Division by a nonzero integer gives exact rational coefficients.
The same polynomial can also be constructed from a dictionary::

    >>> p = ctx.from_dict({
    ...     (1, 1): fmpq(1, 3),
    ...     (1, 0): fmpq(1, 2),
    ...     (0, 1): fmpq(3, 4),
    ... })

For integers :math:`a` and :math:`b \ne 0`, ``fmpq(a, b)`` represents
:math:`a/b` exactly. Integer coefficients are also allowed.

.. autoclass :: flint.fmpq_mpoly_ctx
  :members:
  :inherited-members:
  :undoc-members:

.. autoclass :: flint.fmpq_mpoly
  :members:
  :inherited-members:
  :undoc-members:

.. autoclass :: flint.fmpq_mpoly_vec
  :members:
  :inherited-members:
  :undoc-members:

