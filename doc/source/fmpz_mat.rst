**fmpz_mat** -- matrices over integers
===============================================================================

Rational nullspaces and integer kernels
---------------------------------------

For an integer matrix :math:`A` with :math:`m` rows and :math:`n` columns,
:meth:`~flint.fmpz_mat.nullspace` computes a basis of the right nullspace

.. math::

    \ker_{\mathbb{Q}}(A) = \{x \in \mathbb{Q}^n : Ax = 0\}.

It returns ``(X, nullity)``, where ``X`` is an :math:`n \times n` integer
matrix whose first ``nullity`` columns form this basis. The remaining
columns are zero. The basis vectors span the nullspace over
:math:`\mathbb{Q}`; they need not generate every integer solution by
integer linear combinations.

For example, the matrix in the :meth:`~flint.fmpz_mat.nullspace` example is

.. math::

    A = \begin{pmatrix}
        1 & 2 & 3 & 4 & 5 \\
        6 & 7 & 8 & 9 & 10 \\
        11 & 12 & 13 & 14 & 15
    \end{pmatrix}.

The vector :math:`v = (1,-2,1,0,0)^T` satisfies :math:`Av = 0`, and the
first displayed column of ``X`` is :math:`5v`. Every entry of that
``X`` is divisible by five, so :math:`v` is not an integer linear
combination of its columns.

To obtain a basis of the integer kernel
:math:`\{x \in \mathbb{Z}^n : Ax = 0\}`, use the Hermite normal form of
:math:`A^T` with a transformation matrix:

.. code-block:: python

    H, U = A.transpose().hnf(transform=True)

Here :meth:`~flint.fmpz_mat.hnf` gives :math:`H = U A^T`, where
:math:`U` is unimodular. If :math:`r` is the rank of :math:`A`, the last
:math:`n-r` rows of :math:`H` are zero. The transposes of the corresponding
rows of :math:`U` form a basis of the integer kernel.

Indeed, since :math:`U` is unimodular, every integer vector can be written
uniquely as :math:`x = U^T y` with :math:`y` integral. Then :math:`Ax = 0`
is equivalent to :math:`H^T y = 0`. The first :math:`r` rows of :math:`H`
are linearly independent, so this holds exactly when the first :math:`r`
coordinates of :math:`y` are zero. This proves both generation and
independence of the stated integer kernel basis.

LLL reduction preserves the integer span of a basis; applying it to a
rational nullspace basis cannot recover integer vectors missing from that
span.

API reference
-------------

.. autoclass :: flint.fmpz_mat
  :members:
  :inherited-members:
  :undoc-members:

